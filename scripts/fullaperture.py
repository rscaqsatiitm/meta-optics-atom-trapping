"""Full-aperture simulation of a centred square-array meta-hologram (v1: 48 x 48, v2: ~20,000 sites).

Usage: python scripts/fullaperture.py configs/fullaperture_v1.yaml

Steps: (1) D4-symmetric exact-kernel WGS in the pupil domain; (2) symmetric fine render
grid; (3) the full metasurface (pixel = meta-atom pitch, ~1.4e9 pixels at 18.4 mm) is
generated block by block, optionally with library pillars and width errors, and folded
onto an n_fold^2 focal window (propagation/folded.py); (4) spot metrics for every site;
(5) optional feedback rounds: target weights rescaled by (orbit-mean peak / peak)^gain
and the WGS warm-started. Writes metrics.json and per-site npz.
"""

import gc
import json
import sys
import time

import numpy as np

from metatrap import materials
from metatrap.atoms.power_budget import PowerBudget
from metatrap.hologram import evaluate as E
from metatrap.hologram import symmetric as S
from metatrap.hologram import targets
from metatrap.hologram.exact import OPLTable
from metatrap.hologram.render import fourier_upsample, metasurface_blocks
from metatrap.hologram.wgs import PupilGrid
from metatrap.propagation.folded import fold_blocks, propagate_folded
from metatrap.propagation.layers import LayerStack
from metatrap.rcwa import lookup as L
from metatrap.rcwa.library import Library
from metatrap.runio import REPO_ROOT, load_config, run_dir, write_json, write_params


def build_lookup(cfg, lam):
    lc = cfg.get("library")
    if not lc:
        return None
    lib = Library.load(REPO_ROOT / lc["path"])
    w, t, _ = lib.curve(lc["material"], lc["pitch_um"], lc["height_um"], lam)
    im, ip, ih, il = lib.index(lc["material"], lc["pitch_um"], lc["height_um"], lam)
    n_in = float(materials.get(lib.attrs["substrate"]).n(lam))
    lim = L.FabLimits()
    ok = L.feasible(w, lc["pitch_um"], lc["height_um"], lim)
    if lib.t_cross is not None:
        ok &= lib.t_cross[im, ip, ih, il] <= lim.max_cross
    return L.build(w[ok], t[ok], lc["pitch_um"], lc["height_um"], n_in, 1.0, L.FabLimits(min_feature_um=0, min_gap_um=0, max_aspect=1e9))


def simulate(cfg, stack, lam, na, ill, sites, g, fs, fstab, inv, tab, d, lookup, errors, nf, n, dx):
    holo = S.render_grid_symmetric(g, stack, na, ill, fs, fstab, d.coeffs, tab, lam, upsample=cfg.get("render_upsample", 8))
    t0 = time.time()
    blocks = metasurface_blocks(n, dx, stack, lam, na, None, None, ill, holo=holo, lookup=lookup, pillar_errors=errors, block=cfg.get("block_rows", 128), with_incident=True)
    F, p_pix, p_inc_pix = fold_blocks(blocks, n, nf)
    del holo
    gc.collect()
    foc = propagate_folded(F, dx, lam, stack)
    del F
    half = np.abs(sites).max() + 4 * cfg["site_pitch_um"]
    w = int(np.ceil(2 * half / dx / 2) * 2)
    c = nf // 2 - w // 2
    win = foc[c : c + w, c : c + w]
    I = np.abs(fourier_upsample(win, cfg.get("upsample", 4))) ** 2
    dxu = dx / cfg.get("upsample", 4)
    p_in = p_inc_pix * dx * dx  # efficiencies relative to the light incident on the metasurface
    m = E.spot_metrics(I, dxu, sites, cfg["site_pitch_um"], lam, na, p_in)
    m["metasurface_transmission"] = float(p_pix / p_inc_pix)
    per = m.pop("per_site")
    m["axis_intensity_over_mean_site_peak"] = float(I[I.shape[0] // 2, I.shape[1] // 2] / per["peak"].mean())
    # Stray light: power in the array region outside the site cells, relative to input.
    m["seconds"] = time.time() - t0
    return m, per, p_in


def main(cfg_path):
    cfg = load_config(cfg_path)
    out = run_dir(cfg["run_id"])
    write_params(out, cfg)
    lam, na, pitch, dx = cfg["wavelength_um"], cfg["na"], cfg["site_pitch_um"], cfg["pixel_um"]
    stack = LayerStack.from_list(cfg["stack"])
    radius = stack.aperture_for_na(na, lam)
    ill = dict(cfg["illumination"])
    if "edge_frac" in ill:
        ill["edge_um"] = ill.pop("edge_frac") * radius
    if ill["kind"] == "gaussian" and "waist_over_radius" in ill:
        ill["waist_um"] = ill.pop("waist_over_radius") * radius
    sites = targets.make(**cfg["sites"])
    fs, fstab, inv = S.wedge_sites(sites)
    rmax = np.abs(sites).max() * np.sqrt(2)
    g = PupilGrid.for_sites(lam, pitch, 4, na, 2 * rmax, guard=cfg.get("pupil_guard", 1.25))
    tab = OPLTable(stack, lam, radius + 2 * rmax + 20)
    nf = cfg["n_fold"]
    n = int(np.ceil((2 * radius / dx + 2) / nf) * nf)
    print(f"aperture {2*radius/1000:.2f} mm, metasurface {n}^2 = {n*n/1e9:.2f} Gpx, {len(sites)} sites ({len(fs)} fundamental), "
          f"pupil n={g.n} ds={g.ds:.2e} s_max={g.s_max:.3f}", flush=True)

    # Checkpointing: design state after every feedback round, metrics after every step,
    # so a run longer than one background session resumes where it stopped.
    state_path, metrics_path = out / "state.npz", out / "metrics.json"
    results = json.load(open(metrics_path)) if metrics_path.exists() else {}
    results.setdefault("runs", {})
    results.setdefault("feedback_history", [])
    n_rounds = cfg.get("feedback_rounds", 0)
    state = dict(np.load(state_path)) if state_path.exists() else None
    feedback_done = bool(state is not None and int(state["next_round"]) > n_rounds)

    wp = K = None
    if not feedback_done:
        t0 = time.time()
        wp = S.WedgePupil.build(g, stack, na, ill)
        K = S.kernels_G(wp, fs, tab, lam)
        print(f"K_G {K.shape} {K.nbytes/1e9:.2f} GB in {time.time()-t0:.0f}s", flush=True)
    if state is None:
        d = S.wgs_symmetric(wp, K, fstab, n_iter=cfg.get("wgs_iter", 60), n_fix=cfg.get("wgs_fix", 30))
        mult = 8.0 / fstab
        eta = np.abs(d.mode_amp) ** 2
        mean = np.sum(mult * eta) / mult.sum()
        results["pupil_model"] = {"mode_efficiency_sum": float(np.sum(mult * eta)),
                                  "nonuniformity_rms": float(np.sqrt(np.sum(mult * (eta - mean) ** 2) / mult.sum()) / mean)}
        print("pupil model", results["pupil_model"], flush=True)
        tgt, rnd0 = np.ones(len(fs)), 0
    else:
        d = S.SymWGSResult(psi=state["psi"], coeffs=state["coeffs"], mode_amp=np.zeros(len(fs)), history_nonuniformity=[], history_efficiency=[])
        tgt, rnd0 = state["tgt"], int(state["next_round"])
        print(f"resumed from checkpoint at round {rnd0}", flush=True)
    results.update({"aperture_mm": 2 * radius / 1000, "n_pixels": n * n, "n_sites": len(sites)})

    def save_state(next_round):
        np.savez_compressed(state_path, psi=d.psi, coeffs=d.coeffs, tgt=tgt, next_round=next_round)

    if not feedback_done:
        for rnd in range(rnd0, n_rounds + 1):
            m, per, p_in = simulate(cfg, stack, lam, na, ill, sites, g, fs, fstab, inv, tab, d, None, None, nf, n, dx)
            h = {"round": rnd, **{k: m[k] for k in ("depth_nonuniformity_rms", "depth_min_over_mean", "array_efficiency", "core_efficiency", "axis_intensity_over_mean_site_peak")}}
            results["feedback_history"] = [x for x in results["feedback_history"] if x["round"] != rnd] + [h]
            print("ideal pillars, round", rnd, h, f"{m['seconds']:.0f}s", flush=True)
            if rnd == n_rounds:
                results["runs"]["ideal"] = m
                np.savez_compressed(out / "per_site_ideal.npz", sites=sites, **per)
                save_state(n_rounds + 1)
                write_json(metrics_path, results)
                break
            pk = per["peak"]
            orbit_mean = np.bincount(inv, weights=pk) / np.bincount(inv)
            tgt = tgt * (pk.mean() / orbit_mean) ** cfg.get("feedback_gain", 0.5)
            d = S.wgs_symmetric(wp, K, fstab, n_iter=cfg.get("feedback_wgs_iter", 30), n_fix=None, psi0=d.psi, target=tgt)
            save_state(rnd + 1)
            write_json(metrics_path, results)
    del K
    gc.collect()

    lookup = build_lookup(cfg, lam)
    for case in cfg.get("pillar_cases", []):
        if lookup is None:
            break
        if case["name"] in results["runs"]:
            print("skip (done)", case["name"], flush=True)
            continue
        m, per, _ = simulate(cfg, stack, lam, na, ill, sites, g, fs, fstab, inv, tab, d, lookup, case.get("errors"), nf, n, dx)
        results["runs"][case["name"]] = m
        print(case["name"], {k: round(v, 4) for k, v in m.items() if isinstance(v, float)}, flush=True)
        np.savez_compressed(out / f"per_site_{case['name']}.npz", sites=sites, **per)
        write_json(metrics_path, results)

    # Power budget from the simulated efficiency (ideal or library pillars), with stated extra factors.
    pb_cfg = cfg.get("power_budget")
    if pb_cfg:
        rows = {}
        for name, m in results["runs"].items():
            util = m["core_efficiency"] / pb_cfg["core_fraction_ideal"] * np.prod(list(pb_cfg["extra_factors"].values()))
            pb = PowerBudget(lam_nm=lam * 1e3, utilisation=util, na_eff=pb_cfg["na_eff"])
            rows[name] = {"utilisation": util, **{f"depth_mK_at_{p}W": float(pb.depth_K(p, len(sites)) * 1e3) for p in pb_cfg["powers_W"]}}
        results["power_budget"] = rows
        print("power budget", rows, flush=True)
    write_json(metrics_path, results)


if __name__ == "__main__":
    main(sys.argv[1])
