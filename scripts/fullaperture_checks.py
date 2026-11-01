"""Pre-fabrication checks on a full-aperture design (loads the design state of a fullaperture run).

Usage: python scripts/fullaperture_checks.py configs/checks_v2.yaml

Sections of the config (each optional, results checkpointed in results/<run_id>/checks.json):
  axis_light:   for each listed base run, render with ideal pillars and record the intensity
                on the axis and at every interstitial point (centre of four sites), relative to
                the mean site peak.
  cell:         wall-thickness errors and wall tilts, evaluated by re-propagating one folded
                render of the base run (ideal pillars), each at its best focus (scan) and,
                for tilts, re-centred on the laterally shifted array.
  fabrication:  library pillars with correlated errors (write-field stitching, per-field width
                bias, smooth width noise) and/or uncorrelated noise; one render per case.
"""

import gc
import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np

from metatrap.hologram import symmetric as S
from metatrap.hologram import targets
from metatrap.hologram.exact import OPLTable
from metatrap.hologram.wgs import PupilGrid
from metatrap.propagation.angular_spectrum import transfer_function
from metatrap.propagation.folded import transfer_function_tilted
from metatrap.propagation.layers import LayerStack
from metatrap.runio import REPO_ROOT, load_config, run_dir, write_json, write_params

_spec = importlib.util.spec_from_file_location("fa", Path(__file__).with_name("fullaperture.py"))
FA = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(FA)


def setup(base_cfg_path):
    cfg = load_config(REPO_ROOT / base_cfg_path)
    lam, na, pitch, dx = cfg["wavelength_um"], cfg["na"], cfg["site_pitch_um"], cfg["pixel_um"]
    stack = LayerStack.from_list(cfg["stack"])
    radius = stack.aperture_for_na(na, lam)
    ill = dict(cfg["illumination"])
    if "edge_frac" in ill:
        ill["edge_um"] = ill.pop("edge_frac") * radius
    dz = float(cfg.get("array_defocus_um", 0.0))
    layers = [list(l) for l in cfg["stack"]]
    layers[-1][0] = float(layers[-1][0]) + dz
    stack_sites = LayerStack.from_list(layers)
    sites = targets.make(**cfg["sites"])
    fs, fstab, inv = S.wedge_sites(sites)
    rmax = np.abs(sites).max() * np.sqrt(2)
    g = PupilGrid.for_sites(lam, pitch, 4, na, 2 * rmax, guard=cfg.get("pupil_guard", 1.25))
    tab = OPLTable(stack_sites, lam, radius + 2 * rmax + 20)
    nf = cfg["n_fold"]
    n = int(np.ceil((2 * radius / dx + 2) / nf) * nf)
    st = np.load(REPO_ROOT / "results" / cfg["run_id"] / "state.npz")
    d = S.SymWGSResult(psi=st["psi"], coeffs=st["coeffs"], mode_amp=np.zeros(len(fs)), history_nonuniformity=[], history_efficiency=[])
    return dict(cfg=cfg, lam=lam, na=na, dx=dx, stack=stack, stack_sites=stack_sites, ill=ill, sites=sites, fs=fs,
                fstab=fstab, g=g, tab=tab, nf=nf, n=n, d=d, layers=layers)


def summary(m):
    keys = ("depth_nonuniformity_rms", "depth_min_over_mean", "array_efficiency", "core_efficiency", "strehl_mean",
            "waist_mean_um", "position_error_rms_um", "axis_intensity_over_mean_site_peak", "interstitial_median",
            "interstitial_p99", "interstitial_max", "interstitial_max_radius_um", "mean_site_offset_um")
    return {k: m[k] for k in keys if k in m}


def focus_scan(G, F, p_in, H0, site_shift, centre_um, coarse=(-3, -2, -1, 0, 1, 2, 3), fine=(-0.5, -0.25, 0.25, 0.5)):
    """Best plane (max mean site peak) near centre_um beyond the design plane; returns (dz, metrics)."""
    lam, dx, nf = G["lam"], G["dx"], G["nf"]
    k0 = 2 * np.pi / lam
    f = np.fft.fftfreq(nf, d=dx)
    SX, SY = np.meshgrid(lam * f, lam * f, indexing="ij")
    sz = np.sqrt(np.clip(1 - SX**2 - SY**2, 0, None))

    def ev(dz):
        H = H0 * np.exp(1j * k0 * dz * sz)
        m, per, _ = FA.focal_metrics(G["cfg"], F, dx, lam, G["na"], G["sites"], p_in, H=H, site_shift=site_shift)
        return float(per["peak"].mean()), m

    best = None
    for dz in [centre_um + c for c in coarse]:
        val, m = ev(dz)
        if best is None or val > best[0]:
            best = (val, dz, m)
    c0 = best[1]
    for dz in [c0 + c for c in fine]:
        val, m = ev(dz)
        if val > best[0]:
            best = (val, dz, m)
    return best[1], best[2], best[0]


def main(cfg_path):
    cfg = load_config(cfg_path)
    out = run_dir(cfg["run_id"])
    write_params(out, cfg)
    res_path = out / "checks.json"
    res = json.load(open(res_path)) if res_path.exists() else {}

    def save():
        write_json(res_path, res)

    # 1. Axis light against interstitial background.
    for base in cfg.get("axis_light", {}).get("bases", []):
        key = f"axis_light:{base}"
        if key in res:
            continue
        G = setup(base)
        t0 = time.time()
        F, p_pix, p_inc = FA.render_fold(G["cfg"], G["stack"], G["lam"], G["na"], G["ill"], G["g"], G["fs"], G["fstab"], G["tab"], G["d"], None, None, G["nf"], G["n"], G["dx"])
        m, per, I = FA.focal_metrics(G["cfg"], F, G["dx"], G["lam"], G["na"], G["sites"], p_inc * G["dx"] ** 2, stack_eval=G["stack_sites"])
        c = I.shape[0] // 2
        r = int(round(6.0 / (G["dx"] / G["cfg"].get("upsample", 4))))
        np.save(out / f"axis_patch_{Path(base).stem}.npy", I[c - r : c + r + 1, c - r : c + r + 1] / per["peak"].mean())
        res[key] = summary(m) | {"seconds": time.time() - t0}
        print(key, res[key], flush=True)
        save()
        del F, I
        gc.collect()

    # 2. Cell-wall tolerances (one render, many propagations).
    cc = cfg.get("cell")
    if cc and not all(f"cell:{c['name']}" in res for c in cc["cases"]):
        G = setup(cc["base"])
        F, p_pix, p_inc = FA.render_fold(G["cfg"], G["stack"], G["lam"], G["na"], G["ill"], G["g"], G["fs"], G["fstab"], G["tab"], G["d"], None, None, G["nf"], G["n"], G["dx"])
        p_in = p_inc * G["dx"] ** 2
        wall = cc.get("wall_layer", 1)
        n_wall = G["stack"].layers[wall].n(G["lam"])
        for case in cc["cases"]:
            key = f"cell:{case['name']}"
            if key in res:
                continue
            t0 = time.time()
            layers = [list(l) for l in G["layers"]]
            dt = case.get("dthick_um", 0.0)
            layers[wall][0] = float(layers[wall][0]) + dt
            st_eval = LayerStack.from_list(layers)
            beta = np.deg2rad(case.get("tilt_deg", 0.0))
            if beta == 0:
                H0 = transfer_function(G["nf"], G["dx"], G["lam"], st_eval).astype(np.complex128)
                H0 *= np.exp(-1j * np.angle(H0[0, 0]))
                shift = (0.0, 0.0)
            else:
                H0 = transfer_function_tilted(G["nf"], G["dx"], G["lam"], st_eval, wall, beta)
                d_w = layers[wall][0]
                sx = d_w * np.sin(beta) * (1 - np.cos(beta) / np.sqrt(n_wall**2 - np.sin(beta) ** 2))
                shift = (float(sx), 0.0)
            centre = dt * (1 - 1 / n_wall)  # paraxial focal shift of a thicker plate
            dz, m, _ = focus_scan(G, F, p_in, H0, shift, centre)
            if beta != 0:  # refine the lateral re-centring once
                off = m["mean_site_offset_um"]
                shift = (shift[0] + off[0], shift[1] + off[1])
                dz, m, _ = focus_scan(G, F, p_in, H0, shift, dz, coarse=(-0.5, 0, 0.5), fine=(-0.25, 0.25))
            res[key] = summary(m) | {"best_focus_um": dz, "paraxial_focus_um": centre, "array_shift_um": shift, "case": case, "seconds": time.time() - t0}
            print(key, {k: (round(v, 4) if isinstance(v, float) else v) for k, v in res[key].items()}, flush=True)
            save()
        del F
        gc.collect()

    # 3. Fabrication errors with library pillars.
    fc = cfg.get("fabrication")
    if fc:
        G = None
        for case in fc["cases"]:
            key = f"fab:{case['name']}"
            if key in res:
                continue
            if G is None:
                G = setup(fc["base"])
                lookup = FA.build_lookup(G["cfg"], G["lam"])
            t0 = time.time()
            F, p_pix, p_inc = FA.render_fold(G["cfg"], G["stack"], G["lam"], G["na"], G["ill"], G["g"], G["fs"], G["fstab"], G["tab"], G["d"], lookup, case.get("errors"), G["nf"], G["n"], G["dx"])
            m, per, I = FA.focal_metrics(G["cfg"], F, G["dx"], G["lam"], G["na"], G["sites"], p_inc * G["dx"] ** 2, stack_eval=G["stack_sites"])
            res[key] = summary(m) | {"metasurface_transmission": float(p_pix / p_inc), "case": case, "seconds": time.time() - t0}
            print(key, {k: (round(v, 4) if isinstance(v, float) else v) for k, v in res[key].items()}, flush=True)
            save()
            del F, I
            gc.collect()


if __name__ == "__main__":
    main(sys.argv[1])
