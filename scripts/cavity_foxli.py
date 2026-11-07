"""S3: flat-top sheet mode of the self-intersecting ring cavity, its spectrum and tolerances.
Usage:
    python scripts/cavity_foxli.py configs/cavity_foxli.yaml

Geometry variants (config 'variants'):
  plane_fold: X -> plane fold at d1 -> graded mirror G at path p (symmetric half trip).
  telescope:  X -> d0 -> cylindrical mirror f1 -> f1 + f2 -> cylindrical mirror f2 -> L -> G,
              f2 = mag x f1 (Keplerian) or -mag x f1 (Galilean, separation f1 - |f2|).
For each case the laser is locked to the fundamental (maximum circulating power) and the
steady-state circulating field at the crossing is computed for the given input. Metrics:
enhancement G relative to the analytic value for the lumped losses, flatness of the circulating
intensity over the array's transverse range, and per-site barrier depth and site displacement
over the 10^6-site array (both crossing beams, sheet thickness w_z).
Writes results/<run_id>/foxli.json (resumable) and one npz per design.
"""

import copy
import json
import sys
import time

import numpy as np

from metatrap.cavity import budget as B
from metatrap.cavity import foxli as FL
from metatrap.cavity import lattice as LA
from metatrap.runio import load_config, run_dir, write_json, write_params


def build(cfg, var, n=None, width=None):
    """Grid, target and ideal ring of a variant."""
    lam = cfg["wavelength_um"]
    theta = B.crossing_half_angle(cfg["period_um"], lam)
    h = LA.array_halfrange(cfg["array_side_um"], theta)
    w = LA.supergauss_w_for_flat(h, var["order"], cfg["flat_tol"])
    g = FL.Grid1D(n or cfg["grid"]["n"], width or cfg["grid"]["width_um"], lam)
    F = FL.flat_target(g, w, var["order"])
    z = np.zeros(g.n)
    if var["geometry"] == "plane_fold":
        p = var["p_um"]
        ring = FL.RingStrip(g, p, var["d1_um"], z, z, ap_g_um=var["ap_g_um"], ap_fold_um=var["ap_fold_um"])
    elif var["geometry"] == "telescope":
        f1, mag = var["f1_um"], var["mag"]
        f2 = -mag * f1 if var.get("galilean") else mag * f1
        els = FL.telescope_elements(g, f1, f2, var["d0_um"], var["ap1_um"], var["ap2_um"])
        p = var["d0_um"] + f1 + f2 + var["L_um"]
        ring = FL.RingStrip(g, p, 0.0, z, z, ap_g_um=var["ap_g_um"], elements=els)
    else:
        raise ValueError(var["geometry"])
    ph, ok = FL.graded_phase(g, F, p, ring=ring)
    ring.phase_g1, ring.phase_g2 = ph, ph
    return {"lam": lam, "theta": theta, "h": h, "w": w, "g": g, "F": F, "ring": ring, "ph": ph, "ok": ok, "p": p}


def with_phases(ring, ph1, ph2):
    r = copy.copy(ring)
    r.phase_g1, r.phase_g2 = ph1, ph2
    return r


def with_element_phase(ring, index, extra):
    """Copy of the ring with an extra phase on the incoming pass of element `index` in BOTH half
    trips, i.e. the same error on two distinct physical mirrors (the semantics of the committed S3
    runs). Independent errors on every mirror: RingStrip.perturbed (scripts/cavity_correction.py)."""
    zc = ring.elements[index][0]
    i = [j for j, z in enumerate(ring.stop_positions()) if abs(z - zc) < 1e-6][0]
    return ring.perturbed([(1, i, extra), (2, i, extra)])


def evaluate(cfg, s, ring, e_in=None, g_target=None):
    loss_other = cfg["loss_other"]
    g_target = cfg["g_target"] if g_target is None else g_target
    F, x = s["F"], s["g"].x
    e_in = F if e_in is None else e_in
    t_in = B.coupler_for_enhancement(g_target, loss_other)
    lam0 = np.vdot(F, ring.roundtrip(F)) / np.vdot(F, F)
    res = FL.lock_and_solve(ring, e_in, t_in, loss_other, -np.angle(lam0))
    e = res["field"]
    eb = ring.half(e, 1)
    m = LA.site_map(e, eb, x, s["lam"], s["theta"], cfg["array_side_um"], cfg["map_sites_per_side"], w_z=cfg["w_z_um"], flip_b=True)
    g_an = float(B.enhancement(t_in, loss_other))
    out = {
        "t_in": t_in, "enhancement": res["enhancement"], "enhancement_analytic": g_an,
        "g_ratio": res["enhancement"] / g_an, "linewidth_rad": res["linewidth_rad"],
        "flat_A": FL.flatness(e, x, s["h"]), "flat_B": FL.flatness(eb, x, s["h"]),
        "array": LA.map_metrics(m), "roundtrip_loss_of_target": float(1 - abs(lam0) ** 2),
    }
    return out, e, m


def main(cfg_path):
    cfg = load_config(cfg_path)
    out_dir = run_dir(cfg["run_id"])
    write_params(out_dir, cfg)
    path = out_dir / "foxli.json"
    results = json.loads(path.read_text()) if path.exists() else {}

    def todo(key):
        return key not in results

    def save(key, val):
        results[key] = val
        write_json(path, results)
        a = val.get("array", {})
        print(f"{key}: G/G0={val.get('g_ratio', float('nan')):.3f} depth_rms={a.get('depth_rms', float('nan')):.4f} "
              f"dy_rms={a.get('dy_rms_um', float('nan')) * 1e3:.1f} nm", flush=True)

    V = cfg["variants"]
    lam = cfg["wavelength_um"]
    k = 2 * np.pi / lam

    # 1. Designs: ideal mode, dense spectrum, mirror phase, intensity on the mirrors.
    for name in cfg["designs"]:
        key = f"design_{name}"
        if not todo(key):
            continue
        t0 = time.time()
        s = build(cfg, V[name])
        ev, e, m = evaluate(cfg, s, s["ring"])
        g, F = s["g"], s["F"]
        uG = s["ring"].to_mirror(F)
        iG = np.abs(uG) ** 2
        sel = iG > 1e-3 * iG.max()
        sd = build(cfg, V[name], n=cfg["spectrum_grid"]["n"], width=cfg["spectrum_grid"]["width_um"])
        spec = FL.dense_spectrum(sd["ring"], sd["F"])
        mirrors = {}
        for i, (z, hw, ph) in enumerate(s["ring"].elements):
            ui = np.abs(s["ring"]._segment(F, 0.0, z - 1e-3)) ** 2
            mirrors[f"element{i}_z{z/1e3:.0f}mm"] = {"peak_over_flat_intensity": float(ui.max()), "rms_halfwidth_um": float(np.sqrt(np.sum(g.x**2 * ui) / np.sum(ui))), "half_width_um": hw}
        ev.update({
            "variant": V[name], "p_um": s["p"], "w_um": s["w"], "flat_halfrange_um": s["h"],
            "equivalent_width_um": LA.supergauss_equivalent_width(s["w"], V[name]["order"]),
            "fresnel_number_plain": s["w"] ** 2 / (lam * s["p"]),
            "graded_mirror": FL.phase_gradient_stats(s["ph"], g.x, sel, lam),
            "graded_mirror_peak_over_flat_intensity": float(iG.max() / np.abs(F[g.n // 2]) ** 2),
            "graded_mirror_rms_halfwidth_um": float(np.sqrt(np.sum(g.x**2 * iG) / np.sum(iG))),
            "elements": mirrors, "spectrum": spec, "seconds": time.time() - t0,
        })
        save(key, ev)
        np.savez_compressed(out_dir / f"{key}.npz", x=g.x, field=e, phase_g=s["ph"], iG=iG, depth_map=m["depth"], dy_map=m["dy"])

    # Resolution check of the base design.
    base = cfg["base_variants"][0]
    key = f"resolution_check_{base}"
    if todo(key):
        rc = {}
        for n in cfg["resolution_n"]:
            s = build(cfg, V[base], n=n)
            ev, _, _ = evaluate(cfg, s, s["ring"])
            rc[str(n)] = {"g_ratio": ev["g_ratio"], "flat_A_rms": ev["flat_A"]["rms"], "depth_rms": ev["array"]["depth_rms"], "dy_rms_um": ev["array"]["dy_rms_um"]}
        save(key, rc)

    # 2. Tolerances per variant: figure error on both graded mirrors, on the first telescope
    #    mirror (first half trip), x^2 error and tilt of G1.
    fe = cfg["figure_error"]
    for name in cfg["tolerance_variants"]:
        s = build(cfg, V[name])
        g, w, ph, ring = s["g"], s["w"], s["ph"], s["ring"]
        for rms in fe["rms_waves"]:
            for corr in fe["corr_um"]:
                for seed in fe["seeds"]:
                    key = f"figure_{name}_rms{rms}_corr{corr/1e3:.1f}mm_seed{seed}"
                    if todo(key):
                        rng = np.random.default_rng(seed)
                        hw = s["ring"].to_mirror(s["F"])
                        iw = np.abs(hw) ** 2
                        wG = 2 * np.sqrt(np.sum(g.x**2 * iw) / np.sum(iw))
                        s1 = FL.smooth_screen(g, 4 * np.pi * rms, corr * wG / w, rng, wG)
                        s2 = FL.smooth_screen(g, 4 * np.pi * rms, corr * wG / w, rng, wG)
                        ev, _, _ = evaluate(cfg, s, with_phases(ring, ph + s1, ph + s2))
                        ev.update({"variant": name, "rms_waves": rms, "corr_um_scaled_to_X": corr, "beam_halfwidth_on_G_um": wG})
                        save(key, ev)
                    key = f"figure_primary_{name}_rms{rms}_corr{corr/1e3:.1f}mm_seed{seed}"
                    if V[name]["geometry"] == "telescope" and todo(key):
                        rng = np.random.default_rng(seed + 100)
                        sc = FL.smooth_screen(g, 4 * np.pi * rms, corr, rng, w)
                        ev, _, _ = evaluate(cfg, s, with_element_phase(ring, 0, sc))
                        ev.update({"variant": name, "rms_waves": rms, "corr_um": corr})
                        save(key, ev)
        for sag in cfg["curvature_error_waves"]:
            key = f"curvature_{name}_sag{sag}"
            if todo(key):
                iw = np.abs(ring.to_mirror(s["F"])) ** 2
                wG = 2 * np.sqrt(np.sum(g.x**2 * iw) / np.sum(iw))
                ev, _, _ = evaluate(cfg, s, with_phases(ring, ph + 2 * k * sag * lam * (g.x / wG) ** 2, ph))
                ev.update({"variant": name, "sag_waves_at_beam_halfwidth": sag, "beam_halfwidth_on_G_um": wG})
                save(key, ev)
        for tilt in cfg["tilt_urad"]:
            key = f"tilt_{name}_{tilt}urad"
            if todo(key):
                ev, e, _ = evaluate(cfg, s, with_phases(ring, ph + 2 * k * tilt * 1e-6 * g.x, ph))
                i = np.abs(e) ** 2
                ev.update({"variant": name, "tilt_urad": tilt, "centroid_A_um": float(np.sum(g.x * i) / np.sum(i))})
                save(key, ev)
            key = f"tilt_primary_{name}_{tilt}urad"
            if V[name]["geometry"] == "telescope" and todo(key):
                ev, e, _ = evaluate(cfg, s, with_element_phase(ring, 0, 2 * k * tilt * 1e-6 * g.x))
                ev.update({"variant": name, "tilt_urad": tilt})
                save(key, ev)

    # 3. Single-variant studies on the base variants.
    for base in cfg["base_variants"]:
        single_variant_studies(cfg, V, base, todo, save, k)
    print("done", flush=True)


def single_variant_studies(cfg, V, base, todo, save, k):
    s = build(cfg, V[base])
    g, w, ph, ring, F = s["g"], s["w"], s["ph"], s["ring"], s["F"]
    order = V[base]["order"]
    im = cfg["input_mismatch"]
    cases = {f"input_{base}_width{sc}": LA.supergauss(g.x, w * sc, order).astype(complex) for sc in im["width_scale"]}
    cases.update({f"input_{base}_offset{o}um": LA.supergauss(g.x - o, w, order).astype(complex) for o in im["offset_um"]})
    ws = np.linspace(0.5, 1.5, 201) * w
    ovs = [abs(np.vdot(F, np.exp(-(g.x / q) ** 2))) ** 2 / (np.vdot(F, F).real * np.sum(np.exp(-2 * (g.x / q) ** 2))) for q in ws]
    cases[f"input_{base}_gaussian_best"] = np.exp(-(g.x / ws[int(np.argmax(ovs))]) ** 2).astype(complex)
    for key, e_in in cases.items():
        if todo(key):
            ev, _, _ = evaluate(cfg, s, ring, e_in=e_in)
            ev["mode_overlap"] = float(abs(np.vdot(F, e_in)) ** 2 / (np.vdot(F, F).real * np.vdot(e_in, e_in).real))
            save(key, ev)
    ms = cfg["metasurface"]
    for lev in ms["levels"]:
        key = f"meta_{base}_levels{lev}"
        if todo(key):
            q = FL.quantise(ph, lev)
            ev, _, _ = evaluate(cfg, s, with_phases(ring, q, q))
            ev["levels"] = lev
            save(key, ev)
    for sig in ms["pixel_noise_rad"]:
        key = f"meta_{base}_noise{sig}rad"
        if todo(key):
            rng = np.random.default_rng(7)
            ev, _, _ = evaluate(cfg, s, with_phases(ring, ph + sig * rng.standard_normal(g.n), ph + sig * rng.standard_normal(g.n)))
            ev.update({"pixel_noise_rad": sig, "dx_um": g.dx})
            save(key, ev)
    fs = cfg["finesse_scan"]
    rng = np.random.default_rng(1)
    iw = np.abs(ring.to_mirror(F)) ** 2
    wG = 2 * np.sqrt(np.sum(g.x**2 * iw) / np.sum(iw))
    s1 = FL.smooth_screen(g, 4 * np.pi * fs["rms_waves"], fs["corr_um"] * wG / w, rng, wG)
    s2 = FL.smooth_screen(g, 4 * np.pi * fs["rms_waves"], fs["corr_um"] * wG / w, rng, wG)
    for gt in fs["g_target"]:
        for label, (a1, a2) in {"ideal": (ph, ph), "figure": (ph + s1, ph + s2)}.items():
            key = f"finesse_{base}_G{gt:.0f}_{label}"
            if todo(key):
                ev, _, _ = evaluate(cfg, s, with_phases(ring, a1, a2), g_target=gt)
                ev["finesse"] = float(B.finesse(ev["t_in"], cfg["loss_other"]))
                save(key, ev)
    if cfg.get("reference_plane_mirrors"):
        key = f"reference_{base}_plane_mirrors"
        if todo(key):
            z = np.zeros(g.n)
            ev, _, _ = evaluate(cfg, s, with_phases(ring, z, z))
            save(key, ev)


if __name__ == "__main__":
    main(sys.argv[1])
