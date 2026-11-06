"""Correction of the lattice depth non-uniformity by re-shaping the graded mirrors G1, G2.
Usage:
    python scripts/cavity_correction.py configs/cavity_correction.yaml

Errors. Independent smooth figure errors (and optionally random tilts) on every mirror of the
ring: the 2 x (number of folds) fold/telescope reflections of the two half trips and the two
graded mirrors. Correlation lengths are referred to the crossing and scaled by the beam size on
each mirror; rms values are surface heights (reflected phase 4 pi h / lambda).

Corrections. G1 shapes half trip 1 (arm A -> arm B), G2 half trip 2 (B -> A).
  open_loop: G phases redesigned through the half trips as measured, phi_w = arg(K_out^+ T) -
      arg(K_in T) (RingStrip.design_phase); the fabricated surface is phi_w minus the measured G
      error, so the true G phase is phi_w + (true - measured G error). Metrology noise is a smooth
      random error map of the given rms on every mirror (and a tilt error, when tilts are present).
  closed_loop: no knowledge of the errors. Each round the per-site barrier depths (101 x 101
      sampled sites, optional multiplicative measurement noise) are fitted as a product
      a(xi_A) b(xi_B) of the two beams' transverse profiles (log-linear least squares with a
      smoothness penalty); the target profiles T_A, T_B are multiplied by (a, b)^-gain and the
      G phases redesigned with the nominal model.
Writes results/<run_id>/correction.json (resumable).
"""

import json
import sys
import time

import numpy as np

from metatrap.cavity import budget as B
from metatrap.cavity import foxli as FL
from metatrap.cavity import lattice as LA
from metatrap.runio import REPO_ROOT, load_config, run_dir, write_json, write_params

sys.path.insert(0, str(REPO_ROOT / "scripts"))
import cavity_foxli as S  # noqa: E402


def beam_halfwidths(ring, F):
    """2 x rms half-width of the unperturbed beam on every stop of half trip 1, and on G."""
    out = []
    eG = ring.to_mirror(F, 1)
    after = eG * ring._ap_g * np.exp(1j * ring.phase_g1)
    x = ring.grid.x
    for z in ring.stop_positions():
        u = ring._segment(F, 0.0, z - 1e-6, 1) if z < ring.p_um else ring._segment(after, ring.p_um, z - 1e-6, 1)
        i = np.abs(u) ** 2
        out.append(2 * np.sqrt(np.sum(x**2 * i) / np.sum(i)))
    i = np.abs(eG) ** 2
    return out, 2 * np.sqrt(np.sum(x**2 * i) / np.sum(i))


def make_errors(ring, F, w_x, rng, fig_waves, corr_um, tilt_urad):
    """Independent figure + tilt phase maps for every mirror: (changes, g_err)."""
    g = ring.grid
    k = 2 * np.pi / g.lam_um
    sizes, wG = beam_halfwidths(ring, F)
    changes, g_err = [], {}
    for which in (1, 2):
        for i, wi in enumerate(sizes):
            ph = np.zeros(g.n)
            if fig_waves:
                ph += FL.smooth_screen(g, 4 * np.pi * fig_waves, corr_um * wi / w_x, rng, wi)
            if tilt_urad:
                ph += 2 * k * rng.normal(0.0, tilt_urad) * 1e-6 * g.x
            changes.append((which, i, ph))
        ph = np.zeros(g.n)
        if fig_waves:
            ph += FL.smooth_screen(g, 4 * np.pi * fig_waves, corr_um * wG / w_x, rng, wG)
        if tilt_urad:
            ph += 2 * k * rng.normal(0.0, tilt_urad) * 1e-6 * g.x
        g_err[which] = ph
    return changes, g_err


def stroke(phi_new, phi_old, ring, F, lam):
    """Correction relative to the nominal design, over the beam on G: rms and PV of surface (nm)."""
    i = np.abs(ring.to_mirror(F, 1)) ** 2
    sel = i > 1e-3 * i.max()
    d = np.angle(np.exp(1j * (phi_new - phi_old)))[sel]
    d = d - np.average(d, weights=i[sel])
    h = d * lam * 1e3 / (4 * np.pi)
    return {"rms_nm": float(np.sqrt(np.average(h**2, weights=i[sel]))), "pv_nm": float(h.max() - h.min())}


def fit_profiles(m, theta, half, nb=48, smooth=1.0):
    """Fit log depth = alpha(xi_A) + beta(xi_B), xi_A = u_A, xi_B = -u_B (flipped arm B)."""
    X, Y, D = m["x"].ravel(), m["y"].ravel(), m["depth"].ravel()
    sn, ct = np.sin(theta), np.cos(theta)
    xa, xb = -X * sn + Y * ct, -(X * sn + Y * ct)
    knots = np.linspace(-half, half, nb)

    def hat(u):
        f = np.clip(np.interp(u, knots, np.arange(nb)), 0, nb - 1 - 1e-9)
        i0 = np.floor(f).astype(int)
        a = f - i0
        M = np.zeros((len(u), nb))
        M[np.arange(len(u)), i0] = 1 - a
        M[np.arange(len(u)), i0 + 1] += a
        return M

    A = np.hstack([hat(xa), hat(xb)])
    y = np.log(np.clip(D, 1e-6, None))
    L = np.diff(np.eye(nb), 2, axis=0)
    R = np.block([[L, np.zeros_like(L)], [np.zeros_like(L), L]]) * smooth
    gauge = np.concatenate([np.ones(nb), -np.ones(nb)])[None, :] * 10.0
    coef, *_ = np.linalg.lstsq(np.vstack([A, R, gauge]), np.concatenate([y, np.zeros(len(R) + 1)]), rcond=None)
    return knots, coef[:nb], coef[nb:]


def profile_on_grid(x, knots, vals, half):
    v = np.interp(np.clip(x, -half, half), knots, vals)
    return np.exp(v - vals.mean())


def main(cfg_path):
    cfg = load_config(cfg_path)
    out_dir = run_dir(cfg["run_id"])
    write_params(out_dir, cfg)
    path = out_dir / "correction.json"
    results = json.loads(path.read_text()) if path.exists() else {}
    fcfg = load_config(REPO_ROOT / cfg["foxli_config"])
    fcfg.update({k: cfg[k] for k in ("g_target", "map_sites_per_side") if k in cfg})
    lam = fcfg["wavelength_um"]

    def save(key, val):
        results[key] = val
        write_json(path, results)
        print(key, json.dumps(val, default=float)[:300], flush=True)

    for vname in cfg["variants"]:
        s = S.build(fcfg, fcfg["variants"][vname])
        ring0, F, ph0, w_x, h = s["ring"], s["F"], s["ph"], s["w"], s["h"]
        for case in cfg["cases"]:
            for seed in case.get("seeds", [1]):
                tag = f"{vname}_{case['name']}_seed{seed}"
                rng = np.random.default_rng(seed)
                changes, g_err = make_errors(ring0, F, w_x, rng, case.get("fig_waves", 0.0), case.get("corr_um", 2000.0), case.get("tilt_urad", 0.0))
                true = ring0.perturbed(changes)

                def with_g(ring, p1, p2):
                    r = ring.perturbed()
                    r.phase_g1, r.phase_g2 = p1 + g_err[1], p2 + g_err[2]
                    return r

                key = f"{tag}_uncorrected"
                if key not in results:
                    ev, _, _ = S.evaluate(fcfg, s, with_g(true, ph0, ph0))
                    save(key, {"depth_rms": ev["array"]["depth_rms"], "min_over_mean": ev["array"]["depth_min_over_mean"], "g_ratio": ev["g_ratio"], "dy_affine_res_nm": ev["array"]["dy_affine_residual_rms_um"] * 1e3})
                # Open loop: perfect metrology and noisy metrology.
                for noise_nm in cfg["metrology_noise_nm"]:
                    key = f"{tag}_openloop_noise{noise_nm}nm"
                    if key in results:
                        continue
                    nrng = np.random.default_rng(1000 + seed)
                    nwave = noise_nm * 1e-3 / lam
                    n_changes, n_g = make_errors(ring0, F, w_x, nrng, nwave, cfg["metrology_corr_um"], cfg["metrology_tilt_urad"] if case.get("tilt_urad") else 0.0)
                    meas = ring0.perturbed([(w, i, e + n) for (w, i, e), (_, _, n) in zip(changes, n_changes)])
                    p1, p2 = meas.design_phase(F, F, 1), meas.design_phase(F, F, 2)
                    # Fabricated surface: p - measured G error (g_err + n_g); with the true error g_err
                    # on top, the true G phase is p - n_g.
                    r = true.perturbed()
                    r.phase_g1, r.phase_g2 = p1 - n_g[1], p2 - n_g[2]
                    ev, _, _ = S.evaluate(fcfg, s, r)
                    save(key, {"depth_rms": ev["array"]["depth_rms"], "min_over_mean": ev["array"]["depth_min_over_mean"], "g_ratio": ev["g_ratio"],
                               "dy_affine_res_nm": ev["array"]["dy_affine_residual_rms_um"] * 1e3,
                               "stroke_G1": stroke(p1 - g_err[1] - n_g[1], ph0, ring0, F, lam), "stroke_G2": stroke(p2 - g_err[2] - n_g[2], ph0, ring0, F, lam)})
                # Closed loop from atom-measured depths (nominal model only).
                if not case.get("closed_loop", False):
                    continue
                for mnoise in cfg["closed_loop"]["measurement_noise"]:
                    key = f"{tag}_closedloop_noise{mnoise}"
                    if key in results:
                        continue
                    t0 = time.time()
                    mrng = np.random.default_rng(2000 + seed)
                    TA, TB = F.copy(), F.copy()
                    p1, p2 = ph0.copy(), ph0.copy()
                    hist = []
                    for rnd in range(cfg["closed_loop"]["rounds"] + 1):
                        ev, _, m = S.evaluate(fcfg, s, with_g(true, p1, p2))
                        hist.append({"round": rnd, "depth_rms": ev["array"]["depth_rms"], "min_over_mean": ev["array"]["depth_min_over_mean"], "g_ratio": ev["g_ratio"],
                                     "dy_affine_res_nm": ev["array"]["dy_affine_residual_rms_um"] * 1e3})
                        print(f"  {key} round {rnd}: {hist[-1]['depth_rms']*100:.3f} % G/G0 {hist[-1]['g_ratio']:.3f}", flush=True)
                        if rnd == cfg["closed_loop"]["rounds"]:
                            break
                        mm = dict(m)
                        mm["depth"] = m["depth"] * (1 + mnoise * mrng.standard_normal(m["depth"].shape))
                        knots, al, be = fit_profiles(mm, s["theta"], h * 1.02, cfg["closed_loop"]["n_bins"], cfg["closed_loop"]["smooth"])
                        gain = cfg["closed_loop"]["gain"]
                        TA = TA * profile_on_grid(ring0.grid.x, knots, al, h * 1.02) ** (-gain)
                        TB = TB * profile_on_grid(ring0.grid.x, knots, be, h * 1.02) ** (-gain)
                        p1, p2 = ring0.design_phase(TA, TB, 1), ring0.design_phase(TB, TA, 2)
                    save(key, {"history": hist, "final": hist[-1], "stroke_G1": stroke(p1, ph0, ring0, F, lam), "stroke_G2": stroke(p2, ph0, ring0, F, lam), "seconds": time.time() - t0})
                # Closed loop from cameras that image the transmitted circulating beams (arms A and B)
                # onto planes conjugate to the crossing; no atoms, no mirror metrology.
                cam = cfg.get("camera")
                if not cam or vname not in cam["variants"]:
                    continue
                for cn in cam["noise"]:
                    key = f"{tag}_closedloop_camera_{cn['name']}"
                    if key in results:
                        continue
                    t0 = time.time()
                    crng = np.random.default_rng(3000 + seed)
                    px = np.arange(-h * 1.02, h * 1.02, cam["pixel_um"])
                    pats = [fixed_pattern(px, cn, crng, ring0.grid) for _ in range(2)]  # one camera per arm
                    TA, TB = F.copy(), F.copy()
                    p1, p2 = ph0.copy(), ph0.copy()
                    hist = []
                    for rnd in range(cam["rounds"] + 1):
                        r = with_g(true, p1, p2)
                        ev, e, _ = S.evaluate(fcfg, s, r)
                        hist.append({"round": rnd, "depth_rms": ev["array"]["depth_rms"], "min_over_mean": ev["array"]["depth_min_over_mean"], "g_ratio": ev["g_ratio"],
                                     "dy_affine_res_nm": ev["array"]["dy_affine_residual_rms_um"] * 1e3})
                        print(f"  {key} round {rnd}: {hist[-1]['depth_rms']*100:.3f} % G/G0 {hist[-1]['g_ratio']:.3f}", flush=True)
                        if rnd == cam["rounds"]:
                            break
                        eb = r.half(e, 1)
                        knots = np.linspace(-h * 1.02, h * 1.02, cam["n_bins"])
                        va = camera_profile(e, ring0.grid.x, px, knots, pats[0], cn["random"], crng)
                        vb = camera_profile(eb, ring0.grid.x, px, knots, pats[1], cn["random"], crng)
                        TA = TA * profile_on_grid(ring0.grid.x, knots, va, h * 1.02) ** (-cam["gain"])
                        TB = TB * profile_on_grid(ring0.grid.x, knots, vb, h * 1.02) ** (-cam["gain"])
                        p1, p2 = ring0.design_phase(TA, TB, 1), ring0.design_phase(TB, TA, 2)
                    save(key, {"history": hist, "final": hist[-1], "camera": cn, "stroke_G1": stroke(p1, ph0, ring0, F, lam), "seconds": time.time() - t0})
    print("done", flush=True)


def fixed_pattern(px, cn, rng, grid):
    """Static camera response error on the pixel grid: white per-pixel (PRNU) and smooth (flat field)."""
    pat = 1.0 + cn.get("fixed_white", 0.0) * rng.standard_normal(len(px))
    if cn.get("fixed_smooth", 0.0):
        sm = FL.smooth_screen(grid, cn["fixed_smooth"], cn.get("smooth_corr_um", 1000.0), rng, px.max())
        pat = pat * (1.0 + np.interp(px, grid.x, sm))
    return pat


def camera_profile(e, x, px, knots, pattern, random, rng):
    """Log-amplitude profile (on the knots) of a beam's intensity seen by a camera: sampled on pixels,
    multiplied by the static pattern and by fresh random noise, then averaged into bins."""
    ip = np.interp(px, x, np.abs(e) ** 2) * pattern * (1.0 + random * rng.standard_normal(len(px)))
    idx = np.clip(np.searchsorted(knots, px) - 1, 0, len(knots) - 1)
    mean = np.array([ip[idx == i].mean() if np.any(idx == i) else np.nan for i in range(len(knots))])
    mean = np.where(np.isfinite(mean), mean, np.nanmean(mean))
    return 0.5 * np.log(np.clip(mean, 1e-12, None))


if __name__ == "__main__":
    main(sys.argv[1])
