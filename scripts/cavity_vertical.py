"""Uniformity and plane flatness of the vertical lattice ring (flat-top in both directions).
Usage:
    python scripts/cavity_vertical.py configs/cavity_vertical.yaml

The vertical ring's beams propagate along the array diagonal and cross at +-theta_z in a vertical
plane (period a_z). Its mirrors are separable (cylindrical telescope and graded phase in the
horizontal direction x, graded phase only in the vertical direction z), so its mode is modelled
as the product of two 1D Fox-Li modes, each driven and locked separately (valid when one mode
dominates in each direction). A site at (x, y, 0) sees beam A at the vertical coordinate
zeta_A = +s sin(theta_z) of that beam and beam B at -s sin(theta_z), with s the coordinate along the
diagonal; reflections in the vertical ring plane flip zeta between the arms, not x.
Barrier depth ~ |Ex_A(u, s)| |Ex_B(u, s)| |Ez_A(zeta_A, s)| |Ez_B(zeta_B, s)|; vertical plane
displacement dz = -a_z dphi / (2 pi). Writes vertical.json.
"""

import sys

import numpy as np

from metatrap.cavity import budget as B
from metatrap.cavity import foxli as FL
from metatrap.cavity import lattice as LA
from metatrap.runio import REPO_ROOT, load_config, run_dir, write_json, write_params

sys.path.insert(0, str(REPO_ROOT / "scripts"))
import cavity_correction as C  # noqa: E402


def build(cfg, d):
    lam = cfg["wavelength_um"]
    g = FL.Grid1D(d["grid_n"], d["grid_width_um"], lam)
    w = LA.supergauss_w_for_flat(d["flat_half_um"], d["order"], cfg["flat_tol"])
    F = FL.flat_target(g, w, d["order"])
    t = cfg["telescope"]
    if isinstance(d["telescope"], dict):  # own focal lengths at the same mirror positions, f1' + f2' = f1 + f2
        f1z = (1 + t["mag"]) * t["f1_um"] / (1 + d["telescope"]["mag"])
        els = FL.telescope_elements(g, f1z, d["telescope"]["mag"] * f1z, t["d0_um"], d["ap1_um"], d["ap2_um"])
    elif d["telescope"]:
        f2 = t["mag"] * t["f1_um"]
        els = FL.telescope_elements(g, t["f1_um"], f2, t["d0_um"], d["ap1_um"], d["ap2_um"])
    else:  # same mirror positions, flat in this direction
        els = [(t["d0_um"], d["ap1_um"], None), (t["d0_um"] + t["f1_um"] + t["mag"] * t["f1_um"], d["ap2_um"], None)]
    p = t["d0_um"] + t["f1_um"] + t["mag"] * t["f1_um"] + t["L_um"]
    ring = FL.RingStrip(g, p, 0.0, np.zeros(g.n), np.zeros(g.n), ap_g_um=d["ap_g_um"], elements=els)
    ph = ring.design_phase(F, F, 1)
    ring.phase_g1 = ring.phase_g2 = ph
    return {"g": g, "F": F, "w": w, "ring": ring, "ph": ph}


def solve(cfg, s, ring):
    t_in = B.coupler_for_enhancement(cfg["g_target"], cfg["loss_other"])
    lam0 = np.vdot(s["F"], ring.roundtrip(s["F"])) / np.vdot(s["F"], s["F"])
    res = FL.lock_and_solve(ring, s["F"], t_in, cfg["loss_other"], -np.angle(lam0))
    e = res["field"]
    return e, ring.half(e, 1), res["enhancement"] / float(B.enhancement(t_in, cfg["loss_other"]))


def depth_map(cfg, sx, ex, sz, ez, return_map=False):
    lam = cfg["wavelength_um"]
    side = cfg["array_side_um"]
    thz = B.crossing_half_angle(cfg["vertical_period_um"], lam)
    hs = side / np.sqrt(2) * 1.02
    s_planes = np.linspace(-hs, hs, 81)
    stx = [LA.field_stack(e, sx["g"].x, lam, s_planes) for e in ex]
    stz = [LA.field_stack(e, sz["g"].x, lam, s_planes) for e in ez]
    n = cfg["map_sites_per_side"]
    gg = np.linspace(-side / 2, side / 2, n)
    X, Y = np.meshgrid(gg, gg, indexing="xy")
    s = (X + Y) / np.sqrt(2)
    u = (Y - X) / np.sqrt(2)
    ea = LA._interp_field(stx[0], sx["g"].x, s_planes, u, s)
    eb = LA._interp_field(stx[1], sx["g"].x, s_planes, u, s)
    za = LA._interp_field(stz[0], sz["g"].x, s_planes, s * np.sin(thz), s)
    zb = LA._interp_field(stz[1], sz["g"].x, s_planes, s * np.sin(thz), s)  # flipped arm: xi = -zeta_B = s sin
    amp = np.abs(ea) * np.abs(eb) * np.abs(za) * np.abs(zb)
    dphi = np.angle(ea * za * np.conj(eb * zb))
    c = n // 2
    dphi = np.angle(np.exp(1j * (dphi - dphi[c, c])))
    dz = -cfg["vertical_period_um"] * dphi / (2 * np.pi)
    m = {"x": X, "y": Y, "depth": amp / amp[c, c], "dy": dz}
    out = LA.map_metrics(m)
    res = {"depth_rms": out["depth_rms"], "min_over_mean": out["depth_min_over_mean"], "plane_dz_rms_um": out["dy_rms_um"],
           "plane_dz_pv_um": out["dy_pv_um"], "plane_dz_affine_res_rms_um": out["dy_affine_residual_rms_um"]}
    if return_map:
        m.update({"u": u, "s": s})
        return res, m
    return res


def fit_us(m, half_u, half_s, nb, smooth=1.0):
    """Fit log depth = alpha(u) + beta(s) (u: horizontal transverse, s: along the diagonal)."""
    def hat(v, half):
        knots = np.linspace(-half, half, nb)
        f = np.clip(np.interp(v, knots, np.arange(nb)), 0, nb - 1 - 1e-9)
        i0 = np.floor(f).astype(int)
        a = f - i0
        M = np.zeros((len(v), nb))
        M[np.arange(len(v)), i0] = 1 - a
        M[np.arange(len(v)), i0 + 1] += a
        return knots, M
    ku, Mu = hat(m["u"].ravel(), half_u)
    ks, Ms = hat(m["s"].ravel(), half_s)
    A = np.hstack([Mu, Ms])
    Lp = np.diff(np.eye(nb), 2, axis=0) * smooth
    R = np.block([[Lp, np.zeros_like(Lp)], [np.zeros_like(Lp), Lp]])
    gauge = np.concatenate([np.zeros(nb), np.ones(nb)])[None, :] * 10.0
    y = np.log(np.clip(m["depth"].ravel(), 1e-6, None))
    coef, *_ = np.linalg.lstsq(np.vstack([A, R, gauge]), np.concatenate([y, np.zeros(len(R) + 1)]), rcond=None)
    return ku, coef[:nb], ks, coef[nb:]


def closed_loop(cfg, sx, sz, seed):
    """Atom-data feedback with errors in both directions: the u-profile correction is shared by the
    two beams (only their product is seen), the s-profile maps onto the vertical flat-tops at
    zeta = s sin(theta_z) of both beams."""
    cl = cfg["closed_loop"]
    thz = B.crossing_half_angle(cfg["vertical_period_um"], cfg["wavelength_um"])
    rng = np.random.default_rng(seed)
    ch_x, ge_x = C.make_errors(sx["ring"], sx["F"], sx["w"], rng, cfg["fig_waves"], cfg["corr_um"], 0.0)
    ch_z, ge_z = C.make_errors(sz["ring"], sz["F"], sz["w"], rng, cfg["fig_waves"], cfg["corr_um"] * sz["w"] / sx["w"], 0.0)
    tx, tz = sx["ring"].perturbed(ch_x), sz["ring"].perturbed(ch_z)
    TAx, TBx, TAz, TBz = sx["F"].copy(), sx["F"].copy(), sz["F"].copy(), sz["F"].copy()
    px = [sx["ph"].copy(), sx["ph"].copy()]
    pz = [sz["ph"].copy(), sz["ph"].copy()]
    half_u = cfg["array_side_um"] / np.sqrt(2) * 1.02
    hist = []
    for rnd in range(cl["rounds"] + 1):
        rx, rz = tx.perturbed(), tz.perturbed()
        rx.phase_g1, rx.phase_g2 = px[0] + ge_x[1], px[1] + ge_x[2]
        rz.phase_g1, rz.phase_g2 = pz[0] + ge_z[1], pz[1] + ge_z[2]
        ex, ez = solve(cfg, sx, rx), solve(cfg, sz, rz)
        res, m = depth_map(cfg, sx, ex[:2], sz, ez[:2], return_map=True)
        hist.append(dict(res, round=rnd, g_ratio_x=ex[2], g_ratio_z=ez[2]))
        print(f"  closed loop seed {seed} round {rnd}: {res['depth_rms'] * 100:.3f} %  G {ex[2]:.3f}/{ez[2]:.3f}", flush=True)
        if rnd == cl["rounds"]:
            break
        ku, al, ks, be = fit_us(m, half_u, half_u, cl["n_bins"])
        g = cl["gain"]
        ax = np.exp(np.interp(np.clip(sx["g"].x, -half_u, half_u), ku, al - al.mean()))
        zeta = sz["g"].x
        bz = np.exp(np.interp(np.clip(zeta / np.sin(thz), -half_u, half_u), ks, be - be.mean()))
        TAx, TBx = TAx * ax ** (-g / 2), TBx * ax ** (-g / 2)
        TAz, TBz = TAz * bz ** (-g / 2), TBz * bz ** (-g / 2)
        px = [sx["ring"].design_phase(TAx, TBx, 1), sx["ring"].design_phase(TBx, TAx, 2)]
        pz = [sz["ring"].design_phase(TAz, TBz, 1), sz["ring"].design_phase(TBz, TAz, 2)]
    return hist


def main(cfg_path):
    cfg = load_config(cfg_path)
    out = run_dir(cfg["run_id"])
    write_params(out, cfg)
    sx, sz = build(cfg, cfg["x"]), build(cfg, cfg["z"])
    rep = {"x_design": {"w_um": sx["w"], "W_eq_um": LA.supergauss_equivalent_width(sx["w"], cfg["x"]["order"])},
           "z_design": {"w_um": sz["w"], "H_eq_um": LA.supergauss_equivalent_width(sz["w"], cfg["z"]["order"])},
           "spectrum_z": FL.dense_spectrum(build(cfg, dict(cfg["z"], grid_n=2048))["ring"], build(cfg, dict(cfg["z"], grid_n=2048))["F"], n_modes=4),
           "spectrum_x": FL.dense_spectrum(build(cfg, dict(cfg["x"], grid_n=4096))["ring"], build(cfg, dict(cfg["x"], grid_n=4096))["F"], n_modes=4),
           "cases": {}}
    ex0, ez0 = solve(cfg, sx, sx["ring"]), solve(cfg, sz, sz["ring"])
    rep["cases"]["ideal"] = dict(depth_map(cfg, sx, ex0[:2], sz, ez0[:2]), g_ratio_x=ex0[2], g_ratio_z=ez0[2])
    print("ideal", rep["cases"]["ideal"], flush=True)
    for direction in ("x", "z"):
        s_ = sx if direction == "x" else sz
        for seed in cfg["seeds"]:
            rng = np.random.default_rng(seed)
            changes, g_err = C.make_errors(s_["ring"], s_["F"], s_["w"], rng, cfg["fig_waves"], cfg["corr_um"] * s_["w"] / sx["w"], 0.0)
            true = s_["ring"].perturbed(changes)
            for label, (p1, p2) in {"uncorrected": (s_["ph"], s_["ph"]),
                                    "openloop_perfect": (true.design_phase(s_["F"], s_["F"], 1), true.design_phase(s_["F"], s_["F"], 2))}.items():
                r = true.perturbed()
                r.phase_g1, r.phase_g2 = p1 + g_err[1] * (label == "uncorrected"), p2 + g_err[2] * (label == "uncorrected")
                e = solve(cfg, s_, r)
                ex, ez = (e, ez0) if direction == "x" else (ex0, e)
                key = f"fig{cfg['fig_waves']}_{direction}_seed{seed}_{label}"
                rep["cases"][key] = dict(depth_map(cfg, sx, ex[:2], sz, ez[:2]), g_ratio=e[2])
                print(key, rep["cases"][key], flush=True)
    for seed in cfg.get("closed_loop", {}).get("seeds", []):
        rep["cases"][f"fig{cfg['fig_waves']}_both_seed{seed}_closedloop"] = {"history": closed_loop(cfg, sx, sz, seed)}
    write_json(out / "vertical.json", rep)


if __name__ == "__main__":
    main(sys.argv[1])
