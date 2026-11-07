"""Heat budget of the cavity lattice: where the laser power goes, mirror temperatures in vacuum,
and the temperature stability the alignment needs.
Usage:
    python scripts/cavity_heat_budget.py configs/cavity_heat_budget.yaml

Power flow per ring (DERIVED). With input P_in, mode matching eta and the over-coupled coupler T = 3L:
  reflected to an external dump:  (1 - eta) P_in + eta P_in ((T - L)/(T + L))^2 = (1 - 0.75 eta) P_in
  dissipated in the ring:         P_circ L = 0.75 eta P_in
The ring's loss L (10 reflections of m ppm plus 2 graded mirrors of g each) is split into coating
absorption (alpha per reflection), graded-mirror loss absorbed (fraction f_g of g) and the rest,
scattered or transmitted out of the mode. Circulating powers: results/cavity_scaling/scaling.json
(adopted design at each size). Mirror temperature: steady-state axisymmetric conduction
(src/metatrap/thermal/axisym_fd.py), absorbed power on an equal-area disc of the spot, faces cooled
by radiation only (vacuum, linearised h = 4 eps sigma T^3), rim either clamped to a mount at the
chamber temperature or radiating too (worst case, free mirror). Surface bulge (1 + nu) alpha int T dz,
as in scripts/cavity_thermal.py. Writes heat.json.
"""

import json
import sys

import numpy as np

from metatrap.runio import REPO_ROOT, load_config, run_dir, write_json, write_params
from metatrap.thermal.axisym_fd import SUBSTRATES, ThermalProblem, solve, surface_temperature

SIGMA = 5.670374419e-8


def mirror_thermal(m, p_abs, h_rad, lam_m):
    sub = SUBSTRATES[m["substrate"]]
    out = {}
    for rim in ("clamped", "radiating"):
        prob = ThermalProblem(radius_m=m["radius_m"], thickness_m=m["thickness_m"], substrate=sub, beam="top_hat",
                              beam_radius_m=m["spot_radius_m"], h_faces=h_rad, rim="fixed" if rim == "clamped" else "convective",
                              h_rim=h_rad, nr=200, nz=50, z_power=2.0)
        rc, zc, T = solve(prob, p_abs, 1.0)
        ts = surface_temperature(zc, T)
        dz = np.diff(prob.thickness_m * np.linspace(0, 1, len(zc) + 1) ** prob.z_power)
        u = (1 + sub.poisson) * sub.alpha_exp * np.sum(T * dz[None, :], axis=1)
        a = m["spot_radius_m"]
        u_rel = u - np.interp(a, rc, u)
        out[rim] = {"dT_peak_K": float(ts.max()), "dT_rim_K": float(T[-1].mean()), "bulge_centre_minus_edge_nm": float(u_rel[0] * 1e9),
                    "bulge_waves": float(u_rel[0] / lam_m)}
    return out


def main(cfg_path):
    cfg = load_config(cfg_path)
    out = run_dir(cfg["run_id"])
    write_params(out, cfg)
    sc = json.loads((REPO_ROOT / cfg["scaling_json"]).read_text())["sizes"]
    lam_m = cfg["wavelength_nm"] * 1e-9
    h_rad = 4 * cfg["emissivity"] * SIGMA * cfg["chamber_K"] ** 3
    eta = cfg["mode_matching"]
    rep = {"h_rad_W_m2K": h_rad, "cases": []}
    for key, row in sc.items():
        p = row["power"]
        rings = {"in-plane (each of two)": p["circulating_inplane_per_axis_W"], "vertical": p["circulating_vertical_W"]}
        for g_case in cfg["enhancements"]:
            G, nref, m_loss, g_loss = g_case["G"], cfg["reflections"], cfg["mirror_loss"], g_case["graded_extra_loss"]
            L = nref * m_loss + 2 * g_loss
            for absb in cfg["absorption_cases"]:
                case = {"n_sites": row["n_sites"], "G": G, "L": L, "absorption": absb["name"], "rings": {}}
                total = {"input_W": 0.0, "reflected_W": 0.0, "absorbed_in_vacuum_W": 0.0, "scattered_or_transmitted_W": 0.0}
                for rname, pc in rings.items():
                    p_in = pc / (G * eta)
                    refl = (1 - 0.75 * eta) * p_in
                    diss = pc * L
                    a_coat = pc * nref * absb["alpha_per_reflection"]
                    a_grad = pc * 2 * g_loss * absb["graded_loss_absorbed_fraction"]
                    scat = diss - a_coat - a_grad
                    per_mirror = pc * absb["alpha_per_reflection"]
                    hottest = pc * (absb["alpha_per_reflection"] + g_loss * absb["graded_loss_absorbed_fraction"])
                    mk = "inplane_graded" if rname.startswith("in-plane") else "vertical_graded"
                    th = mirror_thermal(cfg["mirrors"][mk], hottest, h_rad, lam_m)
                    case["rings"][rname] = {"P_circ_W": pc, "P_in_W": p_in, "reflected_W": refl, "dissipated_W": diss,
                                            "absorbed_coatings_W": a_coat, "absorbed_graded_extra_W": a_grad, "scattered_or_transmitted_W": scat,
                                            "absorbed_per_plain_mirror_W": per_mirror, "absorbed_graded_mirror_W": hottest, "graded_mirror": th}
                    n = 2 if rname.startswith("in-plane") else 1
                    total["input_W"] += n * p_in
                    total["reflected_W"] += n * refl
                    total["absorbed_in_vacuum_W"] += n * (a_coat + a_grad)
                    total["scattered_or_transmitted_W"] += n * scat
                case["total"] = total
                rep["cases"].append(case)
    # Temperature stability for the alignment and the lock (DERIVED from placeholder expansion coefficients).
    st = cfg["stability"]
    nu_laser = 299792458.0 / lam_m
    rep["stability"] = {}
    for mat, alpha in st["cte_per_K"].items():
        rep["stability"][mat] = {
            "cte_per_K": alpha,
            "resonance_shift_MHz_per_K": nu_laser * alpha / 1e6,
            "gradient_K_for_tilt_tolerance": st["tilt_tolerance_rad"] * st["plate_thickness_m"] / (alpha * st["lever_m"]),
            "uniform_K_for_100MHz_detuning_change": 100e6 / (nu_laser * alpha) if alpha > 0 else None,
        }
    # Mirror warm-up time constant, radiation only.
    for mk, m in cfg["mirrors"].items():
        vol = np.pi * m["radius_m"] ** 2 * m["thickness_m"]
        area = 2 * np.pi * m["radius_m"] ** 2 + 2 * np.pi * m["radius_m"] * m["thickness_m"]
        rep.setdefault("time_constants_s", {})[mk] = float(vol * m["density_kg_m3"] * m["heat_capacity_J_kgK"] / (h_rad * area))
    write_json(out / "heat.json", rep)
    for c in rep["cases"]:
        t = c["total"]
        print(f"N={c['n_sites']:.0e} G={c['G']} {c['absorption']:>22s}: in {t['input_W']:.3f} W, dumped {t['reflected_W']:.3f}, "
              f"absorbed in vacuum {t['absorbed_in_vacuum_W'] * 1e3:.1f} mW, scattered {t['scattered_or_transmitted_W']:.3f} W")
        for rn, r in c["rings"].items():
            g = r["graded_mirror"]
            print(f"    {rn:24s} hottest mirror {r['absorbed_graded_mirror_W'] * 1e3:7.2f} mW: clamped dT {g['clamped']['dT_peak_K']:.3f} K, "
                  f"bulge {g['clamped']['bulge_centre_minus_edge_nm']:.3f} nm; free dT {g['radiating']['dT_peak_K']:.2f} K")
    print(json.dumps(rep["stability"], indent=1))
    print(rep["time_constants_s"])


if __name__ == "__main__":
    main(sys.argv[1])
