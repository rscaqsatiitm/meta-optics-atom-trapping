"""Thermal-mechanical damage assessment and power limits. Usage:
    python scripts/damage_assessment.py configs/damage_v2.yaml

For each substrate / film-absorption case the steady-state temperature field is solved at a
reference power and every criterion is scaled linearly (all are linear in absorbed power):
substrate tensile stress against strength / safety factor, peak temperature rise against a
limit, focal shift against the 1 um budget, plus film mismatch stress, the largest tolerable
absorbing defect and a two-photon-absorption bound. Writes damage.json.
"""

import sys

import numpy as np

from metatrap.propagation.layers import LayerStack
from metatrap.runio import load_config, run_dir, write_json, write_params
from metatrap.thermal import damage as D
from metatrap.thermal.axisym_fd import SUBSTRATES, ThermalProblem, opl_change, solve, surface_temperature
from metatrap.thermal.thermo_optic import wavefront_effect


def main(cfg_path):
    cfg = load_config(cfg_path)
    out = run_dir(cfg["run_id"])
    write_params(out, cfg)
    stack = LayerStack.from_list(cfg["stack"])
    lam, na = cfg["wavelength_um"], cfg["na"]
    R_ap = stack.aperture_for_na(na, lam) * 1e-6
    area_cm2 = np.pi * (R_ap * 100) ** 2
    p_ref = cfg["reference_power_W"]
    lim = cfg["limits"]
    rows = []
    for case in cfg["cases"]:
        sub = SUBSTRATES[case["substrate"]]
        prob = ThermalProblem(radius_m=cfg["die_radius_mm"] * 1e-3, thickness_m=case["thickness_mm"] * 1e-3, substrate=sub,
                              beam="top_hat", beam_radius_m=R_ap, h_faces=10.0, rim=cfg.get("rim", "fixed"), nr=220, nz=40, z_power=1.5)
        rc, zc, T = solve(prob, p_ref, case["a_film"], case.get("alpha_bulk_per_m", 0.0))
        ts = surface_temperature(zc, T)
        dz = np.diff(prob.thickness_m * np.linspace(0, 1, len(zc) + 1) ** prob.z_power)
        T_bar = np.sum(T * dz[None, :], axis=1) / prob.thickness_m
        sr, st, st_edge = D.disc_stresses(rc, T_bar, prob.radius_m, sub)
        sig_mem = max(float(sr.max()), float(st.max()), float(st_edge))
        sig_z = D.through_thickness_stress(zc, T, sub)
        sig_max = sig_mem + max(sig_z, 0.0)
        wf = wavefront_effect(stack, lam, na, {"kind": "top_hat"}, rc, opl_change(prob, rc, zc, T))
        film = cfg["films"][case["film"]]
        dsig_film = D.film_mismatch_stress_per_K(film["biaxial_modulus_GPa"] * (1 - film.get("nu", 0.0)), film.get("nu", 0.0), film["alpha"], sub) * ts.max()
        allow = sub.strength_MPa * 1e6 / lim["stress_safety_factor"]
        p_stress = p_ref * allow / sig_max if sig_max > 0 else np.inf
        p_temp = p_ref * lim["dT_K"] / ts.max()
        p_focus = p_ref * lim["focal_shift_um"] / abs(wf["dz_um"]) if wf["dz_um"] else np.inf
        row = {
            "name": case["name"], "case": case,
            "at_ref_power": {"P_W": p_ref, "intensity_W_cm2": p_ref / area_cm2, "dT_max_K": float(ts.max()),
                             "sigma_membrane_MPa": sig_mem / 1e6, "sigma_through_thickness_MPa": sig_z / 1e6,
                             "sigma_tensile_max_MPa": sig_max / 1e6, "film_mismatch_stress_MPa": dsig_film / 1e6,
                             "focal_shift_um": wf["dz_um"], "strehl_refocused": wf["strehl_refocused"]},
            "power_limit_W": {"stress": p_stress, "temperature": p_temp, "focal_shift": p_focus},
            "allowable_stress_MPa": allow / 1e6,
        }
        row["power_limit_W"]["binding"] = min(row["power_limit_W"], key=lambda k: row["power_limit_W"][k])
        rows.append(row)
        pl = row["power_limit_W"]
        print(f"{case['name']}: dT={ts.max():.2f} K, sigma={sig_max/1e6:.3f} MPa (allow {allow/1e6:.1f}), film {dsig_film/1e6:.2f} MPa, dz={wf['dz_um']:.3f} um at {p_ref} W"
              f" -> limits stress {pl['stress']:.0f} W, temp {pl['temperature']:.0f} W, focus {pl['focal_shift']:.0f} W (binding: {pl['binding']})", flush=True)
    # Defect hot spots and TPA bound, per substrate, over the stated powers.
    defects = {}
    for name in ("fused_silica", "sapphire"):
        kappa = SUBSTRATES[name].kappa
        defects[name] = {f"{p}W": {"intensity_W_cm2": p / area_cm2,
                                   "max_absorbing_defect_radius_um_for_dT": float(lim["defect_dT_K"] * kappa / (p / area_cm2 * 1e4) * 1e6),
                                   "dT_for_100um_radius_defect_K": float(D.defect_rise(p / area_cm2 * 1e4, 100e-6, kappa))}
                         for p in cfg["powers_W"]}
    tpa = {f"{p}W": D.tpa_over_linear(p / area_cm2 * 1e4, cfg["tpa"]["enhancement"], cfg["tpa"]["beta_cm_per_GW"], cfg["tpa"]["alpha_linear_per_cm"]) for p in cfg["powers_W"]}
    print("defects", defects, "\nTPA/linear", tpa, flush=True)
    write_json(out / "damage.json", {"aperture_area_cm2": area_cm2, "cases": rows, "defects": defects, "tpa_over_linear": tpa, "limits": lim})


if __name__ == "__main__":
    main(sys.argv[1])
