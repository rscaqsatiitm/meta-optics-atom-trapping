"""High-power thermal campaign (hypothesis H2a). Usage: python scripts/run_thermal.py configs/thermal_h2a.yaml

For each case (substrate, thickness, beam, power, film absorption) solves the
steady-state temperature field, converts it to an optical-path change and evaluates
the focal shift and Strehl ratio of the tweezers through the stated layer stack.
"""

import sys

import numpy as np
from scipy.optimize import brentq

from metatrap.propagation.layers import LayerStack
from metatrap.runio import load_config, run_dir, write_json, write_params
from metatrap.thermal.axisym_fd import SUBSTRATES, ThermalProblem, opl_change, solve, surface_temperature
from metatrap.thermal.thermo_optic import wavefront_effect


def evaluate(case, p_W, a_film, alpha_bulk=0.0):
    stack = LayerStack.from_list(case["stack"])
    lam, na = case["wavelength_um"], case["na"]
    ap_um = stack.aperture_for_na(na, lam)
    beam = case["beam"]
    waist = beam.get("waist_over_radius", None)
    prob = ThermalProblem(
        radius_m=case["substrate_radius_mm"] * 1e-3, thickness_m=case["thickness_mm"] * 1e-3,
        substrate=SUBSTRATES[case["substrate"]], beam=beam["kind"], beam_radius_m=ap_um * 1e-6,
        beam_waist_m=(waist * ap_um * 1e-6) if waist else None, h_faces=case.get("h_faces", 10.0),
        rim=case.get("rim", "fixed"), nr=200, nz=30, z_power=1.5,
    )
    rc, zc, T = solve(prob, p_W, a_film, alpha_bulk)
    dopl = opl_change(prob, rc, zc, T)
    ill = {"kind": "top_hat"} if beam["kind"] == "top_hat" else {"kind": "gaussian", "waist_um": waist * ap_um}
    wf = wavefront_effect(stack, lam, na, ill, rc, dopl)
    ts = surface_temperature(zc, T)
    inside = rc <= ap_um * 1e-6
    wf.update({
        "P_W": p_W, "a_film": a_film, "alpha_bulk_per_m": alpha_bulk,
        "dT_max_K": float(ts.max()), "dT_mean_beam_K": float(np.average(ts[inside], weights=rc[inside])),
        "aperture_diameter_mm": 2 * ap_um * 1e-3,
        "pitch_scale_ppm": float(SUBSTRATES[case["substrate"]].alpha_exp * np.average(ts[inside], weights=rc[inside]) * 1e6),
    })
    return wf


def main(cfg_path):
    cfg = load_config(cfg_path)
    out = run_dir(cfg["run_id"])
    write_params(out, cfg)
    res = {"sweeps": [], "inversions": []}
    for case in cfg["cases"]:
        rows = []
        for p in case["powers_W"]:
            for a in case["a_film"]:
                r = evaluate(case, p, a, case.get("alpha_bulk_per_m", 0.0))
                rows.append(r)
                print(case["name"], f"P={p} W a={a:g}: dTmax={r['dT_max_K']:.2f} K dz={r['dz_um']:.3f} um S_ref={r['strehl_refocused']:.4f} S_fix={r['strehl_fixed_plane']:.4f}", flush=True)
        # Absorbed fraction that gives a 1 um focal shift at the largest power.
        pmax = max(case["powers_W"])
        f = lambda la: abs(evaluate(case, pmax, 10**la)["dz_um"]) - 1.0  # noqa: E731
        try:
            a1 = 10 ** brentq(f, -8, -1, xtol=1e-3)
        except ValueError:
            a1 = float("nan")
        res["sweeps"].append({"case": case, "rows": rows, "a_film_for_1um_shift_at_Pmax": a1})
        print(case["name"], f"a_film for |dz| = 1 um at {pmax} W: {a1:.3g}", flush=True)
    for inv in cfg.get("inversions", []):
        case = next(c for c in cfg["cases"] if c["name"] == inv["case"])
        g = lambda la: evaluate(case, inv["P_W"], 10**la)[inv["metric"]] - inv["target"]  # noqa: E731
        try:
            a = 10 ** brentq(g, -7, -0.5, xtol=1e-3)
            r = evaluate(case, inv["P_W"], a)
            other = [evaluate(next(c for c in cfg["cases"] if c["name"] == o["case"]), o["P_W"], a) | {"case": o["case"]} for o in inv.get("then_evaluate", [])]
        except ValueError:
            a, r, other = float("nan"), None, []
        res["inversions"].append({"spec": inv, "a_film": a, "at_target": r, "consequences": other})
        print("inversion", inv["case"], f"a_film={a:.3g}", r, other, flush=True)
    write_json(out / "thermal.json", res)


if __name__ == "__main__":
    main(sys.argv[1])
