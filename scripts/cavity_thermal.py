"""S5: intensities, absorbing-defect limits and thermal deformation inside the cavity lattice.
Usage:
    python scripts/cavity_thermal.py configs/cavity_thermal.yaml

Intensity on each optic: every crossing beam has the peak intensity I1 at the atoms (fixed by
the barrier depth); on a mirror it is I1 x (1D horizontal peak ratio from S3) x (w_z at the
atoms / w_z on the mirror) for the in-plane rings (vertical Gaussian), and I1 / (area ratio)
for the vertical ring. Defect criterion of thermal/damage.py: a fully absorbing disc of radius a
on a half-space heats by I a / kappa; the limit is 50 K. Thermal deformation: the axisymmetric
solver of thermal/axisym_fd.py with an equal-area top-hat source; surface bulge of a thick
mirror u(r) = (1 + nu) alpha int T dz (DERIVED, constrained-expansion estimate), compared with
the figure tolerances of S3. Writes thermal.json.
"""

import json
import sys

import numpy as np

from metatrap.cavity import budget as B
from metatrap.runio import REPO_ROOT, load_config, run_dir, write_json, write_params
from metatrap.thermal.axisym_fd import SUBSTRATES, ThermalProblem, solve, surface_temperature


def main(cfg_path):
    cfg = load_config(cfg_path)
    out = run_dir(cfg["run_id"])
    write_params(out, cfg)
    lam = cfg["wavelength_nm"]
    lam_m = lam * 1e-9
    fox = json.loads((REPO_ROOT / cfg["foxli_json"]).read_text())
    lay = json.loads((REPO_ROOT / cfg["layout_json"]).read_text())
    bud = json.loads((REPO_ROOT / cfg["budget_json"]).read_text())
    des = fox[cfg["design_key"]]
    vd = [r for r in lay["vertical_designs"] if r.get("f_z_m") == cfg["vertical_f_z_m"]][0]
    w_x = des["equivalent_width_um"] * 1e-6
    wz0 = cfg["w_z_m"]
    zr = np.pi * wz0**2 / lam_m
    ratio_G = des["graded_mirror_peak_over_flat_intensity"]
    ratio_S = [v for k, v in des["elements"].items() if k.startswith("element1")][0]["peak_over_flat_intensity"]
    rep = {"cases": [], "defects": {}, "deformation": []}
    for c in bud["cases"]:
        if c["n_sites"] != cfg["n_sites"] or c["depth_mK"] not in cfg["depths_mK"]:
            continue
        pc = c["circulating"]
        i1 = pc["beam_intensity_W_cm2"] * 1e4  # W/m^2 at the atoms, each beam
        wz_wall = wz0 * np.sqrt(1 + (cfg["wall_distance_m"] / zr) ** 2)
        optics = {
            "cell wall (through-cell variant)": i1 * wz0 / wz_wall,
            "primary P1": i1 * wz0 / vd["w_z_at_primary_m"],
            "secondary P2": i1 * ratio_S * wz0 / vd["w_z_at_secondary_m"],
            "graded mirror G": i1 * ratio_G * wz0 / vd["w_z_at_graded_m"],
        }
        geo_v = c["geometry"]["vertical"]
        # placeholder for a mirror with 3x expansion in both directions; the adopted vertical ring is
        # computed in scripts/cavity_vertical_mirrors.py (docs/cavity_concept.md, Section 5)
        optics_v = {"vertical ring at the atoms": i1, "vertical ring, mirror with 3x expansion": i1 / 9.0}
        row = {"match": c["match"], "depth_mK": c["depth_mK"], "I_atoms_W_cm2": i1 * 1e-4,
               "P_circ_inplane_W": pc["inplane_per_axis_W"], "P_circ_vertical_W": pc["vertical_W"],
               "intensity_W_cm2": {k: v * 1e-4 for k, v in {**optics, **optics_v}.items()},
               "max_defect_radius_um": {}}
        for sub in ("fused_silica", "sapphire"):
            kap = SUBSTRATES[sub].kappa
            row["max_defect_radius_um"][sub] = {k: cfg["defect_dT_K"] * kap / v * 1e6 for k, v in {**optics, **optics_v}.items()}
        rep["cases"].append(row)
    # Thermal deformation of the graded mirror at the reference (equal-depth, 0.2 mK) power.
    ref = [r for r in rep["cases"] if r["match"] == "equal_depth" and r["depth_mK"] == cfg["depths_mK"][0]][0]
    p_c = ref["P_circ_inplane_W"]
    # Equal-area disc of the spot on G: horizontal equivalent width M W_x (flat level / peak ratio), vertical Gaussian.
    w_x_G = w_x / ratio_G
    area = w_x_G * np.sqrt(np.pi / 2) * vd["w_z_at_graded_m"]
    a_eq = np.sqrt(area / np.pi)
    for case in cfg["deformation_cases"]:
        sub = SUBSTRATES[case["substrate"]]
        prob = ThermalProblem(radius_m=cfg["mirror_radius_m"], thickness_m=cfg["mirror_thickness_m"], substrate=sub,
                              beam="top_hat", beam_radius_m=a_eq, h_faces=10.0, rim="fixed", nr=240, nz=60, z_power=2.0)
        rc, zc, T = solve(prob, p_c, case["absorption"])
        ts = surface_temperature(zc, T)
        dz = np.diff(prob.thickness_m * np.linspace(0, 1, len(zc) + 1) ** prob.z_power)
        u = (1 + sub.poisson) * sub.alpha_exp * np.sum(T * dz[None, :], axis=1)
        u_rel = u - np.interp(a_eq, rc, u)
        sel = rc <= a_eq
        w_r = rc[sel] * np.ones_like(rc[sel])
        u_beam = u[sel] - np.sum(u[sel] * w_r) / np.sum(w_r)
        rep["deformation"].append({
            **case, "P_circ_W": p_c, "P_abs_W": p_c * case["absorption"], "spot_equivalent_radius_um": a_eq * 1e6,
            "dT_surface_max_K": float(ts.max()), "bulge_centre_minus_edge_nm": float(u_rel[0] * 1e9),
            "bulge_rms_over_spot_waves": float(np.sqrt(np.sum(u_beam**2 * w_r) / np.sum(w_r)) / lam_m),
        })
    write_json(out / "thermal.json", rep)
    for r in rep["cases"]:
        print(r["match"], r["depth_mK"], {k: round(v) for k, v in r["intensity_W_cm2"].items()})
        print("   defect FS", {k: round(v, 2) for k, v in r["max_defect_radius_um"]["fused_silica"].items()})
        print("   defect sapphire", {k: round(v, 1) for k, v in r["max_defect_radius_um"]["sapphire"].items()})
    for d in rep["deformation"]:
        print(d)


if __name__ == "__main__":
    main(sys.argv[1])
