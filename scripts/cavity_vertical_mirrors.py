"""Peak intensities on the mirrors of the vertical-lattice ring, and a z-expanding telescope.
Usage:
    python scripts/cavity_vertical_mirrors.py configs/cavity_vertical_mirrors.yaml

The vertical ring's mode is the product of two 1D modes (scripts/cavity_vertical.py): x (6.9 mm
flat-top, 0.15x cylindrical telescope) and z (order-10 flat-top of 0.19 mm half-height, no telescope
in the adopted design). The peak intensity of the circulating beam on a mirror is
P_circ max|ex|^2 / int|ex|^2 max|ez|^2 / int|ez|^2 (beam cross-section; the surface value is lower
by cos AOI). Variants add a Keplerian telescope in z at the same mirror positions (magnification
mag > 1, f1z + f2z = f1 + f2) to spread the beam vertically on P2 and G. For each variant: ideal
depth uniformity and enhancement, the z spectrum, lambda/1000 z errors (uncorrected and
perfect-metrology redesign), peak intensities and absorbing-defect limits (criterion of
scripts/cavity_thermal.py). Writes mirrors.json.
"""

import json
import sys

import numpy as np

from metatrap.cavity import foxli as FL
from metatrap.runio import REPO_ROOT, load_config, run_dir, write_json, write_params
from metatrap.thermal.axisym_fd import SUBSTRATES

sys.path.insert(0, str(REPO_ROOT / "scripts"))
import cavity_correction as C  # noqa: E402
import cavity_vertical as V  # noqa: E402


def stop_fields(s):
    """Normalised peak intensity (1/um) of the ideal mode at X, at every stop of half trip 1 and at G."""
    r, g = s["ring"], s["g"]
    F = s["F"].astype(complex)
    out = {"X": 0.0}
    for zs, _ in r._stops_half[1]:
        out[zs] = 0.0
    out[r.p_um] = 0.0
    for z in out:
        u = F if z == "X" else r._segment(F, 0.0, z - 1e-3, 1)
        i = np.abs(u) ** 2
        out[z] = float(i.max() / (i.sum() * g.dx))
    return out


def main(cfg_path):
    cfg = load_config(cfg_path)
    out = run_dir(cfg["run_id"])
    write_params(out, cfg)
    vcfg = load_config(REPO_ROOT / cfg["vertical_config"])
    tcfg = load_config(REPO_ROOT / cfg["thermal_config"])
    lay = json.loads((REPO_ROOT / cfg["layout_json"]).read_text())["chosen"]
    bud = json.loads((REPO_ROOT / cfg["budget_json"]).read_text())
    p_circ = {c["depth_mK"]: c["circulating"]["vertical_W"] for c in bud["cases"]
              if c["n_sites"] == cfg["n_sites"] and c["match"] == "equal_depth" and c["depth_mK"] in cfg["depths_mK"]}
    t = vcfg["telescope"]
    p = t["d0_um"] + t["f1_um"] * (1 + t["mag"]) + t["L_um"]
    names = {"X": "atoms (X)", t["d0_um"]: "P1", t["d0_um"] + t["f1_um"] * (1 + t["mag"]): "P2", p: "G",
             p + t["L_um"]: "P2 (return)", 2 * p - t["d0_um"]: "P1 (return)"}
    aoi = {"P1": max(lay["aoi_deg"]["P1"], lay["aoi_deg"]["P3"]), "P2": min(lay["aoi_deg"]["P2"], lay["aoi_deg"]["P4"]),
           "G": min(lay["aoi_deg"]["G1"], lay["aoi_deg"]["G2"])}
    sx = V.build(vcfg, vcfg["x"])
    ex0 = V.solve(vcfg, sx, sx["ring"])
    ix = stop_fields(sx)
    rep = {"P_circ_W": p_circ, "aoi_deg_smallest": aoi, "variants": {}}
    for var in cfg["z_variants"]:
        dz = dict(vcfg["z"], **{k: v for k, v in var.items() if k != "name"})
        sz = V.build(vcfg, dz)
        ez0 = V.solve(vcfg, sz, sz["ring"])
        res = {"ideal": dict(V.depth_map(vcfg, sx, ex0[:2], sz, ez0[:2]), g_ratio_x=ex0[2], g_ratio_z=ez0[2])}
        sp = V.build(vcfg, dict(dz, grid_n=max(cfg["spectrum_grid_n"], dz["grid_n"] // 4)))
        res["spectrum_z"] = FL.dense_spectrum(sp["ring"], sp["F"], n_modes=4)
        iz = stop_fields(sz)
        res["peak_intensity_kW_cm2"] = {}
        for z, nm in names.items():
            if z in ix and z in iz:
                res["peak_intensity_kW_cm2"][nm] = {str(d): P * ix[z] * iz[z] * 1e8 / 1e3 for d, P in p_circ.items()}
        res["max_defect_radius_um"] = {
            sub: {nm: {d: tcfg["defect_dT_K"] * SUBSTRATES[sub].kappa / (v * 1e7) * 1e6 for d, v in vals.items()}
                  for nm, vals in res["peak_intensity_kW_cm2"].items()} for sub in ("fused_silica", "sapphire")}
        for seed in cfg["error_seeds"]:
            rng = np.random.default_rng(seed)
            changes, g_err = C.make_errors(sz["ring"], sz["F"], sz["w"], rng, vcfg["fig_waves"], vcfg["corr_um"] * sz["w"] / sx["w"], 0.0)
            true = sz["ring"].perturbed(changes)
            for label in ("uncorrected", "openloop_perfect"):
                r = true.perturbed()
                if label == "uncorrected":
                    r.phase_g1, r.phase_g2 = sz["ph"] + g_err[1], sz["ph"] + g_err[2]
                else:
                    r.phase_g1, r.phase_g2 = true.design_phase(sz["F"], sz["F"], 1), true.design_phase(sz["F"], sz["F"], 2)
                e = V.solve(vcfg, sz, r)
                res[f"fig{vcfg['fig_waves']}_z_seed{seed}_{label}"] = dict(V.depth_map(vcfg, sx, ex0[:2], sz, e[:2]), g_ratio=e[2])
        rep["variants"][var["name"]] = res
        print(var["name"], json.dumps({k: v for k, v in res.items() if k != "max_defect_radius_um"}, default=float)[:1500], flush=True)
    write_json(out / "mirrors.json", rep)


if __name__ == "__main__":
    main(sys.argv[1])
