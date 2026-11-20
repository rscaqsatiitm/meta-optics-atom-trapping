"""Absorbing-defect limits of the vertical-ring mirrors from stored peak intensities.

Usage: python scripts/cavity_defect_radii.py configs/cavity_defect_radii.yaml

For each mirrors.json listed (scripts/cavity_vertical_mirrors.py), the peak intensities on the
mirrors (kW/cm^2, from the Fox-Li modes, independent of the substrate) are converted into the
largest fully absorbing defect radius for a defect_dT_K rise on a half-space,
a = defect_dT_K x kappa / I, with kappa from metatrap.thermal.axisym_fd (DERIVED from SIM).
Output: results/cavity_defect_radii/defect_radii.json (per file, variant, surface and depth; and
the smallest radius per variant with the surface that sets it).
"""

import json
import sys

from metatrap.runio import REPO_ROOT, load_config, run_dir, write_json, write_params
from metatrap.thermal.axisym_fd import SUBSTRATES


def main(cfg_path):
    cfg = load_config(cfg_path)
    out = run_dir(cfg["run_id"])
    res = {"kappa_W_mK": {s: SUBSTRATES[s].kappa for s in ("fused_silica", "sapphire")}, "files": {}}
    for f in cfg["mirrors_json"]:
        m = json.load(open(REPO_ROOT / f))
        fr = {}
        for vn, v in m["variants"].items():
            radii = {sub: {nm: {d: cfg["defect_dT_K"] * SUBSTRATES[sub].kappa / (val * 1e7) * 1e6 for d, val in vals.items()}
                           for nm, vals in v["peak_intensity_kW_cm2"].items()} for sub in ("fused_silica", "sapphire")}
            smallest = {}
            for sub, per in radii.items():
                for nm, vals in per.items():
                    for d, a in vals.items():
                        key = f"{sub}|{d}"
                        if key not in smallest or a < smallest[key][0]:
                            smallest[key] = (a, nm)
            fr[vn] = {"max_defect_radius_um": radii, "smallest": {k: {"radius_um": a, "surface": nm} for k, (a, nm) in smallest.items()}}
            print(f, vn, {k: (round(a, 4), nm) for k, (a, nm) in smallest.items()})
        res["files"][f] = fr
    write_json(out / "defect_radii.json", res)
    write_params(out, cfg)


if __name__ == "__main__":
    main(sys.argv[1])
