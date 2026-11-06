"""Mechanical tolerances of the 42 mm imaging meta-objective over the 10^6-site field.

Usage: python scripts/million_tolerances.py configs/million.yaml
Air-spaced aplanatic doublet, NA 0.5, 852.3 nm, air gap 25 mm (f = 42 mm), construction to 1.6 R1.
Field points -2828 / 0 / +2828 um; one common plane re-optimised per case.
Output: results/million/tolerances_f42.json.
"""

import json
import sys

import numpy as np

from metatrap.metalens import doublet as D
from metatrap.propagation.layers import LayerStack
from metatrap.runio import REPO_ROOT, load_config, write_json


def main(cfg_path):
    cfg = load_config(cfg_path)
    tc = cfg["tolerances_f42"]
    lam, na = 0.8523, 0.5
    cell = [[tc["air_gap_um"], "air"], [4000, "fused_silica"], [9000, "air"]]
    f0 = LayerStack.from_list(cell).aperture_for_na(na, lam) / na
    d = D.make(D.DoubletSpec(na, cell, [[f0, "air"]], n_pts=2001, ext=1.6), lam)
    R1, f = d.info["R1_um"], d.info["f_um"]
    sn = D.system(d, lam)
    hs = tc["heights_um"]
    th = [D.theta_for_height(sn, R1, lam, h, f) for h in hs]
    path = REPO_ROOT / "results" / "million" / "tolerances_f42.json"
    res = {}
    for case in tc["cases"]:
        pert = D.Perturbation(**{k: v for k, v in case.items() if k != "name"})
        dz, q = D.best_plane(D.system(d, lam, pert), R1, lam, th, n_rays=32, span_um=8.0)
        res[case["name"]] = {"strehl": [x["strehl"] for x in q], "min": min(x["strehl"] for x in q), "dz_um": dz, "case": case}
        print(case["name"], round(dz, 2), np.round([x["strehl"] for x in q], 3), flush=True)
        write_json(path, res)


if __name__ == "__main__":
    main(sys.argv[1])
