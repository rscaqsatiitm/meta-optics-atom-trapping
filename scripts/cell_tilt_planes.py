"""Cell-tilt tolerance with the atom plane fixed to the cell or perpendicular to the objective.

Usage: python scripts/cell_tilt_planes.py

The cell-tilt rows of the metalens report (16 mm doublet, +/-400 um) and of the 10^6 report
(42 mm objective, +/-2.83 mm) were computed with the evaluation plane tilted with the cell
(doublet.Perturbation, cell_tilt_plane="cell"): at field height h the atom plane then lies
h x tilt away from the objective's focal surface. When the atom plane is set by the optics
(traps made through the same objective), the evaluation plane is perpendicular to the
objective axis (cell_tilt_plane="objective") and only the aberration of the tilted wall
remains. Both are computed here, with one common plane re-optimised per case (best_plane), as
in the original tables. Output: results/cell_tilt_planes.json.
"""

import importlib.util
from pathlib import Path

import numpy as np

from metatrap.metalens import doublet as D
from metatrap.propagation.layers import LayerStack
from metatrap.runio import REPO_ROOT, load_config, write_json


def doublet_16mm():
    spec_mod = importlib.util.spec_from_file_location("md", Path(__file__).with_name("metalens_doublet.py"))
    md = importlib.util.module_from_spec(spec_mod)
    spec_mod.loader.exec_module(md)
    cfg = load_config(REPO_ROOT / "configs" / "metalens_doublet.yaml")
    lam = cfg["wavelengths_um"][0]
    d = D.make(md.spec_for(cfg), lam)
    return d, lam, cfg["tolerance_heights_um"], {}


def objective_42mm():
    lam, na = 0.8523, 0.5
    cell = [[25000, "air"], [4000, "fused_silica"], [9000, "air"]]
    f0 = LayerStack.from_list(cell).aperture_for_na(na, lam) / na
    d = D.make(D.DoubletSpec(na, cell, [[f0, "air"]], n_pts=2001, ext=1.6), lam)
    return d, lam, [-2828, 0, 2828], {"n_rays": 32, "span_um": 8.0}


def main():
    res = {}
    for name, build, tilts_deg in (("doublet_16mm", doublet_16mm, [0.05, 0.1, 0.25, 0.5]), ("objective_42mm", objective_42mm, [0.05, 0.1])):
        d, lam, hs, kw = build()
        R1, f = d.info["R1_um"], d.info["f_um"]
        sn = D.system(d, lam)
        th = [D.theta_for_height(sn, R1, lam, h, f) for h in hs]
        rows = {"heights_um": hs}
        for t in tilts_deg:
            for plane in ("cell", "objective"):
                pert = D.Perturbation(cell_tilt_rad=float(np.deg2rad(t)), cell_tilt_plane=plane)
                dz, q = D.best_plane(D.system(d, lam, pert), R1, lam, th, **kw)
                rows[f"{t:g}deg:{plane}"] = {"strehl": [x["strehl"] for x in q], "min": min(x["strehl"] for x in q), "dz_um": dz}
                print(name, t, plane, round(dz, 2), np.round([x["strehl"] for x in q], 3), flush=True)
        res[name] = rows
    write_json(REPO_ROOT / "results" / "cell_tilt_planes.json", res)


if __name__ == "__main__":
    main()
