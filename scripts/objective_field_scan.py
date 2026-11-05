"""Field radius of the air-spaced aplanatic meta-objective versus focal length and NA.

Usage: python scripts/objective_field_scan.py configs/million.yaml

For each NA and air gap between the focusing plate and the cell (which sets the focal length:
R1 = R2 = H(NA) through the stack, f = R1 / NA), the plate spacing is scanned around f, and for
the best spacing the largest field radius with a minimum common-plane Strehl ratio >= 0.8 and
>= 0.5 (all field points from 0 to r in one focal plane) is found by bisection. Output:
results/million/objective_field.json.
"""

import json
import sys

import numpy as np

from metatrap.metalens import doublet as D
from metatrap.runio import REPO_ROOT, load_config, write_json


def min_strehl(s, R1, lam, f, r, n=6):
    H = list(np.linspace(0, r, n))
    dz, q = D.common_plane(s, R1, lam, H, f, n_rays=28)
    return min(x["strehl"] for x in q)


def field_radius(s, R1, lam, f, thr, r_hi):
    lo, hi = 0.0, r_hi
    if min_strehl(s, R1, lam, f, hi) >= thr:
        return hi
    for _ in range(9):
        mid = 0.5 * (lo + hi)
        if min_strehl(s, R1, lam, f, mid) >= thr:
            lo = mid
        else:
            hi = mid
    return lo


def main(cfg_path):
    cfg = load_config(cfg_path)
    oc = cfg["objective_scan"]
    path = REPO_ROOT / "results" / "million" / "objective_field.json"
    res = json.load(open(path)) if path.exists() else {}
    lam = oc["wavelength_um"]
    for na in oc["na"]:
        for gap in oc["air_gap_um"]:
            key = f"NA{na}_gap{gap}"
            if key in res:
                continue
            cell = [[gap, "air"]] + [list(l) for l in oc["cell_rest"]]
            ext = min(oc.get("ext", 1.1), 0.95 / na)  # keep the construction below grazing directions
            from metatrap.propagation.layers import LayerStack

            f0 = LayerStack.from_list(cell).aperture_for_na(na, lam) / na
            best = None
            for frac in oc["spacing_fracs"]:
                try:
                    d = D.make(D.DoubletSpec(na, cell, [[frac * f0, "air"]], n_pts=2001, ext=ext), lam)
                except Exception as e:
                    print(key, frac, "failed:", e, flush=True)
                    continue
                s = D.system(d, lam)
                m = min_strehl(s, d.info["R1_um"], lam, d.info["f_um"], 0.06 * f0)
                if best is None or m > best[0]:
                    best = (m, frac, d, s)
            _, frac, d, s = best
            R1, f = d.info["R1_um"], d.info["f_um"]
            r08 = field_radius(s, R1, lam, f, 0.8, 0.3 * f)
            r05 = field_radius(s, R1, lam, f, 0.5, 0.3 * f)
            res[key] = {"na": na, "air_gap_um": gap, "f_um": f, "aperture_diameter_mm": 2 * R1 / 1000, "spacing_um": frac * f0,
                        "r_strehl08_um": r08, "r_strehl05_um": r05,
                        "sites_4um_08": 2 * r08**2 / 16, "sites_3um_08": 2 * r08**2 / 9, "sites_4um_05": 2 * r05**2 / 16}
            print(key, {k: (round(v, 1) if isinstance(v, float) else v) for k, v in res[key].items()}, flush=True)
            write_json(path, res)


if __name__ == "__main__":
    main(sys.argv[1])
