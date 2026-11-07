"""Exact-kernel focus of single hologram sites far off axis at the present cell geometry.

Usage: python scripts/hologram_large_field.py configs/scaling_100k.yaml

Same model as scripts/scaling_100k.py (aperture of NA 0.55 focusing through a 15.5 mm vacuum
equivalent of the cell, exact per-site kernel, so no field aberration), for site radii out to
6.4 mm. Reports the peak relative to the on-axis site (D), the 1/e^2 half-widths (radial,
tangential), the largest local deflection s_max, and the meta-atom efficiency factors of the
one-family supercell curve (eta: peak factor; eta_power: transmitted fraction).
Output: results/scaling_100k/hologram_large_field.json. (This file was first written by an
uncommitted inline run with a two-family efficiency curve; this script reproduces it.)
"""

import sys

import numpy as np

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from scaling_100k import eta_curve, spot  # noqa: E402

from metatrap.runio import REPO_ROOT, load_config, write_json  # noqa: E402

WD = 15500.0
RADII = [0, 1600, 2300, 3200, 4500, 6400]


def main(cfg_path):
    cfg = load_config(cfg_path)
    lam, na = cfg["wavelength_um"], cfg["na"]
    R = WD * na / np.sqrt(1 - na**2)
    eta_s = eta_curve(cfg)
    out, p0 = {}, None
    for r in RADII:
        q = spot(WD, R, float(r), lam, cfg["n_aperture"], eta_s)
        p0 = q["peak"] if r == 0 else p0
        out[str(r)] = {"D": q["peak"] / p0, "w": [q["w_radial_um"], q["w_tangential_um"]], "s_max": q["s_max"],
                       "eta": q["eta_site"], "eta_power": q["eta_site_power"]}
        print(r, {k: (np.round(v, 4).tolist() if isinstance(v, list) else round(v, 4)) for k, v in out[str(r)].items()}, flush=True)
    write_json(REPO_ROOT / "results" / "scaling_100k" / "hologram_large_field.json", out)


if __name__ == "__main__":
    main(sys.argv[1])
