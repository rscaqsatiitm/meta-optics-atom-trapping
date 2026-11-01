"""Aperture-averaged supercell efficiency relative to the local-periodic approximation.

Usage: python scripts/edge_factor.py results/library_analysis_1075/supercell.json configs/hologram_b1_cell20mm.yaml

The ratio eta_RCWA / eta_LPA (mean of x and y polarisation) is known at the discrete
deflections sin(theta) = lambda / (N p) of the blazed supercells. It is interpolated
linearly in sin(theta), with ratio = 1 at sin(theta) = 0 (a uniform array is exact in the
LPA) and held constant beyond the largest computed angle, then averaged over the
metasurface with the incident intensity as weight, using the ray mapping r(s) of the
layer stack. The x/y average approximates the azimuthal average of a linearly polarised
beam, for which the local polarisation relative to the deflection plane rotates around
the aperture.
"""

import json
import sys

import numpy as np

from metatrap.propagation.illumination import profile
from metatrap.propagation.layers import LayerStack
from metatrap.runio import load_config, write_json, REPO_ROOT


def main(sc_path, cfg_path):
    sc = json.load(open(sc_path))
    cfg = load_config(cfg_path)
    lam, na = cfg["wavelength_um"], cfg["na"]
    stack = LayerStack.from_list(cfg["stack"])
    R = stack.aperture_for_na(na, lam)
    r = np.linspace(0, R, 4001)
    s = stack.s_of_r(lam)(r)
    out = {}
    for mat, rows in sc.items():
        by_s = {}
        for row in rows:
            by_s.setdefault(round(row["sin_theta"], 6), []).append(row["eta_rcwa"] / row["eta_lpa"])
        ss = np.array(sorted(by_s))
        ratio = np.array([np.mean(by_s[k]) for k in ss])
        f = np.interp(s, np.concatenate([[0.0], ss]), np.concatenate([[1.0], ratio]))
        res = {"sin_theta": ss, "ratio_xy_mean": ratio}
        for ill in cfg["illuminations"]:
            if ill["kind"] == "top_hat":
                tag = "top_hat" if not ill.get("edge_frac") else f"top_hat_edge{ill['edge_frac']:g}"
                amp = profile("top_hat", r, R, edge_um=ill.get("edge_frac", 0.0) * R)
            else:
                tag = f"gauss_w{ill['waist_over_radius']}"
                amp = profile("gaussian", r, R, waist_um=ill["waist_over_radius"] * R)
            w = amp**2 * r
            res[tag] = float(np.sum(w * f) / np.sum(w))
        out[mat] = res
        print(mat, {k: v for k, v in res.items() if not isinstance(v, np.ndarray)}, "ratios", np.round(ratio, 3), "at s", ss)
    write_json(REPO_ROOT / "results" / "library_analysis_1075" / "edge_factor.json", out)


if __name__ == "__main__":
    main(*sys.argv[1:3])
