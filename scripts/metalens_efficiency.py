"""Focusing efficiency estimate of the doublet plates from blazed-supercell RCWA.

Usage: python scripts/metalens_efficiency.py configs/metalens_doublet.yaml 0.8523:results/library_analysis_852 1.0644:results/library_analysis_1064

For each wavelength and film, the first-order efficiency eta_RCWA(s) of the library's chosen
pillar family (x and y polarisation averaged, approximating the azimuthal average of a
linearly polarised beam) is interpolated linearly in the deflection s = |dphi/dr| / k0, with
eta(0) = the library's LPA efficiency (a uniform array is exact in the LPA) and held constant
beyond the largest computed deflection. It is averaged over each plate with the incident
intensity (uniform over the stop, mapped to the focusing plate along the design rays). The
doublet efficiency is the product of the two plate averages; this neglects correlation
between the zones of the two plates (a ray at h1 meets zone h2(h1), so the product of
area averages slightly overestimates the loss where both are large; the corrector is
weak, so the error is small). Output: results/<run_id>/efficiency.json.
"""

import json
import sys

import numpy as np

from metatrap.metalens import doublet as D
from metatrap.propagation.layers import LayerStack
from metatrap.runio import REPO_ROOT, load_config, write_json


def eta_curve(rows, eta0):
    by_s = {}
    for r in rows:
        by_s.setdefault(round(r["sin_theta"], 6), []).append(r["eta_rcwa"])
    ss = np.array(sorted(by_s))
    eta = np.array([np.mean(by_s[k]) for k in ss])
    return np.concatenate([[0.0], ss]), np.concatenate([[eta0], eta])


def main(cfg_path, *pairs):
    cfg = load_config(cfg_path)
    spec = D.DoubletSpec(cfg["na"], cfg["cell"], [list(l) for l in cfg["spacer"]])
    d0 = D.make(spec, cfg["wavelengths_um"][0])
    R1 = d0.info["R1_um"]
    out = {}
    for pr in pairs:
        lam_s, path = pr.split(":", 1)
        lam = float(lam_s)
        d = D.make(spec, lam, R1_um=R1)
        k0 = 2 * np.pi / lam
        h1 = np.linspace(0, R1, 4001)
        h2 = LayerStack.from_list(cfg["cell"]).ray_height(h1 / d.info["f_um"], lam)
        s1 = np.abs(d.m1.dphase_dr(h1, lam)) / k0
        s2 = np.abs(d.m2.dphase_dr(h2, lam)) / k0
        w = h1  # uniform intensity over the stop, area weight 2 pi h1 dh1 (followed along rays)
        sc = json.load(open(REPO_ROOT / path / "supercell.json"))
        sel = json.load(open(REPO_ROOT / path / "selection.json"))
        res = {}
        for mat, rows in sc.items():
            ss, eta = eta_curve(rows, sel[mat]["eta"])
            e1 = float(np.sum(w * np.interp(s1, ss, eta)) / np.sum(w))
            e2 = float(np.sum(w * np.interp(s2, ss, eta)) / np.sum(w))
            # Ring-resolved efficiency of the focusing plate (outer rings carry the high-NA light).
            rings = {f"{a:.2f}-{b:.2f} R1": float(np.sum((w * np.interp(s2, ss, eta))[(h1 >= a * R1) & (h1 < b * R1)]) / np.sum(w[(h1 >= a * R1) & (h1 < b * R1)]))
                     for a, b in ((0, 0.5), (0.5, 0.8), (0.8, 0.9), (0.9, 1.0))}
            res[mat] = {"family": sel[mat], "sin_theta": ss.tolist(), "eta_xy_mean": eta.tolist(), "corrector": e1, "focusing": e2,
                        "doublet": e1 * e2, "focusing_rings": rings}
            print(lam, mat, {k: (round(v, 4) if isinstance(v, float) else v) for k, v in res[mat].items() if k in ("corrector", "focusing", "doublet")},
                  {k: round(v, 3) for k, v in rings.items()}, "eta(s):", np.round(eta, 3), "at", np.round(ss, 3), flush=True)
        out[str(lam)] = res
    write_json(REPO_ROOT / "results" / cfg["run_id"] / "efficiency.json", out)


if __name__ == "__main__":
    main(*sys.argv[1:])
