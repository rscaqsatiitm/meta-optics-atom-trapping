"""Blazed-supercell efficiency at each operating wavelength for atoms chosen at lambda_d.

Usage: python scripts/multiwave_supercells.py configs/multiwave.yaml

For the family in cfg["supercell"], N-pillar supercells are built with widths chosen for the
phases 2 pi k / N at lambda_d (from the curve computed as in multiwave_atoms.py). They are
solved by RCWA at each operating wavelength, where they deflect to sin(theta) = lambda / (N p),
for x and y polarisation. The first-order efficiency versus deflection at each wavelength feeds
the doublet efficiency estimate (scripts/multiwave_efficiency in this file, --doublet).
Output: results/<run_id>/supercells.json (checkpointed per case).
"""

import json
import sys

import numpy as np

from metatrap import materials
from metatrap.rcwa.cell import PillarCell
from metatrap.rcwa.lookup import FabLimits, densify, feasible
from metatrap.rcwa.solver import solve_grating_orders, solve_pillars
from metatrap.rcwa.supercell import supercell_eps
from metatrap.runio import load_config, run_dir, write_json


def main(cfg_path):
    cfg = load_config(cfg_path)
    sc = cfg["supercell"]
    out = run_dir(cfg["run_id"])
    path = out / "supercells.json"
    res = json.load(open(path)) if path.exists() else {}
    mat, p, h, lam_d = sc["material"], sc["pitch_um"], sc["height_um"], sc["design_um"]
    lim = FabLimits(**cfg["fab_limits"])
    wf = np.linspace(cfg["width_frac"]["min"], cfg["width_frac"]["max"], cfg["width_frac"]["n"])
    w = wf * p
    ns_d = float(materials.get("fused_silica").n(lam_d))
    r = solve_pillars(PillarCell(p, h), w, lam_d, complex(materials.get(mat).eps(lam_d)), ns_d**2, 1.0, 1.0, cfg["n_terms"], cfg["n_grid"])
    ok = feasible(w, p, h, lim)
    wd, td = densify(w[ok], r.t_x[ok])
    for lam in cfg["operating_um"]:
        eps_p = complex(materials.get(mat).eps(lam))
        ns = float(materials.get("fused_silica").n(lam))
        for N in sc["n_cells"]:
            for pol in ("x", "y"):
                key = f"{lam}:{N}:{pol}"
                if key in res:
                    continue
                ph = 2 * np.pi * np.arange(N) / N
                k = [int(np.argmax(np.real(td * np.exp(-1j * q)))) for q in ph]
                widths = wd[k]
                eps = supercell_eps(widths, p, eps_p, 1.0, sc["grid_per_cell"])
                orders, eff, R = solve_grating_orders(eps, N * p, p, h, lam, ns**2, 1.0, n_terms=sc["terms_per_cell"] * N, polarisation=pol)
                e = {tuple(int(v) for v in o): float(x) for o, x in zip(orders, eff)}
                res[key] = {"lambda_um": lam, "n_cells": N, "polarisation": pol, "sin_theta": lam / (N * p), "eta_plus1": e.get((1, 0), 0.0),
                            "eta_minus1": e.get((-1, 0), 0.0), "eta_zero": e.get((0, 0), 0.0), "R": float(R), "widths_um": widths.tolist()}
                print(key, {k2: round(v, 4) for k2, v in res[key].items() if isinstance(v, float)}, flush=True)
                write_json(path, res)


if __name__ == "__main__":
    main(sys.argv[1])
