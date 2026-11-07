"""Film absorption of the selected hole and pillar families (lossy film, k from the material model).

Usage: python scripts/hole_vs_pillar_absorption.py configs/hole_vs_pillar.yaml

For the selected family of each kind (results/hole_vs_pillar/ranking.json, baseline and relaxed
fabrication limits), the 16 meta-atoms of the 16-level blazed grating are recomputed with the
complex film index n + ik and the absorbed fraction A = 1 - T - R is averaged over them. The
global heating of the metasurface scales with this fraction (docs/campaign1_report.md, Section 5).
Output: results/hole_vs_pillar/absorption.json.
"""

import json
import sys

import numpy as np

from metatrap import materials as M
from metatrap.rcwa import lookup as L
from metatrap.rcwa.cell import PillarCell
from metatrap.rcwa.solver import solve_eps
from metatrap.runio import REPO_ROOT, load_config, write_json


def main(cfg_path):
    cfg = load_config(cfg_path)
    d = REPO_ROOT / "results" / cfg["run_id"]
    sweep, rank = json.load(open(d / "sweep.json")), json.load(open(d / "ranking.json"))
    lam = cfg["wavelength_um"]
    mat = M.get(cfg["material"])
    n_c = complex(mat.n(lam), mat.k(lam))
    n_s = float(M.get(cfg["substrate"]).n(lam))
    out = {"k": float(mat.k(lam))}
    for lim_name, lim_cfg in cfg["fab_limits"].items():
        lim = L.FabLimits(**lim_cfg)
        for fam in ("pillar", "hole"):
            sel = rank[lim_name][fam]["selected_shortest_within_0.3pct"]
            p, h = sel["pitch_um"], sel["height_um"]
            r = sweep[f"{fam}:{p:.2f}:{h:.2f}"]
            w = np.array(r["widths_um"])
            t = np.array(r["t_re"]) + 1j * np.array(r["t_im"])
            ok = L.feasible(w, p, h, lim)
            wd, td = L.densify(w[ok], t[ok])
            off = L.eta_n(wd, td, 1.0 / n_s, 16)[1]
            w16 = wd[L.best_widths(wd, td, 2 * np.pi * np.arange(16) / 16, off)]
            if fam == "hole":
                eps = n_c**2 + (1.0 - n_c**2) * PillarCell(pitch=p, height=h, shape="circle").mask(w16, cfg["n_grid"])
            else:
                eps = 1.0 + (n_c**2 - 1.0) * PillarCell(pitch=p, height=h, shape="square").mask(w16, cfg["n_grid"])
            res = solve_eps(eps.astype(complex), p, h, lam, n_s**2, 1.0, n_terms=cfg["n_terms"], formulation=cfg["formulation"])
            A = 1 - np.real(res.T) - np.real(res.R)
            out[f"{lim_name}:{fam}"] = {"pitch_um": p, "height_um": h, "widths_um": w16.tolist(), "absorbed_mean": float(np.mean(A)),
                                       "absorbed_max": float(np.max(A)), "bulk_film_equivalent": float(4 * np.pi * mat.k(lam) * h / lam)}
            print(lim_name, fam, f"p {p} h {h}: absorbed mean {np.mean(A):.2e}, max {np.max(A):.2e}; unpatterned film of this thickness {4*np.pi*mat.k(lam)*h/lam:.2e}", flush=True)
    write_json(d / "absorption.json", out)


if __name__ == "__main__":
    main(sys.argv[1])
