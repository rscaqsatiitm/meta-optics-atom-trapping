"""Single-pillar responses at a shared design wavelength and the cross-wavelength efficiency.

Usage: python scripts/multiwave_atoms.py configs/multiwave.yaml

For each candidate pillar family (material, pitch, height) the square-pillar transmission is
solved at the design wavelength(s) and at every operating wavelength (same numerics as
library_main: 197 terms, 256 grid). Atoms are chosen for the target phase at the design
wavelength; at each operating wavelength the first-order efficiency, the zero-order power and
the transmitted power follow (multiwave.crosswave_efficiency, here on the computed curves).
Output: results/<run_id>/atoms.json.
"""

import sys

import numpy as np

from metatrap import materials
from metatrap.rcwa.cell import PillarCell
from metatrap.rcwa.lookup import FabLimits, densify, feasible
from metatrap.rcwa.solver import solve_pillars
from metatrap.runio import load_config, run_dir, write_json, write_params


def curve(mat, p, h, lam, wf, n_terms, n_grid):
    w = wf * p
    r = solve_pillars(PillarCell(p, h), w, lam, complex(materials.get(mat).eps(lam)), float(materials.get("fused_silica").n(lam)) ** 2,
                      1.0, 1.0, n_terms, n_grid)
    return w, r.t_x, r.T0


def cross(curves, lam_d, lams, ok, n_theta=720):
    w, td, T0d = curves[lam_d]
    th = 2 * np.pi * (np.arange(n_theta) + 0.5) / n_theta
    wd, tdd = densify(w[ok], td[ok])
    k = np.argmax(np.real(tdd[None, :] * np.exp(-1j * th)[:, None]), axis=1)
    out = {}
    for lam in lams:
        w2, t2, T02 = curves[lam]
        _, td2 = densify(w2[ok], t2[ok])
        nr = float(np.median(T02[ok] / np.maximum(np.abs(t2[ok]) ** 2, 1e-12)))
        tt = td2[k]
        out[str(lam)] = {"eta_first": nr * float(abs(np.mean(tt * np.exp(-1j * th))) ** 2), "zero_order": nr * float(abs(np.mean(tt)) ** 2),
                         "transmitted": nr * float(np.mean(np.abs(tt) ** 2))}
    return out


def main(cfg_path):
    cfg = load_config(cfg_path)
    out = run_dir(cfg["run_id"])
    write_params(out, cfg)
    lim = FabLimits(**cfg["fab_limits"])
    wf = np.linspace(cfg["width_frac"]["min"], cfg["width_frac"]["max"], cfg["width_frac"]["n"])
    lams = cfg["operating_um"]
    res = []
    for fam in cfg["families"]:
        mat, p, h = fam["material"], fam["pitch_um"], fam["height_um"]
        curves = {lam: curve(mat, p, h, lam, wf, cfg["n_terms"], cfg["n_grid"]) for lam in sorted(set(lams + cfg["design_um"]))}
        ok = feasible(wf * p, p, h, lim)
        for lam_d in cfg["design_um"]:
            row = {"material": mat, "pitch_um": p, "height_um": h, "design_um": lam_d, "eff": cross(curves, lam_d, lams, ok)}
            res.append(row)
            print(mat, p, h, lam_d, {k: (round(v["eta_first"], 3), round(v["zero_order"], 3)) for k, v in row["eff"].items()}, flush=True)
        write_json(out / "atoms.json", res)


if __name__ == "__main__":
    main(sys.argv[1])
