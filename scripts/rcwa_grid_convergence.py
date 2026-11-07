"""Permittivity-grid convergence of the pillar RCWA (fmmax, jones_direct_fourier, 193 terms).

Usage: python scripts/rcwa_grid_convergence.py

For the two selected 1075 nm families (LPCVD p 0.60 / h 1.6 um; Si-rich p 0.50 / h 1.5 um) and
widths spanning their fabricable range down to a 60 nm gap, computes t at n_grid = 256, 512
and 1024 and reports the phase and amplitude change relative to 1024. (campaign1_report
Section 3 quoted 0.01 to 0.03 rad between 256 and 1024 points from a run that was not stored;
this stores the check.) Output: results/rcwa_grid_convergence.json.
"""

import numpy as np

from metatrap import materials
from metatrap.rcwa.cell import PillarCell
from metatrap.rcwa.solver import solve_pillars
from metatrap.runio import REPO_ROOT, write_json

LAM = 1.075
FAMILIES = {"sin_lpcvd": (0.60, 1.6), "sin_si_rich": (0.50, 1.5)}


def main():
    n_in = float(materials.get("fused_silica").n(LAM))
    out = {}
    for mat, (p, h) in FAMILIES.items():
        n_p = complex(materials.get(mat).n(LAM))
        widths = np.linspace(0.15, p - 0.06, 8)
        cell = PillarCell(pitch=p, height=h)
        t = {g: solve_pillars(cell, widths, LAM, n_p**2, n_in**2, n_terms=197, n_grid=g).t_x for g in (256, 512, 1024)}
        dphi = {g: np.abs(np.angle(t[g] / t[1024])) for g in (256, 512)}
        out[mat] = {"widths_um": widths.tolist(), "pitch_um": p, "height_um": h,
                    **{f"max_dphase_rad_{g}_vs_1024": float(dphi[g].max()) for g in (256, 512)},
                    **{f"dphase_rad_{g}_vs_1024": dphi[g].tolist() for g in (256, 512)},
                    "max_damp_256_vs_1024": float(np.max(np.abs(np.abs(t[256]) - np.abs(t[1024]))))}
        print(mat, {k: round(v, 4) for k, v in out[mat].items() if k.startswith("max")}, flush=True)
    write_json(REPO_ROOT / "results" / "rcwa_grid_convergence.json", out)


if __name__ == "__main__":
    main()
