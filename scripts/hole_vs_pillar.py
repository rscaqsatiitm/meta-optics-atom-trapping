"""Hole (inverse, "Mie-void") meta-atoms in a Si-rich SiN film versus pillars, at 1075 nm.

Usage: python scripts/hole_vs_pillar.py configs/hole_vs_pillar.yaml

Transmissive, polarisation-independent propagation-phase meta-atoms for the B1 hologram:
(i) circular air holes etched through a film of thickness h on fused silica (film continuous
around the holes), and (ii) square pillars of the same material and height (the present
design family), both on a square lattice of pitch p, light incident from the substrate, exit
into air. Both are computed with the same solver settings (fmmax, pol_fourier formulation,
which has no spurious cross-polarisation for C4v cells; see tests/test_rcwa.py), so the
comparison is like for like. Output: results/hole_vs_pillar/sweep.json (complex t per width,
per (family, pitch, height)); ranking by scripts/hole_vs_pillar_rank.py.
"""

import json
import sys
import time

import numpy as np

from metatrap import materials as M
from metatrap.rcwa.cell import PillarCell
from metatrap.rcwa.solver import solve_eps
from metatrap.runio import REPO_ROOT, load_config, run_dir, write_json, write_params


def main(cfg_path):
    cfg = load_config(cfg_path)
    out = run_dir(cfg["run_id"])
    write_params(out, cfg)
    path = out / "sweep.json"
    res = json.load(open(path)) if path.exists() else {}
    lam = cfg["wavelength_um"]
    n_f = complex(M.get(cfg["material"]).n(lam))
    n_s = float(M.get(cfg["substrate"]).n(lam))
    for fam in cfg["families"]:
        for p in cfg["pitches_um"]:
            w = np.linspace(cfg["w_min_um"], p - cfg["w_min_um"], cfg["n_widths"])
            for h in cfg["heights_um"]:
                key = f"{fam}:{p:.2f}:{h:.2f}"
                if key in res:
                    continue
                t0 = time.time()
                if fam == "hole":
                    m = PillarCell(pitch=p, height=h, shape="circle").mask(w, cfg["n_grid"])
                    eps = n_f**2 + (1.0 - n_f**2) * m  # air inside the hole
                else:
                    m = PillarCell(pitch=p, height=h, shape="square").mask(w, cfg["n_grid"])
                    eps = 1.0 + (n_f**2 - 1.0) * m  # film material inside the pillar
                r = solve_eps(eps.astype(complex), p, h, lam, n_s**2, 1.0, n_terms=cfg["n_terms"], formulation=cfg["formulation"])
                res[key] = {"family": fam, "pitch_um": p, "height_um": h, "widths_um": w.tolist(),
                            "t_re": np.real(r.t_x).tolist(), "t_im": np.imag(r.t_x).tolist(), "T": np.real(r.T).tolist(),
                            "R": np.real(r.R).tolist(), "t_cross_max": float(np.max(r.t_cross)), "n_terms": int(r.n_terms),
                            "seconds": time.time() - t0}
                print(key, f"{time.time() - t0:.0f}s, phase range {np.ptp(np.unwrap(np.angle(r.t_x))):.2f} rad, max t_cross {np.max(r.t_cross):.1e}", flush=True)
                write_json(path, res)


if __name__ == "__main__":
    main(sys.argv[1])
