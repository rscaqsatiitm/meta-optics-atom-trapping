"""Two-wavelength meta-atom screen (RCWA). Usage: python scripts/run_dual_library.py configs/dual_library.yaml

One checkpoint file per (family, material, pitch, height, wavelength) in output_dir; existing
files are skipped, so the run can be restarted. A convergence spot check (reduced against
full numerics) is written to results/<run_id>/convergence.json first.
"""

import json
import sys
import time

import numpy as np

from metatrap import materials
from metatrap.rcwa.molecules import Rules, candidates, masks, material_groups
from metatrap.rcwa.solver import solve_eps
from metatrap.runio import REPO_ROOT, load_config, run_dir, write_json, write_params


def solve_group(fam, params, mat, sub, p, h, lam, n_terms, n_grid, batch):
    eps_p = complex(materials.get(mat).eps(lam))
    ns = float(materials.get(sub).n(lam))
    tx, T, R, tc = [], [], [], []
    for i in range(0, len(params), batch):
        m = masks(fam, params[i : i + batch], p, n_grid)
        eps = 1 + (eps_p - 1) * m
        r = solve_eps(eps, p, h, lam, ns**2, 1.0, n_terms)
        tx.append(r.t_x)
        T.append(r.T)
        R.append(r.R)
        tc.append(r.t_cross)
    return np.concatenate(tx), np.concatenate(T), np.concatenate(R), np.concatenate(tc), ns


def main(cfg_path):
    cfg = load_config(cfg_path)
    out = run_dir(cfg["run_id"])
    write_params(out, cfg)
    od = REPO_ROOT / cfg["output_dir"]
    od.mkdir(parents=True, exist_ok=True)
    rules0 = cfg["rules"]
    cv = cfg.get("convergence")
    cpath = out / "convergence.json"
    groups = material_groups(cfg)
    if cv and not cpath.exists():
        rows = []
        mat, sub, p, h = groups[-1]
        for fam in cfg["families"]:
            lam = cfg["wavelengths_um"][0]
            c = candidates(fam, p, h, Rules(step_um=cfg["families"][fam]["step_um"], **rules0))
            pick = [c[i] for i in np.linspace(0, len(c) - 1, cv["n_atoms"]).astype(int)]
            a = solve_group(fam, pick, mat, sub, p, h, lam, cfg["n_terms"], cfg["n_grid"], cfg["batch"])[0]
            b = solve_group(fam, pick, mat, sub, p, h, lam, cv["n_terms"], cv["n_grid"], cfg["batch"])[0]
            dphi = np.angle(a * np.conj(b))
            rows.append({"family": fam, "material": mat, "params": pick, "dphase_rad": dphi.tolist(), "dabs_rel": (np.abs(a) / np.abs(b) - 1).tolist()})
            print("convergence", fam, mat, np.round(dphi, 3), np.round(np.abs(a) / np.abs(b) - 1, 4), flush=True)
        write_json(cpath, rows)
    for fam, fc in cfg["families"].items():
        for mat, sub, p, h in groups:
            params = candidates(fam, p, h, Rules(step_um=fc["step_um"], **rules0))
            for lam in cfg["wavelengths_um"]:
                f = od / f"{fam}_{mat}_p{p:.2f}_h{h:.1f}_l{lam:.4f}.npz"
                if f.exists():
                    continue
                t0 = time.time()
                tx, T, R, tc, ns = solve_group(fam, params, mat, sub, p, h, lam, cfg["n_terms"], cfg["n_grid"], cfg["batch"])
                np.savez(f, params=np.array(params), t=tx, T=T, R=R, t_cross=tc, n_in=ns, pitch=p, height=h, lam=lam, substrate=sub)
                print(f"{f.name}: {len(params)} atoms, {time.time()-t0:.0f} s", flush=True)
    print("dual library done", flush=True)


if __name__ == "__main__":
    main(sys.argv[1])
