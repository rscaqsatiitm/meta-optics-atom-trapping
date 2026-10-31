"""Convergence of the pillar response with Fourier order, formulation and grid.

Usage: python scripts/rcwa_convergence.py configs/rcwa_convergence.yaml
Writes results/<run_id>/convergence.json and convergence.png.
"""

import sys
import time

import numpy as np

from metatrap import materials
from metatrap.rcwa.cell import PillarCell
from metatrap.rcwa.solver import solve_pillars, solve_pillars_grcwa
from metatrap.runio import load_config, run_dir, write_json, write_params


def main(cfg_path):
    cfg = load_config(cfg_path)
    out = run_dir(cfg["run_id"])
    write_params(out, cfg)
    lam = cfg["wavelength_um"]
    n_sub = materials.get(cfg["substrate"]).n(lam)
    res = {"cases": []}
    for case in cfg["cases"]:
        n_pil = materials.get(case["material"]).n(lam)
        cell = PillarCell(case["pitch_um"], case["height_um"], cfg.get("shape", "square"))
        widths = np.array(case["widths_um"])
        w_ref = case["ref_width_um"]
        allw = np.concatenate([[w_ref], widths])
        entry = {"case": case, "fmmax": {}, "grcwa": {}}
        for form in cfg["formulations"]:
            rows = []
            for nt in cfg["n_terms"]:
                t0 = time.time()
                r = solve_pillars(cell, allw, lam, n_pil**2, n_sub**2, n_terms=nt, n_grid=cfg["n_grid"], formulation=form)
                rel = r.t_x[1:] / r.t_x[0]
                rows.append({"n_terms": int(r.n_terms), "phase_rel": np.angle(rel), "T0": r.T0[1:], "T": r.T[1:], "sec": time.time() - t0})
                print(form, nt, r.n_terms, np.round(np.angle(rel), 4), np.round(r.T0[1:], 4), f"{time.time()-t0:.1f}s", flush=True)
            entry["fmmax"][form] = rows
        rows = []
        for nt in cfg["n_terms_grcwa"]:
            t0 = time.time()
            vals = [solve_pillars_grcwa(cell, w, lam, n_pil**2, n_sub**2, n_terms=nt, n_grid=cfg["n_grid"]) for w in allw]
            t = np.array([v[0] for v in vals])
            rows.append({"n_terms": nt, "phase_rel": np.angle(t[1:] / t[0]), "T": np.array([v[1] for v in vals])[1:], "sec": time.time() - t0})
            print("grcwa", nt, np.round(np.angle(t[1:] / t[0]), 4), f"{time.time()-t0:.1f}s", flush=True)
        entry["grcwa"] = rows
        res["cases"].append(entry)
    write_json(out / "convergence.json", res)
    plot(res, out)


def plot(res, out):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n = len(res["cases"])
    fig, axes = plt.subplots(1, n, figsize=(5 * n, 4), squeeze=False)
    for ax, entry in zip(axes[0], res["cases"]):
        for form, rows in entry["fmmax"].items():
            nts = [r["n_terms"] for r in rows]
            ph = np.array([np.unwrap(r["phase_rel"]) for r in rows])
            for j in range(ph.shape[1]):
                ax.plot(nts, ph[:, j] - ph[-1, j], "o-", label=f"fmmax {form} w={entry['case']['widths_um'][j]}")
        rows = entry["grcwa"]
        nts = [r["n_terms"] for r in rows]
        ref = np.array([np.unwrap(r) for r in [entry["fmmax"]["jones_direct_fourier"][-1]["phase_rel"]]])[0]
        ph = np.array([np.unwrap(r["phase_rel"]) for r in rows])
        for j in range(ph.shape[1]):
            ax.plot(nts, np.angle(np.exp(1j * (ph[:, j] - ref[j]))), "s--", label=f"grcwa w={entry['case']['widths_um'][j]}")
        ax.axhline(0.02, color="k", lw=0.5)
        ax.axhline(-0.02, color="k", lw=0.5)
        c = entry["case"]
        ax.set_title(f"{c['material']} p={c['pitch_um']} h={c['height_um']}")
        ax.set_xlabel("number of Fourier terms")
        ax.set_ylabel("relative phase minus best fmmax JDF value (rad)")
        ax.set_ylim(-0.3, 0.3)
        ax.legend(fontsize=6)
    fig.tight_layout()
    fig.savefig(out / "convergence.png", dpi=130)


if __name__ == "__main__":
    main(sys.argv[1])
