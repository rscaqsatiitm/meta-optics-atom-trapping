"""Rank hole and pillar meta-atom families from results/hole_vs_pillar/sweep.json.

Usage: python scripts/hole_vs_pillar_rank.py configs/hole_vs_pillar.yaml

Per (family, pitch, height) and fabrication-limit set (min feature, min gap = wall between holes,
max aspect = depth / smallest hole diameter or height / smallest pillar width), on the feasible
widths: eta_16 (16-level blazed grating in the local-periodic approximation, Fresnel loss
included, n_out / n_in), its average over a +/-10 nm width bias (eta_robust, the library
selection metric), phase coverage, minimum zeroth-order transmittance, and the material fill
fraction averaged over the 16 selected meta-atoms (the film absorption scales with it).
Output: results/hole_vs_pillar/ranking.json.
"""

import json
import sys

import numpy as np

from metatrap import materials as M
from metatrap.rcwa import lookup as L
from metatrap.runio import REPO_ROOT, load_config, write_json


def fill_fraction(fam, w, p):
    return 1 - np.pi * w**2 / (4 * p**2) if fam == "hole" else (w / p) ** 2


def main(cfg_path):
    cfg = load_config(cfg_path)
    res = json.load(open(REPO_ROOT / "results" / cfg["run_id"] / "sweep.json"))
    n_s = float(M.get(cfg["substrate"]).n(cfg["wavelength_um"]))
    out = {}
    for lim_name, lim_cfg in cfg["fab_limits"].items():
        lim = L.FabLimits(**lim_cfg)
        rows = []
        for key, r in res.items():
            w = np.array(r["widths_um"])
            t = np.array(r["t_re"]) + 1j * np.array(r["t_im"])
            T = np.array(r["T"])
            ok = L.feasible(w, r["pitch_um"], r["height_um"], lim)
            if ok.sum() < 4:
                continue
            wd, td = L.densify(w[ok], t[ok])
            eta_r, eta, _ = L.eta_robust(wd, td, 1.0 / n_s, 16)
            off = L.eta_n(wd, td, 1.0 / n_s, 16)[1]
            idx = L.best_widths(wd, td, 2 * np.pi * np.arange(16) / 16, off)
            rows.append({"family": r["family"], "pitch_um": r["pitch_um"], "height_um": r["height_um"],
                         "eta16": eta, "eta_robust": eta_r, "phase_coverage_rad": L.phase_coverage(wd, td),
                         "T_min": float(T[ok].min()), "w_min_um": float(w[ok].min()), "w_max_um": float(w[ok].max()),
                         "fill_mean_selected": float(np.mean(fill_fraction(r["family"], wd[idx], r["pitch_um"]))),
                         "t_cross_max": r["t_cross_max"]})
        out[lim_name] = {}
        for fam in ("hole", "pillar"):
            fr = sorted([x for x in rows if x["family"] == fam], key=lambda x: -x["eta_robust"])
            if not fr:
                continue
            best = fr[0]["eta_robust"]
            near = [x for x in fr if x["eta_robust"] >= best - 0.003]
            pick = min(near, key=lambda x: x["height_um"])  # shortest within 0.3 % of the best, as for the pillar library
            out[lim_name][fam] = {"best": fr[0], "selected_shortest_within_0.3pct": pick, "top10": fr[:10], "n_feasible": len(fr)}
            print(lim_name, fam, "best", {k: (round(v, 3) if isinstance(v, float) else v) for k, v in fr[0].items()})
            print(lim_name, fam, "selected", {k: (round(v, 3) if isinstance(v, float) else v) for k, v in pick.items()}, flush=True)
    write_json(REPO_ROOT / "results" / cfg["run_id"] / "ranking.json", out)


if __name__ == "__main__":
    main(sys.argv[1])
