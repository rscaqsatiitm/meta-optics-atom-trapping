"""Site statistics and feedback-convergence diagnostics of full-aperture runs.

Usage: python scripts/fullaperture_diagnostics.py RUN_ID [RUN_ID ...]

From the per-site files of results/<run_id>/ (per_site_*.npz, git-ignored, written by
scripts/fullaperture.py) computes and stores in results/<run_id>/diagnostics.json:
(i) for each final per-site file: rms depth non-uniformity, its split into orbit-to-orbit and
within-orbit parts (D4 orbits), the radially binned mean depth (100 um bins), the fraction of
the variance explained by a smooth radial fit (quartic in r), and the correlation between
nearest-neighbour sites; (ii) for successive feedback rounds (per_site_round<k>.npz): the
regression slope and correlation of the per-orbit depth errors of round k+1 on round k.
"""

import json
import re
import sys

import numpy as np

from metatrap.hologram import symmetric as S
from metatrap.runio import REPO_ROOT, write_json


def site_stats(sites, pk, pitch):
    pk = pk / pk.mean()
    _, _, inv = S.wedge_sites(sites)
    om = np.bincount(inv, weights=pk) / np.bincount(inv)
    r = np.hypot(sites[:, 0], sites[:, 1])
    bins = np.arange(0, r.max() + 100, 100)
    idx = np.digitize(r, bins) - 1
    radial = [float(pk[idx == k].mean()) for k in range(len(bins) - 1) if np.any(idx == k)]
    c = np.polyfit(r**2, pk, 2)
    fit = np.polyval(c, r**2)
    # nearest neighbours along x on the square grid
    # keys in half-pitch units (even arrays have sites at half-integer multiples of the pitch)
    key = {(round(2 * x / pitch), round(2 * y / pitch)): i for i, (x, y) in enumerate(sites)}
    a, b = [], []
    for (ix, iy), i in key.items():
        j = key.get((ix + 2, iy))
        if j is not None:
            a.append(pk[i]); b.append(pk[j])
    return {"rms": float(pk.std()), "orbit_rms": float(np.std(om[inv])), "within_orbit_rms": float(np.std(pk / om[inv])),
            "radial_bin_mean_min": float(min(radial)), "radial_bin_mean_max": float(max(radial)), "radial_bin_um": 100,
            "radial_fit_variance_fraction": float(np.var(fit) / np.var(pk)), "neighbour_correlation": float(np.corrcoef(a, b)[0, 1])}


def main(run_ids):
    for run in run_ids:
        d = REPO_ROOT / "results" / run
        pitch = json.load(open(d / "params.json"))["config"]["site_pitch_um"]
        out = {"final": {}, "feedback": []}
        for f in sorted(d.glob("per_site_*.npz")):
            if f.stem.startswith("per_site_round"):
                continue
            z = np.load(f)
            out["final"][f.stem.replace("per_site_", "")] = site_stats(z["sites"], z["peak"], pitch)
        rounds = sorted((int(re.findall(r"\d+", f.stem)[-1]), f) for f in d.glob("per_site_round*.npz"))
        for (k0, f0), (k1, f1) in zip(rounds, rounds[1:]):
            a, b = np.load(f0), np.load(f1)
            _, _, inv = S.wedge_sites(a["sites"])
            ea = np.bincount(inv, weights=a["peak"]) / np.bincount(inv)
            eb = np.bincount(inv, weights=b["peak"]) / np.bincount(inv)
            ea, eb = ea / ea.mean() - 1, eb / eb.mean() - 1
            out["feedback"].append({"from_round": k0, "to_round": k1, "slope": float(np.polyfit(ea, eb, 1)[0]),
                                    "correlation": float(np.corrcoef(ea, eb)[0, 1])})
        write_json(d / "diagnostics.json", out)
        print(run, json.dumps(out)[:600], flush=True)


if __name__ == "__main__":
    main(sys.argv[1:])
