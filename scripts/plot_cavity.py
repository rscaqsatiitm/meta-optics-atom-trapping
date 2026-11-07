"""Figures for docs/cavity_concept.md from the S1 and S3 results.
Usage:
    python scripts/plot_cavity.py
Writes results/cavity_budget/power_scaling.png and results/cavity_foxli/tolerances.png.
"""

import json
from collections import defaultdict

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from metatrap.runio import REPO_ROOT  # noqa: E402

SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6"]
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"


def style(ax):
    ax.set_facecolor(SURFACE)
    ax.grid(True, color=GRID, lw=0.6)
    for s in ax.spines.values():
        s.set_color(INK2)
        s.set_linewidth(0.6)
    ax.tick_params(colors=INK2, labelsize=8)
    ax.xaxis.label.set_color(INK)
    ax.yaxis.label.set_color(INK)


def power_scaling():
    b = json.loads((REPO_ROOT / "results/cavity_budget/budget.json").read_text())
    rows = sorted([c for c in b["cases"] if c["depth_mK"] == 0.2 and c["match"] == "equal_depth"], key=lambda c: c["n_sites"])
    n = np.array([c["n_sites"] for c in rows])
    series = [
        ("Meta-hologram tweezers", [c["tweezer_hologram_W"] for c in rows]),
        ("Lenslet-array tweezers", [c["tweezer_lenslet_W"] for c in rows]),
        ("Lattice, free space", [c["free_space_lattice_W"] for c in rows]),
        ("Cavity lattice, through cell (AR 0.1 %)", [c["cavity_input_W"]["cell_AR"]["P_in_W"] for c in rows]),
        ("Cavity lattice, Brewster windows", [c["cavity_input_W"]["brewster"]["P_in_W"] for c in rows]),
        ("Cavity lattice, in vacuum", [c["cavity_input_W"]["vacuum"]["P_in_W"] for c in rows]),
    ]
    fig, ax = plt.subplots(figsize=(6.4, 4.2), facecolor=SURFACE)
    style(ax)
    for (lab, y), col in zip(series, SERIES):
        ax.loglog(n, y, "-", color=col, lw=2, marker="o", ms=4, label=lab)
    ax.set_xlabel("Number of sites")
    ax.set_ylabel("Laser power at the optics input (W)")
    ax.legend(fontsize=7, frameon=False, labelcolor=INK)
    ax.set_title("0.2 mK barrier, 4 um pitch, 1075 nm, Cs", fontsize=9, color=INK)
    fig.tight_layout()
    fig.savefig(REPO_ROOT / "results/cavity_budget/power_scaling.png", dpi=150, facecolor=SURFACE)


def tolerances():
    r = json.loads((REPO_ROOT / "results/cavity_foxli/foxli.json").read_text())
    variants = [("plane150", "Plane folds, p 150 mm"), ("plane300", "Plane folds, p 300 mm"),
                ("telK_M0.2", "Telescope 0.2x (Keplerian)"), ("telK_M0.1", "Telescope 0.1x (Keplerian)"),
                ("telG_M0.2", "Telescope 0.2x (Galilean)")]
    fig, axs = plt.subplots(1, 2, figsize=(8.6, 3.6), facecolor=SURFACE)
    for ax in axs:
        style(ax)
    for (v, lab), col in zip(variants, SERIES):
        acc = defaultdict(list)
        for k, e in r.items():
            if k.startswith(f"figure_{v}_rms"):
                acc[e["rms_waves"]].append(e["array"]["depth_rms"])
        if acc:
            xs = sorted(acc)
            axs[0].loglog(xs, [np.median(acc[x]) for x in xs], "-o", color=col, lw=2, ms=4, label=lab)
        tl = sorted((e["tilt_urad"], e["array"]["depth_rms"]) for k, e in r.items() if k.startswith(f"tilt_{v}_"))
        if tl:
            axs[1].loglog(*zip(*tl), "-o", color=col, lw=2, ms=4, label=lab)
    axs[0].set_xlabel("Graded-mirror figure error, rms (waves of surface)")
    axs[0].set_ylabel("Depth non-uniformity over 10^6 sites (rms)")
    axs[1].set_xlabel("Graded-mirror tilt (urad)")
    axs[1].set_ylabel("Depth non-uniformity over 10^6 sites (rms)")
    axs[0].legend(fontsize=7, frameon=False, labelcolor=INK)
    fig.tight_layout()
    fig.savefig(REPO_ROOT / "results/cavity_foxli/tolerances.png", dpi=150, facecolor=SURFACE)


if __name__ == "__main__":
    power_scaling()
    tolerances()
