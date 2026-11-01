"""Figure for the doublet study: Strehl ratio versus field height (three architectures, 852 nm)
and versus wavelength detuning of one phase map. Usage: python scripts/plot_metalens_doublet.py"""

import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from metatrap.runio import REPO_ROOT

C = ["#2a78d6", "#eb6834", "#1baf7a"]  # categorical slots 1-3 (validated)
res = json.load(open(REPO_ROOT / "results/metalens_doublet/doublet.json"))
fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
lam0 = "0.8523"
rows = [("Hyperbolic singlet", res["baselines"]["hyperbolic_singlet"]["free_focus"]),
        ("Aplanat on 2 mm silica", res["baselines"]["silica_spaced_aplanat"]["free_focus"]),
        ("Air-spaced aplanat (16 mm)", res["headroom"][lam0]["free_focus"])]
for (lbl, r), c in zip(rows, C):
    h = [x["height_um"] for x in r]
    s = [x["strehl"] for x in r]
    ax[0].plot(h, s, "-o", color=c, lw=2, ms=6, label=lbl)
ax[0].axvline(399, color="0.6", lw=1, ls=":")
ax[0].text(410, 0.30, "v2 array corner", color="0.35", fontsize=8)
ax[0].set_xlabel("image height in the atom plane (um)")
ax[0].set_ylabel("Strehl ratio (ray trace, best focus)")
ax[0].set_ylim(0, 1.03)
ax[0].legend(fontsize=8, frameon=False, loc="center right")
ch = [r for r in res["chromatic"]["detuning"] if r["detuning_nm"] > 0]
dl = np.array([r["detuning_nm"] for r in ch])
ax[1].semilogx(dl, [r["strehl_fixed_plane"] for r in ch], "-o", color=C[0], lw=2, ms=6, label="fixed input, nominal plane")
ax[1].semilogx(dl, [r["strehl_vergence_corrected"][0] for r in ch], "-o", color=C[1], lw=2, ms=6, label="input vergence re-focused (on axis)")
ax[1].axhline(0.8, color="0.6", lw=1, ls=":")
for r in ch:
    if r["detuning_nm"] in (42.3, 212.1):
        ax[1].annotate({42.3: "894.6 nm", 212.1: "1064.4 nm"}[r["detuning_nm"]], (r["detuning_nm"], 0.13), fontsize=8, color="0.35", ha="center")
ax[1].set_xlabel("detuning from the 852.3 nm design (nm)")
ax[1].set_ylabel("on-axis Strehl ratio")
ax[1].set_ylim(0, 1.03)
ax[1].legend(fontsize=8, frameon=False, loc="upper right")
for a in ax:
    a.grid(color="0.9", lw=0.6)
    for s in ("top", "right"):
        a.spines[s].set_visible(False)
fig.tight_layout()
fig.savefig(REPO_ROOT / "results/metalens_doublet/doublet_summary.png", dpi=130)
