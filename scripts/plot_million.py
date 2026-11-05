"""Figure for docs/million_sites.md: objective field radius versus focal length and NA, and the
imaging ensquared energy across the 10^6-site field. Usage: python scripts/plot_million.py"""

import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from metatrap.runio import REPO_ROOT

C = ["#2a78d6", "#eb6834", "#1baf7a"]
f = json.load(open(REPO_ROOT / "results/million/objective_field.json"))
e = json.load(open(REPO_ROOT / "results/million/imaging_ee.json"))
fig, ax = plt.subplots(1, 2, figsize=(10, 3.9))
for na, c in zip((0.5, 0.6, 0.7), C):
    rows = sorted((v["f_um"] / 1000, v["r_strehl08_um"] / 1000) for v in f.values() if v["na"] == na)
    ax[0].plot([r[0] for r in rows], [r[1] for r in rows], "-o", color=c, lw=2, ms=6, label=f"NA {na}")
for y, lbl in ((2.828, "10^6 sites, 4 um pitch"), (2.121, "10^6 sites, 3 um pitch")):
    ax[0].axhline(y, color="0.6", lw=1, ls=":")
    ax[0].text(17, y + 0.05, lbl, color="0.35", fontsize=8, ha="left")
ax[0].set_xlabel("objective focal length f (mm)")
ax[0].set_ylabel("field radius, Strehl >= 0.8 (mm)")
ax[0].legend(fontsize=8, frameon=False, loc="lower right")
for n, c, lbl in (("f16", C[0], "f = 16 mm"), ("f25", C[1], "f = 25 mm"), ("f42", C[2], "f = 42 mm")):
    pts = sorted((v["field_um"] / 1000, v["ensquared_energy"]) for k, v in e.items() if k.startswith(n + ":"))
    ax[1].plot([p[0] for p in pts], [p[1] for p in pts], "-o", color=c, lw=2, ms=6, label=lbl)
ax[1].axvline(2.828, color="0.6", lw=1, ls=":")
ax[1].text(2.78, 0.42, "4 um-pitch corner", color="0.35", fontsize=8, ha="right")
ax[1].set_xlabel("field radius (mm)")
ax[1].set_ylabel("ensquared energy in the 4 um box")
ax[1].set_ylim(0.3, 1.0)
ax[1].legend(fontsize=8, frameon=False, loc="lower left")
for a in ax:
    a.grid(color="0.9", lw=0.6)
    for sp in ("top", "right"):
        a.spines[sp].set_visible(False)
fig.tight_layout()
fig.savefig(REPO_ROOT / "results/million/million_summary.png", dpi=130)
