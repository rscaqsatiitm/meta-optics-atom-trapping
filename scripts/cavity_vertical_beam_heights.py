"""Beam extent on every mirror of the vertical-lattice ring, in both directions, for the z-telescope variants.
Usage:
    python scripts/cavity_vertical_beam_heights.py configs/cavity_vertical_beam_heights.yaml

For the ideal design mode (scripts/cavity_vertical.py:build) the field arriving at each stop of half
trip 1 (P1, P2, G, P2', P1') is propagated from the crossing; we report the rms half-width and the
half-width that encloses a fraction (1 - clip) of the power, in the horizontal direction u (strip
width) and in the ring plane z (strip height; projected on the mirror face as / cos AOI is NOT applied
here). These set the mirror strip sizes of the fold layout (scripts/cavity_vertical_layout.py).
Writes heights.json.
"""

import sys

import numpy as np

from metatrap.runio import REPO_ROOT, load_config, run_dir, write_json, write_params

sys.path.insert(0, str(REPO_ROOT / "scripts"))
import cavity_vertical as V  # noqa: E402


def extents(s, clip):
    r, g = s["ring"], s["g"]
    F = s["F"].astype(complex)
    t = sorted([zs for zs, _ in r._stops_half[1]] + [r.p_um])
    out = {}
    for z in t:
        i = np.abs(r._segment(F, 0.0, z - 1e-3, 1)) ** 2
        c = np.cumsum(i) / i.sum()
        lo, hi = g.x[np.searchsorted(c, clip / 2)], g.x[np.searchsorted(c, 1 - clip / 2)]
        out[f"{z / 1e3:.0f}mm"] = {"rms_half_width_um": float(np.sqrt(np.sum(g.x**2 * i) / i.sum())),
                                    "enclosed_half_width_um": float(max(abs(lo), abs(hi)))}
    return out


def main(cfg_path):
    cfg = load_config(cfg_path)
    out = run_dir(cfg["run_id"])
    write_params(out, cfg)
    vcfg = load_config(REPO_ROOT / cfg["vertical_config"])
    rep = {"clip": cfg["clip"], "x": extents(V.build(vcfg, vcfg["x"]), cfg["clip"]), "z": {}}
    for var in cfg["z_variants"]:
        dz = dict(vcfg["z"], **{k: v for k, v in var.items() if k != "name"})
        rep["z"][var["name"]] = extents(V.build(vcfg, dz), cfg["clip"])
        print(var["name"], {k: round(v["enclosed_half_width_um"]) for k, v in rep["z"][var["name"]].items()}, flush=True)
    print("x", {k: round(v["enclosed_half_width_um"]) for k, v in rep["x"].items()})
    write_json(out / "heights.json", rep)


if __name__ == "__main__":
    main(sys.argv[1])
