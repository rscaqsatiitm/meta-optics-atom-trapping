"""Grid convergence of the vertical-mode (z) spectra of the vertical-lattice ring, with and without a
vertical telescope.
Usage:
    python scripts/cavity_vertical_spectra_check.py configs/cavity_vertical_spectra_check.yaml

scripts/cavity_vertical.py and scripts/cavity_vertical_mirrors.py computed these spectra on different,
coarse grids (2048 or 4096 points over 40 mm) and quoted different numbers for the same mode. Here the
ideal design is rebuilt on a series of grids and the k largest round-trip eigenvalues are found with
ARPACK (src/metatrap/cavity/foxli.py: eigenmodes; the same method as sparse_spectrum, now used by both
scripts); we report the fundamental's loss and the competing
modes (lowest loss, smallest phase separation) per grid. Writes spectra.json.
"""

import sys

import numpy as np

from metatrap.cavity import foxli as FL
from metatrap.runio import REPO_ROOT, load_config, run_dir, write_json, write_params

sys.path.insert(0, str(REPO_ROOT / "scripts"))
import cavity_vertical as V  # noqa: E402


def spectrum(s, k):
    vals, vecs = FL.eigenmodes(s["ring"], k=k, v0=s["F"].astype(complex))
    F = s["F"].astype(complex)
    ov = np.abs(vecs.conj().T @ F) ** 2 / (np.vdot(F, F).real * np.sum(np.abs(vecs) ** 2, axis=0))
    i0 = int(np.argmax(ov))
    d = np.angle(vals / vals[i0])
    loss = 1 - np.abs(vals) ** 2
    others = [j for j in range(len(vals)) if j != i0]
    lowest = min(others, key=lambda j: loss[j])
    nearest = min((j for j in others if loss[j] < 0.01), key=lambda j: abs(d[j]), default=None)
    out = {"fundamental_loss": float(loss[i0]), "fundamental_overlap": float(ov[i0]),
           "lowest_loss_mode": {"dphase_rad": float(d[lowest]), "loss": float(loss[lowest])}}
    if nearest is not None:
        out["nearest_mode_below_1pct_loss"] = {"dphase_rad": float(d[nearest]), "loss": float(loss[nearest])}
    return out


def main(cfg_path):
    cfg = load_config(cfg_path)
    out = run_dir(cfg["run_id"])
    write_params(out, cfg)
    vcfg = load_config(REPO_ROOT / cfg["vertical_config"])
    rep = {}
    for var in cfg["variants"]:
        rep[var["name"]] = {}
        for n, width in cfg["grids"]:
            dz = dict(vcfg[var["direction"]], **{k: v for k, v in var.items() if k not in ("name", "direction")}, grid_n=n, grid_width_um=width)
            s = V.build(vcfg, dz)
            res = spectrum(s, cfg["k"])
            res["dx_um"] = width / n
            rep[var["name"]][f"{n}_{width / 1e3:.0f}mm"] = res
            print(var["name"], n, width, {k: (round(v, 9) if isinstance(v, float) else v) for k, v in res.items()}, flush=True)
            write_json(out / "spectra.json", rep)


if __name__ == "__main__":
    main(sys.argv[1])
