"""Validate folded (block-wise) propagation against the full-grid angular spectrum.

Usage: python scripts/validate_folded.py configs/folded_validation_5mm.yaml
Designs the exact-kernel hologram for the 5 mm bench, propagates the same rendered
metasurface (i) on the full 10,240^2 grid and (ii) folded onto n_fold^2, and compares
the focal field and the spot metrics in the array window.
"""

import gc
import sys
import time

import numpy as np

from metatrap.hologram import evaluate as E
from metatrap.hologram import exact as X
from metatrap.hologram import targets
from metatrap.hologram.render import fourier_upsample, metasurface_blocks, metasurface_field
from metatrap.hologram.wgs import PupilGrid
from metatrap.propagation.angular_spectrum import propagate
from metatrap.propagation.folded import fold_blocks, propagate_folded
from metatrap.propagation.layers import LayerStack
from metatrap.runio import load_config, run_dir, write_json, write_params


def main(cfg_path):
    cfg = load_config(cfg_path)
    out = run_dir(cfg["run_id"])
    write_params(out, cfg)
    lam, na, pitch, dx = cfg["wavelength_um"], cfg["na"], cfg["site_pitch_um"], cfg["pixel_um"]
    stack = LayerStack.from_list(cfg["stack"])
    ill = cfg["illumination"]
    sites = targets.make(**cfg["sites"])
    radius = stack.aperture_for_na(na, lam)
    nf = cfg["n_fold"]
    n = int(np.ceil(2 * radius / dx / nf + 0.02) * nf)
    g = PupilGrid.for_sites(lam, pitch, 4, na, 2 * np.abs(sites).max() + pitch, guard=1.6)
    ps = X.PupilSampling.build(g, stack, na, ill)
    tab = X.OPLTable(stack, lam, radius + 2 * np.abs(sites).max() + 20)
    K = X.kernels(ps, sites, tab, lam)
    d = X.wgs_exact(ps, K, n_iter=60, n_fix=30)
    del K
    gc.collect()
    holo = X.render_grid(ps, stack, sites, d.coeffs, lam, na, kind="exact", table=tab, upsample=8)
    half = np.abs(sites).max() + 4 * pitch
    w = int(np.ceil(2 * half / dx / 2) * 2)

    t0 = time.time()
    f = metasurface_field(n, dx, stack, lam, na, ps, d.psi, ill, holo=holo)
    p_in = float(np.sum(np.abs(f) ** 2) * dx * dx)
    full = propagate(f, dx, lam, stack)
    del f
    c0 = n // 2 - w // 2
    win_full = np.array(full[c0 : c0 + w, c0 : c0 + w])
    del full
    gc.collect()
    t_full = time.time() - t0

    t0 = time.time()
    F, _ = fold_blocks(metasurface_blocks(n, dx, stack, lam, na, ps, d.psi, ill, holo=holo), n, nf)
    foc = propagate_folded(F, dx, lam, stack)
    c1 = nf // 2 - w // 2
    win_fold = foc[c1 : c1 + w, c1 : c1 + w]
    t_fold = time.time() - t0

    rel = float(np.linalg.norm(win_fold - win_full) / np.linalg.norm(win_full))
    res = {"n_full": n, "n_fold": nf, "window_um": nf * dx, "rel_l2_error_window": rel, "t_full_s": t_full, "t_fold_s": t_fold}
    ref_peak = None
    for name, win in (("full", win_full), ("folded", win_fold)):
        I = np.abs(fourier_upsample(win, 4)) ** 2
        m = E.spot_metrics(I, dx / 4, sites, pitch, lam, na, p_in)
        pk = m.pop("per_site")["peak"]
        res[name] = {k: m[k] for k in ("depth_nonuniformity_rms", "depth_min_over_mean", "array_efficiency", "core_efficiency")}
        if ref_peak is None:
            ref_peak = pk
        else:
            res["peak_rel_diff_max"] = float(np.max(np.abs(pk / ref_peak - 1)))
    print(res, flush=True)
    write_json(out / "metrics.json", res)


if __name__ == "__main__":
    main(sys.argv[1])
