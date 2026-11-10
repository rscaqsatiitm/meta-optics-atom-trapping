"""Render accuracy check: windowed-sinc (tiled) and FFT-upsampled hologram renders against an exact
per-pixel kernel sum, on a small design.

Usage: python scripts/tiled_render_check.py configs/tiled_render_check.yaml

The metasurface phase at each pixel is the argument of the site sum q(s) = sum_j K_j(s) c_j at
that pixel's direction cosine s. Three ways to obtain it are compared, for the same stored design
and the same propagation and metrics (unfolded grid, so nothing is aliased):
(i)   exact: q evaluated at every pixel with the exact kernels (hologram/lazy.LazyKG);
(ii)  tiled: windowed-sinc interpolation of q from the coarse pupil grid (hologram/tiled.py);
(iii) fft: the FFT-upsampled fine render grid of the untiled lazy pipeline (lazy.render_grid_symmetric_lazy).
Reported: per-site depth rms difference of (ii) and (iii) against (i), and each one's own
non-uniformity. Output: results/<run_id>/metrics.json. Supports the statement in
docs/scaling_100k.md on the render.
"""

import importlib.util
import sys
import types
from pathlib import Path

import numpy as np

from metatrap.hologram import lazy as LZK
from metatrap.hologram.lazy import LazyKG
from metatrap.hologram.tiled import focal_metrics_tiled, tiled_render_grid
from metatrap.propagation.folded import fold_blocks
from metatrap.runio import REPO_ROOT, load_config, run_dir, write_json, write_params

_spec = importlib.util.spec_from_file_location("fws", Path(__file__).with_name("fold_window_study.py"))
FWS = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(FWS)
FA = FWS.FA


def exact_grid(G, coeffs, upsample):
    tg = tiled_render_grid(G["g"], G["stack"], G["na"], G["ill"], G["fs"], G["fstab"], coeffs, G["tab"], G["lam"], upsample=upsample)
    x = coeffs / G["fstab"]
    stack, lam = G["stack"], G["lam"]

    def sample(fi, fj):
        sx, sy = (np.ravel(fi) - tg.c0) * tg.ds, (np.ravel(fj) - tg.c0) * tg.ds
        s = np.hypot(sx, sy)
        s2 = np.minimum(s, 0.95)
        rho = stack.ray_height(s2, lam)
        ux = np.where(s > 0, sx / np.maximum(s, 1e-30), 0.0)
        uy = np.where(s > 0, sy / np.maximum(s, 1e-30), 0.0)
        wp = types.SimpleNamespace(r=np.column_stack([-ux * rho, -uy * rho]), opl0=stack.opl(s2, lam))
        q = LazyKG(wp, G["fs"], G["tab"], lam, float(rho.max()) + 1.5 * float(np.abs(G["fs"]).max()) + 50) @ x
        z = np.where(s**2 <= tg.lim, q / np.maximum(np.abs(q), 1e-30), 0)
        return z.astype(np.complex64).reshape(np.shape(fi))

    tg.sample = sample
    return tg


def main(cfg_path):
    cfg = load_config(cfg_path)
    out = run_dir(cfg["run_id"])
    write_params(out, cfg)
    G = FWS.geometry(cfg)
    coeffs = np.load(REPO_ROOT / cfg["state"])["coeffs"]
    up = cfg.get("render_upsample", 5)
    renders = {
        "exact": lambda: exact_grid(G, coeffs, up),
        "tiled": lambda: tiled_render_grid(G["g"], G["stack"], G["na"], G["ill"], G["fs"], G["fstab"], coeffs, G["tab"], G["lam"], upsample=up),
        "fft": lambda: LZK.render_grid_symmetric_lazy(G["g"], G["stack"], G["na"], G["ill"], G["fs"], G["fstab"], coeffs, G["tab"], G["lam"], upsample=up),
    }
    res, peaks = {"n_metasurface": G["n"], "n_fold": G["n"], "state": cfg["state"]}, {}
    for name, make in renders.items():
        holo = make()
        blocks = FA.metasurface_blocks(G["n"], G["dx"], G["stack"], G["lam"], G["na"], None, None, G["ill"], holo=holo,
                                       block=cfg.get("block_rows", 64), with_incident=True)
        F, p_pix, p_inc = fold_blocks(blocks, G["n"], G["n"], dtype=np.complex64)
        m, per, _ = focal_metrics_tiled(cfg, F, G["dx"], G["lam"], G["na"], G["sites"], p_inc * G["dx"] ** 2,
                                        stack_eval=G["stack_sites"], tile_px=cfg.get("metrics_tile_px", 256))
        peaks[name] = per["peak"]
        res[name] = {k: m[k] for k in ("depth_nonuniformity_rms", "array_efficiency", "core_efficiency")}
        print(name, res[name], flush=True)
    for name in ("tiled", "fft"):
        d = peaks[name] / peaks[name].mean() - peaks["exact"] / peaks["exact"].mean()
        res[f"{name}_vs_exact_depth_rms"] = float(np.std(d))
        res[f"{name}_vs_exact_depth_max"] = float(np.max(np.abs(d)))
    print({k: v for k, v in res.items() if "vs_exact" in k}, flush=True)
    write_json(out / "metrics.json", res)


if __name__ == "__main__":
    main(sys.argv[1])
