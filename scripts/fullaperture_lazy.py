"""Full-aperture design and simulation with matrix-free exact kernels (very large arrays).

Usage: python scripts/fullaperture_lazy.py CONFIG [--gpu] [--tiled]

Identical to scripts/fullaperture.py (design by D4-symmetric weighted GS with exact kernels,
feedback rounds, folded full-aperture propagation, library pillars), except that the kernel
matrix and the render-grid site sum are evaluated on the fly (hologram/lazy.py), so memory no
longer scales with (pupil samples x sites).

--gpu    kernel products on a CUDA GPU (hologram/lazy_gpu.py, numba.cuda, single precision).
--tiled  no fine render grid and no full upsampled focal window in memory
         (hologram/tiled.py); required for 10^6 sites on a 16 GB machine. The fine
         hologram values are evaluated by a windowed-sinc interpolant of the coarse pupil
         grid (closer to the exact field than the FFT-upsampled grid, tests/test_tiled.py),
         and the focal field is upsampled and measured tile by tile. With fold_dtype:
         complex64 and n_fold equal to the metasurface size (36,864 at 18.4 mm), the field
         is not folded at all (unfolded grid, ~11 GB, propagated in place).
"""

import importlib.util
import sys
from pathlib import Path

import numpy as np

from metatrap.hologram import lazy as L
from metatrap.hologram import symmetric as S

_spec = importlib.util.spec_from_file_location("fa", Path(__file__).with_name("fullaperture.py"))
FA = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(FA)


def _kernel_class(gpu):
    if gpu:
        from metatrap.hologram.lazy_gpu import LazyKGGPU

        return LazyKGGPU
    return L.LazyKG


def install(gpu=False, tiled=False):
    KC = _kernel_class(gpu)

    def _kernels(wp, fs, tab, lam, block=16):
        rho_max = float(np.max(np.hypot(wp.r[:, 0], wp.r[:, 1]))) + float(np.max(np.abs(fs))) * 1.5 + 50.0
        return KC(wp, fs, tab, lam, rho_max)

    S.kernels_G = _kernels
    FA.S = S
    if not tiled:
        S.render_grid_symmetric = L.render_grid_symmetric_lazy
        return
    from metatrap.hologram.tiled import focal_metrics_tiled, tiled_render_grid

    def render_fold(cfg, stack, lam, na, ill, g, fs, fstab, tab, d, lookup, errors, nf, n, dx):
        holo = tiled_render_grid(g, stack, na, ill, fs, fstab, d.coeffs, tab, lam, upsample=cfg.get("render_upsample", 8),
                                 kernel=lambda *a: KC(*a))
        blocks = FA.metasurface_blocks(n, dx, stack, lam, na, None, None, ill, holo=holo, lookup=lookup, pillar_errors=errors,
                                       block=cfg.get("block_rows", 128), with_incident=True)
        F, p_pix, p_inc_pix = FA.fold_blocks(blocks, n, nf, dtype=np.dtype(cfg.get("fold_dtype", "complex128")))
        del holo
        FA.gc.collect()
        return F, p_pix, p_inc_pix

    def focal_metrics(cfg, F, dx, lam, na, sites, p_in, H=None, stack_eval=None, site_shift=(0.0, 0.0)):
        return focal_metrics_tiled(cfg, F, dx, lam, na, sites, p_in, H=H, stack_eval=stack_eval, site_shift=site_shift,
                                   tile_px=cfg.get("metrics_tile_px", 1024))

    FA.render_fold = render_fold
    FA.focal_metrics = focal_metrics


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    install(gpu="--gpu" in sys.argv, tiled="--tiled" in sys.argv)
    FA.main(args[0])
