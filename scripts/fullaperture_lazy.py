"""Full-aperture design and simulation with matrix-free exact kernels (very large arrays).

Usage: python scripts/fullaperture_lazy.py configs/fullaperture_100k.yaml

Identical to scripts/fullaperture.py (design by D4-symmetric weighted GS with exact kernels,
feedback rounds, folded full-aperture propagation, library pillars), except that the kernel
matrix and the render-grid site sum are evaluated on the fly (hologram/lazy.py), so memory no
longer scales with (pupil samples x sites).
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


def _kernels(wp, fs, tab, lam, block=16):
    rho_max = float(np.max(np.hypot(wp.r[:, 0], wp.r[:, 1]))) + float(np.max(np.abs(fs))) * 1.5 + 50.0
    return L.LazyKG(wp, fs, tab, lam, rho_max)


S.kernels_G = _kernels
S.render_grid_symmetric = L.render_grid_symmetric_lazy
FA.S = S

if __name__ == "__main__":
    FA.main(sys.argv[1])
