"""LazyKGGPU against the dense symmetric kernel (runs on the CUDA simulator without a GPU:
NUMBA_ENABLE_CUDASIM=1 python -m pytest tests/test_lazy_gpu.py)."""

import os

import numpy as np
import pytest

from metatrap.hologram import symmetric as S
from metatrap.hologram import targets
from metatrap.hologram.exact import OPLTable
from metatrap.hologram.wgs import PupilGrid
from metatrap.propagation.layers import LayerStack


@pytest.mark.skipif(os.environ.get("NUMBA_ENABLE_CUDASIM") != "1" and not __import__("numba").cuda.is_available(),
                    reason="needs a CUDA GPU or NUMBA_ENABLE_CUDASIM=1")
def test_gpu_kernels_match_dense():
    from metatrap.hologram.lazy_gpu import LazyKGGPU

    stack = LayerStack.from_list([[2500, "air"], [4000, "fused_silica"], [9000, "air"]])
    sst = LayerStack.from_list([[2500, "air"], [4000, "fused_silica"], [9050, "air"]])
    R = stack.aperture_for_na(0.55, 1.075)
    sites = targets.make(kind="square", nx=6, ny=6, pitch_um=4.0)
    fs, fstab, inv = S.wedge_sites(sites)
    rmax = np.abs(sites).max() * np.sqrt(2)
    g = PupilGrid.for_sites(1.075, 4.0, 4, 0.55, 2 * rmax, guard=1.25)
    wp = S.WedgePupil.build(g, stack, 0.55, {"kind": "top_hat"})
    # A subset of pupil samples keeps the simulator fast.
    keep = np.linspace(0, len(wp.opl0) - 1, 300).astype(int)
    for name in ("r", "opl0", "amp", "stab", "sigma", "ij"):
        setattr(wp, name, getattr(wp, name)[keep])
    tab = OPLTable(sst, 1.075, R + 2 * rmax + 20)
    K = S.kernels_G(wp, fs, tab, 1.075)
    KG = LazyKGGPU(wp, fs, tab, 1.075, R + 2 * rmax + 20, threads=32, pairs_per_launch=8 * 300 * 2)  # forces several launches
    rng = np.random.default_rng(1)
    x = (rng.normal(size=len(fs)) + 1j * rng.normal(size=len(fs))).astype(np.complex64)
    v = (rng.normal(size=len(keep)) + 1j * rng.normal(size=len(keep))).astype(np.complex64)
    a, b = K @ x, KG @ x
    assert np.abs(a - b).max() / np.abs(a).max() < 5e-3
    a, b = np.conj(v) @ K, np.conj(v) @ KG
    assert np.abs(a - b).max() / np.abs(a).max() < 5e-3
