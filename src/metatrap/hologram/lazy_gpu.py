"""GPU (CUDA, numba.cuda) version of the matrix-free exact kernels of hologram/lazy.py.

LazyKGGPU has the interface of lazy.LazyKG (K @ x, conj(v) @ K, .shape) and can replace it in
the D4-symmetric weighted GS (scripts/fullaperture_lazy.py --gpu). Design choices for consumer
GPUs, whose double precision runs at 1/32 to 1/64 of single precision:

* All per-pair arithmetic is single precision. The optical-path table and the reference paths
  OPL0 are stored relative to a common offset C = OPL(0), so the values are <= ~10^4 um and the
  float32 rounding of a path difference is <~ 1e-3 um (phase error <~ 6e-3 rad at 1.075 um,
  DERIVED). The phase is reduced to [-pi, pi] before __sinf/__cosf-class evaluation (here
  math.sin / math.cos, which numba maps to fast single-precision intrinsics with fastmath).
* Sums run per thread in float32 with Kahan compensation; the 10^6-term sums of the
  10^6-site design then carry relative errors ~1e-7 instead of ~1e-4. Partial sums of
  successive launches are accumulated in float64.
* Each product is split into launches of about `pairs_per_launch` (default 2e9) kernel
  evaluations, i.e. well under a second on a consumer GPU. A GPU that also drives a display
  under Windows (including WSL2) is reset by the driver (TDR) if one launch exceeds ~2 s.

Test without a GPU: NUMBA_ENABLE_CUDASIM=1 runs the same kernels on the CPU (slowly) and
tests/test_lazy_gpu.py compares them with the CPU kernels.
"""

from __future__ import annotations

import math

import numpy as np
from numba import cuda

from .lazy import _images

TWO_PI = np.float32(2 * np.pi)
INV_TWO_PI = np.float32(1 / (2 * np.pi))


@cuda.jit(fastmath=True)
def _matvec_kernel(rx, ry, opl0, gx, gy, xr, xi, tab, inv_h, k0, j0, j1, outr, outi):
    i = cuda.grid(1)
    if i >= rx.shape[0]:
        return
    ntab = tab.shape[0]
    sr = np.float32(0.0)
    si = np.float32(0.0)
    cr_ = np.float32(0.0)  # Kahan compensations
    ci_ = np.float32(0.0)
    xi0 = rx[i]
    yi0 = ry[i]
    o0 = opl0[i]
    for g in range(8):
        for j in range(j0, j1):
            dx = xi0 - gx[g, j]
            dy = yi0 - gy[g, j]
            u = math.sqrt(dx * dx + dy * dy) * inv_h
            t = int(u)
            if t > ntab - 2:
                t = ntab - 2
            f = u - t
            opl = tab[t] + f * (tab[t + 1] - tab[t])
            ph = -k0 * (opl - o0)
            ph = ph - TWO_PI * math.floor(ph * INV_TWO_PI + np.float32(0.5))
            c = math.cos(ph)
            s = math.sin(ph)
            yr = c * xr[j] - s * xi[j] - cr_
            tr = sr + yr
            cr_ = (tr - sr) - yr
            sr = tr
            yim = c * xi[j] + s * xr[j] - ci_
            ti = si + yim
            ci_ = (ti - si) - yim
            si = ti
    outr[i] += sr  # float64 accumulators across launches
    outi[i] += si


@cuda.jit(fastmath=True)
def _rmatvec_kernel(rx, ry, opl0, gx, gy, ur, ui, tab, inv_h, k0, i0, i1, outr, outi):
    j = cuda.grid(1)
    if j >= gx.shape[1]:
        return
    ntab = tab.shape[0]
    sr = np.float32(0.0)
    si = np.float32(0.0)
    cr_ = np.float32(0.0)
    ci_ = np.float32(0.0)
    for g in range(8):
        a = gx[g, j]
        b = gy[g, j]
        for i in range(i0, i1):
            dx = rx[i] - a
            dy = ry[i] - b
            u = math.sqrt(dx * dx + dy * dy) * inv_h
            t = int(u)
            if t > ntab - 2:
                t = ntab - 2
            f = u - t
            opl = tab[t] + f * (tab[t + 1] - tab[t])
            ph = -k0 * (opl - opl0[i])
            ph = ph - TWO_PI * math.floor(ph * INV_TWO_PI + np.float32(0.5))
            c = math.cos(ph)
            s = math.sin(ph)
            yr = ur[i] * c + ui[i] * s - cr_
            tr = sr + yr
            cr_ = (tr - sr) - yr
            sr = tr
            yim = ur[i] * s - ui[i] * c - ci_
            ti = si + yim
            ci_ = (ti - si) - yim
            si = ti
    outr[j] += sr
    outi[j] += si


class LazyKGGPU:
    """GPU counterpart of lazy.LazyKG (same products and normalisation)."""

    __array_ufunc__ = None

    def __init__(self, wp, fsites, table, lam_um: float, rho_max: float, h_um: float = 0.5, threads: int = 128,
                 pairs_per_launch: float = 2e9):
        rho = np.arange(0.0, rho_max + 2 * h_um, h_um)
        tab = np.asarray(table(rho), np.float64)
        C = tab[0]
        self.d_tab = cuda.to_device((tab - C).astype(np.float32))
        self.d_rx = cuda.to_device(np.asarray(wp.r[:, 0], np.float32))
        self.d_ry = cuda.to_device(np.asarray(wp.r[:, 1], np.float32))
        self.d_opl0 = cuda.to_device((np.asarray(wp.opl0, np.float64) - C).astype(np.float32))
        gx, gy = _images(np.asarray(fsites[:, 0], np.float64), np.asarray(fsites[:, 1], np.float64))
        self.d_gx = cuda.to_device(gx.astype(np.float32))
        self.d_gy = cuda.to_device(gy.astype(np.float32))
        self.inv_h = np.float32(1.0 / h_um)
        self.k0 = np.float32(2 * np.pi / lam_um)
        self.m, self.nf = len(wp.opl0), gx.shape[1]
        self.shape = (self.m, self.nf)
        self.dtype = np.complex64
        self.nbytes = 0
        self.threads = threads
        self.pairs_per_launch = pairs_per_launch

    def __matmul__(self, x):
        x = np.asarray(x, np.complex64)
        outr = cuda.to_device(np.zeros(self.m, np.float64))
        outi = cuda.to_device(np.zeros(self.m, np.float64))
        d_xr, d_xi = cuda.to_device(np.ascontiguousarray(x.real)), cuda.to_device(np.ascontiguousarray(x.imag))
        blocks = (self.m + self.threads - 1) // self.threads
        step = max(1, int(self.pairs_per_launch // (8 * self.m)))
        for j0 in range(0, self.nf, step):
            _matvec_kernel[blocks, self.threads](self.d_rx, self.d_ry, self.d_opl0, self.d_gx, self.d_gy, d_xr, d_xi,
                                                 self.d_tab, self.inv_h, self.k0, j0, min(j0 + step, self.nf), outr, outi)
            cuda.synchronize()
        return (outr.copy_to_host() + 1j * outi.copy_to_host()).astype(np.complex64)

    def __rmatmul__(self, v):
        u = np.conj(np.asarray(v, np.complex64))
        outr = cuda.to_device(np.zeros(self.nf, np.float64))
        outi = cuda.to_device(np.zeros(self.nf, np.float64))
        d_ur, d_ui = cuda.to_device(np.ascontiguousarray(u.real)), cuda.to_device(np.ascontiguousarray(u.imag))
        blocks = (self.nf + self.threads - 1) // self.threads
        step = max(1, int(self.pairs_per_launch // (8 * self.nf)))
        for i0 in range(0, self.m, step):
            _rmatvec_kernel[blocks, self.threads](self.d_rx, self.d_ry, self.d_opl0, self.d_gx, self.d_gy, d_ur, d_ui,
                                                  self.d_tab, self.inv_h, self.k0, i0, min(i0 + step, self.m), outr, outi)
            cuda.synchronize()
        return (outr.copy_to_host() + 1j * outi.copy_to_host()).astype(np.complex64)
