"""Matrix-free exact kernels for very large symmetric designs (numba).

LazyKG behaves like the dense K_G of symmetric.kernels_G (M wedge pupil samples x N_f
fundamental sites, each column summed over the 8 D4 images of its site) for the two products
the weighted GS loop needs, K @ x and conj(v) @ K, but evaluates the kernel entries on the fly:

    K_G[i, j] = sum_g exp(-i k0 [OPL(|r_i - g x_j|) - OPL0_i]).

OPL(rho) is read from a uniform table with linear interpolation. The interpolation error is
h^2 / 8 x max |d^2 OPL / d rho^2| <= h^2 / (8 z_eff); with h = 0.5 um and z_eff ~ 14 mm it is
2e-6 um, or 1e-5 rad (DERIVED). Memory is O(M + N_f) instead of O(M N_f): the 100,000-site
kernel would need ~50 GB dense.

render_grid_symmetric_lazy is symmetric.render_grid_symmetric with the site loop replaced by
one lazy product.
"""

from __future__ import annotations

import numba as nb
import numpy as np

from .symmetric import RenderGrid, WedgePupil, unfold

# D4 images (x, y) -> (sx * (x or y), sy * (y or x)), as in symmetric.D4.
_D4 = np.array([[1, 1, 0], [-1, 1, 0], [1, -1, 0], [-1, -1, 0], [1, 1, 1], [-1, 1, 1], [1, -1, 1], [-1, -1, 1]], np.int64)


@nb.njit(inline="always", fastmath=True)
def _sincos(ph):
    # Range reduction to [-pi, pi] and Taylor/Horner polynomials (error <= 1e-6, DERIVED:
    # pi^17 / 17! = 8e-7 for sin, pi^18 / 18! = 1.4e-7 for cos); vectorisable, unlike libm calls.
    x = ph - 6.283185307179586 * np.floor(ph * 0.15915494309189535 + 0.5)
    x2 = x * x
    sn = x * (1.0 + x2 * (-1.0 / 6 + x2 * (1.0 / 120 + x2 * (-1.0 / 5040 + x2 * (1.0 / 362880 + x2 * (-1.0 / 39916800
         + x2 * (1.0 / 6227020800 + x2 * (-1.0 / 1307674368000))))))))
    cs = 1.0 + x2 * (-0.5 + x2 * (1.0 / 24 + x2 * (-1.0 / 720 + x2 * (1.0 / 40320 + x2 * (-1.0 / 3628800 + x2 * (1.0 / 479001600
         + x2 * (-1.0 / 87178291200 + x2 * (1.0 / 20922789888000))))))))
    return sn, cs


@nb.njit(parallel=True, fastmath=True, cache=True)
def _matvec(rx, ry, opl0, gx, gy, xr, xi, tab, h, k0):
    # gx, gy: (8, N_f) image coordinates; returns K @ x.
    m = rx.shape[0]
    nfs = gx.shape[1]
    outr = np.zeros(m)
    outi = np.zeros(m)
    ntab = tab.shape[0]
    inv_h = 1.0 / h
    for i in nb.prange(m):
        sr = 0.0
        si = 0.0
        xi_ = rx[i]
        yi_ = ry[i]
        o0 = opl0[i]
        for g in range(8):
            for j in range(nfs):
                dx = xi_ - gx[g, j]
                dy = yi_ - gy[g, j]
                u = np.sqrt(dx * dx + dy * dy) * inv_h
                t = min(int(u), ntab - 2)
                f = u - t
                opl = tab[t] + f * (tab[t + 1] - tab[t])
                sn, cs = _sincos(-k0 * (opl - o0))
                sr += cs * xr[j] - sn * xi[j]
                si += cs * xi[j] + sn * xr[j]
        outr[i] = sr
        outi[i] = si
    return outr, outi


@nb.njit(parallel=True, fastmath=True, cache=True)
def _rmatvec(rx, ry, opl0, gx, gy, ur, ui, tab, h, k0):
    # returns sum_i conj(u_i) K_ij
    m = rx.shape[0]
    nfs = gx.shape[1]
    outr = np.zeros(nfs)
    outi = np.zeros(nfs)
    ntab = tab.shape[0]
    inv_h = 1.0 / h
    for j in nb.prange(nfs):
        sr = 0.0
        si = 0.0
        for g in range(8):
            a = gx[g, j]
            b = gy[g, j]
            for i in range(m):
                dx = rx[i] - a
                dy = ry[i] - b
                u = np.sqrt(dx * dx + dy * dy) * inv_h
                t = min(int(u), ntab - 2)
                f = u - t
                opl = tab[t] + f * (tab[t + 1] - tab[t])
                sn, cs = _sincos(-k0 * (opl - opl0[i]))
                sr += ur[i] * cs + ui[i] * sn
                si += ur[i] * sn - ui[i] * cs
        outr[j] = sr
        outi[j] = si
    return outr, outi


def _images(fx, fy):
    gx = np.empty((8, len(fx)))
    gy = np.empty((8, len(fx)))
    for g, (a, b, sw) in enumerate(_D4):
        x, y = (fy, fx) if sw else (fx, fy)
        gx[g], gy[g] = a * x, b * y
    return np.ascontiguousarray(gx), np.ascontiguousarray(gy)


class LazyKG:
    __array_ufunc__ = None  # make numpy defer `ndarray @ LazyKG` to __rmatmul__

    def __init__(self, wp: WedgePupil, fsites: np.ndarray, table, lam_um: float, rho_max: float, h_um: float = 0.5):
        self.rx = np.ascontiguousarray(wp.r[:, 0], np.float64)
        self.ry = np.ascontiguousarray(wp.r[:, 1], np.float64)
        self.opl0 = np.ascontiguousarray(wp.opl0, np.float64)
        self.gx, self.gy = _images(np.asarray(fsites[:, 0], np.float64), np.asarray(fsites[:, 1], np.float64))
        rho = np.arange(0.0, rho_max + 2 * h_um, h_um)
        self.tab = np.ascontiguousarray(table(rho), np.float64)
        self.h = float(h_um)
        self.k0 = 2 * np.pi / lam_um
        self.shape = (len(self.rx), self.gx.shape[1])
        self.dtype = np.complex64
        self.nbytes = 0  # matrix-free

    def __matmul__(self, x):
        x = np.asarray(x, np.complex128)
        r, i = _matvec(self.rx, self.ry, self.opl0, self.gx, self.gy, np.ascontiguousarray(x.real), np.ascontiguousarray(x.imag),
                       self.tab, self.h, self.k0)
        return (r + 1j * i).astype(np.complex64)

    def __rmatmul__(self, v):
        # numpy expression np.conj(u) @ K arrives here with v = conj(u); return v @ K = sum_i v_i K_ij.
        v = np.asarray(v, np.complex128)
        u = np.conj(v)
        r, i = _rmatvec(self.rx, self.ry, self.opl0, self.gx, self.gy, np.ascontiguousarray(u.real), np.ascontiguousarray(u.imag),
                        self.tab, self.h, self.k0)
        return (r + 1j * i).astype(np.complex64)


def kernels_G_lazy(wp, fsites, table, lam_um, rho_max):
    return LazyKG(wp, fsites, table, lam_um, rho_max)


def render_grid_symmetric_lazy(grid, stack, na, illumination, fsites, fstab, coeffs, table, lam_um, margin: int = 8, upsample: int = 8,
                               rho_max: float | None = None):
    """symmetric.render_grid_symmetric with the site loop done by one lazy product."""
    wp = WedgePupil.build(grid, stack, na, illumination, s_max=na + margin * grid.ds)
    rmax = float(np.max(np.hypot(wp.r[:, 0], wp.r[:, 1]))) + float(np.max(np.abs(fsites))) * 1.5 + 50 if rho_max is None else rho_max
    K = LazyKG(wp, fsites, table, lam_um, rmax)
    q = K @ (coeffs / fstab)
    full = unfold(grid.n, wp.ij, q.astype(np.complex64))
    n, m = grid.n, grid.n * upsample
    import scipy.fft as sfft

    spec = sfft.fftshift(sfft.fft2(sfft.ifftshift(full), workers=-1))
    del full
    big = np.zeros((m, m), np.complex64)
    c = m // 2 - n // 2
    big[c : c + n, c : c + n] = spec
    del spec
    fine = sfft.fftshift(sfft.ifft2(sfft.ifftshift(big), workers=-1)).astype(np.complex64)
    del big
    s_f = (np.arange(m) - m // 2).astype(np.float32) * np.float32(grid.ds / upsample)
    lim = np.float32((na + 0.5 * margin * grid.ds) ** 2)
    # Row blocks keep the temporaries small at m ~ 1e4.
    for a in range(0, m, 2048):
        r2 = s_f[a : a + 2048, None] ** 2 + s_f[None, :] ** 2
        blk = fine[a : a + 2048]
        mag = np.abs(blk)
        fine[a : a + 2048] = np.where(r2 <= lim, blk / np.maximum(mag, 1e-30), 0).astype(np.complex64)
    return RenderGrid(field=fine, ds=grid.ds / upsample, c0=m // 2)
