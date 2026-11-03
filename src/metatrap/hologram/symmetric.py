"""Exact-kernel hologram design for D4-symmetric arrays (square arrays centred on the axis).

The layer stack is rotationally symmetric and a centred square array is invariant under
the dihedral group D4 (8 elements: rotations by 90 deg and reflections). A D4-symmetric
pupil phase produces a D4-symmetric focal field, so the design needs only the pupil
samples in the wedge 0 <= sigma_y <= sigma_x and the sites with 0 < y <= x.

With K(r, x) = exp(-i k0 [OPL(r -> x) - OPL0(r)]) and the group sum
K_G(i, j) = sum_{g in D4} K(r_i, g x_j), the full-pupil WGS steps reduce to

    E_j = sum_{i in wedge} (w_i / t_i) P_i conj(K_G(i, j)) / W,
    Q_i = sum_{j in wedge} (c_j / s_j) K_G(i, j),

where w_i = A_i^2, W = sum over the full pupil of A^2, t_i is the size of the stabiliser of
pupil sample i (1 generic, 2 on the sigma_x axis or the diagonal, 8 at the origin) and
s_j that of site j (1 generic, 2 on the diagonal). One complex64 matrix K_G of size
(pupil samples / 8) x (sites / 8) serves both steps; for 20,164 sites it is ~2 GB instead
of ~140 GB for the unreduced problem.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..propagation.layers import LayerStack
from .exact import OPLTable, RenderGrid
from .wgs import PupilGrid, pupil_amplitude

# D4 acting on (x, y): identity, rotations, reflections.
D4 = [
    lambda x, y: (x, y), lambda x, y: (-y, x), lambda x, y: (-x, -y), lambda x, y: (y, -x),
    lambda x, y: (x, -y), lambda x, y: (-x, y), lambda x, y: (y, x), lambda x, y: (-y, -x),
]


def wedge_sites(sites: np.ndarray, tol: float = 1e-9):
    """Fundamental sites (0 < y <= x), their stabiliser sizes, and the map full -> fundamental index."""
    ax, ay = np.abs(sites[:, 0]), np.abs(sites[:, 1])
    hi, lo = np.maximum(ax, ay), np.minimum(ax, ay)
    if np.any(lo < tol):
        raise ValueError("sites on the x or y axis are not supported (use an even array)")
    key = np.round(np.column_stack([hi, lo]), 9)
    uniq, inv = np.unique(key, axis=0, return_inverse=True)
    stab = np.where(np.abs(uniq[:, 0] - uniq[:, 1]) < tol, 2, 1)
    return uniq, stab, inv.ravel()


@dataclass
class WedgePupil:
    grid: PupilGrid
    ij: np.ndarray  # (M, 2) grid indices of the wedge samples
    sigma: np.ndarray  # (M, 2)
    r: np.ndarray  # (M, 2) metasurface positions
    opl0: np.ndarray
    amp: np.ndarray
    stab: np.ndarray  # stabiliser sizes t_i
    w_total: float  # sum of A^2 over the full pupil

    @classmethod
    def build(cls, grid: PupilGrid, stack: LayerStack, na: float, illumination: dict, s_max: float | None = None) -> "WedgePupil":
        """Wedge samples with A > 0 (s_max None) or all wedge samples with |sigma| <= s_max (for rendering)."""
        a_full = pupil_amplitude(grid, stack, na, illumination)
        n, c = grid.n, grid.n // 2
        k = np.arange(n) - c
        I, J = np.meshgrid(k, k, indexing="ij")
        wedge = (I >= 0) & (J >= 0) & (J <= I) & (I <= c - 1)
        S = np.hypot(I, J) * grid.ds
        sel = wedge & ((a_full > 0) if s_max is None else (S <= s_max))
        ii, jj = I[sel], J[sel]
        stab = np.ones(len(ii), int)
        stab[(jj == 0) | (jj == ii)] = 2
        stab[(ii == 0) & (jj == 0)] = 8
        sx, sy = ii * grid.ds, jj * grid.ds
        # Design: partially covered edge cells use the ray mapping at the rim, as in
        # exact.PupilSampling. Rendering (s_max given): true positions out to s_max.
        s = np.minimum(np.hypot(sx, sy), na if s_max is None else 0.95)
        rho = stack.ray_height(s, grid.lam_um)
        with np.errstate(invalid="ignore", divide="ignore"):
            ux = np.where(s > 0, sx / np.where(s > 0, np.hypot(sx, sy), 1), 0.0)
            uy = np.where(s > 0, sy / np.where(s > 0, np.hypot(sx, sy), 1), 0.0)
        r = np.column_stack([-ux * rho, -uy * rho])
        return cls(grid=grid, ij=np.column_stack([ii + c, jj + c]), sigma=np.column_stack([sx, sy]), r=r,
                   opl0=stack.opl(s, grid.lam_um), amp=a_full[sel], stab=stab, w_total=float(np.sum(a_full**2)))


def kernels_G(wp: WedgePupil, fsites: np.ndarray, table: OPLTable, lam_um: float, block: int = 16) -> np.ndarray:
    """K_G (M x N_f) complex64, summing K over the 8 images of each fundamental site."""
    k0 = 2 * np.pi / lam_um
    m, nf = len(wp.opl0), len(fsites)
    K = np.zeros((m, nf), np.complex64)
    rx, ry = wp.r[:, 0:1], wp.r[:, 1:2]
    for j0 in range(0, nf, block):
        xs = fsites[j0 : j0 + block]
        acc = np.zeros((m, len(xs)), np.complex64)
        for g in D4:
            gx, gy = g(xs[:, 0], xs[:, 1])
            rho = np.hypot(rx - gx[None, :], ry - gy[None, :])
            acc += np.exp(-1j * k0 * (table(rho) - wp.opl0[:, None])).astype(np.complex64)
        K[:, j0 : j0 + block] = acc
    return K


@dataclass
class SymWGSResult:
    psi: np.ndarray
    coeffs: np.ndarray  # c_j on fundamental sites
    mode_amp: np.ndarray  # E_j on fundamental sites
    history_nonuniformity: list
    history_efficiency: list


def mode_amplitudes(wp: WedgePupil, psi: np.ndarray, K: np.ndarray) -> np.ndarray:
    v = ((wp.amp**2 / wp.stab) * np.exp(1j * psi)).astype(np.complex64)
    return np.conj(np.conj(v) @ K) / wp.w_total


def wgs_symmetric(wp: WedgePupil, K: np.ndarray, fstab: np.ndarray, n_iter: int = 60, n_fix: int | None = 30,
                  seed: int = 0, psi0: np.ndarray | None = None, target: np.ndarray | None = None,
                  ckpt_path=None, ckpt_every: int = 10) -> SymWGSResult:
    """Weighted GS on the wedge. With ckpt_path, the loop state (psi, g, fixed phases, iteration)
    is saved every ckpt_every iterations and restored on the next call with the same path, so a
    long design survives a restart; the file is removed when the loop finishes."""
    import os

    rng = np.random.default_rng(seed)
    nf = K.shape[1]
    tgt = np.ones(nf) if target is None else np.asarray(target, float)
    tgt = tgt / np.sqrt(np.mean(tgt**2))
    mult = 8.0 / fstab  # orbit sizes, weights for site statistics over the full array
    if psi0 is None:
        psi = np.angle(K @ ((tgt / fstab) * np.exp(2j * np.pi * rng.random(nf))).astype(np.complex64))
    else:
        psi = psi0.copy()
    g = np.ones(nf)
    fixed = None
    hu, he = [], []
    c = None
    it0 = 0
    if ckpt_path is not None and os.path.exists(str(ckpt_path)):
        z = np.load(str(ckpt_path))
        psi, g, it0 = z["psi"], z["g"], int(z["it"])
        fixed = z["fixed"] if bool(z["has_fixed"]) else None
        hu, he = list(z["hu"]), list(z["he"])
    for it in range(it0, n_iter):
        e = mode_amplitudes(wp, psi, K)
        a = np.abs(e)
        q = a**2 / tgt**2
        mean_q = np.sum(mult * q) / np.sum(mult)
        hu.append(float(np.sqrt(np.sum(mult * (q - mean_q) ** 2) / np.sum(mult)) / mean_q))
        he.append(float(np.sum(mult * a**2)))
        ratio = a / tgt
        g = g * (np.sum(mult * ratio) / np.sum(mult)) / ratio
        if n_fix is not None and it == n_fix:
            fixed = np.angle(e)
        th = fixed if fixed is not None else np.angle(e)
        c = (g * tgt * np.exp(1j * th)).astype(np.complex64)
        psi = np.angle(K @ (c / fstab).astype(np.complex64))
        if ckpt_path is not None and (it + 1) % ckpt_every == 0 and it + 1 < n_iter:
            tmp = str(ckpt_path) + ".tmp.npz"
            np.savez(tmp, psi=psi, g=g, it=it + 1, fixed=(fixed if fixed is not None else np.zeros(1)), has_fixed=fixed is not None,
                     hu=np.array(hu), he=np.array(he))
            os.replace(tmp, str(ckpt_path))
    if ckpt_path is not None and os.path.exists(str(ckpt_path)):
        os.remove(str(ckpt_path))
    e = mode_amplitudes(wp, psi, K)
    return SymWGSResult(psi=psi.astype(np.float64), coeffs=c, mode_amp=e, history_nonuniformity=hu, history_efficiency=he)


def unfold(grid_n: int, ij: np.ndarray, values: np.ndarray, dtype=np.complex64) -> np.ndarray:
    """Fill the full n x n grid from wedge values using D4 (grid centre at n//2)."""
    c = grid_n // 2
    out = np.zeros((grid_n, grid_n), dtype)
    a, b = ij[:, 0] - c, ij[:, 1] - c
    for g in D4:
        ga, gb = g(a, b)
        ok = (ga + c >= 0) & (ga + c < grid_n) & (gb + c >= 0) & (gb + c < grid_n)
        out[ga[ok] + c, gb[ok] + c] = values[ok]
    return out


def render_grid_symmetric(grid: PupilGrid, stack: LayerStack, na: float, illumination: dict, fsites: np.ndarray, fstab: np.ndarray,
                          coeffs: np.ndarray, table: OPLTable, lam_um: float, margin: int = 8, upsample: int = 8) -> RenderGrid:
    """Fine render grid (see exact.render_grid) built from the wedge by symmetry, in complex64."""
    wp = WedgePupil.build(grid, stack, na, illumination, s_max=na + margin * grid.ds)
    k0 = 2 * np.pi / lam_um
    q = np.zeros(len(wp.opl0), np.complex128)
    cw = coeffs / fstab
    rx, ry = wp.r[:, 0:1], wp.r[:, 1:2]
    for j0 in range(0, len(fsites), 16):
        xs = fsites[j0 : j0 + 16]
        for g in D4:
            gx, gy = g(xs[:, 0], xs[:, 1])
            rho = np.hypot(rx - gx[None, :], ry - gy[None, :])
            q += np.exp(-1j * k0 * (table(rho) - wp.opl0[:, None])) @ cw[j0 : j0 + 16]
    full = unfold(grid.n, wp.ij, q.astype(np.complex64))
    n, m = grid.n, grid.n * upsample
    import scipy.fft as sfft  # keeps complex64 (numpy would upcast to complex128)

    spec = sfft.fftshift(sfft.fft2(sfft.ifftshift(full), workers=-1))
    del full
    big = np.zeros((m, m), np.complex64)
    c = m // 2 - n // 2
    big[c : c + n, c : c + n] = spec
    del spec
    fine = sfft.fftshift(sfft.ifft2(sfft.ifftshift(big), workers=-1)).astype(np.complex64)
    del big
    s_f = (np.arange(m) - m // 2).astype(np.float32) * np.float32(grid.ds / upsample)
    r2 = s_f[:, None] ** 2 + s_f[None, :] ** 2
    lim = np.float32((na + 0.5 * margin * grid.ds) ** 2)
    mag = np.abs(fine)
    fine = np.where(r2 <= lim, fine / np.maximum(mag, 1e-30), 0).astype(np.complex64)
    return RenderGrid(field=fine, ds=grid.ds / upsample, c0=m // 2)
