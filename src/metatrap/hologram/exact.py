"""Hologram design and evaluation with exact per-site optical paths.

The Debye (FFT) design assumes that the metasurface phase which focuses to site x_j
differs from the on-axis lens phase by a tilt in direction-cosine space. The exact
phase differs by the full path difference, so here each site j has the kernel

    K_ij = exp(-i k0 [OPL(r_i -> x_j) - OPL0(r_i)]),

where r_i are metasurface sample points (the images of the pupil samples sigma_i under
the ray mapping), OPL(r -> x) is the optical path through the layer stack to atom-plane
point x, and OPL0(r) = OPL(r -> 0). Here, k0 = 2 pi / lambda. The weighted GS loop is

    E_j = sum_i A_i^2 exp(i psi_i) conj(K_ij) / sum_i A_i^2     (mode amplitude of site j)
    psi_i = arg sum_j K_ij g_j |t_j| exp(i theta_j)

with A_i the pupil amplitude. |E_j|^2 is the fraction of the incident power coupled
into an ideal, aberration-free focus at x_j (the 'mode efficiency' of site j).

Validity: scalar, high-Fresnel-number (stationary-phase) model; validated against
full-pixel angular-spectrum propagation in scripts/validate_fullpixel.py.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.interpolate import CubicSpline

from ..propagation.layers import LayerStack
from .field_model import ExactOPL
from .wgs import PupilGrid, pupil_amplitude


@dataclass
class PupilSampling:
    grid: PupilGrid
    mask: np.ndarray  # (n, n) bool, samples inside the NA
    sigma: np.ndarray  # (M, 2) direction cosines towards the axis focus
    r: np.ndarray  # (M, 2) metasurface positions (um)
    amp: np.ndarray  # (M,) pupil amplitude A_i
    opl0: np.ndarray  # (M,) on-axis optical path (um)

    @classmethod
    def build(cls, grid: PupilGrid, stack: LayerStack, na: float, illumination: dict) -> "PupilSampling":
        a_full = pupil_amplitude(grid, stack, na, illumination)
        s1 = grid.s_coords()
        SX, SY = np.meshgrid(s1, s1, indexing="ij")
        S = np.hypot(SX, SY)
        mask = a_full > 0  # includes partially covered edge cells
        sx, sy = SX[mask], SY[mask]
        s = np.minimum(S[mask], na)  # edge cells: evaluate the ray mapping at the rim
        scale = np.where(S[mask] > 0, s / np.where(S[mask] > 0, S[mask], 1), 1.0)
        sx, sy = sx * scale, sy * scale
        rho = stack.ray_height(s, grid.lam_um)
        with np.errstate(invalid="ignore", divide="ignore"):
            ux = np.where(s > 0, sx / np.where(s > 0, s, 1), 0.0)
            uy = np.where(s > 0, sy / np.where(s > 0, s, 1), 0.0)
        # sigma = -r_hat s, so r_hat = -sigma_hat.
        r = np.column_stack([-ux * rho, -uy * rho])
        return cls(grid=grid, mask=mask, sigma=np.column_stack([sx, sy]), r=r, amp=a_full[mask], opl0=stack.opl(s, grid.lam_um))

    def to_grid(self, values: np.ndarray, fill=0.0) -> np.ndarray:
        out = np.full(self.mask.shape, fill, dtype=values.dtype)
        out[self.mask] = values
        return out


class OPLTable:
    """Cubic-spline table of OPL(rho) for fast kernel evaluation (error << 1 nm)."""

    def __init__(self, stack: LayerStack, lam_um: float, rho_max: float, n: int = 400_001):
        ex = ExactOPL(stack, lam_um)
        rho = np.linspace(0.0, rho_max, n)
        opl, _ = ex(rho)
        self._cs = CubicSpline(rho, opl)

    def __call__(self, rho):
        return self._cs(rho)


def kernels(ps: PupilSampling, sites_um: np.ndarray, table: OPLTable, lam_um: float, block: int = 64) -> np.ndarray:
    """K (M x N_sites) as complex64."""
    k0 = 2 * np.pi / lam_um
    m, nsite = len(ps.amp), len(sites_um)
    K = np.empty((m, nsite), np.complex64)
    for j0 in range(0, nsite, block):
        xs = sites_um[j0 : j0 + block]
        rho = np.hypot(ps.r[:, 0:1] - xs[None, :, 0], ps.r[:, 1:2] - xs[None, :, 1])
        d = table(rho) - ps.opl0[:, None]
        K[:, j0 : j0 + block] = np.exp(-1j * k0 * d).astype(np.complex64)
    return K


def debye_kernels(ps: PupilSampling, sites_um: np.ndarray, lam_um: float, block: int = 64) -> np.ndarray:
    """Kernels of the Debye approximation, exp(-i k0 sigma . x_j), for comparison."""
    k0 = 2 * np.pi / lam_um
    K = np.empty((len(ps.amp), len(sites_um)), np.complex64)
    for j0 in range(0, len(sites_um), block):
        K[:, j0 : j0 + block] = np.exp(-1j * k0 * (ps.sigma @ sites_um[j0 : j0 + block].T)).astype(np.complex64)
    return K


@dataclass
class ExactWGSResult:
    psi: np.ndarray  # (M,) pupil phase on the samples
    mode_amp: np.ndarray  # (N,) complex site mode amplitudes E_j
    coeffs: np.ndarray  # (N,) final site coefficients c_j, psi = arg(K c)
    history_nonuniformity: list
    history_efficiency: list


def mode_amplitudes(ps: PupilSampling, psi: np.ndarray, K: np.ndarray) -> np.ndarray:
    w = (ps.amp**2).astype(np.float32)
    v = (w * np.exp(1j * psi)).astype(np.complex64)
    # sum_i v_i conj(K_ij) = conj(sum_i conj(v_i) K_ij); avoids copying conj(K).
    return np.conj(np.conj(v) @ K) / w.sum()


def wgs_exact(ps: PupilSampling, K: np.ndarray, n_iter: int = 60, n_fix: int | None = 30, seed: int = 0, psi0: np.ndarray | None = None, target=None) -> ExactWGSResult:
    rng = np.random.default_rng(seed)
    nsite = K.shape[1]
    tgt = np.ones(nsite) if target is None else np.asarray(target, float)
    tgt = tgt / np.sqrt(np.mean(tgt**2))
    if psi0 is None:
        q = K @ (tgt * np.exp(2j * np.pi * rng.random(nsite))).astype(np.complex64)
        psi = np.angle(q)
    else:
        psi = psi0.copy()
    g = np.ones(nsite)
    fixed = None
    hu, he = [], []
    for it in range(n_iter):
        e = mode_amplitudes(ps, psi, K)
        a = np.abs(e)
        hu.append(float(np.std(a**2 / tgt**2) / np.mean(a**2 / tgt**2)))
        he.append(float(np.sum(a**2)))
        g = g * np.mean(a / tgt) / (a / tgt)
        if n_fix is not None and it == n_fix:
            fixed = np.angle(e)
        th = fixed if fixed is not None else np.angle(e)
        c = (g * tgt * np.exp(1j * th)).astype(np.complex64)
        q = K @ c
        psi = np.angle(q)
    e = mode_amplitudes(ps, psi, K)
    return ExactWGSResult(psi=psi.astype(np.float64), mode_amp=e, coeffs=c, history_nonuniformity=hu, history_efficiency=he)


def focal_field_at(ps: PupilSampling, pupil_field: np.ndarray, points_um: np.ndarray, table: OPLTable, lam_um: float, block: int = 256) -> np.ndarray:
    """Exact-path pupil integral E(x) = sum_i P_i exp(i k0 [OPL(r_i->x) - OPL0_i]) / sum A_i.

    Normalised so that an aberration-free focus of unit-power pupil has |E(x_j)|^2 equal to
    the mode efficiency (1 for the plain lens at its focus). ``pupil_field`` is P_i
    on the samples (complex, amplitude included).
    """
    k0 = 2 * np.pi / lam_um
    out = np.empty(len(points_um), complex)
    w = (pupil_field * ps.amp / np.sum(ps.amp**2)).astype(np.complex64)
    for j0 in range(0, len(points_um), block):
        xs = points_um[j0 : j0 + block]
        rho = np.hypot(ps.r[:, 0:1] - xs[None, :, 0], ps.r[:, 1:2] - xs[None, :, 1])
        d = table(rho) - ps.opl0[:, None]
        out[j0 : j0 + block] = w @ np.exp(1j * k0 * d).astype(np.complex64)
    return out


@dataclass
class RenderGrid:
    """Unit-modulus hologram field exp(i psi) on a fine pupil grid, for rendering."""

    field: np.ndarray  # (m, m) complex, zero outside the computed disc
    ds: float
    c0: int  # index of sigma = 0


def render_grid(ps: PupilSampling, stack: LayerStack, sites_um: np.ndarray, coeffs: np.ndarray, lam_um: float, na: float,
                kind: str = "exact", table: OPLTable | None = None, margin: int = 8, upsample: int = 8, block: int = 64) -> RenderGrid:
    """Back-propagated pupil field Q = sum_j K(sigma, x_j) c_j on a finer grid.

    Q is evaluated on the design grid out to NA + margin samples (so that the Fourier
    interpolation does not ring inside the aperture), Fourier-upsampled by ``upsample``
    and normalised to unit modulus. Bilinear interpolation of this fine grid has a
    focal-plane envelope sinc^2(x ds_fine / lambda), i.e. > 0.98 in intensity over the
    array for upsample = 8 (DERIVED), against ~0.5 at the array corner without it.
    """
    g = ps.grid
    s1 = g.s_coords()
    SX, SY = np.meshgrid(s1, s1, indexing="ij")
    S = np.hypot(SX, SY)
    mask = S <= na + margin * g.ds
    sx, sy, s = SX[mask], SY[mask], S[mask]
    k0 = 2 * np.pi / lam_um
    q = np.zeros(len(s), np.complex128)
    if kind == "exact":
        rho_ms = stack.ray_height(s, lam_um)
        with np.errstate(invalid="ignore", divide="ignore"):
            ux = np.where(s > 0, sx / np.where(s > 0, s, 1), 0.0)
            uy = np.where(s > 0, sy / np.where(s > 0, s, 1), 0.0)
        rx, ry = -ux * rho_ms, -uy * rho_ms
        opl0 = stack.opl(s, lam_um)
        for j0 in range(0, len(sites_um), block):
            xs = sites_um[j0 : j0 + block]
            rho = np.hypot(rx[:, None] - xs[None, :, 0], ry[:, None] - xs[None, :, 1])
            q += np.exp(-1j * k0 * (table(rho) - opl0[:, None])) @ coeffs[j0 : j0 + block]
    else:
        sig = np.column_stack([sx, sy])
        for j0 in range(0, len(sites_um), block):
            q += np.exp(-1j * k0 * (sig @ sites_um[j0 : j0 + block].T)) @ coeffs[j0 : j0 + block]
    full = np.zeros(SX.shape, np.complex128)
    full[mask] = q
    n = g.n
    m = n * upsample
    spec = np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(full)))
    big = np.zeros((m, m), np.complex128)
    c = m // 2 - n // 2
    big[c : c + n, c : c + n] = spec
    fine = np.fft.fftshift(np.fft.ifft2(np.fft.ifftshift(big)))
    s_f = (np.arange(m) - m // 2) * (g.ds / upsample)
    SXf, SYf = np.meshgrid(s_f, s_f, indexing="ij")
    keep = np.hypot(SXf, SYf) <= na + 0.5 * margin * g.ds
    fine = np.where(keep, fine / np.maximum(np.abs(fine), 1e-30), 0.0)
    return RenderGrid(field=fine, ds=g.ds / upsample, c0=m // 2)
