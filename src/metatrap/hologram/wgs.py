"""Weighted Gerchberg-Saxton design of the array hologram in the pupil domain.

Debye picture: behind an ideal lens the focal field is the Fourier transform of the
pupil function P(s) = A(s) exp(i psi(s)), with s = (sx, sy) the direction cosines in
the atom-plane medium. The pupil is sampled with spacing ds on an N x N grid, which
fixes the focal sampling dx = lambda / (N ds) and period L = lambda / ds. Sites are
placed on focal grid points (pitch = m dx), so the WGS update acts on exact site
samples, as in Kim et al., Opt. Lett. 44, 3178 (2019) [REF?: verify citation] and
the metasurface designs of proposal refs. 4 and 14.

The pupil amplitude follows from energy conservation along rays (see
LayerStack.pupil_jacobian): |A(s)| = |E_ms(r(s))| sqrt(J(s) / cos(theta)).

The Debye approximation neglects the field-dependent (off-axis) aberration of a
direct-projection hologram; metatrap.hologram.field_model quantifies it.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import scipy.fft as sfft

from ..propagation.angular_spectrum import WORKERS
from ..propagation.illumination import profile
from ..propagation.layers import LayerStack


@dataclass
class PupilGrid:
    n: int  # FFT size
    ds: float  # direction-cosine spacing
    lam_um: float

    @property
    def dx_focal(self) -> float:
        return self.lam_um / (self.n * self.ds)

    @property
    def period_um(self) -> float:
        return self.lam_um / self.ds

    def s_coords(self) -> np.ndarray:
        return (np.arange(self.n) - self.n // 2) * self.ds

    @classmethod
    def for_sites(cls, lam_um: float, pitch_um: float, oversample: int, na: float, field_um: float, guard: float = 2.0):
        """Grid with dx = pitch/oversample and a focal period >= guard x field extent.

        The direction-cosine half-extent of the grid is sigma_max = lambda / (2 dx), fixed by
        the focal sampling alone. It must exceed the NA with a margin, otherwise the pupil
        is clipped to a square of half-width sigma_max and, at rendering, metasurface pixels
        outside it receive no hologram phase and focus onto the axis as a plain lens.
        ``oversample`` is therefore raised as needed so that sigma_max >= 1.08 NA.
        (Before this check the default oversample of 4 at 4.5 um pitch gave
        sigma_max = 0.478 < NA; see docs/campaign1_report.md, erratum.)
        """
        oversample = max(int(oversample), int(np.ceil(2 * 1.08 * na * pitch_um / lam_um)))
        dx = pitch_um / oversample
        n_min_field = int(np.ceil(guard * field_um / dx))
        n = max(n_min_field, 64)
        n = int(sfft.next_fast_len(n))
        return cls(n=n, ds=lam_um / (n * dx), lam_um=lam_um)

    @property
    def s_max(self) -> float:
        return (self.n // 2 - 1) * self.ds


def pupil_amplitude(grid: PupilGrid, stack: LayerStack, na: float, illumination: dict) -> np.ndarray:
    """Amplitude A(s) on the pupil grid for a given illumination of the metasurface."""
    if grid.s_max < na + 2 * grid.ds:
        raise ValueError(f"pupil grid covers |sigma| <= {grid.s_max:.3f} < NA {na}; refine the focal sampling (oversample)")
    s1 = grid.s_coords()
    SX, SY = np.meshgrid(s1, s1, indexing="ij")
    S = np.hypot(SX, SY)
    inside = S <= na
    r = np.where(inside, stack.ray_height(np.where(inside, S, 0.0), grid.lam_um), 0.0)
    radius = stack.ray_height(na, grid.lam_um)
    kw = {k: v for k, v in illumination.items() if k != "kind"}
    e_ms = profile(illumination["kind"], r, radius, **kw)
    jac = stack.pupil_jacobian(np.where(inside, S, 0.0), grid.lam_um)
    cos_t = np.sqrt(np.clip(1 - S**2, 1e-12, None))
    a = np.where(inside, e_ms * np.sqrt(jac / cos_t), 0.0)
    # Fractional coverage of each pupil cell by the aperture disc (8 x 8 supersampling).
    # Without it the disc edge is a staircase of whole samples, and the design cancels
    # the zero order on the wrong edge geometry: a full-pixel test then shows a zero-order
    # focus of ~1e-3 of the power from the outermost few per cent of the aperture radius.
    cov = disc_coverage(SX, SY, grid.ds, na)
    edge = (cov > 0) & ~inside
    if np.any(edge):
        s_e = np.minimum(S[edge], na)
        r_e = stack.ray_height(s_e, grid.lam_um)
        e_e = profile(illumination["kind"], r_e, radius, **kw)
        a[edge] = e_e * np.sqrt(stack.pupil_jacobian(s_e, grid.lam_um) / np.sqrt(1 - s_e**2))
    a = a * np.sqrt(cov)
    return a / np.sqrt(np.sum(a**2))


def disc_coverage(SX, SY, ds, radius, sub: int = 8):
    """Fraction of each square cell (side ds, centred on SX, SY) inside the disc |s| <= radius."""
    S = np.hypot(SX, SY)
    cov = (S <= radius).astype(float)
    band = np.abs(S - radius) <= ds
    if np.any(band):
        o = (np.arange(sub) + 0.5) / sub - 0.5
        ox, oy = np.meshgrid(o, o, indexing="ij")
        bx = SX[band][:, None] + ds * ox.ravel()[None, :]
        by = SY[band][:, None] + ds * oy.ravel()[None, :]
        cov[band] = np.mean(np.hypot(bx, by) <= radius, axis=1)
    return cov


def site_indices(grid: PupilGrid, sites_um: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """FFT-array indices of the sites (origin at index 0, wrapped)."""
    dx = grid.dx_focal
    q = sites_um / dx
    if not np.allclose(q, np.round(q), atol=1e-6):
        raise ValueError("sites must lie on the focal grid (pitch must be a multiple of dx)")
    q = np.round(q).astype(int) % grid.n
    return q[:, 0], q[:, 1]


@dataclass
class WGSResult:
    phase: np.ndarray  # pupil phase psi(s) on the grid (rad)
    amplitude: np.ndarray
    history_nonuniformity: list
    history_efficiency: list
    weights: np.ndarray


def _to_focal(p):
    return sfft.ifft2(sfft.ifftshift(p), workers=WORKERS, norm="ortho")


def _to_pupil(f):
    return sfft.fftshift(sfft.fft2(f, workers=WORKERS, norm="ortho"))


def wgs(
    amplitude: np.ndarray,
    grid: PupilGrid,
    sites_um: np.ndarray,
    target_weights: np.ndarray | None = None,
    n_iter: int = 60,
    n_fix: int | None = 30,
    seed: int = 0,
    initial_phase: np.ndarray | None = None,
) -> WGSResult:
    """Weighted GS with phase fixing after ``n_fix`` iterations.

    The focal field is exp(+i k s.x), i.e. an inverse DFT of the centred pupil.
    'Efficiency' in the history is the fraction of power at the site samples relative
    to the ideal (all power in the site samples); it tracks convergence only, the
    physical efficiency is computed by hologram.evaluate on a finer grid.
    """
    rng = np.random.default_rng(seed)
    ix, iy = site_indices(grid, sites_um)
    nsite = len(ix)
    tgt = np.ones(nsite) if target_weights is None else np.asarray(target_weights, float)
    tgt = tgt / np.sqrt(np.sum(tgt**2))
    support = amplitude > 0
    if initial_phase is None:
        # Random site phases give a well-conditioned start for periodic arrays.
        f0 = np.zeros((grid.n, grid.n), complex)
        f0[ix, iy] = tgt * np.exp(2j * np.pi * rng.random(nsite))
        psi = np.angle(_to_pupil(f0))
    else:
        psi = initial_phase
    g = np.ones(nsite)
    fixed_phase = None
    hist_u, hist_e = [], []
    for it in range(n_iter):
        p = amplitude * np.exp(1j * psi)
        f = _to_focal(p)
        fs = f[ix, iy]
        a = np.abs(fs)
        hist_u.append(float(np.std(a**2 / tgt**2) / np.mean(a**2 / tgt**2)))
        hist_e.append(float(np.sum(a**2)))
        g = g * np.mean(a / tgt) / (a / tgt)
        if n_fix is not None and it == n_fix:
            fixed_phase = np.angle(fs)
        ph = fixed_phase if fixed_phase is not None else np.angle(fs)
        f_new = np.zeros_like(f)
        f_new[ix, iy] = g * tgt * np.exp(1j * ph)
        psi = np.where(support, np.angle(_to_pupil(f_new)), 0.0)
    return WGSResult(phase=psi, amplitude=amplitude, history_nonuniformity=hist_u, history_efficiency=hist_e, weights=g)


def focal_field(amplitude: np.ndarray, phase: np.ndarray, grid: PupilGrid, pad: int = 1) -> tuple[np.ndarray, float]:
    """Focal field on a grid refined ``pad`` times; returns (field, dx). Field normalised to unit power."""
    n = grid.n
    p = amplitude * np.exp(1j * phase)
    if pad > 1:
        big = np.zeros((n * pad, n * pad), complex)
        c = (n * pad) // 2 - n // 2
        big[c : c + n, c : c + n] = p
        p = big
    f = _to_focal(p)
    f = sfft.fftshift(f)  # optical axis at index N//2
    return f, grid.dx_focal / pad
