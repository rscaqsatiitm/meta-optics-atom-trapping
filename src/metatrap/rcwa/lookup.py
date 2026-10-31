"""Phase-to-width lookup and figures of merit for a meta-atom family.

For one (material, pitch, height, wavelength) the library gives t(w) over widths w.
Fabrication limits remove widths below ``min_feature_um``, gaps p - w below
``min_gap_um`` and aspect ratios h / w above ``max_aspect``.

Figure of merit: the first-order efficiency of an N-level blazed grating built from the
best width for each level, in the local-periodic approximation (LPA),

    eta_N = (n_out / n_in) |(1/N) sum_k t(w_k) exp(-i 2 pi k / N)|^2,

maximised over a global phase offset. Here, n_in and n_out are the incidence and exit
indices (substrate and air), so eta_N includes Fresnel loss, amplitude variation and
residual phase error. eta_N for N -> infinity is the efficiency of a continuous-phase
hologram rendered with this family in the LPA.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.interpolate import CubicSpline


@dataclass
class FabLimits:
    min_feature_um: float = 0.10
    min_gap_um: float = 0.10
    max_aspect: float = 10.0


@dataclass
class Lookup:
    widths_um: np.ndarray  # dense feasible widths
    t: np.ndarray  # complex transmission at those widths
    n_ratio: float  # n_out / n_in
    phase_offset: float  # global offset used by the table
    table_phase: np.ndarray  # target phases of the table bins (rad)
    table_width: np.ndarray  # chosen width for each bin
    table_t: np.ndarray  # complex t of the chosen width

    def width_for_phase(self, phi):
        idx = self._bin(phi)
        return self.table_width[idx]

    def t_for_phase(self, phi):
        idx = self._bin(phi)
        return self.table_t[idx]

    def _bin(self, phi):
        m = len(self.table_phase)
        return np.round(np.mod(np.asarray(phi), 2 * np.pi) / (2 * np.pi) * m).astype(int) % m

    def eta_levels(self, n_levels: int) -> float:
        return eta_n(self.widths_um, self.t, self.n_ratio, n_levels)[0]


def feasible(widths_um, pitch_um, height_um, lim: FabLimits) -> np.ndarray:
    w = np.asarray(widths_um)
    return (w >= lim.min_feature_um) & (pitch_um - w >= lim.min_gap_um) & (height_um / w <= lim.max_aspect)


def densify(widths_um, t, n_dense: int = 2000):
    """Cubic interpolation of Re t and Im t on a dense width grid."""
    cs_r, cs_i = CubicSpline(widths_um, t.real), CubicSpline(widths_um, t.imag)
    wd = np.linspace(widths_um[0], widths_um[-1], n_dense)
    return wd, cs_r(wd) + 1j * cs_i(wd)


def best_widths(wd, td, target_phases, offset=0.0):
    """For each target phase pick the width maximising Re(t exp(-i(phi + offset)))."""
    proj = np.real(td[None, :] * np.exp(-1j * (np.asarray(target_phases)[:, None] + offset)))
    idx = np.argmax(proj, axis=1)
    return idx


def eta_n(wd, td, n_ratio: float, n_levels: int, n_offsets: int = 64):
    """eta_N maximised over the global offset; returns (eta, offset)."""
    phases = 2 * np.pi * np.arange(n_levels) / n_levels
    best = (-1.0, 0.0)
    for off in np.linspace(0, 2 * np.pi, n_offsets, endpoint=False):
        idx = best_widths(wd, td, phases, off)
        amp = np.mean(td[idx] * np.exp(-1j * phases))
        eta = n_ratio * np.abs(amp) ** 2
        if eta > best[0]:
            best = (float(eta), float(off))
    return best


def phase_coverage(wd, td) -> float:
    """Range of the unwrapped transmission phase over the feasible widths (rad)."""
    ph = np.unwrap(np.angle(td))
    return float(ph.max() - ph.min())


def build(widths_um, t, pitch_um, height_um, n_in: float, n_out: float = 1.0, lim: FabLimits | None = None, n_bins: int = 256) -> Lookup | None:
    lim = lim or FabLimits()
    ok = feasible(widths_um, pitch_um, height_um, lim)
    if ok.sum() < 4:
        return None
    wd, td = densify(np.asarray(widths_um)[ok], np.asarray(t)[ok])
    n_ratio = n_out / n_in
    _, off = eta_n(wd, td, n_ratio, 16)
    phases = 2 * np.pi * np.arange(n_bins) / n_bins
    idx = best_widths(wd, td, phases, off)
    return Lookup(widths_um=wd, t=td, n_ratio=n_ratio, phase_offset=off, table_phase=phases, table_width=wd[idx], table_t=td[idx])


def scan_geometries(lib, material: str, lam_um: float, n_in: float, lim: FabLimits | None = None, n_levels: int = 16):
    """Rank every (pitch, height) of a library for one material and wavelength."""
    lim = lim or FabLimits()
    im = lib.materials.index(material)
    il = int(np.argmin(np.abs(lib.wavelength_um - lam_um)))
    rows = []
    for ip, p in enumerate(lib.pitch_um):
        w = lib.width_frac * p
        for ih, h in enumerate(lib.height_um):
            t = lib.t[im, ip, ih, il]
            if np.any(~np.isfinite(t)):
                continue
            ok = feasible(w, p, h, lim)
            if ok.sum() < 4:
                continue
            wd, td = densify(w[ok], t[ok])
            eta, _ = eta_n(wd, td, 1.0 / n_in, n_levels)
            rows.append({
                "pitch_um": float(p), "height_um": float(h), "eta": eta,
                "coverage_rad": phase_coverage(wd, td),
                "T0_min": float((1.0 / n_in) * np.min(np.abs(td) ** 2)),
                "T0_mean": float((1.0 / n_in) * np.mean(np.abs(td) ** 2)),
                "w_min_um": float(w[ok].min()), "w_max_um": float(w[ok].max()),
                "aspect_max": float(h / w[ok].min()),
            })
    rows.sort(key=lambda r: -r["eta"])
    return rows
