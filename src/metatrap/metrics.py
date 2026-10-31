"""Small helpers for focal-spot and array metrics."""

from __future__ import annotations

import numpy as np


def rms_nonuniformity(values) -> float:
    """Relative rms spread std/mean (the 'rms non-uniformity' of proposal Table 1)."""
    v = np.asarray(values, dtype=float)
    return float(np.std(v) / np.mean(v))


def peak_to_peak(values) -> float:
    v = np.asarray(values, dtype=float)
    return float((v.max() - v.min()) / (v.max() + v.min()))


def rayleigh_range(waist_m: float, lam_m: float, n: float = 1.0) -> float:
    return np.pi * waist_m**2 * n / lam_m


def airy_first_zero_radius(lam_m: float, na: float) -> float:
    return 0.61 * lam_m / na


def gaussian_waist_from_na(lam_m: float, na: float) -> float:
    """1/e^2 waist of a Gaussian beam of far-field divergence half-angle asin(NA)."""
    return lam_m / (np.pi * np.tan(np.arcsin(na)))


def strehl_from_rms_phase(rms_rad):
    """Marechal approximation S ~ exp(-sigma^2)."""
    return np.exp(-np.asarray(rms_rad) ** 2)
