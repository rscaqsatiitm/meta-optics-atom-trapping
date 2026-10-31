"""Incident beam profiles on the metasurface (amplitude, radius in um)."""

from __future__ import annotations

import numpy as np


def top_hat(r_um, radius_um: float, edge_um: float = 0.0):
    """Flat-top amplitude with an optional raised-cosine edge of width edge_um."""
    r = np.asarray(r_um, dtype=float)
    if edge_um <= 0:
        return (r <= radius_um).astype(float)
    a = np.ones_like(r)
    zone = (r > radius_um - edge_um) & (r <= radius_um)
    a[zone] = 0.5 * (1 + np.cos(np.pi * (r[zone] - (radius_um - edge_um)) / edge_um))
    a[r > radius_um] = 0.0
    return a


def truncated_gaussian(r_um, radius_um: float, waist_um: float):
    """Gaussian amplitude exp(-r^2/w^2) (1/e^2 intensity radius w) clipped at radius_um."""
    r = np.asarray(r_um, dtype=float)
    return np.exp(-(r**2) / waist_um**2) * (r <= radius_um)


def gaussian_clip_loss(radius_um: float, waist_um: float) -> float:
    """Fraction of a Gaussian beam's power outside the aperture: exp(-2 R^2 / w^2)."""
    return float(np.exp(-2 * radius_um**2 / waist_um**2))


def profile(kind: str, r_um, radius_um: float, **kw):
    if kind == "top_hat":
        return top_hat(r_um, radius_um, kw.get("edge_um", 0.0))
    if kind == "gaussian":
        return truncated_gaussian(r_um, radius_um, kw["waist_um"])
    raise ValueError(kind)
