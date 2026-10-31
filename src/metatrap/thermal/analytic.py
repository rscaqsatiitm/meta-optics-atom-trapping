"""Analytic steady-state heating of a substrate by absorbed laser power.

Surface heat source on a semi-infinite substrate of thermal conductivity kappa
(W/m/K), all heat entering through the illuminated face, adiabatic elsewhere:

* Gaussian source q(r) = (2 P / pi w^2) exp(-2 r^2 / w^2): centre temperature rise
  dT0 = P / (sqrt(2 pi) kappa w)  [standard result; e.g. Lax, J. Appl. Phys. 48, 3919
  (1977), REF?: verify].
* Uniform disc of radius a: dT0 = P / (pi kappa a).

Both scale as P_abs / (kappa w), the scaling used in proposal footnote 1.
"""

from __future__ import annotations

import numpy as np


def gaussian_centre_rise(p_abs_W: float, waist_m: float, kappa: float) -> float:
    return p_abs_W / (np.sqrt(2 * np.pi) * kappa * waist_m)


def disc_centre_rise(p_abs_W: float, radius_m: float, kappa: float) -> float:
    return p_abs_W / (np.pi * kappa * radius_m)


def disc_rise_profile(r_m, p_abs_W: float, radius_m: float, kappa: float, n_k: int = 4000):
    """Surface temperature rise of a uniform disc source (Hankel integral, numerical)."""
    from scipy.special import j0, j1

    r = np.atleast_1d(np.asarray(r_m, dtype=float))
    q0 = p_abs_W / (np.pi * radius_m**2)
    k = np.linspace(1e-6, 60, n_k) / radius_m
    dk = k[1] - k[0]
    integrand = j1(k[None, :] * radius_m) * j0(k[None, :] * r[:, None]) / k[None, :]
    return q0 * radius_m / kappa * np.sum(integrand, axis=1) * dk
