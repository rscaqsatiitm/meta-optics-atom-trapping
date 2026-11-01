"""Damage model for the metasurface under high-power CW illumination.

At the CW intensities of interest (tens of W/cm^2), laser damage of dielectric optics is
thermal: absorbed power heats the part until it cracks, delaminates or degrades (vendor CW
damage ratings of dielectric mirrors near 1.07 um are ~2e4 W/cm^2 for mm-size beams; Thorlabs
catalogue, retrieved 2026-11-01). The criteria evaluated here, all from the steady-state
temperature field of axisym_fd.solve:

1. Substrate thermo-elastic stress. Thin free-edged disc in plane stress with the
   thickness-averaged temperature rise T(r) (Timoshenko & Goodier):
       sigma_r(r) = alpha E [ (1/R^2) int_0^R T r dr - (1/r^2) int_0^r T r dr ]
       sigma_t(r) = alpha E [ (1/R^2) int_0^R T r dr + (1/r^2) int_0^r T r dr - T(r) ],
   plus the self-equilibrated through-thickness part
       sigma_z-part(r, z) = -alpha E / (1 - nu) [T(r, z) - T_lin(r, z)],
   where T_lin is the linear fit of T across the thickness (the linear part only bends a free
   plate). The largest tensile stress is compared with strength / safety factor.
2. Film thermal-mismatch stress increment: d sigma_f = E_f / (1 - nu_f) (alpha_s - alpha_f) dT.
3. Local hot spots from fully absorbing defects or contamination of radius a: the centre
   temperature rise of a uniform disc source q = I on a half-space, dT = I a / kappa.
4. Two-photon absorption bound: alpha_2 = beta I_local; with the bulk-silicon value
   beta = 1.5 cm/GW at 1.06 um (an upper bound for the wider-gap nitrides) and a local field
   enhancement factor, compared with the linear absorption coefficient.
"""

from __future__ import annotations

import numpy as np

from .axisym_fd import Substrate


def disc_stresses(rc, T_bar, R, sub: Substrate):
    """Membrane stresses (Pa) of a free-edged thin disc for thickness-averaged temperature T_bar(r)."""
    aE = sub.alpha_exp * sub.youngs_GPa * 1e9
    edges = np.concatenate([[0.0], 0.5 * (rc[1:] + rc[:-1]), [R]])
    dA = 0.5 * (edges[1:] ** 2 - edges[:-1] ** 2)  # int r dr per cell
    cum = np.cumsum(T_bar * dA)  # int_0^r T r dr at outer cell edges
    prev = np.concatenate([[0.0], cum[:-1]])
    inner = prev + T_bar * 0.5 * (rc**2 - edges[:-1] ** 2)  # exact to the cell centre for piecewise-constant T
    total = cum[-1]
    with np.errstate(invalid="ignore", divide="ignore"):
        avg_r = np.where(rc > 0, inner / rc**2, 0.5 * T_bar[0])
    sr = aE * (total / R**2 - avg_r)
    st = aE * (total / R**2 + avg_r - T_bar)
    st_edge = aE * (2 * total / R**2 - T_bar[-1])
    return sr, st, st_edge


def through_thickness_stress(zc, T, sub: Substrate):
    """Max tensile self-equilibrated stress (Pa) from the non-linear part of T across the thickness."""
    aE = sub.alpha_exp * sub.youngs_GPa * 1e9 / (1 - sub.poisson)
    A = np.vstack([np.ones_like(zc), zc]).T
    coef, *_ = np.linalg.lstsq(A, T.T, rcond=None)
    T_lin = (A @ coef).T
    sig = -aE * (T - T_lin)
    return float(sig.max())


def film_mismatch_stress_per_K(E_f_GPa, nu_f, alpha_f, sub: Substrate):
    """Film stress increment per kelvin of heating (Pa/K); negative = compressive."""
    return E_f_GPa * 1e9 / (1 - nu_f) * (sub.alpha_exp - alpha_f)


def defect_rise(intensity_W_m2, radius_m, kappa):
    """Centre temperature rise (K) of a fully absorbing disc defect on a half-space."""
    return intensity_W_m2 * radius_m / kappa


def tpa_over_linear(intensity_W_m2, enhancement, beta_cm_per_GW, alpha_lin_per_cm):
    """Ratio of two-photon to linear absorption coefficient."""
    i_gw_cm2 = intensity_W_m2 * enhancement / 1e13
    return beta_cm_per_GW * i_gw_cm2 / alpha_lin_per_cm
