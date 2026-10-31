"""Ground-state dynamic polarisability and dipole-trap depth for alkali atoms.

Model: the two D lines are summed explicitly with Steck's matrix elements, and all
other contributions (core, higher np states) are lumped into a frequency-independent
remainder fixed by Steck's measured static polarisability,

    alpha(w) = (2/3) sum_{D1,D2} |d_k|^2 w_k / (w_k^2 - w^2) + alpha_rest,
    alpha_rest = alpha0_measured - (2/3) sum |d_k|^2 / w_k.

Here, d_k is the reduced matrix element in Steck's convention (units e a0), w_k and
w are angular frequencies in atomic units, and alpha is in atomic units (a0^3).
The constant remainder neglects the dispersion of the 7p (Cs) or 6p (Rb) terms;
for Cs these lie near 455-459 nm, so at 1064 nm their ~2 a.u. static contribution
grows by ~1.2x, an error below 0.1% of alpha (DERIVED, order of magnitude).
Vector and tensor parts vanish for linearly polarised light on a J = 1/2 state.

Trap depth for peak intensity I (W/m^2): U = -(2 pi a0^3 / c) alpha_au I.
"""

from __future__ import annotations

import numpy as np

from ..constants import A0, C, CESIUM, HARTREE_CM, HZ_PER_VCM2_TO_AU, KB, AlkaliLines

# U/k_B in kelvin per (atomic unit of polarisability x W/m^2).
DEPTH_K_PER_AU_PER_WM2 = 2.0 * np.pi * A0**3 / C / KB


def _omega_au(lam_nm):
    """Angular frequency in atomic units for a vacuum wavelength in nm."""
    return (1.0e7 / np.asarray(lam_nm, dtype=float)) / HARTREE_CM


def alpha_static_measured_au(atom: AlkaliLines = CESIUM) -> float:
    return atom.alpha0_hz_per_vcm2 * HZ_PER_VCM2_TO_AU


def alpha_rest_au(atom: AlkaliLines = CESIUM) -> float:
    """Static remainder (core + higher states) implied by the measured alpha0."""
    w1, w2 = _omega_au(atom.lambda_d1_nm), _omega_au(atom.lambda_d2_nm)
    d_lines = (2.0 / 3.0) * (atom.d1_au**2 / w1 + atom.d2_au**2 / w2)
    return alpha_static_measured_au(atom) - d_lines


def alpha_au(lam_nm, atom: AlkaliLines = CESIUM):
    """Scalar ground-state polarisability (atomic units) at vacuum wavelength lam_nm."""
    w = _omega_au(lam_nm)
    w1, w2 = _omega_au(atom.lambda_d1_nm), _omega_au(atom.lambda_d2_nm)
    a = (2.0 / 3.0) * (
        atom.d1_au**2 * w1 / (w1**2 - w**2) + atom.d2_au**2 * w2 / (w2**2 - w**2)
    )
    return a + alpha_rest_au(atom)


def depth_K(intensity_W_m2, lam_nm, atom: AlkaliLines = CESIUM):
    """Trap depth |U|/k_B in kelvin for a red-detuned trap at peak intensity I."""
    return DEPTH_K_PER_AU_PER_WM2 * alpha_au(lam_nm, atom) * np.asarray(intensity_W_m2)


def intensity_for_depth(depth_K_, lam_nm, atom: AlkaliLines = CESIUM):
    """Peak intensity (W/m^2) needed for a trap depth (K)."""
    return np.asarray(depth_K_) / (DEPTH_K_PER_AU_PER_WM2 * alpha_au(lam_nm, atom))


def airy_peak_intensity(power_W, lam_m, na):
    """Peak intensity of a uniformly filled, aberration-free focus (paraxial).

    I0 = pi NA^2 P / lambda^2, with P the total power in the focal pattern.
    This is the relation used for the proposal's site-count model (Section 4 (ii)).
    """
    return np.pi * na**2 * np.asarray(power_W) / lam_m**2


def gaussian_peak_intensity(power_W, waist_m):
    """I0 = 2 P / (pi w0^2) for a Gaussian focus of 1/e^2 radius w0."""
    return 2.0 * np.asarray(power_W) / (np.pi * waist_m**2)


def gaussian_trap_frequencies_Hz(depth_K_, waist_m, lam_m, atom: AlkaliLines = CESIUM):
    """Radial and axial trap frequencies (Hz) in the harmonic Gaussian-beam approximation."""
    from ..constants import AMU

    m = atom.mass_amu * AMU
    u = np.asarray(depth_K_) * KB
    zr = np.pi * waist_m**2 / lam_m
    wr = np.sqrt(4.0 * u / (m * waist_m**2))
    wz = np.sqrt(2.0 * u / (m * zr**2))
    return wr / (2 * np.pi), wz / (2 * np.pi)
