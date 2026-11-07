"""Power and loss budget for cavity-enhanced lattices and the tweezer baselines.

Cavity. A ring cavity with input-coupler power transmission T and all other round-trip
power losses L (mirror transmission, scatter and absorption, AR surfaces, diffraction) has,
on resonance and for a perfectly mode-matched input, the circulating-power enhancement

    G = T / (1 - r)^2,      r = sqrt((1 - T)(1 - L)),

with r the round-trip amplitude factor. G is maximal for the impedance-matched coupler
T = L, where G = 1 / L exactly. The finesse is F = pi sqrt(r) / (1 - r) (~ pi G for T = L).
Any element inside the cavity counts with every pass per round trip.

Lattice (sheet beams). One build-up cavity per lattice axis whose mode crosses the atom
plane twice at a full angle 2 theta (a self-intersecting ring, travelling wave), so the two
crossing beams carry the same circulating power P_c. Each beam is a sheet: flat-top of
equivalent width W in the atom plane, Gaussian of 1/e^2 half-thickness w_z perpendicular to
it, peak intensity I = P_c / (W sqrt(pi/2) w_z). The two beams interfere with contrast c
(c = 1 for polarisation perpendicular to the crossing plane), so the intensity modulation,
i.e. the barrier between neighbouring sites along this axis, is 4 c sqrt(I_a I_b). The
period is a = lambda / (2 sin theta).

Tweezer baselines (SIM, this repository). Direct meta-hologram at NA 0.55, 1075 nm:
143 mK x site / W including its array efficiency (docs/scaling_100k.md, from the 20,164-site
v2 design; the 99,856-site design reproduces it to 2 %). Lenslet array relayed by the
meta-objective: 80 (square) to 95 (hexagonal) mK x site / W (docs/million_sites.md, Section 4).
"""

from __future__ import annotations

import numpy as np

from ..atoms.alkali import depth_K
from ..constants import AMU, CESIUM, KB, AlkaliLines

HOLOGRAM_MK_SITE_PER_W = 143.0  # SIM, docs/scaling_100k.md (Method), results/scaling_100k
LENSLET_SQUARE_MK_SITE_PER_W = 80.0  # SIM/DERIVED, docs/million_sites.md Section 4
LENSLET_HEX_MK_SITE_PER_W = 95.0  # SIM/DERIVED, docs/million_sites.md Section 4
TWEEZER_WAIST_UM = 0.866  # SIM, mean curvature waist, library pillars, results/fullaperture_100k


# ---------------------------------------------------------------- cavity
def roundtrip_amplitude(t_in, loss_rt):
    return np.sqrt((1.0 - np.asarray(t_in)) * (1.0 - np.asarray(loss_rt)))


def enhancement(t_in, loss_rt):
    """Circulating / input power on resonance for a mode-matched ring cavity."""
    r = roundtrip_amplitude(t_in, loss_rt)
    return np.asarray(t_in) / (1.0 - r) ** 2


def g_max(loss_rt):
    """Enhancement with the impedance-matched coupler T = L (exactly 1 / L)."""
    return 1.0 / np.asarray(loss_rt)


def finesse(t_in, loss_rt):
    r = roundtrip_amplitude(t_in, loss_rt)
    return np.pi * np.sqrt(r) / (1.0 - r)


def coupler_for_enhancement(g_target, loss_rt):
    """Over-coupled input transmission T >= L that gives the enhancement g_target (< 1/L).

    Over-coupling (T > L) is the robust choice: lower finesse, wider linewidth, and the
    enhancement depends less on the exact loss. Solved by bisection on T in [L, 1).
    """
    if g_target >= 1.0 / loss_rt:
        raise ValueError("target above the impedance-matched maximum 1/L")
    lo, hi = loss_rt, 0.999999
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if enhancement(mid, loss_rt) > g_target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def roundtrip_loss(elements: dict) -> float:
    """Total round-trip power loss from {name: (count_per_roundtrip, loss_each)}.

    Losses combine multiplicatively: L = 1 - prod (1 - l_i)^n_i.
    """
    keep = 1.0
    for _, (count, loss) in elements.items():
        keep *= (1.0 - float(loss)) ** float(count)
    return 1.0 - keep


def intracavity_element_bound(pass_efficiency: float, passes_per_roundtrip: int) -> float:
    """Upper bound on G when an element of per-pass efficiency eta sits in the cavity."""
    return 1.0 / (1.0 - pass_efficiency**passes_per_roundtrip)


# ---------------------------------------------------------------- lattice
def crossing_half_angle(period_um: float, lam_um: float) -> float:
    """Half-angle theta (rad) between two beams that interfere with period a."""
    return float(np.arcsin(lam_um / (2.0 * period_um)))


def sheet_peak_intensity(power_W, width_eq_m, w_z_m):
    """Peak intensity (W/m^2) of a flat-top x Gaussian sheet of power P."""
    return np.asarray(power_W) / (width_eq_m * np.sqrt(np.pi / 2.0) * w_z_m)


def flat_flat_peak_intensity(power_W, width_eq_m, height_eq_m):
    """Peak intensity (W/m^2) of a beam that is flat-top in both transverse directions."""
    return np.asarray(power_W) / (width_eq_m * height_eq_m)


def depth_per_intensity_K(lam_nm: float, atom: AlkaliLines = CESIUM) -> float:
    """Trap depth (K) per W/m^2 of intensity for a red-detuned trap."""
    return float(depth_K(1.0, lam_nm, atom))


def barrier_depth_K(i_a, i_b, lam_nm: float, contrast: float = 1.0, atom: AlkaliLines = CESIUM):
    """Barrier between neighbouring sites of a two-beam lattice: 4 c sqrt(I_a I_b) x (U/I)."""
    return 4.0 * contrast * np.sqrt(np.asarray(i_a) * np.asarray(i_b)) * depth_per_intensity_K(lam_nm, atom)


def beam_intensity_for_barrier(depth_K_, lam_nm: float, contrast: float = 1.0, atom: AlkaliLines = CESIUM):
    """Intensity of each of two equal beams for a given barrier depth."""
    return np.asarray(depth_K_) / (4.0 * contrast * depth_per_intensity_K(lam_nm, atom))


def lattice_omega(depth_K_, period_m, atom: AlkaliLines = CESIUM):
    """Harmonic angular frequency of a cos^2 lattice site, V = -U cos^2(pi x / a)."""
    m = atom.mass_amu * AMU
    return (np.pi / period_m) * np.sqrt(2.0 * np.asarray(depth_K_) * KB / m)


def tweezer_omega_radial(depth_K_, waist_m, atom: AlkaliLines = CESIUM):
    m = atom.mass_amu * AMU
    return np.sqrt(4.0 * np.asarray(depth_K_) * KB / (m * waist_m**2))


def depth_for_equal_omega(depth_tweezer_K, waist_m, period_m):
    """Lattice barrier depth giving the same radial frequency as a tweezer of the given depth."""
    return depth_tweezer_K * 2.0 * period_m**2 / (np.pi**2 * waist_m**2)


# ---------------------------------------------------------------- tweezer baselines
def tweezer_power_W(n_sites, depth_mK, mk_site_per_W=HOLOGRAM_MK_SITE_PER_W):
    return np.asarray(n_sites) * np.asarray(depth_mK) / mk_site_per_W


def retro_tweezer_gain(eta_return: float) -> float:
    """Peak-intensity gain of a tweezer whose light is returned onto itself with power
    efficiency eta_return (perfect mode overlap, antinode at focus): (1 + sqrt(eta))^2."""
    return (1.0 + np.sqrt(eta_return)) ** 2
