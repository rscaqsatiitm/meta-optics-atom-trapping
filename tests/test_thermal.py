import numpy as np
import pytest

from metatrap.thermal import analytic
from metatrap.thermal.axisym_fd import FUSED_SILICA, Substrate, ThermalProblem, solve, surface_temperature


def test_thin_plate_uniform_disc_matches_analytic():
    sub = Substrate("test", kappa=1.0, dndT=0, alpha_exp=0, poisson=0, n=1.5)
    R, t, a, P = 10e-3, 0.2e-3, 3e-3, 0.1
    prob = ThermalProblem(radius_m=R, thickness_m=t, substrate=sub, beam="top_hat", beam_radius_m=a, h_faces=0.0, nr=300, nz=6)
    rc, zc, T = solve(prob, P, a_film=1.0)
    expect = P / (4 * np.pi * sub.kappa * t) * (1 + 2 * np.log(R / a))
    assert T[0].mean() == pytest.approx(expect, rel=0.02)


def test_thick_substrate_gaussian_matches_half_space():
    """Thick, wide substrate approaches the half-space result P/(sqrt(2 pi) kappa w)."""
    sub = Substrate("test", kappa=1.0, dndT=0, alpha_exp=0, poisson=0, n=1.5)
    w = 0.5e-3
    prob = ThermalProblem(radius_m=40e-3, thickness_m=40e-3, substrate=sub, beam="gaussian", beam_radius_m=4 * w, beam_waist_m=w, h_faces=0.0, nr=320, nz=160, z_power=3.0)
    # Refine z near the surface by using many cells; tolerance accounts for the finite domain.
    rc, zc, T = solve(prob, 1.0, a_film=1.0)
    expect = analytic.gaussian_centre_rise(1.0, w, 1.0)
    assert surface_temperature(zc, T)[0] == pytest.approx(expect, rel=0.05)


def test_chi_fused_silica():
    assert FUSED_SILICA.chi == pytest.approx(9.6e-6 + 0.4495 * 1.16 * 0.52e-6)


def test_disc_stress_uniform_and_hot_core():
    from metatrap.thermal.damage import disc_stresses

    R = 10e-3
    rc = (np.arange(4000) + 0.5) / 4000 * R
    sr, st, st_edge = disc_stresses(rc, np.full_like(rc, 7.0), R, FUSED_SILICA)
    assert np.max(np.abs(sr)) < 1e-6 * 7.0 * 37.8e3 + 1.0 and abs(st_edge) < 1.0  # uniform heating: stress free
    a, T0 = 3e-3, 10.0
    T = np.where(rc < a, T0, 0.0)
    _, _, st_edge = disc_stresses(rc, T, R, FUSED_SILICA)
    aE = FUSED_SILICA.alpha_exp * FUSED_SILICA.youngs_GPa * 1e9
    assert st_edge == pytest.approx(aE * T0 * a**2 / R**2, rel=1e-3)
