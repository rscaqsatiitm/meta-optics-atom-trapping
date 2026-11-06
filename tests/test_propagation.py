"""Analytic anchors for the layer stack and the angular-spectrum propagator."""

import numpy as np
import pytest

from metatrap.propagation.angular_spectrum import grid_coords, propagate
from metatrap.propagation.layers import LayerStack

LAM = 1.075


def test_ray_mapping_vacuum_is_tangent():
    st = LayerStack.from_list([(1000.0, "air")])
    s = np.array([0.1, 0.3, 0.5])
    assert np.allclose(st.ray_height(s, LAM), 1000.0 * np.tan(np.arcsin(s)))


def test_lens_phase_obeys_generalised_snell():
    st = LayerStack.from_list([(2500, "air"), (4000, "fused_silica"), (9000, "air")])
    r = np.linspace(100, 9000, 50)
    h = 1e-3
    dphi = (st.lens_phase(r + h, LAM) - st.lens_phase(r - h, LAM)) / (2 * h)
    s = st.s_of_r(LAM)(r)
    assert np.allclose(dphi, -2 * np.pi / LAM * s, rtol=1e-5)


def test_paraxial_plate_shift():
    st = LayerStack.from_list([(4000, "fused_silica")])
    n = st.indices(LAM)[0]
    assert st.paraxial_focal_shift_um(LAM) == pytest.approx(4000 * (1 - 1 / n))


def test_gaussian_beam_waist_and_gouy():
    """Free-space propagation of a Gaussian waist reproduces w(z) and the Gouy phase."""
    n, dx, w0 = 512, 0.5, 6.0
    x = grid_coords(n, dx)
    X, Y = np.meshgrid(x, x, indexing="ij")
    u0 = np.exp(-(X**2 + Y**2) / w0**2).astype(np.complex64)
    z = 200.0
    u = propagate(u0.copy(), dx, LAM, LayerStack.from_list([(z, "air")]), remove_piston=True)
    zr = np.pi * w0**2 / LAM
    wz = w0 * np.sqrt(1 + (z / zr) ** 2)
    I = np.abs(u) ** 2
    # second moment: <x^2> = w^2 / 4 for intensity exp(-2 x^2 / w^2)
    w_meas = 2 * np.sqrt(np.sum(I * X**2) / np.sum(I))
    assert w_meas == pytest.approx(wz, rel=2e-3)
    gouy = np.angle(u[n // 2, n // 2])  # piston of the kz transfer removed -> pure Gouy + paraxial corrections
    assert gouy == pytest.approx(-np.arctan(z / zr), abs=5e-3)


def test_ideal_lens_through_plate_focuses_to_airy():
    """A lens phase from LayerStack.lens_phase focuses through a plate with Strehl ~1."""
    lam, na = LAM, 0.3
    st = LayerStack.from_list([(150.0, "air"), (100.0, "fused_silica"), (100.0, "air")])
    R = st.aperture_for_na(na, lam)
    dx = 0.25
    n = int(2 ** np.ceil(np.log2(2.4 * R / dx)))
    x = grid_coords(n, dx)
    X, Y = np.meshgrid(x, x, indexing="ij")
    r = np.hypot(X, Y)
    inside = r <= R
    u = np.where(inside, np.exp(1j * st.lens_phase(np.where(inside, r, 0), lam)), 0).astype(np.complex64)
    p_in = np.sum(np.abs(u) ** 2) * dx**2
    out = propagate(u, dx, lam, st)
    I = np.abs(out) ** 2
    c = n // 2
    i_peak = I.max()
    # Airy peak for uniform pupil (paraxial): pi NA^2 P / lambda^2; high-NA corrections ~(1 - NA^2/4 ...)
    i_airy = np.pi * na**2 * p_in / lam**2
    assert np.unravel_index(np.argmax(I), I.shape) == (c, c)
    assert i_peak / i_airy == pytest.approx(1.0, abs=0.06)


def test_uncorrected_plate_defocuses():
    """The same lens designed for air only focuses ~d(1-1/n) further behind a plate."""
    lam, na = LAM, 0.3
    plate = LayerStack.from_list([(150.0, "air"), (100.0, "fused_silica"), (100.0, "air")])
    air_only = LayerStack.from_list([(350.0, "air")])
    shift = 100.0 * (1 - 1 / plate.indices(lam)[1])
    assert plate.paraxial_focal_shift_um(lam) == pytest.approx(shift)
    R = air_only.aperture_for_na(na, lam)
    # Marginal-ray focus of the air design behind the plate is beyond the plate focus.
    assert plate.ray_height(na, lam) < R
