import numpy as np
import pytest

from metatrap.cavity import budget as B
from metatrap.cavity import foxli as FL
from metatrap.cavity import lattice as LA


def test_enhancement_limits():
    # Impedance-matched coupler: G = 1/L exactly; finesse ~ pi G.
    for L in (1e-4, 1e-3, 1e-2):
        assert B.enhancement(L, L) == pytest.approx(1 / L, rel=1e-9)
        assert B.finesse(L, L) == pytest.approx(np.pi / L, rel=1e-2)
    # Lossless cavity, small T: G -> 4 / T.
    assert B.enhancement(1e-3, 0.0) == pytest.approx(4e3, rel=1e-3)
    # Over-coupled design solver.
    t = B.coupler_for_enhancement(300.0, 3e-4)
    assert t > 3e-4 and B.enhancement(t, 3e-4) == pytest.approx(300.0, rel=1e-6)
    # Loss bookkeeping and the intracavity-element bound.
    assert B.roundtrip_loss({"m": (4, 1e-4), "ar": (2, 1e-3)}) == pytest.approx(1 - (1 - 1e-4) ** 4 * (1 - 1e-3) ** 2)
    assert B.intracavity_element_bound(0.83, 2) == pytest.approx(1 / (1 - 0.83**2))


def test_two_beam_lattice_period_and_barrier():
    lam, a = 1.075, 4.0
    th = B.crossing_half_angle(a, lam)
    k = 2 * np.pi / lam
    fam = [[(1.0, k * np.array([np.cos(th), np.sin(th), 0.0])), (1.0, k * np.array([np.cos(th), -np.sin(th), 0.0]))]]
    y = np.linspace(0, 3 * a, 3001)
    pts = np.stack([np.zeros_like(y), y, np.zeros_like(y)], axis=-1)
    i = LA.plane_wave_intensity(pts, [fam[0]])
    # Peaks of 4 I1 every a along y; barrier (max - min) = 4 sqrt(I_a I_b).
    assert i.max() == pytest.approx(4.0, rel=1e-6) and i.min() == pytest.approx(0.0, abs=1e-6)
    peaks = y[1:-1][(i[1:-1] > i[:-2]) & (i[1:-1] > i[2:])]
    assert np.diff(peaks).mean() == pytest.approx(a, rel=1e-3)
    assert B.barrier_depth_K(1.0, 1.0, 1075.0) == pytest.approx(4 * B.depth_per_intensity_K(1075.0))


def test_lattice_omega_matches_curvature():
    u, a = 2e-4, 4e-6
    from metatrap.constants import AMU, CESIUM, KB

    m = CESIUM.mass_amu * AMU
    x = np.linspace(-1e-8, 1e-8, 3)
    v = -u * KB * np.cos(np.pi * x / a) ** 2
    curv = (v[0] - 2 * v[1] + v[2]) / (x[1] - x[0]) ** 2
    assert B.lattice_omega(u, a) == pytest.approx(np.sqrt(curv / m), rel=1e-4)
    # Equal-frequency depth: omega_lattice(U') = omega_tweezer(U).
    w0 = 0.866e-6
    up = B.depth_for_equal_omega(u, w0, a)
    assert B.lattice_omega(up, a) == pytest.approx(B.tweezer_omega_radial(u, w0), rel=1e-9)


def test_supergauss_helpers():
    n, h = 20, 2000.0
    w = LA.supergauss_w_for_flat(h, n, 0.01)
    assert np.abs(LA.supergauss(h, w, n)) ** 2 == pytest.approx(0.99, rel=1e-9)
    x = np.linspace(-3 * w, 3 * w, 200001)
    assert np.trapezoid(np.abs(LA.supergauss(x, w, n)) ** 2, x) == pytest.approx(LA.supergauss_equivalent_width(w, n), rel=1e-6)


def test_foxli_gaussian_mode_of_plane_concave_strip():
    """Spherical graded mirror: the fundamental is the analytic Gaussian with waist at the crossing."""
    lam, p, w0 = 1.075, 50e3, 400.0
    g = FL.Grid1D(2048, 8000.0, lam)
    ph = FL.spherical_phase(g, FL.gaussian_mode_radius(w0, p, lam))
    ring = FL.RingStrip(g, p, 0.5 * p, ph, ph, ap_g_um=3000.0, ap_fold_um=3000.0)
    hg0 = np.exp(-(g.x / w0) ** 2).astype(complex)
    hg1 = g.x * hg0

    def rayleigh(u):
        v = ring.roundtrip(u)
        lam_u = np.vdot(u, v) / np.vdot(u, u)
        fid = abs(np.vdot(u, v)) ** 2 / (np.vdot(u, u).real * np.vdot(v, v).real)
        return lam_u, fid

    l0, f0 = rayleigh(hg0)
    l1, f1 = rayleigh(hg1)
    # Both analytic strip modes reproduce themselves, without loss.
    assert f0 > 1 - 1e-8 and f1 > 1 - 1e-8 and abs(l0) ** 2 > 1 - 1e-8
    # Round-trip Gouy phase of TEM1 relative to TEM0 (1D strip): 2 x 2 atan(p / zR).
    zr = np.pi * w0**2 / lam
    assert abs(np.angle(l1 / l0)) == pytest.approx(4 * np.arctan(p / zr), rel=1e-4)


def test_graded_phase_flat_top_is_eigenmode():
    lam, p = 1.075, 100e3
    g = FL.Grid1D(2048, 12000.0, lam)
    F = FL.flat_target(g, 1500.0, 20)
    ph, _ = FL.graded_phase(g, F, p)
    ring = FL.RingStrip(g, p, 0.5 * p, ph, ph, ap_g_um=4000.0, ap_fold_um=4000.0)
    u = ring.roundtrip(F)
    fid = abs(np.vdot(F, u)) ** 2 / (np.vdot(F, F).real * np.vdot(u, u).real)
    assert fid > 1 - 1e-7
    # Driven, mode-matched: circulating gain equals the analytic enhancement.
    res = FL.lock_and_solve(ring, F, 1e-2, 1e-3, -np.angle(np.vdot(F, u)))
    assert res["enhancement"] == pytest.approx(B.enhancement(1e-2, 1e-3), rel=1e-4)


def test_telescope_design_is_eigenmode():
    lam = 1.075
    g = FL.Grid1D(4096, 16000.0, lam)
    F = FL.flat_target(g, 2000.0, 20)
    els = FL.telescope_elements(g, 80e3, 16e3, 40e3, 5000.0, 2500.0)
    p = 40e3 + 80e3 + 16e3 + 60e3
    ring = FL.RingStrip(g, p, 0.0, np.zeros(g.n), np.zeros(g.n), ap_g_um=2500.0, elements=els)
    ph, _ = FL.graded_phase(g, F, p, ring=ring)
    ring.phase_g1 = ring.phase_g2 = ph
    u = ring.roundtrip(F)
    assert abs(np.vdot(F, u)) ** 2 / (np.vdot(F, F).real * np.vdot(u, u).real) > 1 - 1e-6


def test_separable_propagation_matches_2d():
    """Low-angle separable field: 2D angular-spectrum propagation equals the product of 1D ones."""
    from metatrap.propagation.angular_spectrum import propagate
    from metatrap.propagation.layers import LayerStack

    lam, n, W, z = 1.075, 256, 2000.0, 20e3
    g = FL.Grid1D(n, W, lam)
    fx = FL.flat_target(g, 500.0, 10)
    fy = np.exp(-(g.x / 300.0) ** 2).astype(complex)
    u2 = propagate(np.outer(fy, fx).astype(np.complex64), W / n, lam, LayerStack.from_list([[z, "air"]]), remove_piston=False)
    k = 2 * np.pi / lam
    u1 = np.outer(g.prop(fy, z), g.prop(fx, z)) * np.exp(1j * k * z)
    assert np.max(np.abs(u2 - u1)) / np.max(np.abs(u1)) < 2e-3


def test_site_map_displacement_formula():
    lam, a, side = 1.075, 4.0, 400.0
    th = B.crossing_half_angle(a, lam)
    g = FL.Grid1D(4096, 4000.0, lam)
    flat = np.ones(g.n, dtype=complex)
    m0 = LA.site_map(flat, flat, g.x, lam, th, side, 21)
    assert np.allclose(m0["depth"], 1.0, atol=1e-6) and np.allclose(m0["dy"], 0.0, atol=1e-6)
    # Small tilt alpha of beam A: phase k alpha u_A, site shift -a dphi / (2 pi) relative to the centre.
    k, alpha = 2 * np.pi / lam, 2e-5
    m = LA.site_map(np.exp(1j * k * alpha * g.x), flat, g.x, lam, th, side, 21)
    X, Y = m["x"], m["y"]
    u_a = -X * np.sin(th) + Y * np.cos(th)
    expect = -a * (k * alpha * u_a) / (2 * np.pi)
    assert np.max(np.abs(m["dy"] - expect)) < 2e-3
