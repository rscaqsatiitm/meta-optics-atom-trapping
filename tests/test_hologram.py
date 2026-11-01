"""Anchors for the pupil-domain hologram design."""

import numpy as np
import pytest

from metatrap.hologram import exact as X
from metatrap.hologram import targets
from metatrap.hologram import wgs as W
from metatrap.hologram.field_model import spot_aberrations
from metatrap.propagation.layers import LayerStack

LAM = 1.075
STACK = LayerStack.from_list([(2500, "air"), (4000, "fused_silica"), (9000, "air")])


def test_wgs_small_array_converges():
    sites = targets.square_grid(10, 10, 4.5)
    g = W.PupilGrid.for_sites(LAM, 4.5, 8, 0.5, 45.0)
    a = W.pupil_amplitude(g, STACK, 0.5, {"kind": "top_hat"})
    r = W.wgs(a, g, sites, n_iter=60, n_fix=30)
    assert r.history_nonuniformity[-1] < 0.01


def test_exact_kernel_equals_debye_on_axis_and_lens_efficiency_one():
    g = W.PupilGrid.for_sites(LAM, 4.5, 4, 0.5, 20.0)
    ps = X.PupilSampling.build(g, STACK, 0.5, {"kind": "top_hat"})
    tab = X.OPLTable(STACK, LAM, STACK.aperture_for_na(0.5, LAM) + 50)
    site = np.zeros((1, 2))
    K = X.kernels(ps, site, tab, LAM)
    assert np.allclose(K, 1.0, atol=1e-4)  # the on-axis kernel is the lens itself
    e = X.mode_amplitudes(ps, np.zeros(len(ps.amp)), K)
    assert abs(e[0]) ** 2 == pytest.approx(1.0, abs=1e-6)


def test_exact_design_beats_debye_off_axis():
    sites = targets.square_grid(16, 16, 4.5, center=(150.0, 0.0))
    g = W.PupilGrid.for_sites(LAM, 4.5, 4, 0.55, 2 * 190.0, guard=1.2)
    ps = X.PupilSampling.build(g, STACK, 0.55, {"kind": "top_hat"})
    tab = X.OPLTable(STACK, LAM, STACK.aperture_for_na(0.55, LAM) + 500)
    K = X.kernels(ps, sites, tab, LAM)
    Kd = X.debye_kernels(ps, sites, LAM)
    rd = X.wgs_exact(ps, Kd, n_iter=40, n_fix=20)
    re = X.wgs_exact(ps, K, n_iter=40, n_fix=20)
    ed = np.abs(X.mode_amplitudes(ps, rd.psi, K)) ** 2
    ee = np.abs(X.mode_amplitudes(ps, re.psi, K)) ** 2
    assert ee.std() / ee.mean() < 0.01
    assert ed.std() / ed.mean() > ee.std() / ee.mean()
    assert ee.sum() > ed.sum()


def test_aberration_model_on_axis_and_symmetry():
    a = spot_aberrations(STACK, LAM, 0.55, np.array([[0.0, 0.0], [100.0, 0.0], [0.0, 100.0]]), {"kind": "top_hat"}, n_pupil=96)
    assert a["strehl"][0] == pytest.approx(1.0, abs=1e-9)
    assert a["strehl"][1] == pytest.approx(a["strehl"][2], rel=1e-3)
    assert a["dz_um"][1] < 0  # field curvature towards the metasurface


def test_symmetric_design_matches_full_kernels():
    """D4-reduced WGS: unfolding the wedge phase and evaluating with full kernels reproduces E_j."""
    from metatrap.hologram import symmetric as S

    sites = targets.square_grid(12, 12, 4.5)
    g = W.PupilGrid.for_sites(LAM, 4.5, 4, 0.55, 2 * np.abs(sites).max() + 4.5, guard=1.8)
    tab = X.OPLTable(STACK, LAM, STACK.aperture_for_na(0.55, LAM) + 200)
    fs, fstab, inv = S.wedge_sites(sites)
    assert len(fs) == 6 * 7 // 2 and fstab.sum() == len(fs) + 6
    wp = S.WedgePupil.build(g, STACK, 0.55, {"kind": "top_hat"})
    KG = S.kernels_G(wp, fs, tab, LAM)
    r = S.wgs_symmetric(wp, KG, fstab, n_iter=40, n_fix=20)
    ps = X.PupilSampling.build(g, STACK, 0.55, {"kind": "top_hat"})
    full_psi = S.unfold(g.n, wp.ij, r.psi.astype(np.float64), dtype=np.float64)
    psi_samples = full_psi[ps.mask]
    K = X.kernels(ps, sites, tab, LAM)
    e_full = X.mode_amplitudes(ps, psi_samples, K)
    assert np.allclose(np.abs(e_full), np.abs(r.mode_amp[inv]), rtol=2e-3, atol=1e-6)
    ee = np.abs(e_full) ** 2
    assert ee.std() / ee.mean() < 0.01
