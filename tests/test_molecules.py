"""Anchors for the two-wavelength meta-atom tools."""

import numpy as np

from metatrap.rcwa.molecules import Rules, candidates, dual_efficiency, masks


def test_masks_fill_fractions():
    m = masks("ring", [(0.4, 0.2), (0.3, 0.0)], 0.55, 128)
    assert abs(m[0].mean() * 0.55**2 - (0.4**2 - 0.2**2)) < 2e-3
    assert abs(m[1].mean() * 0.55**2 - 0.09) < 2e-3
    d = masks("duo", [(0.2, 0.2)], 0.55, 128)
    assert abs(d[0].mean() * 0.55**2 - 0.08) < 2e-3
    # Square symmetry (polarisation independence at normal incidence).
    assert np.allclose(d[0], d[0].T) and np.allclose(m[0], m[0][::-1])


def test_candidates_respect_rules():
    r = Rules(min_feature_um=0.1, max_aspect=10.0)
    for a, b in candidates("ring", 0.55, 1.6, r):
        assert a >= 0.16 - 1e-9 and a <= 0.45 + 1e-9
        if b > 0:
            assert b >= 0.1 - 1e-9 and (a - b) / 2 >= 0.16 - 1e-9


def test_dual_efficiency_limits():
    g = 2 * np.pi * np.arange(24) / 24
    A, B = np.meshgrid(g, g, indexing="ij")
    full = dual_efficiency(np.exp(1j * A.ravel()), np.exp(1j * B.ravel()))
    assert full[0] > 0.98 and full[1] > 0.98
    # A one-parameter family with tied phases (phi_b = 0.8 phi_a) cannot cover the torus.
    ph = np.linspace(0, 2 * np.pi, 200, endpoint=False)
    tied = dual_efficiency(np.exp(1j * ph), np.exp(1j * 0.8 * ph))
    assert tied[0] + tied[1] < 1.5
