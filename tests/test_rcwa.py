"""RCWA wrapper conventions (rcwa/solver.py) against Fresnel and thin-film results."""

import numpy as np

from metatrap.rcwa.cell import PillarCell
from metatrap.rcwa.solver import solve_eps, solve_pillars


def airy_t(n1, n2, n3, h, lam):
    """Thin-film transmission coefficient of a slab n2 (thickness h) between n1 and n3."""
    t12, t23 = 2 * n1 / (n1 + n2), 2 * n2 / (n2 + n3)
    r12, r23 = (n1 - n2) / (n1 + n2), (n2 - n3) / (n2 + n3)
    d = 2 * np.pi / lam * n2 * h
    return t12 * t23 * np.exp(1j * d) / (1 + r12 * r23 * np.exp(2j * d))


def test_homogeneous_slab_is_airy_thin_film():
    n1, n2, n3, h, lam = 1.45, 2.0, 1.0, 0.9, 1.075
    for n_layer in (n2, n3):  # dense slab, and an "empty" layer (pure Fresnel interface plus propagation)
        eps = np.full((1, 32, 32), n_layer**2, complex)
        r = solve_eps(eps, 0.5, h, lam, n1**2, n3**2, n_terms=9)
        t = airy_t(n1, n_layer, n3, h, lam)
        assert abs(r.t_x[0] - t) < 1e-6
        assert abs(r.T0[0] - (n3 / n1) * abs(t) ** 2) < 1e-6
        assert abs(r.T[0] + r.R[0] - 1) < 1e-6  # lossless


def test_square_pillar_polarisation_independent():
    cell = PillarCell(pitch=0.5, height=1.5)
    r = solve_pillars(cell, [0.2, 0.3], 1.075, 2.3**2, 1.45**2, n_terms=49, n_grid=64)
    assert np.allclose(r.t_x, r.t_y, atol=2e-4)
    # Cross-polarisation of a C4v pillar vanishes by symmetry. pol_fourier reproduces this;
    # the default jones_direct_fourier formulation leaves a spurious t_cross of ~1e-2 here
    # (0.038 / 0.003 for widths 0.2 / 0.3 um at 49 terms; see the solver docstring).
    rp = solve_pillars(cell, [0.2, 0.3], 1.075, 2.3**2, 1.45**2, n_terms=49, n_grid=64, formulation="pol_fourier")
    assert np.all(rp.t_cross < 1e-3)
