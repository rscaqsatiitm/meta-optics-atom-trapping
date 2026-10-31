"""Deflection efficiency of blazed supercells built from library pillars.

A supercell of N pillars at pitch p (period Lambda = N p) with phases 2 pi k / N
deflects normally incident light into the first order at sin(theta) = lambda / Lambda
in air. This is the local behaviour of the meta-hologram at the radius where the lens
deflects by theta, so efficiency versus sin(theta) bounds the efficiency of the outer
aperture, where the local-periodic approximation is least accurate (NA 0.55 at the
rim corresponds to theta = 33 deg).
"""

from __future__ import annotations

import numpy as np

from .cell import PillarCell
from .lookup import Lookup
from .solver import solve_grating_orders


def supercell_eps(widths_um, pitch_um: float, eps_pillar: complex, eps_clad: complex, n_grid: int) -> np.ndarray:
    cell = PillarCell(pitch_um, 1.0)
    blocks = cell.permittivity(widths_um, eps_pillar, eps_clad, n_grid)  # (N, ng, ng)
    return np.concatenate(list(blocks), axis=0)  # (N ng, ng), x along the deflection direction


def deflection_efficiency(lookup: Lookup, n_cells: int, pitch_um: float, height_um: float, lam_um: float, eps_pillar: complex, eps_sub: complex, n_grid: int = 64, terms_per_cell: int = 120, polarisation: str = "x"):
    """First-order transmitted efficiency of an N-pillar blazed supercell (fmmax).

    Returns a dict with the RCWA efficiency, the LPA prediction eta_N of the same
    widths, the deflection sin(theta) and the reflected power.
    """
    phases = 2 * np.pi * np.arange(n_cells) / n_cells
    # Ascending phase along +x deflects towards +x for exp(+i phi) transmission.
    widths = lookup.width_for_phase(phases)
    t = lookup.t_for_phase(phases)
    eta_lpa = float(lookup.n_ratio * np.abs(np.mean(t * np.exp(-1j * phases))) ** 2)
    eps = supercell_eps(widths, pitch_um, eps_pillar, 1.0, n_grid)
    orders, eff, R = solve_grating_orders(eps, n_cells * pitch_um, pitch_um, height_um, lam_um, eps_sub, 1.0, n_terms=terms_per_cell * n_cells, polarisation=polarisation)
    e = {tuple(o): float(v) for o, v in zip(orders, eff)}
    first = max(e.get((1, 0), 0.0), e.get((-1, 0), 0.0))
    return {
        "n_cells": n_cells,
        "sin_theta": lam_um / (n_cells * pitch_um),
        "eta_rcwa": first,
        "eta_plus1": e.get((1, 0), 0.0),
        "eta_minus1": e.get((-1, 0), 0.0),
        "eta_zero": e.get((0, 0), 0.0),
        "eta_lpa": eta_lpa,
        "R": R,
        "polarisation": polarisation,
        "widths_um": widths,
    }
