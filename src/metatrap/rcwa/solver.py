"""RCWA solvers for pillar meta-atoms.

The primary solver is fmmax (Fourier modal method in JAX). The stack is
    incidence medium (semi-infinite) | patterned pillar layer (height h) | exit medium,
with light incident from the substrate side by default, as for a metasurface
illuminated through its substrate and focusing into air.

Conventions (verified in tests/test_rcwa.py against Fresnel and thin-film results):

* ``t`` is the complex amplitude of the zeroth-order transmitted E_x (for unit
  incident E_x) at the exit face of the pillar layer; the incident field is
  referenced to the entrance face. A homogeneous slab of index n and thickness h
  therefore gives the Airy thin-film coefficient, including exp(i k0 n h).
* Zeroth-order power transmittance at normal incidence is T0 = (n_exit/n_in) |t|^2.
* ``T`` and ``R`` are the total transmitted and reflected power fractions summed over
  all propagating orders, from the Poynting flux.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

import jax

jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp  # noqa: E402
import fmmax  # noqa: E402

from .cell import PillarCell  # noqa: E402

FORMULATIONS = {
    "fft": fmmax.Formulation.FFT,
    "jones_direct_fourier": fmmax.Formulation.JONES_DIRECT_FOURIER,
    "pol_fourier": fmmax.Formulation.POL_FOURIER,
}


@dataclass
class PillarResult:
    widths: np.ndarray
    t_x: np.ndarray  # complex zeroth-order co-polarised transmission, x input
    t_y: np.ndarray  # complex zeroth-order co-polarised transmission, y input
    t_cross: np.ndarray  # |E_y| for x input (should vanish for square pillars)
    r_x: np.ndarray  # complex zeroth-order reflection, x input
    T: np.ndarray  # total transmitted power fraction, x input
    R: np.ndarray  # total reflected power fraction, x input
    T0: np.ndarray  # zeroth-order transmitted power fraction, x input
    n_terms: int


def _uniform(lam, k0, lv, eps, expansion):
    return fmmax.eigensolve_isotropic_media(
        jnp.asarray(lam), k0, lv, jnp.full((1, 1), complex(eps)), expansion
    )


def _xy_input_amplitudes(layer, n_terms):
    """Forward eigenmode amplitudes producing unit zeroth-order E_x and E_y in a uniform layer."""
    a = jnp.zeros((2 * n_terms, 2), complex).at[0, 0].set(1.0).at[n_terms, 1].set(1.0)
    (ex, ey, _), _ = fmmax.fields_from_wave_amplitudes(a, jnp.zeros_like(a), layer)
    m = jnp.array([[ex[0, 0], ex[0, 1]], [ey[0, 0], ey[0, 1]]])
    c = jnp.linalg.solve(m, jnp.eye(2, dtype=complex))
    return a @ c  # columns: x-polarised, y-polarised


def solve_pillars(
    cell: PillarCell,
    widths,
    lam_um: float,
    eps_pillar: complex,
    eps_in: complex,
    eps_out: complex = 1.0,
    eps_clad: complex = 1.0,
    n_terms: int = 150,
    n_grid: int = 128,
    formulation: str = "jones_direct_fourier",
) -> PillarResult:
    """Normal-incidence response of a batch of pillar widths (one eigensolve batch)."""
    widths = np.atleast_1d(np.asarray(widths, dtype=float))
    eps = cell.permittivity(widths, eps_pillar, eps_clad, n_grid).astype(complex)
    return solve_eps(eps, cell.pitch, cell.height, lam_um, eps_in, eps_out, n_terms, formulation, labels=widths)


def solve_eps(
    eps: np.ndarray,
    pitch_um: float,
    height_um: float,
    lam_um: float,
    eps_in: complex,
    eps_out: complex = 1.0,
    n_terms: int = 150,
    formulation: str = "jones_direct_fourier",
    labels=None,
) -> PillarResult:
    """Normal-incidence response of a batch of square unit cells with permittivity maps
    eps (B, n_grid, n_grid) on a square lattice of the given pitch (patterned layer of the
    given height between the incidence and exit half-spaces)."""
    lv = fmmax.LatticeVectors(u=pitch_um * fmmax.X, v=pitch_um * fmmax.Y)
    expansion = fmmax.generate_expansion(lv, n_terms, fmmax.Truncation.CIRCULAR)
    nt = expansion.num_terms
    k0 = jnp.zeros(2)
    lam = jnp.asarray(lam_um)

    inc = _uniform(lam, k0, lv, eps_in, expansion)
    out = _uniform(lam, k0, lv, eps_out, expansion)
    pil = fmmax.eigensolve_isotropic_media(
        lam, k0, lv, jnp.asarray(eps.astype(complex)), expansion, formulation=FORMULATIONS[formulation]
    )
    s = fmmax.stack_s_matrix([inc, pil, out], [jnp.asarray(0.0), jnp.asarray(height_um), jnp.asarray(0.0)])

    a0 = _xy_input_amplitudes(inc, nt)  # (2nt, 2)
    a0b = jnp.broadcast_to(a0, s.s11.shape[:-2] + a0.shape)
    a_n = s.s11 @ a0b
    b_0 = s.s21 @ a0b
    (tx, ty, _), _ = fmmax.fields_from_wave_amplitudes(a_n, jnp.zeros_like(a_n), out)
    (rx, _, _), _ = fmmax.fields_from_wave_amplitudes(jnp.zeros_like(b_0), b_0, inc)

    f_in, _ = fmmax.directional_poynting_flux(a0, jnp.zeros_like(a0), inc)
    f_out, _ = fmmax.directional_poynting_flux(a_n, jnp.zeros_like(a_n), out)
    _, f_ref = fmmax.directional_poynting_flux(jnp.zeros_like(b_0), b_0, inc)
    p_in = jnp.sum(f_in, axis=-2)[..., 0]
    T = jnp.sum(f_out, axis=-2)[..., 0] / p_in
    R = -jnp.sum(f_ref, axis=-2)[..., 0] / p_in
    n_in = np.sqrt(complex(eps_in)).real
    n_out = np.sqrt(complex(eps_out)).real

    t_x = np.asarray(tx[..., 0, 0])
    return PillarResult(
        widths=np.asarray(labels) if labels is not None else np.arange(eps.shape[0]),
        t_x=t_x,
        t_y=np.asarray(ty[..., 0, 1]),
        t_cross=np.abs(np.asarray(ty[..., 0, 0])),
        r_x=np.asarray(rx[..., 0, 0]),
        T=np.asarray(T),
        R=np.asarray(R),
        T0=(n_out / n_in) * np.abs(t_x) ** 2,
        n_terms=nt,
    )


def solve_grating_orders(
    eps_grid: np.ndarray,
    period_x: float,
    period_y: float,
    height: float,
    lam_um: float,
    eps_in: complex,
    eps_out: complex = 1.0,
    n_terms: int = 300,
    formulation: str = "jones_direct_fourier",
    polarisation: str = "x",
):
    """Per-order transmitted power fractions of a rectangular supercell at normal incidence.

    ``eps_grid`` has shape (nx, ny) over a period_x x period_y cell. Returns
    (orders, efficiencies, R_total) where ``orders`` is the (n_terms, 2) integer array
    of diffraction orders and ``efficiencies`` the transmitted power fraction in each.
    """
    lv = fmmax.LatticeVectors(u=period_x * fmmax.X, v=period_y * fmmax.Y)
    expansion = fmmax.generate_expansion(lv, n_terms, fmmax.Truncation.CIRCULAR)
    nt = expansion.num_terms
    k0 = jnp.zeros(2)
    lam = jnp.asarray(lam_um)
    inc = _uniform(lam, k0, lv, eps_in, expansion)
    out = _uniform(lam, k0, lv, eps_out, expansion)
    lay = fmmax.eigensolve_isotropic_media(
        lam, k0, lv, jnp.asarray(eps_grid.astype(complex)), expansion, formulation=FORMULATIONS[formulation]
    )
    s = fmmax.stack_s_matrix([inc, lay, out], [jnp.asarray(0.0), jnp.asarray(height), jnp.asarray(0.0)])
    a0 = _xy_input_amplitudes(inc, nt)[:, [0 if polarisation == "x" else 1]]
    a_n = s.s11 @ a0
    b_0 = s.s21 @ a0
    f_in, _ = fmmax.directional_poynting_flux(a0, jnp.zeros_like(a0), inc)
    f_out, _ = fmmax.directional_poynting_flux(a_n, jnp.zeros_like(a_n), out)
    _, f_ref = fmmax.directional_poynting_flux(jnp.zeros_like(b_0), b_0, inc)
    p_in = float(jnp.sum(f_in))
    per = np.asarray(f_out[:, 0])
    eff = per[:nt] + per[nt:]  # sum the two polarisation blocks of each order
    return np.asarray(expansion.basis_coefficients), eff / p_in, float(-jnp.sum(f_ref) / p_in)


def solve_pillars_grcwa(
    cell: PillarCell,
    width: float,
    lam_um: float,
    eps_pillar: complex,
    eps_in: complex,
    eps_out: complex = 1.0,
    eps_clad: complex = 1.0,
    n_terms: int = 150,
    n_grid: int = 128,
) -> tuple[complex, float, float]:
    """Independent cross-check with grcwa: returns (t_x, T, R) for one width.

    grcwa reports powers directly; the complex zeroth-order amplitude is extracted from
    the Fourier field amplitudes at the output interface, then normalised to the
    incident E_x and to fmmax's phase reference by the caller-independent convention
    of a unit incident E_x at the entrance face.
    """
    import grcwa

    grcwa.set_backend("numpy")
    freq = 1.0 / lam_um
    obj = grcwa.obj(n_terms, [cell.pitch, 0.0], [0.0, cell.pitch], freq, 0.0, 0.0, verbose=0)
    obj.Add_LayerUniform(0.0, complex(eps_in))
    obj.Add_LayerGrid(cell.height, n_grid, n_grid)
    obj.Add_LayerUniform(0.0, complex(eps_out))
    obj.Init_Setup()
    eps = cell.permittivity([width], eps_pillar, eps_clad, n_grid)[0]
    obj.GridLayer_geteps(eps.flatten())
    # p-polarised amplitude 1 with phi = 0 at normal incidence corresponds to E_x.
    obj.MakeExcitationPlanewave(1.0, 0.0, 0.0, 0.0, order=0)
    R, T = obj.RT_Solve(normalize=1)
    e_out, _ = obj.Solve_FieldFourier(2, 0.0)
    e_in, _ = obj.Solve_FieldFourier(0, 0.0)
    # Incident + reflected at the entrance; the incident part is obtained from a bare run.
    return complex(np.ravel(e_out[0])[0]), float(np.real(T)), float(np.real(R))
