"""Steady-state axisymmetric heat conduction in a metasurface substrate (finite volume).

Geometry: a disc of radius R and thickness t. The metasurface film lies on the face
z = 0 and absorbs a fraction ``a_film`` of the incident power as a surface source with
the beam's radial profile; the substrate absorbs ``alpha_bulk`` (1/m) per unit length
along z (weak absorption, uniform in z). The faces lose heat to ambient through a
linear coefficient h (W/m^2/K, natural convection plus linearised radiation); the rim
is either held at ambient (heat-sunk mount, default) or convective.

Here, T is the temperature rise above ambient in K and kappa the substrate thermal
conductivity (W/m/K). The film is thin (~1.5 um), so its own conduction is neglected.

Thin-plate check (tests/test_thermal.py): for t << R and h -> 0, a uniform disc source
of radius a in a plate with its rim at ambient has the analytic centre rise
    dT0 = P / (4 pi kappa t) [1 + 2 ln(R / a)].
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla


@dataclass
class Substrate:
    name: str
    kappa: float  # W/m/K
    dndT: float  # 1/K at the design wavelength
    alpha_exp: float  # linear expansion 1/K
    poisson: float
    n: float  # refractive index at the design wavelength
    youngs_GPa: float = 0.0
    strength_MPa: float = 0.0  # characteristic (mean) strength used for the stress margin
    source: str = ""

    @property
    def chi(self) -> float:
        """Thermo-optic path coefficient d(OPL)/dT per unit thickness, free plate (plane stress):
        dn/dT + (n - 1)(1 + nu) alpha. Photo-elastic terms are neglected."""
        return self.dndT + (self.n - 1.0) * (1.0 + self.poisson) * self.alpha_exp


# Values with provenance; see docs/model_notes.md.
FUSED_SILICA = Substrate(
    name="fused_silica", kappa=1.38, dndT=9.6e-6, alpha_exp=0.52e-6, poisson=0.16, n=1.4495,
    youngs_GPa=72.7, strength_MPa=52.4,
    source=(
        "kappa 1.38 W/m/K (Corning HPFS 7980 product information as quoted in search; the 2003 sheet "
        "read this session gives 1.30); dn/dT 9.6 ppm/K at 1064 nm, alpha 0.52 ppm/K (5-35 C), E 72.7 GPa, "
        "nu 0.16, modulus of rupture (abraded) 52.4 MPa: Corning 7980 sheet (uqgoptics copy, read 2026-10-31)"
    ),
)
SAPPHIRE = Substrate(
    name="sapphire", kappa=33.0, dndT=13.0e-6, alpha_exp=5.0e-6, poisson=0.25, n=1.7543,
    youngs_GPa=335.0, strength_MPa=400.0,
    source=(
        "kappa 33.0 W/m/K perpendicular to c at 300 K, alpha 5.0 ppm/K perpendicular, E 335 GPa, nu 0.25: "
        "Alkor sapphire data sheet (read 2026-10-31); dn/dT [VALUE?: approx. 13 ppm/K, quoted at 0.546 um, "
        "verify at 1075 nm]; strength [VALUE?: approx. 400 MPa flexural, not retrieved, verify]"
    ),
)
SUBSTRATES = {s.name: s for s in (FUSED_SILICA, SAPPHIRE)}


@dataclass
class ThermalProblem:
    radius_m: float
    thickness_m: float
    substrate: Substrate
    beam: str = "top_hat"  # 'top_hat' or 'gaussian'
    beam_radius_m: float = 10e-3  # aperture radius (top hat) or clip radius (gaussian)
    beam_waist_m: float | None = None  # 1/e^2 intensity radius for 'gaussian'
    h_faces: float = 10.0
    rim: str = "fixed"  # 'fixed' or 'convective'
    h_rim: float = 10.0
    nr: int = 240
    nz: int = 40
    z_power: float = 1.0  # z edges at t (x^p); p > 1 refines cells near the film face z = 0

    def intensity_profile(self, r):
        """Incident intensity per unit power (1/m^2), normalised over the aperture."""
        r = np.asarray(r, dtype=float)
        a = self.beam_radius_m
        if self.beam == "top_hat":
            prof = (r <= a).astype(float)
        else:
            w = self.beam_waist_m
            prof = np.exp(-2 * r**2 / w**2) * (r <= a)
        rr = np.linspace(0, a, 20001)
        pr = (rr <= a) if self.beam == "top_hat" else np.exp(-2 * rr**2 / self.beam_waist_m**2)
        norm = np.trapezoid(pr * 2 * np.pi * rr, rr)
        return prof / norm


def _edges(R, n, a, refine=0.6):
    """Radial cell edges, uniform inside the beam (finer) and stretched outside."""
    n_in = int(n * refine)
    e_in = np.linspace(0, a, n_in + 1)
    if R <= a * 1.0001:
        return np.linspace(0, R, n + 1)
    n_out = n - n_in
    x = np.linspace(0, 1, n_out + 1)[1:]
    e_out = a + (R - a) * x**1.5
    return np.concatenate([e_in, e_out])


def solve(prob: ThermalProblem, p_in_W: float, a_film: float, alpha_bulk_per_m: float = 0.0):
    """Return (r_centres, z_centres, T[nr, nz]) in metres and kelvin."""
    re = _edges(prob.radius_m, prob.nr, min(prob.beam_radius_m, prob.radius_m))
    ze = prob.thickness_m * np.linspace(0, 1, prob.nz + 1) ** prob.z_power
    rc, zc = 0.5 * (re[1:] + re[:-1]), 0.5 * (ze[1:] + ze[:-1])
    dr, dz = np.diff(re), np.diff(ze)
    nr, nz = len(rc), len(zc)
    k = prob.substrate.kappa
    area = np.pi * (re[1:] ** 2 - re[:-1] ** 2)  # annulus areas
    idx = lambda i, j: i * nz + j  # noqa: E731
    rows, cols, vals = [], [], []
    b = np.zeros(nr * nz)
    diag = np.zeros(nr * nz)

    def add(i1, j1, i2, j2, g):
        p, q = idx(i1, j1), idx(i2, j2)
        diag[p] += g
        diag[q] += g
        rows.extend([p, q])
        cols.extend([q, p])
        vals.extend([-g, -g])

    for i in range(nr):
        for j in range(nz):
            if i + 1 < nr:
                g = k * 2 * np.pi * re[i + 1] * dz[j] / (rc[i + 1] - rc[i])
                add(i, j, i + 1, j, g)
            if j + 1 < nz:
                g = k * area[i] / (zc[j + 1] - zc[j])
                add(i, j, i, j + 1, g)
    # Faces: series of half-cell conduction and surface film coefficient.
    for i in range(nr):
        for j, half in ((0, dz[0] / 2), (nz - 1, dz[-1] / 2)):
            g = area[i] / (1.0 / prob.h_faces + half / k) if prob.h_faces > 0 else 0.0
            diag[idx(i, j)] += g
    # Rim.
    for j in range(nz):
        a_rim = 2 * np.pi * re[-1] * dz[j]
        half = dr[-1] / 2
        g = k * a_rim / half if prob.rim == "fixed" else a_rim / (1.0 / prob.h_rim + half / k)
        diag[idx(nr - 1, j)] += g
    # Sources: film on z = 0 and bulk absorption.
    prof = prob.intensity_profile(rc)
    for i in range(nr):
        b[idx(i, 0)] += a_film * p_in_W * prof[i] * area[i]
        if alpha_bulk_per_m > 0:
            for j in range(nz):
                b[idx(i, j)] += alpha_bulk_per_m * p_in_W * prof[i] * area[i] * dz[j]
    A = sp.csr_matrix((vals, (rows, cols)), shape=(nr * nz, nr * nz)) + sp.diags(diag)
    T = spla.spsolve(A.tocsc(), b).reshape(nr, nz)
    return rc, zc, T


def opl_change(prob: ThermalProblem, rc, zc, T, film_dndT: float = 0.0, film_thickness_m: float = 0.0):
    """Optical path change (m) along z at each radius: chi int T dz (+ film term)."""
    dz = np.diff(prob.thickness_m * np.linspace(0, 1, len(zc) + 1) ** prob.z_power)
    d_opl = prob.substrate.chi * np.sum(T * dz[None, :], axis=1)
    if film_dndT and film_thickness_m:
        d_opl += film_dndT * film_thickness_m * T[:, 0]
    return d_opl


def surface_temperature(zc, T):
    """Film-face temperature by linear extrapolation of the first two cell centres to z = 0."""
    return T[:, 0] - zc[0] * (T[:, 1] - T[:, 0]) / (zc[1] - zc[0])
