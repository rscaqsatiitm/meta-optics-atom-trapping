"""Analytic aplanatic metasurface doublet (stigmatic on axis and satisfying the sine condition).

Two flat metasurfaces on the faces of a substrate (thickness t, index n_s), followed by the
cell stack (air gap, wall, vacuum) to the atom plane. For a collimated input, a ray entering at
height h1 is required to leave the second surface with transverse direction cosine
s_f = h1 / f (Abbe sine condition, f = R / NA) and to reach the on-axis focus through the
stack, which fixes its height at the second surface: h2 = H(s_f) = sum_i d_i s_f /
sqrt(n_i^2 - s_f^2). The internal ray direction follows from h2 = h1 + t tan(a). The phase
gradients follow from the generalised Snell law at each surface:

    dphi1/dr (h1) = k0 n_s sin(a),            dphi2/dr (h2) = -k0 s_f - k0 n_s sin(a),

and are integrated to phase profiles (spline representation). A system with zero spherical
aberration and zero linear coma results; astigmatism and field curvature remain (quadratic in
the field angle). Here, k0 = 2 pi / lambda and all lengths are in um.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.interpolate import CubicSpline

from .. import materials
from ..propagation.layers import LayerStack
from .raytrace import Medium


@dataclass
class SplineMeta:
    """Radial metasurface phase given as a spline phi(r) (rad) at the design wavelength."""

    r: np.ndarray
    phi: np.ndarray
    lam_design_um: float
    name: str = ""
    disp: object = None

    def __post_init__(self):
        self._cs = CubicSpline(self.r, self.phi)

    def scale(self, lam_um):
        return 1.0 if self.disp is None else float(self.disp(lam_um))

    def phase(self, r, lam_um):
        # Beyond the design radius: linear continuation consistent with the edge gradient.
        r = np.asarray(r, float)
        rm = self.r[-1]
        inside = self._cs(np.clip(r, 0, rm))
        return self.scale(lam_um) * np.where(r <= rm, inside, self._cs(rm) + self._cs(rm, 1) * (r - rm))

    def dphase_dr(self, r, lam_um):
        r = np.asarray(r, float)
        return self.scale(lam_um) * self._cs(np.clip(r, 0, self.r[-1]), 1)


def design(lam_um: float, na: float, cell: list, substrate: str = "fused_silica", t_um: float = 2000.0, n_pts: int = 4001):
    """Return (meta1, substrate Medium, meta2, cell Media, info) for the aplanatic doublet."""
    st = LayerStack.from_list(cell)
    R2 = st.aperture_for_na(na, lam_um)  # height at the second surface of the marginal ray
    k0 = 2 * np.pi / lam_um
    ns = float(materials.get(substrate).n(lam_um))
    # Choose the stop radius R1 so that the marginal ray needs no extra bend: the sine condition
    # gives f = R1 / NA; R1 is free. Taking R1 = R2 keeps the corrector weak at the edge.
    R1 = R2
    f = R1 / na
    # The construction is extended to 1.1 R1 so that off-axis bundles, which reach the second
    # surface slightly beyond R2, still see the designed phase (the stop stays at R1).
    h1 = np.linspace(0, 1.1 * R1, n_pts)
    s_f = h1 / f
    h2 = st.ray_height(s_f, lam_um)
    tan_a = (h2 - h1) / t_um
    sin_a = tan_a / np.sqrt(1 + tan_a**2)
    g1 = k0 * ns * sin_a
    g2 = -k0 * s_f - k0 * ns * sin_a
    phi1 = np.concatenate([[0.0], np.cumsum(0.5 * (g1[1:] + g1[:-1]) * np.diff(h1))])
    if np.any(np.diff(h2) <= 0):
        raise ValueError("ray heights at the second surface are not monotonic; increase t")
    phi2 = np.concatenate([[0.0], np.cumsum(0.5 * (g2[1:] + g2[:-1]) * np.diff(h2))])
    m1 = SplineMeta(h1, phi1, lam_um, "corrector")
    m2 = SplineMeta(h2, phi2, lam_um, "focusing")
    sub = Medium(ns, t_um, substrate)
    media = [Medium(float(materials.get(m).n(lam_um)), float(t), m) for t, m in cell]
    info = {"R1_um": float(R1), "R2_um": float(R2), "f_um": float(f), "max_internal_angle_deg": float(np.degrees(np.max(np.abs(np.arctan(tan_a))))),
            "max_deflection_corrector": float(np.max(np.abs(ns * sin_a))), "max_deflection_focusing": float(np.max(np.abs(s_f + ns * sin_a)))}
    return m1, sub, m2, media, info
