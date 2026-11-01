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


def _invert_height(st: LayerStack, dh, lam_um: float):
    """Transverse direction cosine p with st.ray_height(p) = dh (vectorised Newton, |p| < n_min)."""
    dh = np.asarray(dh, float)
    n_min = float(st.indices(lam_um).min())
    d, n = st.thicknesses(), st.indices(lam_um)
    p = dh / max(st.total_um / n_min, 1e-9) * n_min * 0.5
    for _ in range(60):
        h = sum(di * p / np.sqrt(ni**2 - p**2) for di, ni in zip(d, n))
        dhdp = sum(di * ni**2 / (ni**2 - p**2) ** 1.5 for di, ni in zip(d, n))
        step = (h - dh) / dhdp
        p = np.clip(p - step, -0.999 * n_min, 0.999 * n_min)
        if np.max(np.abs(step)) < 1e-15:
            break
    return p


def design(lam_um: float, na: float, cell: list, substrate: str = "fused_silica", t_um: float = 2000.0, n_pts: int = 4001,
           spacer: list | None = None, R1_um: float | None = None):
    """Return (meta1, spacer Media, meta2, cell Media, info) for the aplanatic doublet.

    ``spacer`` is the list [[thickness_um, material], ...] between the two metasurfaces
    (default: one layer of ``substrate`` of thickness ``t_um``). The internal ray of the zone
    h1 has transverse direction cosine p with h2 - h1 = sum_i d_i p / sqrt(n_i^2 - p^2), so
    that dphi1/dr (h1) = k0 p and dphi2/dr (h2) = -k0 (s_f + p). ``R1_um`` fixes the stop radius
    (default: R1 = R2, the aperture radius of an ideal singlet at the second surface, which
    keeps the corrector weak; for a multi-wavelength design pass a common value).
    """
    st = LayerStack.from_list(cell)
    sp = LayerStack.from_list(spacer if spacer is not None else [[t_um, substrate]])
    R2 = st.aperture_for_na(na, lam_um)  # height at the second surface of the marginal ray
    k0 = 2 * np.pi / lam_um
    R1 = R2 if R1_um is None else float(R1_um)
    f = R1 / na
    # The construction is extended to 1.1 R1 so that off-axis bundles, which reach the second
    # surface slightly beyond R2, still see the designed phase (the stop stays at R1).
    h1 = np.linspace(0, 1.1 * R1, n_pts)
    s_f = h1 / f
    h2 = st.ray_height(s_f, lam_um)
    if np.any(np.diff(h2) <= 0):
        raise ValueError("ray heights at the second surface are not monotonic")
    p = _invert_height(sp, h2 - h1, lam_um)
    g1 = k0 * p
    g2 = -k0 * (s_f + p)
    phi1 = np.concatenate([[0.0], np.cumsum(0.5 * (g1[1:] + g1[:-1]) * np.diff(h1))])
    phi2 = np.concatenate([[0.0], np.cumsum(0.5 * (g2[1:] + g2[:-1]) * np.diff(h2))])
    m1 = SplineMeta(h1, phi1, lam_um, "corrector")
    m2 = SplineMeta(h2, phi2, lam_um, "focusing")
    sp_media = [Medium(float(materials.get(m).n(lam_um)), float(t), m) for t, m in (spacer if spacer is not None else [[t_um, substrate]])]
    media = [Medium(float(materials.get(m).n(lam_um)), float(t), m) for t, m in cell]
    i_R = np.searchsorted(h1, R1)
    info = {"R1_um": float(R1), "R2_um": float(R2), "f_um": float(f),
            "max_deflection_corrector": float(np.max(np.abs(p[: i_R + 1]))),
            "max_deflection_focusing": float(np.max(np.abs(s_f + p)[: i_R + 1])),
            "max_deflection_focusing_ext": float(np.max(np.abs(s_f + p)))}
    sub = sp_media[0] if len(sp_media) == 1 else sp_media
    return m1, sub, m2, media, info
