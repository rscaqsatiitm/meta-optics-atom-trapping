"""Plane-parallel optical path between the metasurface and the atom plane.

The path is a stack of homogeneous layers, e.g. air gap | glass-cell wall | vacuum.
A ray leaving the metasurface with direction cosine s (relative to the optical axis,
measured in air) obeys n_i sin(theta_i) = s in every layer, so

    r(s)   = sum_i d_i s / sqrt(n_i^2 - s^2)            (radial ray height at the metasurface)
    OPL(s) = sum_i n_i^2 d_i / sqrt(n_i^2 - s^2)         (optical path to the on-axis focus)

Here, d_i and n_i are the thickness and index of layer i; the first layer starts at
the metasurface and the last ends at the atom plane. The phase of an ideal lens that
focuses to the on-axis point through the stack is phi_lens(r) = -k0 OPL(s(r)); this
satisfies d phi / dr = -k0 s (generalised Snell law), see tests/test_propagation.py.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.interpolate import PchipInterpolator

from .. import materials


@dataclass(frozen=True)
class Layer:
    thickness_um: float
    material: str

    def n(self, lam_um: float) -> float:
        return float(materials.get(self.material).n(lam_um))


@dataclass(frozen=True)
class LayerStack:
    layers: tuple[Layer, ...] = field(default_factory=tuple)

    @classmethod
    def from_list(cls, items) -> "LayerStack":
        return cls(tuple(Layer(float(d), str(m)) for d, m in items))

    def indices(self, lam_um: float) -> np.ndarray:
        return np.array([l.n(lam_um) for l in self.layers])

    def thicknesses(self) -> np.ndarray:
        return np.array([l.thickness_um for l in self.layers])

    @property
    def total_um(self) -> float:
        return float(self.thicknesses().sum())

    def ray_height(self, s, lam_um: float):
        s = np.asarray(s, dtype=float)
        n, d = self.indices(lam_um), self.thicknesses()
        return sum(di * s / np.sqrt(ni**2 - s**2) for di, ni in zip(d, n))

    def opl(self, s, lam_um: float):
        s = np.asarray(s, dtype=float)
        n, d = self.indices(lam_um), self.thicknesses()
        return sum(ni**2 * di / np.sqrt(ni**2 - s**2) for di, ni in zip(d, n))

    def s_of_r(self, lam_um: float, s_max: float = 0.95, n_samples: int = 20001):
        """Monotone interpolant r -> s (the ray mapping is strictly increasing)."""
        s = np.linspace(0.0, s_max, n_samples)
        r = self.ray_height(s, lam_um)
        return PchipInterpolator(r, s, extrapolate=False)

    def na_for_aperture(self, radius_um: float, lam_um: float) -> float:
        return float(self.s_of_r(lam_um)(radius_um))

    def aperture_for_na(self, na: float, lam_um: float) -> float:
        """Metasurface radius (um) whose marginal ray reaches the focus with sin(theta) = na."""
        return float(self.ray_height(na, lam_um))

    def lens_phase(self, r_um, lam_um: float):
        """Ideal focusing phase (rad), referenced to zero on axis."""
        k0 = 2 * np.pi / lam_um
        s = self.s_of_r(lam_um)(np.asarray(r_um, dtype=float))
        return -k0 * (self.opl(s, lam_um) - self.opl(0.0, lam_um))

    def pupil_jacobian(self, s, lam_um: float):
        """Area ratio dA_metasurface / dA_s, with dA_s = s ds dphi (direction-cosine area)."""
        s = np.asarray(s, dtype=float)
        n, d = self.indices(lam_um), self.thicknesses()
        r = self.ray_height(s, lam_um)
        drds = sum(di * ni**2 / (ni**2 - s**2) ** 1.5 for di, ni in zip(d, n))
        with np.errstate(invalid="ignore", divide="ignore"):
            jac = np.where(s > 0, r * drds / np.where(s > 0, s, 1.0), drds**2)
        return jac

    def paraxial_focal_shift_um(self, lam_um: float) -> float:
        """Longitudinal focal shift that a plate adds to an uncorrected air-designed lens.

        For each non-air layer, d (1 - 1/n) in the paraxial limit.
        """
        n, d = self.indices(lam_um), self.thicknesses()
        return float(np.sum(d * (1.0 - 1.0 / n)))
