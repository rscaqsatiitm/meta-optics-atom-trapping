"""Field-dependent aberration of a Debye-designed direct-projection hologram.

The pupil-domain (Debye) design gives spot j the metasurface phase
    phi_j(r) = -k0 [OPL0(r) - s(r) r_hat . x_j],
the first-order expansion in x_j of the exact focusing phase -k0 OPL(r -> x_j).
Here, OPL0(r) is the optical path from metasurface point r to the on-axis focus,
s(r) the direction cosine of that ray, and x_j the site position in the atom plane.
The residual wavefront of spot j is therefore

    W_j(r) = OPL(r -> x_j) - OPL0(r) + s(r) r_hat . x_j,

computed exactly through the layer stack (ray traced with a Newton-refined inverse of
the ray mapping). Tilt and defocus are fitted and removed by weighted least squares:
tilt is reported as a lateral position error, and defocus as an axial focal shift
(dz = -c / (dW/dz coefficient), see code). The remaining wavefront gives the peak
Strehl ratio S_j = |sum A e^{i k W}|^2 / (sum A)^2 on the pupil.

This is a scalar, ray-optical model; tests/test_hologram.py compares it against the
full-pixel angular-spectrum propagation of a 5 mm hologram.
"""

from __future__ import annotations

import numpy as np

from ..propagation.illumination import profile
from ..propagation.layers import LayerStack


class ExactOPL:
    """Optical path from a metasurface point to an atom-plane point at transverse distance rho."""

    def __init__(self, stack: LayerStack, lam_um: float, s_max: float = 0.9):
        self.stack, self.lam = stack, lam_um
        self._inv = stack.s_of_r(lam_um, s_max=s_max, n_samples=40001)
        n, d = stack.indices(lam_um), stack.thicknesses()
        self._n, self._d = n, d

    def _rho(self, s):
        return sum(di * s / np.sqrt(ni**2 - s**2) for di, ni in zip(self._d, self._n))

    def _drho(self, s):
        return sum(di * ni**2 / (ni**2 - s**2) ** 1.5 for di, ni in zip(self._d, self._n))

    def __call__(self, rho):
        rho = np.asarray(rho, dtype=float)
        s = self._inv(rho)
        for _ in range(2):  # Newton refinement of the inverse ray mapping
            s = s + (rho - self._rho(s)) / self._drho(s)
        return self.stack.opl(s, self.lam), s


def pupil_samples(stack: LayerStack, lam_um: float, na: float, n: int, illumination: dict):
    """Metasurface sample points (x, y) on an n x n grid over the aperture, with field weights."""
    radius = stack.aperture_for_na(na, lam_um)
    c = (np.arange(n) + 0.5) / n * 2 * radius - radius
    X, Y = np.meshgrid(c, c, indexing="ij")
    R = np.hypot(X, Y)
    inside = R <= radius
    kw = {k: v for k, v in illumination.items() if k != "kind"}
    amp = profile(illumination["kind"], R, radius, **kw)
    return X[inside], Y[inside], amp[inside], radius


def spot_aberrations(stack: LayerStack, lam_um: float, na: float, sites_um: np.ndarray, illumination: dict, n_pupil: int = 192) -> dict:
    """Residual-wavefront Strehl, lateral shift and axial shift for each site."""
    k0 = 2 * np.pi / lam_um
    x, y, amp, _ = pupil_samples(stack, lam_um, na, n_pupil, illumination)
    exact = ExactOPL(stack, lam_um)
    r = np.hypot(x, y)
    opl0, s = exact(r)
    with np.errstate(invalid="ignore", divide="ignore"):
        rhx = np.where(r > 0, x / np.where(r > 0, r, 1), 0.0)
        rhy = np.where(r > 0, y / np.where(r > 0, r, 1), 0.0)
    # Ray direction cosines towards the on-axis focus: sigma = -r_hat s.
    sx, sy = -rhx * s, -rhy * s
    cz = np.sqrt(1 - s**2)  # cos(theta) in the atom-plane medium (vacuum)
    w = amp**2  # intensity weights for the fit
    # Basis: piston, tilt (sigma_x, sigma_y), defocus (cos theta) in the atom medium.
    B = np.column_stack([np.ones_like(sx), sx, sy, cz])
    BW = B * w[:, None]
    G = B.T @ BW
    out = {k: np.zeros(len(sites_um)) for k in ("strehl", "dx_um", "dy_um", "dz_um", "rms_wave")}
    amp_sum = np.sum(amp)
    for j, (xj, yj) in enumerate(sites_um):
        rho = np.hypot(x - xj, y - yj)
        oplj, _ = exact(rho)
        W = oplj - opl0 - (sx * xj + sy * yj)  # um; Debye design = OPL0 + sigma . x_j... see docstring
        coef = np.linalg.solve(G, BW.T @ W)
        res = W - B @ coef
        # A wave exp(i k (sigma.dx + cos(theta) dz)) focuses at (dx, dz); the field phase is
        # -k W, so a fitted W = -(sigma.dx + cz dz) means the spot sits at (dx, dz).
        out["dx_um"][j], out["dy_um"][j], out["dz_um"][j] = -coef[1], -coef[2], -coef[3]
        out["strehl"][j] = np.abs(np.sum(amp * np.exp(1j * k0 * res))) ** 2 / amp_sum**2
        out["rms_wave"][j] = np.sqrt(np.sum(w * res**2) / np.sum(w)) / lam_um
    return out
