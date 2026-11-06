"""Sequential ray tracing through plane-parallel media and planar metasurfaces.

Coordinates: z along the optical axis, light travelling towards +z, from a collimated input
(field angle theta in air) to the atom plane. A ray is described by its transverse position
(x, y) and its transverse optical direction cosines p = n (sx, sy), which are conserved across
plane interfaces (Snell's law) and change at a metasurface by the phase gradient,

    p_out = p_in + grad(phi) / k0       (generalised Snell law),

where phi(x, y) is the phase the metasurface adds to the transmitted field and k0 = 2 pi / lambda.
In a medium of index n a ray with |p| < n advances transversally by t p / sqrt(n^2 - |p|^2)
over a thickness t, and its optical path grows by t n^2 / sqrt(n^2 - |p|^2). The total phase
of the ray is k0 OPL + sum of phi at the metasurfaces; this is stationary along true rays,
which is what makes the OPD of a focused bundle the wavefront aberration.

Elements of a System are Medium(n, thickness) or Meta(phase), in order. Stops are given
separately (radius at a named element).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable

import numpy as np


@dataclass
class Medium:
    n: float
    thickness_um: float
    name: str = ""


@dataclass
class Meta:
    """Radially symmetric metasurface phase phi(r) (rad) at the design wavelength.

    phi(r) = k_d * sum_m c_m r^(2m), r in um, c_m in um^(1-2m); k_d = 2 pi / lambda_d.
    ``disp`` maps the wavelength to the phase scale relative to the design wavelength
    (None: phase unchanged, i.e. a dispersionless phase response).
    """

    coeffs: np.ndarray
    lam_design_um: float
    name: str = ""
    disp: Callable[[float], float] | None = None

    def scale(self, lam_um: float) -> float:
        return 1.0 if self.disp is None else float(self.disp(lam_um))

    def phase(self, r, lam_um: float):
        kd = 2 * np.pi / self.lam_design_um
        r2 = np.asarray(r, float) ** 2
        return self.scale(lam_um) * kd * sum(c * r2 ** (m + 1) for m, c in enumerate(self.coeffs))

    def dphase_dr(self, r, lam_um: float):
        kd = 2 * np.pi / self.lam_design_um
        r = np.asarray(r, float)
        return self.scale(lam_um) * kd * sum(c * 2 * (m + 1) * r ** (2 * m + 1) for m, c in enumerate(self.coeffs))


@dataclass
class Frame:
    """Change of coordinate frame, for decentre and tilt tolerances.

    The downstream elements are rotated by ``tilt_rad`` about the y axis through the current
    origin and then shifted by (dx_um, dy_um) in the rotated frame. Rays are intersected exactly
    with the rotated plane (in the medium of index ``n`` in which the frame change lies, normally
    air), and their positions and transverse direction cosines are re-expressed in the new
    frame. Frame(b, dx, dy) ... Frame(-b, -dx, -dy) around an element tilts and decentres that
    element alone (to first order in the shift for combined tilt and decentre).
    """

    tilt_rad: float = 0.0
    dx_um: float = 0.0
    dy_um: float = 0.0
    n: float = 1.0


def _apply_frame(fr: Frame, x, y, px, py, opl, alive):
    n = fr.n
    p2 = px**2 + py**2
    ok = p2 < n**2
    cz = np.sqrt(np.where(ok, n**2 - p2, 1.0))
    ux, uy, uz = px / n, py / n, cz / n
    sb, cb = np.sin(fr.tilt_rad), np.cos(fr.tilt_rad)
    # Plane through the origin with normal N = (sin b, 0, cos b); ray X = (x, y, 0) + s u.
    s = -(x * sb) / (ux * sb + uz * cb)
    X, Y, Z = x + s * ux, y + s * uy, s * uz
    opl = opl + n * s
    # In-plane axes e1 = (cos b, 0, -sin b), e2 = (0, 1, 0).
    xn = X * cb - Z * sb - fr.dx_um
    yn = Y - fr.dy_um
    pxn = n * (ux * cb - uz * sb)
    pyn = n * uy
    return xn, yn, pxn, pyn, opl, alive & ok


@dataclass
class System:
    elements: list = field(default_factory=list)

    def trace(self, x, y, px, py, lam_um: float, stop: dict | None = None, opl0=None):
        """Trace rays; returns final (x, y, px, py, opl_um, alive) at the end of the last medium.

        ``stop`` = {"index": i, "radius_um": R} vignettes rays outside R at element i.
        """
        k0 = 2 * np.pi / lam_um
        x, y = np.array(x, float), np.array(y, float)
        px, py = np.array(px, float), np.array(py, float)
        # Phase of the incoming plane wave across the launch plane (so that image positions,
        # not only aberrations, come out right).
        opl = px * x + py * y if opl0 is None else np.array(opl0, float)
        alive = np.ones_like(x, bool)
        for i, el in enumerate(self.elements):
            if stop is not None and stop["index"] == i:
                alive &= np.hypot(x, y) <= stop["radius_um"]
            if isinstance(el, Frame):
                x, y, px, py, opl, alive = _apply_frame(el, x, y, px, py, opl, alive)
            elif isinstance(el, Medium):
                p2 = px**2 + py**2
                ok = p2 < el.n**2
                alive &= ok
                cz = np.sqrt(np.where(ok, el.n**2 - p2, 1.0))
                x = x + el.thickness_um * px / cz
                y = y + el.thickness_um * py / cz
                opl = opl + el.thickness_um * el.n**2 / cz
            else:
                r = np.hypot(x, y)
                g = el.dphase_dr(r, lam_um) / k0
                with np.errstate(invalid="ignore", divide="ignore"):
                    ux = np.where(r > 0, x / np.where(r > 0, r, 1), 0.0)
                    uy = np.where(r > 0, y / np.where(r > 0, r, 1), 0.0)
                px = px + g * ux
                py = py + g * uy
                opl = opl + el.phase(r, lam_um) / k0
        return x, y, px, py, opl, alive


def point_source_rays(x, y, theta: float, L_um: float | None):
    """Directions and launch optical path for rays at (x, y) on the first plane (in air).

    L_um = None: collimated beam at field angle theta (about y). Otherwise a spherical wave
    whose centre lies a distance |L| before (L > 0, diverging) or after (L < 0, converging) the
    plane, laterally at -L tan(theta) so that the chief direction stays theta. The launch
    optical path is the spherical-wave eikonal, so that image positions and wavefronts are
    consistent with the collimated case.
    """
    if L_um is None or not np.isfinite(L_um):
        px = np.full_like(x, np.sin(theta))
        return px, np.zeros_like(x), px * x
    xs = -L_um * np.tan(theta)
    d = np.sqrt((x - xs) ** 2 + y**2 + L_um**2)
    sg = np.sign(L_um)
    return sg * (x - xs) / d, sg * y / d, sg * (d - abs(L_um) / np.cos(theta))


def pupil_grid(radius_um: float, n: int):
    """Uniform square grid of ray positions inside a disc (equal-area weights)."""
    c = (np.arange(n) + 0.5) / n * 2 * radius_um - radius_um
    X, Y = np.meshgrid(c, c, indexing="ij")
    m = np.hypot(X, Y) <= radius_um
    return X[m], Y[m]


def spot_quality(xf, yf, pxf, pyf, opl, alive, lam_um: float, n_last: float = 1.0):
    """Strehl ratio and best image point of a focused bundle at the end plane.

    The OPL to a common image point P = (X, Y, dZ) relative to the end plane is, to first order
    in the distance from each ray's end point, opl_i + p_i . (P - q_i) + cz_i dZ. Piston, tilt
    (image position) and defocus (image plane) are removed by least squares; the residual
    wavefront gives the Strehl ratio |<exp(i k0 W)>|^2.
    """
    k0 = 2 * np.pi / lam_um
    a = alive
    px, py = pxf[a], pyf[a]
    cz = np.sqrt(n_last**2 - px**2 - py**2)
    # Phase referenced to the transverse origin: opl - p . q (eikonal in direction space).
    w = opl[a] - (px * xf[a] + py * yf[a])
    B = np.column_stack([np.ones_like(px), px, py, cz])
    coef, *_ = np.linalg.lstsq(B, w, rcond=None)
    res = w - B @ coef
    strehl = float(np.abs(np.mean(np.exp(1j * k0 * res))) ** 2)
    # The bundle focuses where the eikonal is stationary: P = -(coef_x, coef_y), dZ = -coef_z.
    return {"strehl": strehl, "rms_waves": float(np.std(res) / lam_um), "x_um": float(-coef[1]),
            "y_um": float(-coef[2]), "dz_um": float(-coef[3]), "n_rays": int(a.sum()), "residual": res}
