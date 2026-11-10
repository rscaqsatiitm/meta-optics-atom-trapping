"""Air-spaced aplanatic metasurface doublet: system assembly, field evaluation, chromatic and
tolerance analysis (ray optics).

Layout, from the collimated side: corrector metasurface (aperture stop, radius R1) | spacer
(air, optionally with the substrate of the second plate) | focusing metasurface | cell stack
(air gap, wall, vacuum) | atom plane. Phases come from aplanat.design. A metasurface designed
at lambda_d and used at lambda keeps its phase in radians (fixed 2 pi zones), so its deflection
scales as lambda / lambda_d (diffractive dispersion); this is the SplineMeta default.

Field points are specified by their image height in the atom plane; the collimated input
angle that produces it is found by a secant iteration. Image quality is the Strehl ratio from
the direction-space eikonal (design._quality), with free per-field image position and either a
free focus per field or one common focal plane.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import brentq, minimize_scalar

from .. import materials
from . import aplanat
from .design import _quality
from .raytrace import Frame, Medium, System, point_source_rays, pupil_grid


@dataclass
class DoubletSpec:
    na: float
    cell: list  # [[thickness_um, material], ...] after the focusing plate
    spacer: list  # [[thickness_um, material], ...] between the plates
    n_pts: int = 4001
    ext: float = 1.1  # construction radius / R1 (must cover R2 + field height on the focusing plate)


@dataclass
class Doublet:
    m1: object
    m2: object
    info: dict
    lam_design_um: float
    spec: DoubletSpec


def make(spec: DoubletSpec, lam_design_um: float, R1_um: float | None = None) -> Doublet:
    m1, _, m2, _, info = aplanat.design(lam_design_um, spec.na, spec.cell, spacer=spec.spacer, R1_um=R1_um, n_pts=spec.n_pts, ext=spec.ext)
    return Doublet(m1, m2, info, lam_design_um, spec)


def _media(layers, lam_um):
    return [Medium(float(materials.get(m).n(lam_um)), float(t), m) for t, m in layers]


@dataclass
class Perturbation:
    """Mechanical errors. Lengths in um, angles in rad (tilts about the y axis)."""

    spacer_delta_um: float = 0.0  # added to the first spacer layer
    wall_delta_um: float = 0.0  # added to the cell layer with index wall_index
    wall_index: int = 1
    m2_decentre_um: float = 0.0  # focusing plate, along x
    m2_tilt_rad: float = 0.0  # focusing plate about its centre
    cell_tilt_rad: float = 0.0  # cell (wall and atom plane) relative to the objective
    # Evaluation plane for a cell tilt: "cell" (perpendicular to the tilted cell normal, i.e. the
    # atoms lie in a plane fixed to the cell; the field edge is then defocused by h x tilt) or
    # "objective" (perpendicular to the objective axis: the atom plane is set by the optics, as
    # for traps made through the same objective; only the aberration of the tilted wall remains).
    cell_tilt_plane: str = "cell"


def system(d: Doublet, lam_um: float, pert: Perturbation | None = None, m1=None, m2=None) -> System:
    """Elements for evaluation at lam_um (m1, m2 override the plates, e.g. for a dual design)."""
    pert = pert or Perturbation()
    spacer = [list(l) for l in d.spec.spacer]
    spacer[0][0] = float(spacer[0][0]) + pert.spacer_delta_um
    cell = [list(l) for l in d.spec.cell]
    if pert.wall_delta_um:
        cell[pert.wall_index][0] = float(cell[pert.wall_index][0]) + pert.wall_delta_um
    el = [m1 or d.m1] + _media(spacer, lam_um)
    n_sp = el[-1].n
    cell_media = _media(cell, lam_um)
    n_after = cell_media[0].n  # medium after the focusing plate (first cell layer; air in all present designs)
    if pert.m2_decentre_um or pert.m2_tilt_rad:
        el += [Frame(pert.m2_tilt_rad, pert.m2_decentre_um, 0.0, n_sp), m2 or d.m2, Frame(-pert.m2_tilt_rad, -pert.m2_decentre_um, 0.0, n_after)]
    else:
        el += [m2 or d.m2]
    if pert.cell_tilt_rad:
        el += [Frame(pert.cell_tilt_rad, 0.0, 0.0, n_after)]
    el += cell_media
    if pert.cell_tilt_rad and pert.cell_tilt_plane == "objective":
        el += [Frame(-pert.cell_tilt_rad, 0.0, 0.0, cell_media[-1].n)]
    elif pert.cell_tilt_plane not in ("cell", "objective"):
        raise ValueError("cell_tilt_plane must be 'cell' or 'objective'")
    return System(el)


def trace_field(sysm: System, R1: float, lam_um: float, theta: float, L_um=None, n_rays: int = 40):
    x, y = pupil_grid(R1, n_rays)
    px, py, opl0 = point_source_rays(x, y, theta, L_um)
    return sysm.trace(x, y, px, py, lam_um, stop={"index": 0, "radius_um": R1}, opl0=opl0)


def quality(sysm, R1, lam_um, theta, L_um=None, dz_common=None, n_rays=40):
    xf, yf, pxf, pyf, opl, alive = trace_field(sysm, R1, lam_um, theta, L_um, n_rays)
    q = _quality(xf, yf, pxf, pyf, opl, alive, lam_um, dz_common)
    q.pop("residual", None)
    return q


def theta_for_height(sysm, R1, lam_um, h_um, f_um, L_um=None):
    """Input angle whose image lies at height h_um (secant iteration on the chief image)."""
    if h_um == 0:
        return 0.0
    t0 = np.arctan(h_um / f_um)
    x0 = quality(sysm, R1, lam_um, t0, L_um, n_rays=16)["x_um"]
    t1 = t0 * h_um / x0
    for _ in range(4):
        x1 = quality(sysm, R1, lam_um, t1, L_um, n_rays=16)["x_um"]
        if abs(x1 - h_um) < 1e-3:
            break
        t0, t1, x0 = t1, t1 + (h_um - x1) * (t1 - t0) / (x1 - x0), x1
    return float(t1)


def field_scan(sysm, R1, lam_um, heights_um, f_um, L_um=None, dz_common=None, n_rays=40):
    """Strehl and image data per field height; with dz_common=None each field has its own focus."""
    out = []
    for h in heights_um:
        th = theta_for_height(sysm, R1, lam_um, h, f_um, L_um)
        q = quality(sysm, R1, lam_um, th, L_um, dz_common, n_rays)
        out.append({"height_um": float(h), "theta_deg": float(np.degrees(th))} | q)
    return out


def best_plane(sysm, R1, lam_um, thetas, L_um=None, n_rays=40, span_um=6.0, step_um=0.05):
    """Common focal plane maximising the minimum Strehl over the given input angles.

    The Strehl ratio versus defocus is narrow (depth of focus ~ lambda / NA^2) and has side
    lobes, so a grid scan around the mean free-focus position precedes a bounded refinement.
    Returns (dz, [quality per angle])."""
    free = [quality(sysm, R1, lam_um, t, L_um, None, n_rays)["dz_um"] for t in thetas]
    c = float(np.mean(free))
    span = span_um + 0.5 * (max(free) - min(free))

    def min_s(dz):
        return min(quality(sysm, R1, lam_um, t, L_um, dz, n_rays)["strehl"] for t in thetas)

    grid = np.arange(c - span, c + span + step_um / 2, step_um)
    vals = np.array([min_s(z) for z in grid])
    z0 = grid[np.argmax(vals)]
    r = minimize_scalar(lambda z: -min_s(z), bounds=(z0 - step_um, z0 + step_um), method="bounded", options={"xatol": 1e-3})
    dz = float(r.x) if -r.fun >= vals.max() else float(z0)
    return dz, [quality(sysm, R1, lam_um, t, L_um, dz, n_rays) for t in thetas]


def common_plane(sysm, R1, lam_um, heights_um, f_um, L_um=None, n_rays=40, span_um=6.0):
    """Common focal plane maximising the minimum Strehl over the listed field heights."""
    thetas = [theta_for_height(sysm, R1, lam_um, h, f_um, L_um) for h in heights_um]
    dz, q = best_plane(sysm, R1, lam_um, thetas, L_um, n_rays, span_um)
    return dz, [qq | {"height_um": float(h)} for qq, h in zip(q, heights_um)]


def marginal_na(sysm, R1, lam_um, L_um=None):
    x, y = np.array([0.999 * R1]), np.array([0.0])
    px, py, opl0 = point_source_rays(x, y, 0.0, L_um)
    xf, yf, pxf, pyf, opl, alive = sysm.trace(x, y, px, py, lam_um, opl0=opl0)
    return float(abs(pxf[0])) if alive[0] else float("nan")


def vergence_for_focus(sysm, R1, lam_um, dz_target=0.0, n_rays=32, v_max=1 / 20000.0):
    """Input vergence v = 1/L (1/um; L > 0 diverging source at distance L) that brings the
    best on-axis focus to dz_target, then refined to maximise the on-axis Strehl in that plane."""

    def dz_free(v):
        L = None if v == 0 else 1.0 / v
        return quality(sysm, R1, lam_um, 0.0, L, None, n_rays)["dz_um"] - dz_target

    vs = np.linspace(-v_max, v_max, 41)
    vals = []
    for v in vs:
        try:
            vals.append(dz_free(v))
        except Exception:
            vals.append(np.nan)
    vals = np.array(vals)
    ok = np.isfinite(vals)
    sgn = np.sign(vals)
    idx = [i for i in range(len(vs) - 1) if ok[i] and ok[i + 1] and sgn[i] != sgn[i + 1]]
    if not idx:
        raise RuntimeError("no vergence brings the focus to the target plane within the search range")
    i = min(idx, key=lambda j: abs(vs[j]))
    v0 = brentq(dz_free, vs[i], vs[i + 1], xtol=1e-12)
    dv = 0.02 * max(abs(v0), 1e-7)

    def neg(v):
        L = None if v == 0 else 1.0 / v
        return -quality(sysm, R1, lam_um, 0.0, L, dz_target, n_rays)["strehl"]

    r = minimize_scalar(neg, bounds=(v0 - dv, v0 + dv), method="bounded", options={"xatol": 1e-12})
    return float(r.x if -r.fun > -neg(v0) else v0)


def phase_pair_uniformity(ma, mb, lam_a, lam_b, r_max, n_r=200001, bins=36):
    """Area-weighted histogram of (phi_a mod 2pi, phi_b mod 2pi) over a disc: the pairs a dual-
    wavelength meta-atom library must provide. Returns (min, max) bin occupancy relative to
    uniform and the histogram."""
    r = np.linspace(0, r_max, n_r)
    pa = np.mod(ma.phase(r, lam_a), 2 * np.pi)
    pb = np.mod(mb.phase(r, lam_b), 2 * np.pi)
    H, _, _ = np.histogram2d(pa, pb, bins=bins, range=[[0, 2 * np.pi], [0, 2 * np.pi]], weights=r)
    H = H / H.mean()
    return float(H.min()), float(H.max()), H
