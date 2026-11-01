"""Wide-field metalens architectures and phase-profile optimisation (ray optics).

Each architecture is a function of a parameter vector (dimensionless radial-polynomial
coefficients a_m, phi(r) = k_d R sum_m a_m (r / R)^(2m)) that returns a System, the stop and
the ray launch plane. The merit is evaluated at several field angles with one common focal
plane (atoms and camera lie in flat planes, so field curvature counts as an error) and
free image position per field (distortion can be calibrated out).

Architectures (light from the collimated side towards the atoms; the cell stack follows):
  singlet:       substrate | metasurface (stop) | air gap | wall | vacuum
  stop_shifted:  aperture stop in air at distance d before the metasurface (on its substrate)
  doublet:       metasurface 1 (stop) | substrate t_d | metasurface 2 | air gap | wall | vacuum
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import minimize

from .. import materials
from .raytrace import Medium, Meta, System, pupil_grid, spot_quality


def _c(a, R):
    # phi = k_d R sum a_m (r/R)^(2m) = k_d sum (a_m R^(1-2m)) r^(2m)
    return np.array([am * R ** (1 - 2 * (m + 1)) for m, am in enumerate(a)])


@dataclass
class Geometry:
    lam_design_um: float
    cell: list  # [[thickness_um, material], ...] after the last metasurface (air gap, wall, vacuum)
    substrate: str = "fused_silica"
    substrate_um: float = 500.0


def cell_media(geom: Geometry, lam_um: float):
    return [Medium(float(materials.get(m).n(lam_um)), float(t), m) for t, m in geom.cell]


def build(arch: str, params: np.ndarray, geom: Geometry, lam_um: float, R: float, n_terms: int, extra: dict):
    ns = float(materials.get(geom.substrate).n(lam_um))
    el = []
    if arch == "singlet":
        el += [Medium(ns, geom.substrate_um, "substrate"), Meta(_c(params[:n_terms], R), geom.lam_design_um, "lens")]
        stop = {"index": 1, "radius_um": R}
    elif arch == "stop_shifted":
        d = extra["stop_distance_um"]
        el += [Medium(1.0, d, "air"), Medium(ns, geom.substrate_um, "substrate"), Meta(_c(params[:n_terms], extra["lens_R"]), geom.lam_design_um, "lens")]
        stop = {"index": 0, "radius_um": R}
    elif arch == "doublet":
        t_d = extra["doublet_substrate_um"]
        el += [Meta(_c(params[:n_terms], R), geom.lam_design_um, "corrector"), Medium(ns, t_d, "substrate"),
               Meta(_c(params[n_terms : 2 * n_terms], extra["lens_R"]), geom.lam_design_um, "focusing")]
        stop = {"index": 0, "radius_um": R}
    else:
        raise ValueError(arch)
    el += cell_media(geom, lam_um)
    return System(el), stop


def launch(R: float, theta: float, n_rays: int, stop_index: int, arch: str, extra: dict):
    """Rays filling the stop for a collimated beam at field angle theta (in air, about y)."""
    # Rays fill the stop. For stop_shifted/doublet the stop is the first element; for the
    # singlet the stop is the lens itself and rays are launched in the substrate at the lens
    # plane (a plane substrate conserves p and only shifts the bundle laterally).
    x, y = pupil_grid(R, n_rays)
    px = np.full_like(x, np.sin(theta))
    py = np.zeros_like(x)
    return x, y, px, py


def evaluate(arch, params, geom, lam_um, R, n_terms, extra, thetas, n_rays=36, dz_common=None):
    """Spot quality per field angle; with dz_common the image plane is shared (no per-field defocus)."""
    sysm, stop = build(arch, params, geom, lam_um, R, n_terms, extra)
    out = []
    for th in thetas:
        x, y, px, py = launch(R, th, n_rays, stop["index"], arch, extra)
        if arch == "singlet":
            # Launch inside the substrate at the lens plane: skip the substrate medium.
            sys2 = System(sysm.elements[1:])
            xf, yf, pxf, pyf, opl, alive = sys2.trace(x, y, px, py, lam_um, stop={"index": 0, "radius_um": R})
        else:
            xf, yf, pxf, pyf, opl, alive = sysm.trace(x, y, px, py, lam_um, stop=stop)
        out.append(_quality(xf, yf, pxf, pyf, opl, alive, lam_um, dz_common))
    return out


def _quality(xf, yf, pxf, pyf, opl, alive, lam_um, dz_common):
    if dz_common is None:
        return spot_quality(xf, yf, pxf, pyf, opl, alive, lam_um)
    a = alive
    px, py = pxf[a], pyf[a]
    cz = np.sqrt(1 - px**2 - py**2)
    w = opl[a] - (px * xf[a] + py * yf[a]) + cz * dz_common
    B = np.column_stack([np.ones_like(px), px, py])
    coef, *_ = np.linalg.lstsq(B, w, rcond=None)
    res = w - B @ coef
    k0 = 2 * np.pi / lam_um
    return {"strehl": float(np.abs(np.mean(np.exp(1j * k0 * res))) ** 2), "rms_waves": float(np.std(res) / lam_um),
            "x_um": float(-coef[1]), "y_um": float(-coef[2]), "dz_um": float(dz_common), "n_rays": int(a.sum()), "residual": res}


def marginal_na(arch, params, geom, lam_um, R, n_terms, extra):
    """|p| in vacuum of the on-axis marginal ray (defines the NA)."""
    sysm, stop = build(arch, params, geom, lam_um, R, n_terms, extra)
    x, y = np.array([R * 0.999]), np.array([0.0])
    els = sysm.elements[1:] if arch == "singlet" else sysm.elements
    xf, yf, pxf, pyf, opl, alive = System(els).trace(x, y, np.zeros(1), np.zeros(1), lam_um)
    return float(np.hypot(pxf[0], pyf[0])) if alive[0] else np.nan


def optimise(arch, p0, geom, lam_um, R, n_terms, extra, thetas, na_target, weights=None, n_rays=28, maxiter=400):
    """Minimise a soft maximum of the rms wavefront over field angles (common image plane),
    with the on-axis marginal ray held at the target NA. The last parameter is dz_common."""
    weights = np.ones(len(thetas)) if weights is None else np.asarray(weights)

    def merit(v):
        params, dz = v[:-1], v[-1]
        q = evaluate(arch, params, geom, lam_um, R, n_terms, extra, thetas, n_rays=n_rays, dz_common=dz)
        rms = np.array([qq["rms_waves"] for qq in q])
        if not np.all(np.isfinite(rms)) or any(qq["n_rays"] < 10 for qq in q):
            return 1e3
        na = marginal_na(arch, params, geom, lam_um, R, n_terms, extra)
        pen = 1e3 * (na - na_target) ** 2 if np.isfinite(na) else 1e3
        soft = np.log(np.sum(weights * np.exp(20 * rms))) / 20
        return float(soft + pen)

    v0 = np.concatenate([p0, [0.0]])
    res = minimize(merit, v0, method="Powell", options={"maxiter": maxiter, "xtol": 1e-6, "ftol": 1e-9})
    return res.x[:-1], res.x[-1], res


def hyperbolic_start(geom: Geometry, lam_um: float, R_fit: float, n_terms: int, extra_air_um: float = 0.0):
    """Dimensionless coefficients of the exact on-axis focusing phase through the cell stack
    (plus extra_air_um of air before it), fitted on 0..R_fit."""
    from ..propagation.layers import LayerStack

    st = LayerStack.from_list([[float(t) + (extra_air_um if i == 0 else 0.0), m] for i, (t, m) in enumerate(geom.cell)])
    r = np.linspace(0, R_fit, 3000)
    phi_over_kd = st.lens_phase(r, lam_um) / (2 * np.pi / geom.lam_design_um)
    A = np.column_stack([(r / R_fit) ** (2 * (m + 1)) for m in range(n_terms)])
    a, *_ = np.linalg.lstsq(A, phi_over_kd / R_fit, rcond=None)
    return a, st
