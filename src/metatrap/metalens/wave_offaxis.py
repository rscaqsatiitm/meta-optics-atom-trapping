"""Off-axis focal field of the doublet by superposition from the focusing plate.

The field just after the focusing plate is built from rays traced from a uniformly filled
stop through the corrector and the spacer. Geometrical optics is accurate there: the
corrector is weak and the Fresnel number R1^2 / (lambda L) ~ 4,800. Each pupil cell i
contributes amplitude sqrt(J_i) dA1 and phase k0 (opl_i + phi2(r_i) / k0). Here, J_i is the
local area magnification from the stop to the plate, dA1 the stop cell area, and opl_i
includes the corrector phase. The focal field at an atom-plane point P is the sum

    U(P) = sum_i a_i exp(i k0 [opl_i + OPL_stack(|r_i - P|)]),

where OPL_stack(rho) is the exact (stationary) optical path from a plate point to an
atom-plane point at transverse distance rho through the cell layers (hologram.field_model.
ExactOPL). This is the same kernel as the exact hologram model, which was validated against
full-pixel angular-spectrum propagation (campaign 1, scripts/validate_fullpixel.py). The
Strehl ratio is |U(P_best)|^2 / (sum_i |a_i|)^2, with P_best found by a lateral and axial
search. Optional amplitude weights eta(s) (meta-atom efficiency versus local deflection, power)
apodise the pupil.

What this tests beyond the ray-traced Strehl: the image is found as an intensity maximum of a
diffraction sum, rather than from a least-squares fit of the eikonal, and amplitude
(apodisation) effects on the peak and the waist are included. It does not test the field
between the plates (ray model there) or vector effects.
"""

from __future__ import annotations

import numpy as np

from ..hologram.exact import OPLTable
from ..propagation.layers import LayerStack
from .raytrace import Medium, System, point_source_rays


def plate_field(d, lam_um: float, theta: float, n_pupil: int = 400, extra_opl=None):
    """Samples (x2, y2, amplitude, phase/k0) just after the focusing plate for field angle theta.

    extra_opl(x, y): optional optical path (um) added at the stop, e.g. a per-field
    pre-compensation supplied by a lenslet of a lenslet-array relay."""
    R1 = d.info["R1_um"]
    c = (np.arange(n_pupil) + 0.5) / n_pupil * 2 * R1 - R1
    X, Y = np.meshgrid(c, c, indexing="ij")
    inside = np.hypot(X, Y) <= R1
    px, py, opl0 = point_source_rays(X, Y, theta, None)
    if extra_opl is not None:
        opl0 = opl0 + extra_opl(X, Y)
    from .doublet import _media

    sysm = System([d.m1] + _media(d.spec.spacer, lam_um))
    xf, yf, pxf, pyf, opl, alive = sysm.trace(X.ravel(), Y.ravel(), px.ravel(), py.ravel(), lam_um, opl0=opl0.ravel())
    x2, y2, opl = xf.reshape(X.shape), yf.reshape(X.shape), opl.reshape(X.shape)
    # Area magnification from the stop grid to the plate (central differences on the grid).
    dx2 = np.gradient(x2, c, axis=0), np.gradient(x2, c, axis=1)
    dy2 = np.gradient(y2, c, axis=0), np.gradient(y2, c, axis=1)
    J = np.abs(dx2[0] * dy2[1] - dx2[1] * dy2[0])
    r2 = np.hypot(x2, y2)
    opl = opl + d.m2.phase(r2, lam_um) / (2 * np.pi / lam_um)
    dA1 = (c[1] - c[0]) ** 2
    m = inside & alive.reshape(X.shape)
    s2 = np.abs(d.m2.dphase_dr(r2, lam_um)) / (2 * np.pi / lam_um)
    s1 = np.abs(d.m1.dphase_dr(np.hypot(X, Y), lam_um)) / (2 * np.pi / lam_um)
    return {"x": x2[m], "y": y2[m], "amp": np.sqrt(J[m]) * dA1, "opl": opl[m], "s1": s1[m], "s2": s2[m]}


class FocalSum:
    def __init__(self, cell: list, lam_um: float, rho_max_um: float):
        self.cell, self.lam = cell, lam_um
        self.rho_max = rho_max_um
        self._tables = {}

    def table(self, dz: float):
        key = round(dz, 4)
        if key not in self._tables:
            layers = [list(l) for l in self.cell]
            layers[-1][0] = float(layers[-1][0]) + dz
            self._tables[key] = OPLTable(LayerStack.from_list(layers), self.lam, self.rho_max, n=200_001)
            if len(self._tables) > 12:
                self._tables.pop(next(iter(self._tables)))
        return self._tables[key]

    def field(self, pf, amp, P, dz, block=64):
        """Complex field at atom-plane points P (M, 2) in the plane dz beyond the nominal one."""
        k0 = 2 * np.pi / self.lam
        tab = self.table(dz)
        out = np.empty(len(P), complex)
        # Remove a large common phase to keep the exponent well conditioned.
        for a in range(0, len(P), block):
            Pb = P[a : a + block]
            rho = np.hypot(pf["x"][None, :] - Pb[:, 0:1], pf["y"][None, :] - Pb[:, 1:2])
            ph = k0 * (pf["opl"][None, :] + tab(rho))
            out[a : a + block] = np.sum(amp[None, :] * np.exp(1j * ph), axis=1)
        return out


def strehl_offaxis(d, lam_um: float, theta: float, x_ray: float, dz_ray: float, eta=None, n_pupil: int = 400,
                   span_um: float = 0.5, dz_span_um: float = 1.0):
    """Best-focus Strehl ratio and spot widths of field angle theta (ray image at (x_ray, 0, dz_ray)).

    eta: optional (s_grid, eta_grid) of the meta-atom power efficiency versus deflection,
    applied to both plates as amplitude weights sqrt(eta)."""
    pf = plate_field(d, lam_um, theta, n_pupil)
    amp = pf["amp"].copy()
    if eta is not None:
        sg, eg = eta
        amp = amp * np.sqrt(np.interp(pf["s1"], sg, eg) * np.interp(pf["s2"], sg, eg))
    fs = FocalSum(d.spec.cell, lam_um, d.info["R2_um"] * 1.15 + abs(x_ray) + 50)
    ref = np.sum(np.abs(amp))
    # Axial and lateral search for the intensity maximum.
    best = None
    for dz in dz_ray + np.linspace(-dz_span_um, dz_span_um, 9):
        g = np.linspace(-span_um, span_um, 11)
        GX, GY = np.meshgrid(x_ray + g, g, indexing="ij")
        P = np.column_stack([GX.ravel(), GY.ravel()])
        I = np.abs(fs.field(pf, amp, P, dz)) ** 2
        j = int(np.argmax(I))
        if best is None or I[j] > best[0]:
            best = (I[j], dz, P[j])
    # Refine around the best point.
    _, dz0, p0 = best
    for _ in range(2):
        g = np.linspace(-0.08, 0.08, 9)
        cands = []
        for dz in dz0 + np.linspace(-0.12, 0.12, 7):
            GX, GY = np.meshgrid(p0[0] + g, p0[1] + g, indexing="ij")
            P = np.column_stack([GX.ravel(), GY.ravel()])
            I = np.abs(fs.field(pf, amp, P, dz)) ** 2
            j = int(np.argmax(I))
            cands.append((I[j], dz, P[j]))
        best = max(cands, key=lambda c: c[0])
        _, dz0, p0 = best
    # Line profiles through the peak (along x and y), 1/e^2 half-widths.
    r = np.linspace(-1.5, 1.5, 151)
    Ix = np.abs(fs.field(pf, amp, np.column_stack([p0[0] + r, np.full_like(r, p0[1])]), dz0)) ** 2
    Iy = np.abs(fs.field(pf, amp, np.column_stack([np.full_like(r, p0[0]), p0[1] + r]), dz0)) ** 2

    def half_width(I):
        I = I / I.max()
        c = len(I) // 2 + int(np.argmax(I[len(I) // 2 - 5 : len(I) // 2 + 6]) - 5)
        right = c + int(np.argmax(I[c:] < np.exp(-2)))
        left = c - int(np.argmax(I[:c + 1][::-1] < np.exp(-2)))
        return float(0.5 * (r[right] - r[left]))

    ref0 = np.sum(pf["amp"])
    return {"strehl": float(best[0] / ref**2), "peak_rel_uniform_lossless": float(best[0] / ref0**2), "x_um": float(p0[0]), "y_um": float(p0[1]), "dz_um": float(dz0),
            "w_x_um": half_width(Ix), "w_y_um": half_width(Iy), "n_samples": int(len(amp)),
            "transmitted_weight": float(np.sum(amp**2) / np.sum(pf["amp"] ** 2))}


def plate_field_front(d, lam_um: float, theta: float, front: list, stop_index: int, stop_radius: float, launch_radius: float, n_pupil: int = 300):
    """As plate_field, for a system with extra elements in front of the corrector metasurface
    (e.g. a per-wavelength pre-corrector A | air | C). Rays are launched on the first front
    element over launch_radius and vignetted at element stop_index."""
    c = (np.arange(n_pupil) + 0.5) / n_pupil * 2 * launch_radius - launch_radius
    X, Y = np.meshgrid(c, c, indexing="ij")
    px, py, opl0 = point_source_rays(X, Y, theta, None)
    from .doublet import _media

    sysm = System(list(front) + [d.m1] + _media(d.spec.spacer, lam_um))
    xf, yf, pxf, pyf, opl, alive = sysm.trace(X.ravel(), Y.ravel(), px.ravel(), py.ravel(), lam_um,
                                              stop={"index": stop_index, "radius_um": stop_radius}, opl0=opl0.ravel())
    x2, y2, opl, alive = xf.reshape(X.shape), yf.reshape(X.shape), opl.reshape(X.shape), alive.reshape(X.shape)
    dx2 = np.gradient(x2, c, axis=0), np.gradient(x2, c, axis=1)
    dy2 = np.gradient(y2, c, axis=0), np.gradient(y2, c, axis=1)
    J = np.abs(dx2[0] * dy2[1] - dx2[1] * dy2[0])
    r2 = np.hypot(x2, y2)
    opl = opl + d.m2.phase(r2, lam_um) / (2 * np.pi / lam_um)
    dA1 = (c[1] - c[0]) ** 2
    m = alive
    s2 = np.abs(d.m2.dphase_dr(r2, lam_um)) / (2 * np.pi / lam_um)
    return {"x": x2[m], "y": y2[m], "amp": np.sqrt(J[m]) * dA1, "opl": opl[m], "s1": np.zeros(m.sum()), "s2": s2[m]}


def strehl_from_field(d, lam_um, pf, x_ray, dz_ray, span_um=0.5, dz_span_um=1.0):
    """Best-focus Strehl ratio of a precomputed plate field (see strehl_offaxis)."""
    amp = pf["amp"]
    fs = FocalSum(d.spec.cell, lam_um, d.info["R2_um"] * 1.3 + abs(x_ray) + 50)
    ref = np.sum(np.abs(amp))
    best = None
    for dz in dz_ray + np.linspace(-dz_span_um, dz_span_um, 9):
        g = np.linspace(-span_um, span_um, 11)
        GX, GY = np.meshgrid(x_ray + g, g, indexing="ij")
        P = np.column_stack([GX.ravel(), GY.ravel()])
        I = np.abs(fs.field(pf, amp, P, dz)) ** 2
        j = int(np.argmax(I))
        if best is None or I[j] > best[0]:
            best = (I[j], dz, P[j])
    _, dz0, p0 = best
    for _ in range(2):
        g = np.linspace(-0.08, 0.08, 9)
        cands = []
        for dz in dz0 + np.linspace(-0.12, 0.12, 7):
            GX, GY = np.meshgrid(p0[0] + g, p0[1] + g, indexing="ij")
            P = np.column_stack([GX.ravel(), GY.ravel()])
            I = np.abs(fs.field(pf, amp, P, dz)) ** 2
            j = int(np.argmax(I))
            cands.append((I[j], dz, P[j]))
        best = max(cands, key=lambda c: c[0])
        _, dz0, p0 = best
    return {"strehl": float(best[0] / ref**2), "x_um": float(p0[0]), "dz_um": float(dz0)}
