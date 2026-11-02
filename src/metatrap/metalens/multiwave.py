"""One metasurface doublet for several wavelengths: per-wavelength stop correctors.

A metasurface designed at lambda_d keeps its phase in radians at another wavelength, so its
deflection scales as lambda / lambda_d (doublet.py). Refocusing the input (vergence) removes
the focal shift but leaves spherical aberration of several waves (Section 5 of
docs/metalens_report.md). Here each extra wavelength receives its own aspheric phase corrector
in the aperture stop, i.e. in a pupil plane of that wavelength's beam path before the beams are
combined (a refractive asphere, a fixed phase plate or an SLM relayed onto the stop). The
corrector phase is constructed so that the axial focus at that wavelength is stigmatic in the
nominal atom plane: for each stop height h1 we solve for the input transverse direction cosine
p(h1) whose ray, traced through the two (diffractively dispersed) metasurfaces and the cell,
crosses the axis in the atom plane. Then dphi_c/dr = k0 p(h1). The corrector also contains the
vergence. Off axis, the residual is what the corrector cannot remove: the sine condition is met
only at lambda_d, so coma grows with the field.

The first-order diffraction efficiency of the meta-atoms at lambda != lambda_d (atoms chosen
for the target phase at lambda_d) and the power left in the zero order are evaluated from the
single-pillar library (crosswave_efficiency).
"""

from __future__ import annotations

import numpy as np

from .aplanat import SplineMeta
from .doublet import _media
from .raytrace import System


def plates_system(d, lam_um: float, pre=None) -> System:
    el = ([pre] if pre is not None else []) + [d.m1] + _media(d.spec.spacer, lam_um) + [d.m2] + _media(d.spec.cell, lam_um)
    return System(el)


def stop_corrector(d, lam_um: float, n_pts: int = 801, r_ext: float = 1.1) -> SplineMeta:
    """Radial phase at the stop that makes the axial focus at lam_um stigmatic in the nominal plane."""
    R1 = d.info["R1_um"]
    h = np.linspace(0.0, r_ext * R1, n_pts)
    sysm = plates_system(d, lam_um)

    def x_end(p):
        xf, yf, pxf, pyf, opl, alive = sysm.trace(h, np.zeros_like(h), p, np.zeros_like(h), lam_um)
        return xf, alive

    p = np.zeros_like(h)
    for _ in range(60):
        x0, a0 = x_end(p)
        dp = 1e-7
        x1, _ = x_end(p + dp)
        J = (x1 - x0) / dp
        step = np.where(np.abs(J) > 0, x0 / J, 0.0)
        p = p - np.clip(step, -0.02, 0.02)
        if np.max(np.abs(x0[1:])) < 1e-6:
            break
    k0 = 2 * np.pi / lam_um
    g = k0 * p
    phi = np.concatenate([[0.0], np.cumsum(0.5 * (g[1:] + g[:-1]) * np.diff(h))])
    pc = SplineMeta(h, phi, lam_um, f"stop corrector {lam_um}")
    pc.residual_um = float(np.max(np.abs(x0)))
    return pc


def corrector_summary(pc: SplineMeta, R1: float, lam_um: float):
    """Paraxial focal length of the corrector (its vergence) and its aspheric part in waves."""
    r = np.linspace(0, R1, 2001)
    phi = pc.phase(r, lam_um)
    k0 = 2 * np.pi / lam_um
    # Least-squares fit phi = a r^2 + b r^4 + ... ; quadratic part -> focal length (phi = -k r^2 / 2F).
    A = np.column_stack([r**2, r**4, r**6, r**8])
    c, *_ = np.linalg.lstsq(A, phi, rcond=None)
    F = -k0 / (2 * c[0]) if c[0] != 0 else np.inf
    asph = phi - c[0] * r**2
    asph = asph - np.polyval(np.polyfit(r**2, asph, 1), r**2)  # remove best-fit quadratic again
    return {"focal_length_mm": float(F / 1000), "aspheric_pv_waves": float(np.ptp(asph) / (2 * np.pi)),
            "total_pv_waves": float(np.ptp(phi) / (2 * np.pi))}


def crosswave_efficiency(lib, material: str, pitch: float, height: float, lam_d: float, lams, lim, n_theta: int = 720):
    """First-order efficiency and zero-order power at each wavelength of atoms chosen at lam_d.

    For a target phase theta (uniform on [0, 2 pi)) the atom maximising Re[t_d e^{-i theta}] is
    used; at wavelength lam the first order (the same phase profile in radians) receives
    n_ratio |<t_lam e^{-i theta}>|^2 and the zero order n_ratio |<t_lam>|^2."""
    from ..rcwa.lookup import densify, feasible

    w, td, T0d = lib.curve(material, pitch, height, lam_d)
    ok = feasible(w, pitch, height, lim)
    th = 2 * np.pi * (np.arange(n_theta) + 0.5) / n_theta
    wd, tdd = densify(w[ok], td[ok])
    k = np.argmax(np.real(tdd[None, :] * np.exp(-1j * th)[:, None]), axis=1)
    out = {}
    for lam in lams:
        w2, t2, T02 = lib.curve(material, pitch, height, lam)
        wd2, td2 = densify(w2[ok], t2[ok])
        nr = float(np.median(T02[ok] / np.maximum(np.abs(t2[ok]) ** 2, 1e-12)))
        tt = td2[k]
        out[lam] = {"eta_first": nr * float(np.abs(np.mean(tt * np.exp(-1j * th))) ** 2),
                    "zero_order": nr * float(np.abs(np.mean(tt)) ** 2),
                    "transmitted": nr * float(np.mean(np.abs(tt) ** 2))}
    return out


def aplanatic_precorrector(d, lam_um: float, L_um: float, na: float | None = None, n_pts: int = 4001, ext: float = 1.1, D_um: float = 0.0):
    """Two phase surfaces A and C (A a distance L_um in air before C, C in the stop at the
    corrector metasurface) that make the fixed metasurface doublet stigmatic AND aplanatic at
    lam_um. Rays are back-traced from the focus: a final ray of direction cosine s_f meets the
    focusing plate at h2 = H_cell(s_f); undoing its deflection there and at the corrector
    (diffractive dispersion, phases in radians fixed) gives its height h1 and the direction p_c it
    must have just after C. The sine condition puts its input height at A at hA = s_f f with
    f = R1 / NA, and the straight path from A to C fixes the direction q between them. Then
    dphi_A/dr (hA) = k0 q and dphi_C/dr (h1) = k0 (p_c - q). In practice A and C sit in the
    wavelength's own beam path and are relayed onto the stop (a 1:1 relay is assumed ideal), or,
    with D_um > 0, C sits D_um before the corrector metasurface (room for a dichroic combiner, no
    relay; the stop is then at C)."""
    from ..propagation.layers import LayerStack

    na = d.spec.na if na is None else na
    R1 = d.info["R1_um"]
    f = R1 / na
    k0 = 2 * np.pi / lam_um
    s_f = np.linspace(0.0, ext * na, n_pts)
    cell = LayerStack.from_list(d.spec.cell)
    h2 = cell.ray_height(s_f, lam_um)
    p_in = -s_f - d.m2.dphase_dr(h2, lam_um) / k0  # direction between the plates (signed, radial)
    sp = LayerStack.from_list(d.spec.spacer)
    n_sp, t_sp = sp.indices(lam_um), sp.thicknesses()
    h1 = h2 - sum(t * p_in / np.sqrt(n**2 - p_in**2) for t, n in zip(t_sp, n_sp))
    p_c = p_in - d.m1.dphase_dr(h1, lam_um) / k0
    if D_um:
        # C sits D_um in air before the corrector metasurface (no relay): back-propagate.
        h1 = h1 - D_um * p_c / np.sqrt(1 - p_c**2)
    hA = s_f * f
    dlt = h1 - hA
    q = dlt / np.sqrt(L_um**2 + dlt**2)
    if np.any(np.diff(h1) <= 0):
        raise ValueError("h1(s_f) not monotonic")
    gA, gC = k0 * q, k0 * (p_c - q)
    phiA = np.concatenate([[0.0], np.cumsum(0.5 * (gA[1:] + gA[:-1]) * np.diff(hA))])
    phiC = np.concatenate([[0.0], np.cumsum(0.5 * (gC[1:] + gC[:-1]) * np.diff(h1))])
    A = SplineMeta(hA, phiA, lam_um, f"pre-corrector A {lam_um}")
    C = SplineMeta(h1, phiC, lam_um, f"pre-corrector C {lam_um}")
    info = {"f_um": float(f), "R_A_um": float(R1), "h1_max_um": float(np.interp(na, s_f, h1)),
            "max_q": float(np.max(np.abs(q[s_f <= na]))), "max_defl_C": float(np.max(np.abs((p_c - q)[s_f <= na]))),
            "A_pv_waves": float(np.ptp(phiA[s_f <= na]) / (2 * np.pi)), "C_pv_waves": float(np.ptp(phiC[s_f <= na]) / (2 * np.pi))}
    return A, C, info


def precorrected_system(d, lam_um: float, A, C, L_um: float, D_um: float = 0.0) -> System:
    from .raytrace import Medium

    mid = [Medium(1.0, D_um, "air (dichroic combiner)")] if D_um else []
    return System([A, Medium(1.0, L_um, "air (relayed)"), C] + mid + [d.m1] + _media(d.spec.spacer, lam_um) + [d.m2] + _media(d.spec.cell, lam_um))
