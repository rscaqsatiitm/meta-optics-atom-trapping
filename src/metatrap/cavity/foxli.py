"""Fox-Li model of a self-intersecting ring cavity with graded-phase mirrors (1D strip).

Model. The horizontal transverse profile of the sheet-beam mode is computed in one dimension
(strip resonator); the vertical profile is a Gaussian set by cylindrical curvature and is
separable, so a separable mirror phase and a rectangular aperture give a mode that is the
product of the two 1D modes (checked in tests/test_cavity.py). Unfolded, one round trip from
the crossing plane X on arm A is

    X(A) -> prop d1 -> aperture of fold mirror M2 -> prop (p - d1) -> graded mirror G1
         -> prop (p - d1) -> aperture of M3 -> prop d1 -> X(B) -> (same half again, G2) -> X(A).

Here, p is the path from the crossing to the graded mirror; the cavity is symmetric about
each graded mirror, which is what makes the conjugation design exact (below). Plane fold
mirrors only flip the transverse coordinate; three reflections per half trip map the 1D
coordinate xi of arm B to its physical left-hand transverse coordinate u_B = -xi.

Graded-phase design (Belanger and Pare, Opt. Lett. 16, 1057 (1991), DOI 10.1364/OL.16.001057;
not retrieved in session). For a real target profile F at X, the mirror phase
phi_G = -2 arg(P_p F) returns conj(P_p F) after reflection, and a further propagation over p
gives conj(F) = F. Hence F is an exact eigenmode of the half trip without apertures.

Driven cavity. With input coupler transmission T (amplitude t = sqrt(T)), lumped
non-diffractive round-trip power loss L and the round-trip operator R (propagation, apertures,
mirror phases), the steady-state circulating field at X(A) solves

    (I - r e^{i psi} R) E = t E_in,      r = sqrt((1 - T)(1 - L)),

with psi the laser detuning per round trip; the laser is locked to the psi that maximises the
circulating power. Lengths in um.
"""

from __future__ import annotations

import copy
from dataclasses import dataclass, field

import numpy as np
import scipy.sparse.linalg as spla

from .lattice import supergauss


@dataclass
class Grid1D:
    n: int
    width_um: float
    lam_um: float
    _h: dict = field(default_factory=dict, repr=False)

    @property
    def dx(self):
        return self.width_um / self.n

    @property
    def x(self):
        return (np.arange(self.n) - self.n // 2) * self.dx

    def transfer(self, z):
        key = round(float(z), 6)
        if key not in self._h:
            k = 2 * np.pi / self.lam_um
            kx = 2 * np.pi * np.fft.fftfreq(self.n, d=self.dx)
            kz2 = k**2 - kx**2
            self._h[key] = np.where(kz2 > 0, np.exp(1j * (np.sqrt(np.clip(kz2, 0, None)) - k) * z), 0.0)
        return self._h[key]

    def prop(self, u, z):
        if z == 0:
            return u
        return np.fft.ifft(np.fft.fft(u) * self.transfer(z))


@dataclass
class RingStrip:
    """Unfolded round trip of the self-intersecting ring, as two (possibly different) half trips.

    Each half trip runs from the crossing X to the graded mirror G (path p) and back, and is
    symmetric about G: an element at path z < p from X is met again at 2p - z. Elements are
    given as (z_um, half_width_um, phase) with phase None (plane fold mirror) or an array
    (reflection phase, e.g. cylindrical curvature). Without explicit elements, one plane fold
    mirror at d1 is used.
    """

    grid: Grid1D
    p_um: float  # crossing -> graded mirror path
    d1_um: float  # crossing -> fold mirror path (used when elements is None)
    phase_g1: np.ndarray
    phase_g2: np.ndarray
    ap_g_um: float = 5000.0  # half-width of the graded mirrors
    ap_fold_um: float = 5000.0  # half-width of the fold mirror at d1 (used when elements is None)
    amp_g1: np.ndarray | None = None
    amp_g2: np.ndarray | None = None
    extra_planes: list = field(default_factory=list)  # [(z_um, half_width_um)] apertures, NOT mirrored
    phase_fold: np.ndarray | None = None  # phase of the fold at d1 (used when elements is None)
    elements: list | None = None  # [(z_um, half_width_um, phase or None)], mirrored about G

    def __post_init__(self):
        x = self.grid.x
        self._ap_g = (np.abs(x) <= self.ap_g_um).astype(float)
        if self.elements is None:
            self.elements = [(self.d1_um, self.ap_fold_um, self.phase_fold)]
        stops = []
        for z, hw, ph in self.elements:
            t = (np.abs(x) <= hw).astype(complex)
            if ph is not None:
                t = t * np.exp(1j * ph)
            stops += [(z, t), (2 * self.p_um - z, t)]
        stops += [(z, (np.abs(x) <= hw).astype(complex)) for z, hw in self.extra_planes]
        self._stops = sorted(stops, key=lambda s: s[0])
        # Each half trip meets its own physical mirrors: separate stop lists, initially identical.
        self._stops_half = {1: list(self._stops), 2: list(self._stops)}

    def stop_positions(self):
        """Path positions (from X) of the stops of a half trip; index i refers to a distinct mirror."""
        return [z for z, _ in self._stops]

    def perturbed(self, changes=(), g_extra=None):
        """Copy with extra reflection phases: changes = [(which_half, stop_index, phase_array)];
        g_extra = {which_half: phase_array} added to the graded mirrors."""
        r = copy.copy(self)
        r._stops_half = {w: list(s) for w, s in self._stops_half.items()}
        for which, i, extra in changes:
            z, t = r._stops_half[which][i]
            r._stops_half[which][i] = (z, t * np.exp(1j * extra))
        if g_extra:
            r.phase_g1 = self.phase_g1 + g_extra.get(1, 0.0)
            r.phase_g2 = self.phase_g2 + g_extra.get(2, 0.0)
        return r

    def _segment(self, u, z0, z1, which=1):
        """Propagate from path position z0 to z1 (measured from X along half trip `which`),
        applying every element crossed on the way."""
        z = z0
        for zs, t in self._stops_half[which]:
            if z0 < zs <= z1:
                u = self.grid.prop(u, zs - z) * t
                z = zs
        return self.grid.prop(u, z1 - z)

    def _segment_adjoint(self, v, z0, z1, which=1):
        """Adjoint of _segment(., z0, z1, which): maps a field at z1 back to z0."""
        inside = [(zs, t) for zs, t in self._stops_half[which] if z0 < zs <= z1]
        z = z1
        for zs, t in reversed(inside):
            v = self.grid.prop(v, -(z - zs)) * np.conj(t)
            z = zs
        return self.grid.prop(v, -(z - z0))

    def to_mirror(self, u, which=1):
        """Field arriving at the graded mirror of half trip `which` from X."""
        return self._segment(u, 0.0, self.p_um, which)

    def from_mirror_adjoint(self, v, which=1):
        """Adjoint of the return path G -> X of half trip `which`, applied to a field at X."""
        return self._segment_adjoint(v, self.p_um, 2 * self.p_um, which)

    def design_phase(self, target_in, target_out, which=1, floor=1e-4):
        """Phase-only graded mirror that best maps target_in (at X, start of the half trip) onto
        target_out (at X, end of it): arg(K_out^+ T_out) - arg(K_in T_in). For a symmetric,
        unperturbed half trip and T_in = T_out = F real this is -2 arg(K F), the exact design."""
        a = self.to_mirror(np.asarray(target_in, dtype=complex), which)
        b = self.from_mirror_adjoint(np.asarray(target_out, dtype=complex), which)
        ph = np.angle(b * np.conj(a))
        amp = np.abs(a)
        ok = amp >= floor * amp.max()
        ph = np.unwrap(np.where(ok, ph, 0.0))
        idx = np.where(ok)[0]
        ph[: idx[0]] = ph[idx[0]]
        ph[idx[-1] + 1 :] = ph[idx[-1]]
        return ph - ph[self.grid.n // 2]

    def half(self, u, which=1):
        ph = self.phase_g1 if which == 1 else self.phase_g2
        am = self.amp_g1 if which == 1 else self.amp_g2
        u = self._segment(u, 0.0, self.p_um, which)
        u = u * self._ap_g * np.exp(1j * ph)
        if am is not None:
            u = u * am
        return self._segment(u, self.p_um, 2 * self.p_um, which)

    def roundtrip(self, u):
        return self.half(self.half(u, 1), 2)

    def operator(self):
        n = self.grid.n
        return spla.LinearOperator((n, n), matvec=lambda v: self.roundtrip(np.asarray(v).ravel()), dtype=complex)


def telescope_elements(grid: Grid1D, f1_um, f2_um, d0_um, ap1_um, ap2_um):
    """Two cylindrical curved folds forming a demagnifying telescope (magnification f2/f1):
    first mirror at path d0 from X (focal length f1), second at d0 + f1 + f2 (focal length f2)."""
    return [(d0_um, ap1_um, spherical_phase(grid, 2 * f1_um)), (d0_um + f1_um + f2_um, ap2_um, spherical_phase(grid, 2 * f2_um))]


# ---------------------------------------------------------------- designs
def graded_phase(grid: Grid1D, target, p_um, floor=1e-4, d1_um=None, phase_fold=None, ring=None):
    """Conjugating mirror phase -2 arg(K F), K the path X -> G: free space over p, or with a fold
    phase at d1, or the elements of a given RingStrip; held constant where |K F| < floor x max."""
    if ring is not None:
        u = ring.to_mirror(target.astype(complex))
    elif phase_fold is None:
        u = grid.prop(target.astype(complex), p_um)
    else:
        u = grid.prop(grid.prop(target.astype(complex), d1_um) * np.exp(1j * phase_fold), p_um - d1_um)
    a = np.abs(u)
    ph = -2.0 * np.unwrap(np.angle(u))
    ok = a >= floor * a.max()
    idx = np.where(ok)[0]
    ph[: idx[0]] = ph[idx[0]]
    ph[idx[-1] + 1 :] = ph[idx[-1]]
    return ph - ph[grid.n // 2], ok


def spherical_phase(grid: Grid1D, radius_um):
    """Reflection phase of a concave mirror of radius R (sag x^2 / 2R): -k x^2 / R."""
    k = 2 * np.pi / grid.lam_um
    return -k * grid.x**2 / radius_um


def gaussian_mode_radius(w0_um, p_um, lam_um):
    """Mirror radius of a plane-concave cavity of length p whose waist (at the plane) is w0."""
    zr = np.pi * w0_um**2 / lam_um
    return p_um + zr**2 / p_um


def flat_target(grid: Grid1D, w_um, order):
    return supergauss(grid.x, w_um, order).astype(complex)


# ---------------------------------------------------------------- perturbations
def smooth_screen(grid: Grid1D, rms, corr_um, rng, rms_halfwidth_um):
    """Random phase screen with Gaussian correlation length corr_um, rms over |x| <= rms_halfwidth_um."""
    white = rng.standard_normal(grid.n)
    kx = 2 * np.pi * np.fft.fftfreq(grid.n, d=grid.dx)
    s = np.real(np.fft.ifft(np.fft.fft(white) * np.exp(-((kx * corr_um) ** 2) / 8.0)))
    sel = np.abs(grid.x) <= rms_halfwidth_um
    s = s - s[sel].mean()
    return s * rms / s[sel].std()


def quantise(phase, levels):
    step = 2 * np.pi / levels
    return np.round(phase / step) * step


# ---------------------------------------------------------------- solvers
def eigenmodes(ring: RingStrip, k=12, v0=None, tol=1e-10):
    vals, vecs = spla.eigs(ring.operator(), k=k, which="LM", v0=v0, tol=tol, maxiter=20000)
    order = np.argsort(-np.abs(vals))
    return vals[order], vecs[:, order]


def fundamental(ring: RingStrip, target, k=12):
    """Eigenmode with the largest overlap with the target profile."""
    vals, vecs = eigenmodes(ring, k=k, v0=target)
    ov = np.abs(vecs.conj().T @ target) ** 2 / (np.vdot(target, target).real * np.sum(np.abs(vecs) ** 2, axis=0))
    i = int(np.argmax(ov))
    return vals[i], vecs[:, i], vals, vecs, ov


def driven(ring: RingStrip, e_in, t_in, loss_rt, psi, x0=None, tol=1e-9):
    r = np.sqrt((1 - t_in) * (1 - loss_rt))
    g = r * np.exp(1j * psi)
    n = ring.grid.n
    A = spla.LinearOperator((n, n), matvec=lambda v: np.asarray(v).ravel() - g * ring.roundtrip(np.asarray(v).ravel()), dtype=complex)
    b = np.sqrt(t_in) * e_in
    e, info = spla.gmres(A, b, x0=x0, rtol=tol, restart=200, maxiter=200)
    return e, info


def lock_and_solve(ring: RingStrip, e_in, t_in, loss_rt, psi0, span=None, iters=40, recentre=8):
    """Golden-section search of psi around psi0 maximising the circulating power."""
    r = np.sqrt((1 - t_in) * (1 - loss_rt))
    lw = 2 * (1 - r) / np.sqrt(r)  # full width (rad) of the resonance in psi
    span = 3 * lw if span is None else span
    pin = np.vdot(e_in, e_in).real
    cache = {}

    def power(psi):
        e, _ = driven(ring, e_in, t_in, loss_rt, psi, x0=cache.get("x"))
        cache["x"] = e
        return np.vdot(e, e).real / pin, e

    def golden(centre):
        a, b = centre - span, centre + span
        gr = (np.sqrt(5) - 1) / 2
        c, d = b - gr * (b - a), a + gr * (b - a)
        pc, pd = power(c)[0], power(d)[0]
        for _ in range(iters):
            if pc > pd:
                b, d, pd = d, c, pc
                c = b - gr * (b - a)
                pc = power(c)[0]
            else:
                a, c, pc = c, d, pd
                d = a + gr * (b - a)
                pd = power(d)[0]
            if b - a < 1e-3 * lw:
                break
        return 0.5 * (a + b)

    # Self-consistent re-centring: a perturbation can shift the fundamental's round-trip phase by
    # more than the search span. Lock to the Rayleigh-quotient phase of the circulating field,
    # iterate, then refine by golden section; keep the better of this and the direct search.
    psi = psi0
    for _ in range(recentre):
        _, e = power(psi)
        lam_e = np.vdot(e, ring.roundtrip(e)) / np.vdot(e, e)
        new = -np.angle(lam_e)
        new = psi + np.angle(np.exp(1j * (new - psi)))
        if abs(new - psi) < 0.05 * lw:
            psi = new
            break
        psi = new
    best = None
    for centre in {round(psi0, 12), round(psi, 12)}:
        p_ = golden(centre)
        g_, e_ = power(p_)
        if best is None or g_ > best[1]:
            best = (p_, g_, e_)
    psi, g, e = best
    return {"psi": psi, "enhancement": g, "field": e, "linewidth_rad": lw}


def flatness(e, x, half_um):
    """Intensity statistics of a field over |x| <= half_um, relative to its mean there."""
    sel = np.abs(x) <= half_um
    i = np.abs(e[sel]) ** 2
    return {"rms": float(i.std() / i.mean()), "min_over_mean": float(i.min() / i.mean()), "max_over_mean": float(i.max() / i.mean())}


def phase_gradient_stats(phase, x, sel, lam_um):
    """Largest local deflection |dphi/dx| / k of a mirror phase (reflection deflects by this)."""
    k = 2 * np.pi / lam_um
    g = np.gradient(phase, x)[sel] / k
    return {"pv_rad": float(phase[sel].max() - phase[sel].min()), "max_deflection": float(np.abs(g).max())}


def dense_spectrum(ring: RingStrip, target, max_loss=0.01, n_modes=12):
    """Full eigen-decomposition of the round trip (small grids only, n <= ~4096).

    Returns the fundamental (largest overlap with the target) and the low-loss modes ordered by
    their round-trip phase separation |arg(lambda_m / lambda_0)| from it. Eigenvalues with loss
    above max_loss are discarded (they include the grid's unphysical edge modes).
    """
    n = ring.grid.n
    M = np.empty((n, n), dtype=complex)
    for j in range(n):
        e = np.zeros(n, dtype=complex)
        e[j] = 1.0
        M[:, j] = ring.roundtrip(e)
    vals, vecs = np.linalg.eig(M)
    ov = np.abs(vecs.conj().T @ target) ** 2 / (np.vdot(target, target).real * np.sum(np.abs(vecs) ** 2, axis=0))
    i0 = int(np.argmax(ov))
    d = np.angle(vals / vals[i0])
    loss = 1 - np.abs(vals) ** 2
    order = [j for j in np.argsort(np.abs(d)) if j != i0 and loss[j] < max_loss][:n_modes]
    lowest = [j for j in np.argsort(loss) if j != i0][:n_modes]
    return {
        "fundamental_loss": float(loss[i0]),
        "fundamental_overlap": float(ov[i0]),
        "modes": [{"dphase_rad": float(d[j]), "loss": float(loss[j])} for j in order],
        "lowest_loss_modes": [{"dphase_rad": float(d[j]), "loss": float(loss[j])} for j in lowest],
    }
