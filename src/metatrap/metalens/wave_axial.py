"""Axisymmetric scalar wave propagation through the doublet (on-axis field point).

The field is propagated plane to plane with the quasi-discrete Hankel transform (QDHT;
M. Guizar-Sicairos and J. C. Gutierrez-Vega, J. Opt. Soc. Am. A 21, 53 (2004), recalled,
unverified): a radially symmetric field f(r) on [0, R] with spectrum band-limited to
[0, V] is sampled at r_k = j_k R / S and nu_k = j_k / (2 pi R), with j_k the zeros of J0
and S = j_{N+1} = 2 pi R V. The transform matrix T_ij = 2 J0(j_i j_j / S) / (|J1(j_i)|
|J1(j_j)| S) is symmetric and orthogonal. Each homogeneous layer multiplies the angular
spectrum by exp(i k0 d sqrt(n^2 - (lambda nu)^2)); evanescent components decay. Interface
Fresnel coefficients are not included (scalar model).

Sequence: top-hat beam (soft edge) at the corrector, times exp(i phi1); air spacer; times
exp(i phi2); cell stack; on-axis field at distance z beyond the nominal atom plane from the
plane-wave sum. The Strehl ratio at z is |sum c_i e^{i psi_i}|^2 / (sum |c_i|)^2 over the
spectral samples (the field of an in-phase spectrum of the same amplitudes is the
diffraction-limited reference), maximised over z.

This is independent of the ray tracer: it tests the phase profiles, the sampling and the
stationary-phase (ray) assumptions behind the ray-traced Strehl ratios.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import j0, j1

from .. import materials


def bessel_j0_zeros(n: int) -> np.ndarray:
    """First n positive zeros of J0 (McMahon start, Newton refinement)."""
    k = np.arange(1, n + 1, dtype=float)
    b = (k - 0.25) * np.pi
    x = b + 1 / (8 * b) - 124 / (3 * (8 * b) ** 3)
    for _ in range(6):
        x = x + j0(x) / j1(x)  # d J0/dx = -J1
    return x


class QDHT:
    def __init__(self, R: float, V: float):
        N = int(np.ceil(2 * R * V))
        z = bessel_j0_zeros(N + 1)
        self.S = z[-1]
        self.j = z[:-1]
        self.R, self.V, self.N = R, self.S / (2 * np.pi * R), N
        self.r = self.j * R / self.S
        self.nu = self.j / (2 * np.pi * R)
        self.J1 = np.abs(j1(self.j))
        T = np.empty((N, N))
        step = 2000
        for a in range(0, N, step):
            T[a : a + step] = 2 * j0(np.outer(self.j[a : a + step], self.j) / self.S) / (self.J1[a : a + step, None] * self.J1[None, :] * self.S)
        self.T = T

    def forward(self, f):
        g = f * self.R / self.J1
        return (self.T @ g.real + 1j * (self.T @ g.imag)) * self.J1 / self.V

    def inverse(self, F):
        g = F * self.V / self.J1
        return (self.T @ g.real + 1j * (self.T @ g.imag)) * self.J1 / self.R


def _kz(nu, lam_um, n):
    k0 = 2 * np.pi / lam_um
    return k0 * np.sqrt((n**2 - (lam_um * nu) ** 2).astype(complex))


def on_axis_focus(d, lam_um: float, r_max_um: float = 9000.0, v_max_over_k: float = 0.7, edge_um: float = 5.0, m1=None, m2=None,
                  z_span_um: float = 30.0):
    """Strehl ratio and focal position of the on-axis focus (d: doublet.Doublet)."""
    m1, m2 = m1 or d.m1, m2 or d.m2
    R1 = d.info["R1_um"]
    q = QDHT(r_max_um, v_max_over_k / lam_um)
    r, nu = q.r, q.nu
    amp = 0.5 * (1 - np.tanh((r - R1) / edge_um)) if edge_um > 0 else (r <= R1).astype(float)
    u = amp * np.exp(1j * m1.phase(r, lam_um))
    for t, m in d.spec.spacer:
        u = q.inverse(q.forward(u) * np.exp(1j * _kz(nu, lam_um, float(materials.get(m).n(lam_um))) * float(t)))
    u = u * np.exp(1j * m2.phase(r, lam_um))
    # Energy beyond the computational window after the spacer (aliasing guard).
    spill = float(np.sum(np.abs(u[r > 0.98 * r_max_um]) ** 2 * r[r > 0.98 * r_max_um]) / np.sum(np.abs(u) ** 2 * r))
    F = q.forward(u)
    H = np.ones_like(nu, dtype=complex)
    for t, m in d.spec.cell:
        H = H * np.exp(1j * _kz(nu, lam_um, float(materials.get(m).n(lam_um))) * float(t))
    kz_v = _kz(nu, lam_um, 1.0)
    prop = np.abs(lam_um * nu) < 1.0
    w = 1.0 / (np.pi * q.R**2 * q.J1**2)  # Fourier-Bessel weights: f(r) = sum F(nu_m) w_m J0(2 pi nu_m r)
    c = (F * H * w)[prop]
    kzp = kz_v[prop].real
    # Fraction of the spectrum (power) beyond 1.1 x the design NA (stray, mis-directed light).
    pw = np.abs(F) ** 2 * w
    frac_high = float(pw[np.abs(lam_um * nu) > 1.1 * d.spec.na].sum() / pw.sum())

    def strehl(z):
        return float(np.abs(np.sum(c * np.exp(1j * kzp * z))) ** 2 / np.sum(np.abs(c)) ** 2)

    zs = np.linspace(-z_span_um, z_span_um, 601)
    s = np.array([strehl(z) for z in zs])
    z0 = zs[np.argmax(s)]
    best = minimize_scalar(lambda z: -strehl(z), bounds=(z0 - 0.2, z0 + 0.2), method="bounded", options={"xatol": 1e-4})
    zb = float(best.x)
    # Radial intensity profile in the best plane (direct inverse sum), and its 1/e^2 radius.
    rr = np.linspace(0, 2.0, 201)
    cf = (F * H * np.exp(1j * kz_v * zb))[prop]
    prof = np.array([np.abs(np.sum(cf * w[prop] * j0(2 * np.pi * nu[prop] * x))) ** 2 for x in rr])
    prof = prof / prof[0]
    i_e2 = np.argmax(prof < np.exp(-2))
    return {"N": int(q.N), "strehl": strehl(zb), "best_focus_um": zb, "strehl_at_nominal_plane": strehl(0.0), "spill_fraction": spill,
            "spectrum_fraction_beyond_1.1NA": frac_high, "waist_1e2_um": float(rr[i_e2]), "profile_r_um": rr.tolist(), "profile": prof.tolist()}
