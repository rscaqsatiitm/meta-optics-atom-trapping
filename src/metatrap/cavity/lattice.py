"""Large-period lattice from crossing sheet beams: geometry, flat-top profiles and site maps.

Geometry of one lattice axis (here the axis whose fringes run along y). Two beams cross in
the atom plane z = 0 with propagation directions (cos t, +sin t) (arm A) and (cos t, -sin t)
(arm B), t the crossing half-angle. Their field near the atoms is

    E = E_A(u_A, s_A) exp(i k s_A) + E_B(u_B, s_B) exp(i k s_B),

where s is the coordinate along each beam and u the transverse coordinate in the plane,

    s_A =  x cos t + y sin t,   u_A = -x sin t + y cos t,
    s_B =  x cos t - y sin t,   u_B =  x sin t + y cos t.

The interference term 2 |E_A| |E_B| cos(2 k y sin t + dphi), dphi = arg E_A - arg E_B, makes
fringes of period a = lambda / (2 sin t). The barrier between neighbouring sites is
proportional to |E_A| |E_B|; a site is displaced along y by

    delta_y = -a dphi / (2 pi)

from the ideal lattice. Hence a differential wavefront error dw = dphi / k between the two
beams moves a site by dw / (2 sin t), 3.7 dw at a = 4 um and 1075 nm.

Square array of side S centred on the crossing point: the transverse coordinate of each beam
spans |u| <= (S/2)(cos t + sin t), and so does the propagation coordinate s.
"""

from __future__ import annotations

import numpy as np
from math import gamma

from ..propagation.angular_spectrum import freq_grid


# ---------------------------------------------------------------- flat-top profiles
def supergauss(x, w, n):
    """Super-Gaussian amplitude exp(-|x/w|^n)."""
    return np.exp(-np.abs(np.asarray(x) / w) ** n)


def supergauss_w_for_flat(half_flat, n, tol=0.01):
    """Parameter w such that the intensity exp(-2|x/w|^n) is >= 1 - tol for |x| <= half_flat."""
    return half_flat / (-0.5 * np.log(1.0 - tol)) ** (1.0 / n)


def supergauss_equivalent_width(w, n):
    """Equivalent width int I dx / I_max of the intensity exp(-2|x/w|^n): 2 w Gamma(1 + 1/n) 2^(-1/n)."""
    return 2.0 * w * gamma(1.0 + 1.0 / n) * 2.0 ** (-1.0 / n)


def array_halfrange(side, theta):
    """Half-range of the transverse (and propagation) coordinate of a crossing beam over the array."""
    return 0.5 * side * (np.cos(theta) + np.sin(theta))


# ---------------------------------------------------------------- 1D propagation along a sheet
def propagate_1d(u, dx, lam, z):
    """Scalar angular-spectrum propagation of a 1D field (paraxial carrier removed)."""
    k = 2 * np.pi / lam
    fx = freq_grid(len(u), dx).astype(np.float64)
    kx = 2 * np.pi * fx
    kz2 = k**2 - kx**2
    h = np.where(kz2 > 0, np.exp(1j * (np.sqrt(np.clip(kz2, 0, None)) - k) * z), 0.0)
    return np.fft.ifft(np.fft.fft(u) * h)


def field_stack(u0, x, lam, s_planes):
    """Field of a 1D beam on the planes s (same units as x), from the field u0 at s = 0."""
    dx = x[1] - x[0]
    return np.array([propagate_1d(u0, dx, lam, s) for s in s_planes])


def _interp_field(stack, x, s_planes, u, s):
    """Bilinear interpolation of the complex field stack[s, x] at points (u, s)."""
    fx = np.interp(u, x, np.arange(len(x)))
    fs = np.interp(s, s_planes, np.arange(len(s_planes)))
    i0 = np.clip(np.floor(fx).astype(int), 0, len(x) - 2)
    j0 = np.clip(np.floor(fs).astype(int), 0, len(s_planes) - 2)
    ax, as_ = fx - i0, fs - j0
    f00 = stack[j0, i0]
    f01 = stack[j0, i0 + 1]
    f10 = stack[j0 + 1, i0]
    f11 = stack[j0 + 1, i0 + 1]
    return (1 - as_) * ((1 - ax) * f00 + ax * f01) + as_ * ((1 - ax) * f10 + ax * f11)


def sheet_factor(s, w_z, lam):
    """Amplitude and Gouy phase of a 1D Gaussian sheet (thickness waist w_z at s = 0) in z = 0."""
    zr = np.pi * w_z**2 / lam
    q = 1.0 + (s / zr) ** 2
    return q ** (-0.25), -0.5 * np.arctan(s / zr)


def site_map(field_a, field_b, x, lam, theta, side, n_side, w_z=None, n_planes=81, flip_b=False):
    """Barrier depth and site displacement over a square array from the two beams' fields at s = 0.

    field_a, field_b: complex 1D fields (transverse coordinate x) of arms A and B at the crossing
    plane; units of x, lam, side and w_z must agree (um here). flip_b: if the 1D coordinate of arm
    B is the mirror image of its physical left-hand transverse coordinate (odd number of in-plane
    reflections between the arms). Returns per-site arrays on an n_side x n_side grid spanning the array.
    """
    h = array_halfrange(side, theta) * 1.02
    s_planes = np.linspace(-h, h, n_planes)
    st_a = field_stack(field_a, x, lam, s_planes)
    st_b = field_stack(field_b, x, lam, s_planes)
    g = np.linspace(-side / 2, side / 2, n_side)
    X, Y = np.meshgrid(g, g, indexing="xy")
    ct, sn = np.cos(theta), np.sin(theta)
    s_a, u_a = X * ct + Y * sn, -X * sn + Y * ct
    s_b, u_b = X * ct - Y * sn, X * sn + Y * ct
    ea = _interp_field(st_a, x, s_planes, u_a, s_a)
    eb = _interp_field(st_b, x, s_planes, -u_b if flip_b else u_b, s_b)
    amp = np.abs(ea) * np.abs(eb)
    dphi = np.angle(ea * np.conj(eb))
    if w_z is not None:
        fa, ga = sheet_factor(s_a, w_z, lam)
        fb, gb = sheet_factor(s_b, w_z, lam)
        amp = amp * fa * fb
        dphi = dphi + (ga - gb)
    period = lam / (2 * sn)
    c = n_side // 2
    dphi = np.angle(np.exp(1j * (dphi - dphi[c, c])))
    dy = -period * dphi / (2 * np.pi)
    return {"x": X, "y": Y, "depth": amp / amp[c, c], "dy": dy, "period": period}


def affine_residual(X, Y, d):
    """Residual of a site-displacement map after removing offset and linear terms (pitch, rotation)."""
    A = np.column_stack([np.ones(X.size), X.ravel(), Y.ravel()])
    coef, *_ = np.linalg.lstsq(A, d.ravel(), rcond=None)
    return (d.ravel() - A @ coef).reshape(d.shape)


def map_metrics(m):
    d, dy = m["depth"], m["dy"]
    res = affine_residual(m["x"], m["y"], dy)
    return {
        "depth_rms": float(np.std(d) / np.mean(d)),
        "depth_min_over_mean": float(d.min() / d.mean()),
        "depth_max_over_mean": float(d.max() / d.mean()),
        "dy_rms_um": float(np.std(dy)),
        "dy_pv_um": float(dy.max() - dy.min()),
        "dy_affine_residual_rms_um": float(np.std(res)),
        "dy_affine_residual_pv_um": float(res.max() - res.min()),
    }


# ---------------------------------------------------------------- plane-wave lattice (site shape)
def plane_wave_intensity(points, families):
    """Intensity at points (..., 3) from mutually incoherent families of coherent plane waves.

    families: list of lists of (amplitude, k-vector (3,)) tuples; units: k in rad per unit of points.
    """
    pts = np.asarray(points)
    tot = np.zeros(pts.shape[:-1])
    for fam in families:
        e = np.zeros(pts.shape[:-1], dtype=complex)
        for amp, kv in fam:
            e += amp * np.exp(1j * pts @ np.asarray(kv))
        tot += np.abs(e) ** 2
    return tot
