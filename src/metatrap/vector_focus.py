"""Vector (Richards-Wolf / Debye) focal field of the metasurface pupil.

The focal field is the plane-wave superposition E(x) = sum_sigma A(sigma) p(sigma)
exp(i k0 [sigma . x + cos(theta) z]) over the pupil samples, with A the energy-conserving
pupil amplitude of the scalar model (hologram.wgs.pupil_amplitude) and p(sigma) the unit
polarisation of each plane wave. The scalar model uses p = 1. Two polarisation models
bracket the behaviour of a polarisation-insensitive metasurface illuminated with x
polarisation:

* 'rotation' (aplanatic lens): the field rotates in the meridional plane,
  p = (cos(t) cos^2(f) + sin^2(f), (cos(t) - 1) sin(f) cos(f), -sin(t) cos(f)).
* 'projection': the transmitted x-polarised field is projected onto the plane transverse
  to the ray, p = x_hat - (x_hat . k_hat) k_hat, normalised.

Here, t is the polar angle of the ray, f its azimuth, and k_hat its unit wavevector.

Vector light shift of an alkali ground state (Grimm et al., Adv. At. Mol. Opt. Phys. 42,
95 (2000), eq. 20, extended with counter-rotating terms): for circularity
C = Im(E* x E) . b_hat / |E|^2 along the quantisation axis b_hat,
U = U_s [1 + r_v g_F m_F C], with r_v = (V_2 - V_1) / (2 S_2 + S_1),
S_i = 1/(w_i - w) + 1/(w_i + w), V_i = 1/(w_i - w) - 1/(w_i + w). The vector part acts as
a fictitious magnetic field B_fict = r_v U_s C / mu_B, independent of m_F.
"""

from __future__ import annotations

import numpy as np

from .constants import CESIUM, AlkaliLines

MU_B_HZ_PER_G = 1.399_624_49e6  # Bohr magneton / h in Hz/G (CODATA 2018)
KB_HZ_PER_K = 2.083_661_912e10  # k_B / h in Hz/K


def vector_ratio(lam_nm: float, atom: AlkaliLines = CESIUM) -> float:
    """r_v = (V_2 - V_1) / (2 S_2 + S_1) (signed; negative for red detuning of both lines)."""
    w = 1.0 / lam_nm
    w1, w2 = 1.0 / atom.lambda_d1_nm, 1.0 / atom.lambda_d2_nm
    S1, S2 = 1 / (w1 - w) + 1 / (w1 + w), 1 / (w2 - w) + 1 / (w2 + w)
    V1, V2 = 1 / (w1 - w) - 1 / (w1 + w), 1 / (w2 - w) - 1 / (w2 + w)
    return (V2 - V1) / (2 * S2 + S1)


def polarisation(sx, sy, model: str):
    s = np.hypot(sx, sy)
    ct = np.sqrt(1 - s**2)
    with np.errstate(invalid="ignore", divide="ignore"):
        cf = np.where(s > 0, sx / np.where(s > 0, s, 1), 1.0)
        sf = np.where(s > 0, sy / np.where(s > 0, s, 1), 0.0)
    if model == "rotation":
        px = ct * cf**2 + sf**2
        py = (ct - 1) * sf * cf
        pz = -s * cf
    elif model == "projection":
        kx, ky, kz = sx, sy, ct
        px, py, pz = 1 - kx * kx, -kx * ky, -kx * kz
        nrm = np.sqrt(px**2 + py**2 + pz**2)
        px, py, pz = px / nrm, py / nrm, pz / nrm
    else:
        raise ValueError(model)
    return px, py, pz


def focal_field(sx, sy, amp, lam_um, points, model: str | None, block: int = 2048):
    """Field components at points (P, 3) in um; model None gives the scalar field."""
    k0 = 2 * np.pi / lam_um
    cz = np.sqrt(1 - sx**2 - sy**2)
    if model is None:
        comps = [amp]
    else:
        px, py, pz = polarisation(sx, sy, model)
        comps = [amp * px, amp * py, amp * pz]
    out = np.zeros((len(comps), len(points)), complex)
    W = np.stack(comps).astype(np.complex128)
    for j0 in range(0, len(points), block):
        p = points[j0 : j0 + block]
        ph = np.exp(1j * k0 * (np.outer(sx, p[:, 0]) + np.outer(sy, p[:, 1]) + np.outer(cz, p[:, 2])))
        out[:, j0 : j0 + block] = W @ ph
    return out


def curvature_waist(profile_1d, dx):
    """1/e^2 radius of the Gaussian with the same peak curvature (centre sample at the middle)."""
    c = len(profile_1d) // 2
    i0, im, ip = profile_1d[c], profile_1d[c - 1], profile_1d[c + 1]
    curv = (ip - 2 * i0 + im) / dx**2
    return float(np.sqrt(-4 * i0 / curv))
