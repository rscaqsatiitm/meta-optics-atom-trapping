"""Two-parameter, C4v-symmetric meta-atoms for independent phase control at two wavelengths.

Families (square lattice of pitch p, layer height h, all polarisation independent at normal
incidence because each unit cell has the full square symmetry):

  ring: square ring of outer width wo and square hole wh (wh = 0: solid square pillar).
  duo:  two square pillars per cell, w1 at the cell centre and w2 at the cell corners
        (a centred-square arrangement of two sizes; nearest-neighbour distance p / sqrt(2)).

A target pair (theta_a, theta_b) of wrapped phases at the two wavelengths is realised by the
atom maximising Re[t_a exp(-i theta_a)] + Re[t_b exp(-i theta_b)]. The two-wavelength focusing
efficiency of a library is, at each wavelength, |<t(atom(target)) exp(-i theta)>|^2 (times the
power normalisation n_out / n_in) averaged over target pairs uniformly distributed on the
phase torus, which is what an aplanatic doublet at 852 and 1064 nm requires
(docs/metalens_report.md, Section 5).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


def _grid(pitch: float, n_grid: int, subpixel: int):
    ns = n_grid * subpixel
    x = (np.arange(ns) + 0.5) / ns * pitch - pitch / 2
    return np.meshgrid(x, x, indexing="ij")


def _square(X, Y, cx, cy, w, pitch):
    """Square of side w centred at (cx, cy), periodic in the cell."""
    dx = np.abs((X - cx + pitch / 2) % pitch - pitch / 2)
    dy = np.abs((Y - cy + pitch / 2) % pitch - pitch / 2)
    return (dx <= w / 2) & (dy <= w / 2)


def masks(family: str, params, pitch: float, n_grid: int, subpixel: int = 4) -> np.ndarray:
    """Fill fractions (B, n_grid, n_grid) for a list of (a, b) parameter pairs."""
    X, Y = _grid(pitch, n_grid, subpixel)
    out = np.empty((len(params), n_grid, n_grid))
    for i, (a, b) in enumerate(params):
        if family == "ring":
            m = _square(X, Y, 0, 0, a, pitch) & ~(_square(X, Y, 0, 0, b, pitch) if b > 0 else np.zeros_like(X, bool))
        elif family == "duo":
            m = (_square(X, Y, 0, 0, a, pitch) if a > 0 else np.zeros_like(X, bool)) | (
                _square(X, Y, pitch / 2, pitch / 2, b, pitch) if b > 0 else np.zeros_like(X, bool))
        else:
            raise ValueError(family)
        out[i] = m.reshape(n_grid, subpixel, n_grid, subpixel).mean(axis=(1, 3))
    return out


@dataclass(frozen=True)
class Rules:
    min_feature_um: float = 0.10  # smallest pillar, wall, hole or gap
    max_aspect: float = 10.0  # height / smallest solid feature (pillar width or ring wall)
    step_um: float = 0.01


def candidates(family: str, pitch: float, height: float, rules: Rules):
    """Fabricable (a, b) pairs on the parameter grid."""
    st, mf = rules.step_um, rules.min_feature_um
    solid_min = max(mf, height / rules.max_aspect)
    out = []
    if family == "ring":
        for wo in np.arange(solid_min, pitch - mf + 1e-9, st):
            out.append((round(wo, 4), 0.0))
            for wh in np.arange(mf, wo - 2 * solid_min + 1e-9, st):
                out.append((round(wo, 4), round(wh, 4)))
    elif family == "duo":
        # Squares at the centre (w1) and at the corners (w2) are offset by (p/2, p/2): their
        # facing corners are separated by sqrt(2) g with g = p/2 - (w1 + w2)/2.
        ws = [0.0] + list(np.arange(solid_min, pitch - mf + 1e-9, st))
        for w1 in ws:
            for w2 in ws:
                if w1 == 0 and w2 == 0:
                    continue
                if w1 > 0 and w2 > 0 and np.sqrt(2) * (pitch / 2 - (w1 + w2) / 2) < mf:
                    continue
                out.append((round(w1, 4), round(w2, 4)))
    else:
        raise ValueError(family)
    return out


def dual_efficiency(ta: np.ndarray, tb: np.ndarray, n_ratio_a: float = 1.0, n_ratio_b: float = 1.0, bins: int = 48, wa: float = 1.0, wb: float = 1.0):
    """Two-wavelength focusing efficiency of a library (complex transmissions ta, tb per atom).

    Returns (eta_a, eta_b, rms phase error a, rms phase error b, index of chosen atom per target)."""
    th = 2 * np.pi * (np.arange(bins) + 0.5) / bins
    A, B = np.meshgrid(th, th, indexing="ij")
    ea = np.exp(-1j * A.ravel())
    eb = np.exp(-1j * B.ravel())
    score = wa * np.real(ta[None, :] * ea[:, None]) + wb * np.real(tb[None, :] * eb[:, None])
    k = np.argmax(score, axis=1)
    ca = ta[k] * ea
    cb = tb[k] * eb
    eta_a = n_ratio_a * float(np.abs(ca.mean()) ** 2)
    eta_b = n_ratio_b * float(np.abs(cb.mean()) ** 2)
    # Phase errors relative to the mean (common) offset.
    pa = np.angle(ca * np.conj(ca.mean()))
    pb = np.angle(cb * np.conj(cb.mean()))
    return eta_a, eta_b, float(np.sqrt(np.mean(pa**2))), float(np.sqrt(np.mean(pb**2))), k


def offset_scan(ta, tb, n_ratio_a=1.0, n_ratio_b=1.0, bins=48, n_off=12):
    """dual_efficiency is offset-free by construction (targets cover the whole torus); this
    returns it for a few weightings between the wavelengths, to expose the trade-off."""
    rows = []
    for wb in (0.5, 1.0, 2.0):
        ea, eb, sa, sb, _ = dual_efficiency(ta, tb, n_ratio_a, n_ratio_b, bins, 1.0, wb)
        rows.append({"weight_b": wb, "eta_a": ea, "eta_b": eb, "rms_a": sa, "rms_b": sb})
    return rows
