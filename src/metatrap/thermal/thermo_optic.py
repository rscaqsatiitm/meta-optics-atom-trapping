"""Effect of a thermal optical-path change in the metasurface substrate on the tweezers.

The substrate path change dOPL(r) (from axisym_fd) adds a phase k0 dOPL(r) at
metasurface radius r. Mapped to the pupil through the ray mapping of the layer stack,
it is fitted with piston and defocus (cos theta in the atom medium). The defocus
coefficient gives the axial shift dz of every focus; the residual gives the Strehl
ratio after refocusing. The Strehl ratio in the original focal plane (no refocus)
includes the defocus term. Tilt is absent by symmetry.

Scalar model; lateral (magnification) changes of the array are reported separately
from the thermal expansion of the metasurface: d(pitch)/pitch = alpha_exp dT_mean.
"""

from __future__ import annotations

import numpy as np

from ..hologram.field_model import pupil_samples
from ..propagation.layers import LayerStack


def wavefront_effect(stack: LayerStack, lam_um: float, na: float, illumination: dict, r_m, d_opl_m, n_pupil: int = 256) -> dict:
    k0 = 2 * np.pi / lam_um
    x, y, amp, _ = pupil_samples(stack, lam_um, na, n_pupil, illumination)
    r_um = np.hypot(x, y)
    W = np.interp(r_um * 1e-6, r_m, d_opl_m) * 1e6  # um
    s = stack.s_of_r(lam_um)(r_um)
    cz = np.sqrt(1 - s**2)
    w = amp**2
    B = np.column_stack([np.ones_like(cz), cz])
    coef = np.linalg.solve(B.T @ (B * w[:, None]), (B * w[:, None]).T @ W)
    res = W - B @ coef
    piston_only = W - np.sum(w * W) / np.sum(w)
    norm = np.sum(amp) ** 2
    return {
        "dz_um": float(-coef[1]),
        "strehl_refocused": float(np.abs(np.sum(amp * np.exp(1j * k0 * res))) ** 2 / norm),
        "strehl_fixed_plane": float(np.abs(np.sum(amp * np.exp(1j * k0 * piston_only))) ** 2 / norm),
        "rms_residual_waves": float(np.sqrt(np.sum(w * res**2) / np.sum(w)) / lam_um),
        "pv_opl_um": float(W.max() - W.min()),
    }
