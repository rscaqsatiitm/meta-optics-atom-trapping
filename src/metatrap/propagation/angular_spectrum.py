"""Scalar angular-spectrum propagation through plane-parallel layers.

For a field U(x, y) sampled on an N x N grid of spacing dx (um), propagation through
a layer of index n and thickness d multiplies its spectrum by exp(i kz d), with
kz = sqrt((n k0)^2 - kx^2 - ky^2); evanescent components are discarded. A stack of
layers is a product of such transfer functions, which is exact for the scalar field
in plane-parallel media. Fresnel transmission at the interfaces is optional (scalar
average of the s and p coefficients); it changes the overall power by a few per cent
and slightly apodises the high-angle part of the pupil.

Arrays are complex64 and transforms use scipy.fft with multithreading, so a
12k x 12k field (1.2 GB) fits in the 15 GB of the development container.
"""

from __future__ import annotations

import os

import numpy as np
import scipy.fft as sfft

from .layers import LayerStack

WORKERS = int(os.environ.get("METATRAP_FFT_WORKERS", os.cpu_count() or 1))


def freq_grid(n: int, dx: float):
    f = sfft.fftfreq(n, d=dx).astype(np.float32)
    return f  # cycles per um


def _fresnel_scalar(n1, n2, cos1, cos2):
    ts = 2 * n1 * cos1 / (n1 * cos1 + n2 * cos2)
    tp = 2 * n1 * cos1 / (n2 * cos1 + n1 * cos2)
    # Power transmission factor for each polarisation, averaged (unpolarised-like weighting).
    Ts = (n2 * cos2) / (n1 * cos1) * np.abs(ts) ** 2
    Tp = (n2 * cos2) / (n1 * cos1) * np.abs(tp) ** 2
    return np.sqrt(0.5 * (Ts + Tp))


def transfer_function(n: int, dx: float, lam_um: float, stack: LayerStack, n_start: float = 1.0, fresnel: bool = False, block_rows=None):
    """Transfer function H(fx, fy) of the stack, returned as complex64 (n x n), or a row block."""
    fx = freq_grid(n, dx)
    rows = fx if block_rows is None else fx[block_rows]
    fx2 = (rows[:, None] ** 2 + fx[None, :] ** 2).astype(np.float64)
    s2 = (lam_um**2) * fx2  # (direction cosine)^2 in vacuum units
    phase = np.zeros_like(fx2)
    amp = np.ones_like(fx2)
    valid = s2 < 1.0
    n_prev, cos_prev = n_start, np.sqrt(np.clip(1 - s2 / n_start**2, 0, None))
    k0 = 2 * np.pi / lam_um
    for layer in stack.layers:
        ni = layer.n(lam_um)
        arg = ni**2 - s2
        valid &= arg > 0
        kz = k0 * np.sqrt(np.clip(arg, 0, None))
        phase += kz * layer.thickness_um
        if fresnel:
            cos_i = np.sqrt(np.clip(1 - s2 / ni**2, 1e-12, None))
            amp *= _fresnel_scalar(n_prev, ni, np.clip(cos_prev, 1e-12, None), cos_i)
            n_prev, cos_prev = ni, cos_i
    h = np.where(valid, amp * np.exp(1j * phase), 0.0)
    return h.astype(np.complex64)


def propagate(u: np.ndarray, dx: float, lam_um: float, stack: LayerStack, fresnel: bool = False, remove_piston: bool = True) -> np.ndarray:
    """Propagate the field u (complex, N x N, spacing dx) through the stack. Overwrites u."""
    n = u.shape[0]
    assert u.shape == (n, n)
    u = u.astype(np.complex64, copy=False)
    spec = sfft.fft2(u, workers=WORKERS, overwrite_x=True)
    del u
    piston = 1.0
    if remove_piston:
        h00 = transfer_function(n, dx, lam_um, stack, fresnel=fresnel, block_rows=slice(0, 1))[0, 0]
        piston = np.complex64(np.exp(-1j * np.angle(h00)))
    block = max(1, (64 * 1024 * 1024) // (8 * n))  # ~64 MB of float64 temporaries per block
    for r0 in range(0, n, block):
        sl = slice(r0, min(n, r0 + block))
        spec[sl] *= transfer_function(n, dx, lam_um, stack, fresnel=fresnel, block_rows=sl) * piston
    return sfft.ifft2(spec, workers=WORKERS, overwrite_x=True)


def propagate_to_window(u: np.ndarray, dx: float, lam_um: float, stack: LayerStack, window: int, fresnel: bool = False) -> np.ndarray:
    """Propagate and return only the central window x window region of the output.

    Saves memory for focusing problems, where the focal pattern occupies a small
    region of the full aperture grid.
    """
    out = propagate(u, dx, lam_um, stack, fresnel=fresnel)
    n = out.shape[0]
    c0 = n // 2 - window // 2
    return np.array(out[c0 : c0 + window, c0 : c0 + window])


def grid_coords(n: int, dx: float) -> np.ndarray:
    """Pixel-centre coordinates with the optical axis on pixel n//2."""
    return (np.arange(n) - n // 2).astype(np.float64) * dx
