"""Block-wise propagation of very large apertures to a small focal window (folding).

The focal field inside a window of side W = N_f dx is all that is needed for the
tweezer metrics. Sampling the spectrum of the focal field at df = 1 / W is equivalent
to replacing the focal field by its periodisation with period W, and sampling the
metasurface spectrum at df is equivalent to summing (folding) the metasurface field
modulo W:

    S(m df) = sum_r E(r) exp(-i 2 pi m df . r) = FFT_{N_f}[ sum_k E(u + k W) ].

Hence  E_focal(x) (periodised) = IFFT_{N_f}[ FFT_{N_f}(E_fold) H ],  with H the
transfer function of the layer stack sampled at df. The result is exact for points
whose true focal field has negligible copies from |x + k W| (k != 0), i.e. it aliases
only the stray light that falls more than W/2 from the axis. The metasurface field is
generated and folded block by block, so the memory needed is that of one block plus the
N_f x N_f buffer, independent of the aperture size. Validation against a full-grid FFT:
scripts/validate_folded.py.
"""

from __future__ import annotations

import numpy as np
import scipy.fft as sfft

from .angular_spectrum import WORKERS, transfer_function
from .layers import LayerStack


def fold_blocks(blocks, n: int, n_fold: int):
    """Fold (row_start, block[, incident]) pieces of an n x n field (axis at n//2) modulo n_fold.

    Returns (folded field with the axis at index 0, transmitted power sum |E|^2 per pixel
    area); if the pieces carry an incident-power term, a third item with its sum is returned.
    n must be a multiple of n_fold.
    """
    if n % n_fold:
        raise ValueError("n must be a multiple of n_fold")
    F = np.zeros((n_fold, n_fold), np.complex128)
    power = 0.0
    incident = None
    c = n // 2
    for item in blocks:
        r0, blk = item[0], item[1]
        if len(item) > 2:
            incident = (incident or 0.0) + item[2]
        b = blk.shape[0]
        power += float(np.sum(np.abs(blk) ** 2))
        cols = np.roll(blk, -c, axis=1).astype(np.complex128)  # column c -> 0; n % n_fold == 0
        cols = cols.reshape(b, n // n_fold, n_fold).sum(axis=1)
        rows = (np.arange(r0, r0 + b) - c) % n_fold
        np.add.at(F, rows, cols)
    if incident is not None:
        return F, power, incident
    return F, power


def propagate_folded(F: np.ndarray, dx: float, lam_um: float, stack: LayerStack, fresnel: bool = False) -> np.ndarray:
    """Focal field from a folded metasurface field; returned with the axis at index N//2."""
    n = F.shape[0]
    spec = sfft.fft2(F, workers=WORKERS)
    h = transfer_function(n, dx, lam_um, stack, fresnel=fresnel).astype(np.complex128)
    h *= np.exp(-1j * np.angle(h[0, 0]))
    out = sfft.ifft2(spec * h, workers=WORKERS)
    return np.fft.fftshift(out)


def transfer_function_tilted(n: int, dx: float, lam_um: float, stack: LayerStack, tilt_layer: int, tilt_rad: float,
                             extra_vacuum_um: float = 0.0) -> np.ndarray:
    """Transfer function of the stack with one plane-parallel layer tilted about the y axis.

    A plane wave of direction k = (sx, sy, sz) crossing a slab of thickness d and index n_t
    whose normal N = (sin b, 0, cos b) makes angle t with k gains, relative to the vacuum it
    replaces, the phase k0 d [sqrt(n_t^2 - sin^2 t) - cos t] (exact for plane waves; Fresnel
    losses neglected). All other layers must be vacuum/air (n = 1), as in the cell stacks.
    The overall translation of the slab does not matter (plane-wave transfer functions are
    translation invariant). ``extra_vacuum_um`` adds free propagation for focus scans.
    """
    k0 = 2 * np.pi / lam_um
    f = np.fft.fftfreq(n, d=dx)
    SX, SY = np.meshgrid(lam_um * f, lam_um * f, indexing="ij")
    s2 = SX**2 + SY**2
    valid = s2 < 1.0
    sz = np.sqrt(np.clip(1.0 - s2, 0.0, None))
    total = float(np.sum(stack.thicknesses())) + extra_vacuum_um
    lay = stack.layers[tilt_layer]
    d, nt = lay.thickness_um, lay.n(lam_um)
    for i, l in enumerate(stack.layers):
        if i != tilt_layer and abs(l.n(lam_um) - 1.0) > 1e-6:
            raise ValueError("only one non-vacuum layer is supported")
    cos_t = SX * np.sin(tilt_rad) + sz * np.cos(tilt_rad)
    slab = d * (np.sqrt(np.clip(nt**2 - (1.0 - cos_t**2), 0.0, None)) - cos_t)
    phase = k0 * (total * sz + slab)
    h = np.where(valid, np.exp(1j * phase), 0.0)
    return h * np.exp(-1j * np.angle(h[0, 0]))
