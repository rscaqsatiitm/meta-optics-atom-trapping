"""Render a pupil-domain hologram design onto the metasurface pixel grid.

Each metasurface pixel at (x, y) (pitch = meta-atom pitch) receives the phase

    phi(x, y) = phi_lens(r) + psi(sigma(r)),   sigma = -r_hat s(r),

where phi_lens is the exact focusing phase through the layer stack and psi the
hologram phase designed on the pupil grid (bilinear interpolation of exp(i psi)).
The pixel transmission is then either ideal, exp(i phi), or the library value
t(w(phi)) of the meta-atom chosen by the lookup table (local-periodic approximation).
The incident illumination amplitude multiplies the result. Construction is done in
row blocks so that 10^8-pixel apertures fit in memory.
"""

from __future__ import annotations

import numpy as np
from scipy.ndimage import map_coordinates

from ..propagation.illumination import profile
from ..propagation.layers import LayerStack
from .exact import PupilSampling, RenderGrid


def metasurface_blocks(
    n: int,
    dx: float,
    stack: LayerStack,
    lam_um: float,
    na: float,
    ps: PupilSampling | None,
    psi: np.ndarray | None,
    illumination: dict,
    lookup=None,
    block: int = 256,
    extra_phase=None,
    holo: RenderGrid | None = None,
    pillar_errors: dict | None = None,
    with_incident: bool = False,
):
    """Yield (row_start, complex64 block) of the field just after the metasurface.

    With ``with_incident`` the generator yields (row_start, block, incident_power_sum),
    where the last item is sum |illumination amplitude|^2 over the block (per pixel area),
    so that efficiencies can be quoted relative to the light incident on the metasurface.

    The optical axis is at pixel n//2. ``ps``/``psi`` may be None for a plain lens.
    ``holo`` (a fine RenderGrid from exact.render_grid) takes precedence over
    ``ps``/``psi`` and should be used for arrays wider than a few tens of um.
    ``extra_phase(r_um)`` adds a radial phase if given.

    With ``lookup`` each pixel takes the library transmission of the pillar chosen for
    its phase (local-periodic approximation), scaled by sqrt(n_out/n_in) so that |E|^2
    is a power density. ``pillar_errors`` = {"sigma_um": s, "bias_um": b, "seed": k}
    adds a uniform width bias and independent Gaussian width noise per pixel before the
    transmission is looked up (requires ``lookup``).
    """
    k0 = 2 * np.pi / lam_um
    radius = stack.aperture_for_na(na, lam_um)
    s_of_r = stack.s_of_r(lam_um)
    opl_axis = stack.opl(0.0, lam_um)
    kw = {k: v for k, v in illumination.items() if k != "kind"}
    holo_re = holo_im = None
    if holo is not None:
        holo_re, holo_im = holo.field.real, holo.field.imag
        ds, c0 = holo.ds, holo.c0
    elif ps is not None:
        hg = ps.to_grid(np.exp(1j * psi).astype(np.complex128), fill=0.0)
        holo_re, holo_im = hg.real, hg.imag
        ds, c0 = ps.grid.ds, ps.grid.n // 2
    rng = np.random.default_rng(pillar_errors.get("seed", 0)) if pillar_errors else None
    amp_scale = np.sqrt(lookup.n_ratio) if lookup is not None else 1.0
    coords = (np.arange(n) - n // 2) * dx
    for r0 in range(0, n, block):
        xs = coords[r0 : r0 + block]
        X, Y = np.meshgrid(xs, coords, indexing="ij")
        R = np.hypot(X, Y)
        inside = R <= radius
        if not inside.any():
            continue
        Ri, Xi, Yi = R[inside], X[inside], Y[inside]
        s = s_of_r(Ri)
        phi = -k0 * (stack.opl(s, lam_um) - opl_axis)
        if extra_phase is not None:
            phi = phi + extra_phase(Ri)
        if holo_re is not None:
            with np.errstate(invalid="ignore", divide="ignore"):
                ux = np.where(Ri > 0, Xi / np.where(Ri > 0, Ri, 1), 0.0)
                uy = np.where(Ri > 0, Yi / np.where(Ri > 0, Ri, 1), 0.0)
            fi = (-ux * s) / ds + c0
            fj = (-uy * s) / ds + c0
            hr = map_coordinates(holo_re, [fi, fj], order=1, mode="constant")
            hi = map_coordinates(holo_im, [fi, fj], order=1, mode="constant")
            phi = phi + np.angle(hr + 1j * hi)
        if lookup is None:
            t = np.exp(1j * phi)
        elif pillar_errors:
            w = lookup.width_for_phase(phi) + pillar_errors.get("bias_um", 0.0)
            if pillar_errors.get("sigma_um", 0.0) > 0:
                w = w + pillar_errors["sigma_um"] * rng.standard_normal(w.shape)
            t = lookup.t_of_width(w) * amp_scale
        else:
            t = lookup.t_for_phase(phi) * amp_scale
        amp = profile(illumination["kind"], Ri, radius, **kw)
        blk = np.zeros(R.shape, np.complex64)
        blk[inside] = (amp * t).astype(np.complex64)
        if with_incident:
            yield r0, blk, float(np.sum(amp**2))
        else:
            yield r0, blk


def metasurface_field(n: int, dx: float, stack: LayerStack, lam_um: float, na: float, ps, psi, illumination: dict, **kw) -> np.ndarray:
    """Complex64 field (n x n) just after the metasurface; see metasurface_blocks."""
    out = np.zeros((n, n), np.complex64)
    for r0, blk in metasurface_blocks(n, dx, stack, lam_um, na, ps, psi, illumination, **kw):
        out[r0 : r0 + blk.shape[0]] = blk
    return out


def fourier_upsample(field: np.ndarray, factor: int) -> np.ndarray:
    """Band-limited interpolation of a square field by zero-padding its spectrum."""
    n = field.shape[0]
    spec = np.fft.fftshift(np.fft.fft2(field))
    m = n * factor
    big = np.zeros((m, m), complex)
    c = m // 2 - n // 2
    big[c : c + n, c : c + n] = spec
    return np.fft.ifft2(np.fft.ifftshift(big)) * factor**2
