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


def _hash_normal(i, j, seed: int, salt: int):
    """Deterministic standard-normal values for integer grid indices (i, j) (vectorised).

    Used for spatially correlated fabrication errors, which must be identical for every
    block of the metasurface that touches the same write field or noise cell.
    """
    def mix(x):
        x = (x ^ (x >> np.uint64(30))) * np.uint64(0xBF58476D1CE4E5B9)
        x = (x ^ (x >> np.uint64(27))) * np.uint64(0x94D049BB133111EB)
        return x ^ (x >> np.uint64(31))

    with np.errstate(over="ignore"):
        a = (np.asarray(i, np.int64).astype(np.uint64) * np.uint64(0x9E3779B97F4A7C15)
             + np.asarray(j, np.int64).astype(np.uint64) * np.uint64(0xC2B2AE3D27D4EB4F)
             + np.uint64(seed * 1_000_003 + salt * 7919))
        u1 = (mix(a) >> np.uint64(11)).astype(np.float64) / float(1 << 53)
        u2 = (mix(a + np.uint64(0x632BE59BD9B4E019)) >> np.uint64(11)).astype(np.float64) / float(1 << 53)
    u1 = np.clip(u1, 1e-300, 1.0)
    return np.sqrt(-2.0 * np.log(u1)) * np.cos(2 * np.pi * u2)


def _smooth_noise(x, y, length: float, seed: int, salt: int):
    """Unit-variance value noise with correlation length ~length (bilinear on a grid of that spacing)."""
    gx, gy = x / length, y / length
    i0, j0 = np.floor(gx).astype(np.int64), np.floor(gy).astype(np.int64)
    u, v = gx - i0, gy - j0
    w = [(1 - u) * (1 - v), u * (1 - v), (1 - u) * v, u * v]
    nodes = [(0, 0), (1, 0), (0, 1), (1, 1)]
    val = sum(wk * _hash_normal(i0 + a, j0 + b, seed, salt) for wk, (a, b) in zip(w, nodes))
    return val / np.sqrt(sum(wk**2 for wk in w))


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
    holo_sample = getattr(holo, "sample", None) if holo is not None else None
    if holo_sample is not None:  # tiled render grid (hologram/tiled.py): fine values on demand
        ds, c0 = holo.ds, holo.c0
        holo_re = True
    elif holo is not None:
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
        Ra = Ri  # true radius for the illumination amplitude
        dw_corr = 0.0
        if pillar_errors and pillar_errors.get("field_um"):
            fsz, sd = pillar_errors["field_um"], pillar_errors.get("seed", 0)
            fx = np.floor((Xi + radius) / fsz).astype(np.int64)
            fy = np.floor((Yi + radius) / fsz).astype(np.int64)
            st = pillar_errors.get("stitch_sigma_um", 0.0)
            if st > 0:
                Xi = Xi - st * _hash_normal(fx, fy, sd, 1)
                Yi = Yi - st * _hash_normal(fx, fy, sd, 2)
                Ri = np.hypot(Xi, Yi)
            fb = pillar_errors.get("field_bias_sigma_um", 0.0)
            if fb > 0:
                dw_corr = dw_corr + fb * _hash_normal(fx, fy, sd, 3)
        if pillar_errors and pillar_errors.get("corr_sigma_um", 0.0) > 0:
            dw_corr = dw_corr + pillar_errors["corr_sigma_um"] * _smooth_noise(Xi + radius, Yi + radius, pillar_errors["corr_length_um"], pillar_errors.get("seed", 0), 4)
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
            if holo_sample is not None:
                phi = phi + np.angle(holo_sample(fi, fj))
            else:
                hr = map_coordinates(holo_re, [fi, fj], order=1, mode="constant")
                hi = map_coordinates(holo_im, [fi, fj], order=1, mode="constant")
                phi = phi + np.angle(hr + 1j * hi)
        if lookup is None:
            t = np.exp(1j * phi)
        elif pillar_errors:
            w = lookup.width_for_phase(phi) + pillar_errors.get("bias_um", 0.0) + dw_corr
            if pillar_errors.get("sigma_um", 0.0) > 0:
                w = w + pillar_errors["sigma_um"] * rng.standard_normal(w.shape)
            t = lookup.t_of_width(w) * amp_scale
        else:
            t = lookup.t_for_phase(phi) * amp_scale
        amp = profile(illumination["kind"], Ra, radius, **kw)
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
