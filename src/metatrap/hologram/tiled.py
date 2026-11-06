"""Memory-bounded rendering and focal-plane metrics for very large arrays (10^6 sites).

TiledRenderGrid replaces the fine RenderGrid of symmetric.render_grid_symmetric, which for a
10^6-site design would be a ~44,500^2 complex64 array (16 GB). It keeps only the coarse pupil
grid and, for each block of metasurface rows, Fourier-upsamples the coarse sub-array that the
block needs (plus a margin of `margin` coarse samples, tapered) before the bilinear lookup,
normalising to unit modulus inside the same limit as the full render. A row block of the
metasurface maps to a narrow band of the pupil grid, so the temporary is a band, not the grid.

focal_metrics_tiled replaces fullaperture.focal_metrics: the folded focal field is upsampled
and measured tile by tile (each tile padded by a margin that covers the per-site analysis
window), with the same per-site metrics as hologram.evaluate.spot_metrics.

Both are checked against the untiled versions on small problems in tests/test_tiled.py.
"""

from __future__ import annotations

import numpy as np
import scipy.fft as sfft
from numba import njit, prange

from ..propagation.folded import propagate_folded
from .evaluate import _quad_peak


def fourier_interp_2d(x: np.ndarray, factor: int) -> np.ndarray:
    """Band-limited interpolation with sample alignment: output[k] <-> input coordinate k / factor."""
    a, b = x.shape
    A, B = a * factor, b * factor
    spec = sfft.fft2(x, workers=-1)
    big = np.zeros((A, B), np.complex64)
    ha, hb = (a + 1) // 2, (b + 1) // 2
    big[:ha, :hb] = spec[:ha, :hb]
    big[:ha, B - (b - hb):] = spec[:ha, hb:]
    big[A - (a - ha):, :hb] = spec[ha:, :hb]
    big[A - (a - ha):, B - (b - hb):] = spec[ha:, hb:]
    return (sfft.ifft2(big, workers=-1) * (factor * factor)).astype(np.complex64)


def _kernel_table(W: int, beta: float, res: int = 1024) -> np.ndarray:
    """Kaiser-windowed sinc on [0, W] sampled at 1/res (linear interpolation at use)."""
    t = np.arange(W * res + 2) / res
    win = np.where(t < W, np.i0(beta * np.sqrt(np.clip(1 - (t / W) ** 2, 0, None))) / np.i0(beta), 0.0)
    return (np.sinc(t) * win).astype(np.float32)


def _sinc_upsample_1d(x: np.ndarray, factor: int, axis: int, W: int, beta: float) -> np.ndarray:
    x = np.moveaxis(x, axis, 0)
    a = x.shape[0]
    xp = np.concatenate([np.zeros((W,) + x.shape[1:], x.dtype), x, np.zeros((W,) + x.shape[1:], x.dtype)])
    out = np.zeros((a * factor,) + x.shape[1:], np.complex64)
    k = np.arange(-W + 1, W + 1)
    for r in range(factor):
        t = r / factor - k
        h = np.sinc(t) * np.i0(beta * np.sqrt(np.clip(1 - (t / W) ** 2, 0, None))) / np.i0(beta)
        acc = np.zeros(x.shape, np.complex64)
        for kk, hk in zip(k, h):
            if hk != 0.0:
                acc += np.float32(hk) * xp[W + kk : W + kk + a]
        out[r::factor] = acc
    return np.moveaxis(out, 0, axis)


def sinc_upsample_2d(x: np.ndarray, factor: int, W: int = 8, beta: float = 9.0) -> np.ndarray:
    """Local band-limited interpolation, output[k] <-> input coordinate k / factor.

    Separable Kaiser-windowed sinc of half-width W samples. For fields band-limited well
    inside Nyquist (a focal field sampled at dx <= lambda / (4 NA) is band-limited to half of
    it) the error is far below that of FFT interpolation of a non-periodic tile, whose
    wrap-around discontinuity rings as 1/distance into the tile."""
    return _sinc_upsample_1d(_sinc_upsample_1d(np.asarray(x, np.complex64), factor, 0, W, beta), factor, 1, W, beta)


@njit(parallel=True, fastmath=True, cache=True)
def _sinc_interp(cre, cim, ci, cj, tab, W, res, out_re, out_im):
    n0, n1 = cre.shape
    for p in prange(ci.shape[0]):
        x, y = ci[p], cj[p]
        i0 = int(np.floor(x))
        j0 = int(np.floor(y))
        wx = np.empty(2 * W, np.float32)
        wy = np.empty(2 * W, np.float32)
        for k in range(2 * W):
            u = abs(x - (i0 - W + 1 + k)) * res
            t = int(u)
            f = u - t
            wx[k] = tab[t] + f * (tab[t + 1] - tab[t]) if t < tab.shape[0] - 1 else 0.0
            u = abs(y - (j0 - W + 1 + k)) * res
            t = int(u)
            f = u - t
            wy[k] = tab[t] + f * (tab[t + 1] - tab[t]) if t < tab.shape[0] - 1 else 0.0
        sr = 0.0
        si = 0.0
        for a in range(2 * W):
            ii = (i0 - W + 1 + a) % n0  # periodic, as the FFT upsampling of the reference
            rr = 0.0
            ri = 0.0
            for b in range(2 * W):
                jj = (j0 - W + 1 + b) % n1
                rr += wy[b] * cre[ii, jj]
                ri += wy[b] * cim[ii, jj]
            sr += wx[a] * rr
            si += wx[a] * ri
        out_re[p] = sr
        out_im[p] = si


class TiledRenderGrid:
    """Coarse pupil grid with fine values evaluated on demand (no fine grid in memory).

    The reference (render_grid_symmetric) Fourier-upsamples the coarse grid by `upsample`,
    normalises to unit modulus inside the render limit and interpolates bilinearly. Here the
    band-limited interpolant is evaluated directly at each requested point with a separable
    Kaiser-windowed sinc of half-width W coarse samples, then normalised. The coarse grid is
    oversampled (pupil_guard > 1), so W = 8 gives interpolation errors well below the
    bilinear error of the reference (tests/test_tiled.py)."""

    def __init__(self, coarse: np.ndarray, ds_coarse: float, upsample: int, na: float, margin_render: int, W: int = 8, beta: float = 9.0):
        self.cre = np.ascontiguousarray(coarse.real, np.float32)
        self.cim = np.ascontiguousarray(coarse.imag, np.float32)
        self.n = coarse.shape[0]  # axis at n // 2
        self.U = upsample
        self.ds = ds_coarse / upsample  # fine-grid spacing, as RenderGrid.ds
        self.c0 = (self.n * upsample) // 2  # fine-grid centre index, as RenderGrid.c0
        self.lim = (na + 0.5 * margin_render * ds_coarse) ** 2
        self.W, self.res = W, 1024
        self.tab = _kernel_table(W, beta, self.res)

    def sample(self, fi: np.ndarray, fj: np.ndarray) -> np.ndarray:
        """Unit-modulus hologram field at fine-grid coordinates (fi, fj); zero outside the limit."""
        shape = np.shape(fi)
        fi = np.ravel(fi).astype(np.float64)
        fj = np.ravel(fj).astype(np.float64)
        ci = (fi - self.c0) / self.U + self.n // 2
        cj = (fj - self.c0) / self.U + self.n // 2
        out_re = np.empty(fi.size, np.float32)
        out_im = np.empty(fi.size, np.float32)
        _sinc_interp(self.cre, self.cim, ci, cj, self.tab, self.W, self.res, out_re, out_im)
        z = out_re + 1j * out_im
        r2 = ((fi - self.c0) * self.ds) ** 2 + ((fj - self.c0) * self.ds) ** 2
        z = np.where(r2 <= self.lim, z / np.maximum(np.abs(z), 1e-30), 0)
        return z.astype(np.complex64).reshape(shape)


def tiled_render_grid(grid, stack, na, illumination, fsites, fstab, coeffs, table, lam_um, margin: int = 8, upsample: int = 8, kernel=None):
    """As lazy.render_grid_symmetric_lazy, returning a TiledRenderGrid (no fine grid in memory).

    kernel(wp, fsites, table, lam_um, rho_max) builds the lazy kernel (CPU or GPU)."""
    from .lazy import LazyKG
    from .symmetric import WedgePupil, unfold

    wp = WedgePupil.build(grid, stack, na, illumination, s_max=na + margin * grid.ds)
    rmax = float(np.max(np.hypot(wp.r[:, 0], wp.r[:, 1]))) + float(np.max(np.abs(fsites))) * 1.5 + 50
    K = (kernel or (lambda *a: LazyKG(*a)))(wp, fsites, table, lam_um, rmax)
    q = K @ (coeffs / fstab)
    full = unfold(grid.n, wp.ij, q.astype(np.complex64))
    return TiledRenderGrid(full, grid.ds, upsample, na, margin)


def _site_metrics(I, origin_ij, dxu, sites, pitch, lam, na, out, idx):
    """Per-site metrics on an intensity tile whose element (0, 0) is global index origin_ij."""
    r_search = max(1, int(round(0.35 * pitch / dxu)))
    r_cell, r_core = 0.5 * pitch, 0.61 * lam / na
    rr = int(np.ceil(r_cell / dxu)) + 2
    yy, xx = np.mgrid[-rr : rr + 1, -rr : rr + 1]
    oi, oj, c = origin_ij
    for k, (sx, sy) in zip(idx, sites):
        i0, j0 = c + int(round(sx / dxu)) - oi, c + int(round(sy / dxu)) - oj
        win = I[i0 - r_search : i0 + r_search + 1, j0 - r_search : j0 + r_search + 1]
        a, b = np.unravel_index(np.argmax(win), win.shape)
        pi, pj = i0 - r_search + a, j0 - r_search + b
        peak, ox, oy, dxx, dyy = _quad_peak(I[pi - 1 : pi + 2, pj - 1 : pj + 2])
        out["peak"][k] = peak
        out["x"][k], out["y"][k] = (pi + oi - c + ox) * dxu, (pj + oj - c + oy) * dxu
        sub = I[pi - rr : pi + rr + 1, pj - rr : pj + rr + 1]
        d = np.hypot((xx - ox) * dxu, (yy - oy) * dxu)
        out["p_cell"][k] = np.sum(sub[d <= r_cell]) * dxu * dxu
        out["p_core"][k] = np.sum(sub[d <= r_core]) * dxu * dxu
        curv = 0.5 * (dxx + dyy) / dxu**2
        out["waist"][k] = np.sqrt(-4.0 * peak / curv) if curv < 0 else np.nan


def focal_metrics_tiled(cfg, F, dx, lam, na, sites, p_in, H=None, stack_eval=None, site_shift=(0.0, 0.0), tile_px: int = 1024):
    """Tiled counterpart of fullaperture.focal_metrics (same keys; I is returned as None)."""
    nf = F.shape[0]
    foc = propagate_folded(F, dx, lam, stack_eval, low_memory=True) if H is None else np.fft.fftshift(np.fft.ifft2(np.fft.fft2(F) * H))
    up = cfg.get("upsample", 4)
    dxu = dx / up
    pitch = cfg["site_pitch_um"]
    tsites = sites + np.array(site_shift)
    marg = int(np.ceil((pitch * 1.5) / dx)) + 4  # covers the search and cell windows
    cU = (nf // 2) * up  # global upsampled index of the axis
    out = {k: np.zeros(len(tsites)) for k in ("peak", "x", "y", "p_cell", "p_core", "waist")}
    inter_vals, inter_rad = [], []
    axis_I = None
    gi = np.floor(tsites[:, 0] / dx).astype(int) + nf // 2
    gj = np.floor(tsites[:, 1] / dx).astype(int) + nf // 2
    lo_i, hi_i, lo_j, hi_j = gi.min(), gi.max() + 1, gj.min(), gj.max() + 1
    for a0 in range(lo_i, hi_i, tile_px):
        for b0 in range(lo_j, hi_j, tile_px):
            sel = np.where((gi >= a0) & (gi < a0 + tile_px) & (gj >= b0) & (gj < b0 + tile_px))[0]
            if len(sel) == 0:
                continue
            i0, i1 = a0 - marg, min(a0 + tile_px, hi_i) + marg
            j0, j1 = b0 - marg, min(b0 + tile_px, hi_j) + marg
            tile = foc[i0:i1, j0:j1]
            I = np.abs(sinc_upsample_2d(tile, up)) ** 2
            _site_metrics(I, (i0 * up, j0 * up, cU), dxu, tsites[sel], pitch, lam, na, out, sel)
            if i0 <= nf // 2 < i1 and j0 <= nf // 2 < j1:
                axis_I = I[(nf // 2 - i0) * up, (nf // 2 - j0) * up]
            # Interstitial points (centres of four sites) inside this tile's core.
            for s_ in sel:
                cx, cy = tsites[s_, 0] + 0.5 * pitch, tsites[s_, 1] + 0.5 * pitch
                if cx > tsites[:, 0].max() or cy > tsites[:, 1].max():
                    continue
                pi = int(round(cx / dxu)) + cU - i0 * up
                pj = int(round(cy / dxu)) + cU - j0 * up
                if 0 <= pi < I.shape[0] and 0 <= pj < I.shape[1]:
                    inter_vals.append(I[pi, pj])
                    inter_rad.append(np.hypot(cx, cy))
            del I
    pk = out["peak"]
    pk_mean = pk.mean()
    res = {
        "n_sites": int(len(tsites)),
        "depth_nonuniformity_rms": float(np.std(pk) / pk_mean),
        "depth_min_over_mean": float(pk.min() / pk_mean),
        "array_efficiency": float(out["p_cell"].sum() / p_in),
        "core_efficiency": float(out["p_core"].sum() / p_in),
        "cell_power_nonuniformity_rms": float(np.std(out["p_cell"]) / np.mean(out["p_cell"])),
        "waist_mean_um": float(np.nanmean(out["waist"])),
        "waist_std_um": float(np.nanstd(out["waist"])),
        "position_error_rms_um": float(np.sqrt(np.mean((out["x"] - tsites[:, 0]) ** 2 + (out["y"] - tsites[:, 1]) ** 2))),
        "axis_intensity_over_mean_site_peak": float(axis_I / pk_mean) if axis_I is not None else float("nan"),
    }
    if inter_vals:
        v = np.array(inter_vals) / pk_mean
        res |= {"interstitial_median": float(np.median(v)), "interstitial_p99": float(np.percentile(v, 99)),
                "interstitial_max": float(v.max()), "interstitial_max_radius_um": float(inter_rad[int(np.argmax(v))])}
    res["mean_site_offset_um"] = (float(np.mean(out["x"] - tsites[:, 0])), float(np.mean(out["y"] - tsites[:, 1])))
    return res, out, None
