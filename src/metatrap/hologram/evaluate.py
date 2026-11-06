"""Spot and array metrics from a sampled focal-plane intensity.

Definitions (used consistently in the report):

* peak intensity I_j: maximum of the intensity within 0.35 pitch of site j, refined by
  a 3 x 3 quadratic fit; trap depth is proportional to I_j (red-detuned, linear pol.).
* trap-depth non-uniformity: std(I_j) / mean(I_j) (proposal Table 1 metric).
* cell power P_j: power within a disc of radius pitch/2 about the fitted peak.
* array efficiency: sum_j P_j / P_incident.
* core power: power within the first-dark-ring radius 0.61 lambda / NA.
* Strehl S_j: I_j / (P_core,j x I_ref / P_core,ref), the peak per unit core power
  relative to an aberration-free reference spot computed with the same pupil.
* waist w_j: 1/e^2 radius of the Gaussian with the same peak curvature,
  w^2 = -8 I / (d^2I/dx^2 + d^2I/dy^2), i.e. -4 I over the mean curvature.
* position error: fitted peak minus target site.
"""

from __future__ import annotations

import numpy as np


def _quad_peak(patch: np.ndarray):
    """Sub-pixel peak and curvatures from a 3x3 patch (log-free quadratic fit)."""
    c = patch[1, 1]
    dxx = patch[2, 1] - 2 * c + patch[0, 1]
    dyy = patch[1, 2] - 2 * c + patch[1, 0]
    dx = 0.5 * (patch[2, 1] - patch[0, 1])
    dy = 0.5 * (patch[1, 2] - patch[1, 0])
    ox = -dx / dxx if dxx < 0 else 0.0
    oy = -dy / dyy if dyy < 0 else 0.0
    ox, oy = np.clip(ox, -0.5, 0.5), np.clip(oy, -0.5, 0.5)
    peak = c + 0.5 * (dx * ox + dy * oy)
    return peak, ox, oy, dxx, dyy


def spot_metrics(intensity: np.ndarray, dx: float, sites_um: np.ndarray, pitch_um: float, lam_um: float, na: float, p_incident: float = 1.0, ref: dict | None = None) -> dict:
    """Metrics for every site; intensity has the optical axis at index N//2."""
    n = intensity.shape[0]
    c = n // 2
    r_search = max(1, int(round(0.35 * pitch_um / dx)))
    r_cell = 0.5 * pitch_um
    r_core = 0.61 * lam_um / na
    rr = int(np.ceil(r_cell / dx)) + 2
    yy, xx = np.mgrid[-rr : rr + 1, -rr : rr + 1]
    out = {k: np.zeros(len(sites_um)) for k in ("peak", "x", "y", "p_cell", "p_core", "waist")}
    pix_area = dx * dx
    for j, (sx, sy) in enumerate(sites_um):
        i0, j0 = c + int(round(sx / dx)), c + int(round(sy / dx))
        win = intensity[i0 - r_search : i0 + r_search + 1, j0 - r_search : j0 + r_search + 1]
        a, b = np.unravel_index(np.argmax(win), win.shape)
        pi, pj = i0 - r_search + a, j0 - r_search + b
        peak, ox, oy, dxx, dyy = _quad_peak(intensity[pi - 1 : pi + 2, pj - 1 : pj + 2])
        px, py = (pi - c + ox) * dx, (pj - c + oy) * dx
        sub = intensity[pi - rr : pi + rr + 1, pj - rr : pj + rr + 1]
        d = np.hypot((xx - ox) * dx, (yy - oy) * dx)
        out["peak"][j] = peak
        out["x"][j], out["y"][j] = px, py
        out["p_cell"][j] = np.sum(sub[d <= r_cell]) * pix_area
        out["p_core"][j] = np.sum(sub[d <= r_core]) * pix_area
        curv = 0.5 * (dxx + dyy) / dx**2
        out["waist"][j] = np.sqrt(-4.0 * peak / curv) if curv < 0 else np.nan
    pk = out["peak"]
    res = {
        "n_sites": int(len(sites_um)),
        "depth_nonuniformity_rms": float(np.std(pk) / np.mean(pk)),
        "depth_min_over_mean": float(pk.min() / pk.mean()),
        "array_efficiency": float(out["p_cell"].sum() / p_incident),
        "core_efficiency": float(out["p_core"].sum() / p_incident),
        "cell_power_nonuniformity_rms": float(np.std(out["p_cell"]) / np.mean(out["p_cell"])),
        "waist_mean_um": float(np.nanmean(out["waist"])),
        "waist_std_um": float(np.nanstd(out["waist"])),
        "position_error_rms_um": float(np.sqrt(np.mean((out["x"] - sites_um[:, 0]) ** 2 + (out["y"] - sites_um[:, 1]) ** 2))),
        "per_site": out,
    }
    if ref is not None:
        strehl = (pk / out["p_core"]) / (ref["peak"] / ref["p_core"])
        res["per_site"]["strehl"] = strehl
        res["strehl_mean"] = float(np.mean(strehl))
        res["strehl_min"] = float(np.min(strehl))
    return res


def single_spot_reference(intensity: np.ndarray, dx: float, lam_um: float, na: float, pitch_um: float) -> dict:
    """Peak and core power of an on-axis reference spot (for Strehl normalisation)."""
    m = spot_metrics(intensity, dx, np.zeros((1, 2)), pitch_um, lam_um, na)
    return {"peak": float(m["per_site"]["peak"][0]), "p_core": float(m["per_site"]["p_core"][0]), "waist_um": float(m["per_site"]["waist"][0])}


def zero_order_and_ghosts(intensity: np.ndarray, dx: float, sites_um: np.ndarray, pitch_um: float, p_incident: float = 1.0) -> dict:
    """Largest stray peak outside the site cells relative to the mean site peak."""
    n = intensity.shape[0]
    c = n // 2
    mask = np.ones_like(intensity, dtype=bool)
    rr = int(np.ceil(0.5 * pitch_um / dx))
    for sx, sy in sites_um:
        i0, j0 = c + int(round(sx / dx)), c + int(round(sy / dx))
        mask[i0 - rr : i0 + rr + 1, j0 - rr : j0 + rr + 1] = False
    stray = np.where(mask, intensity, 0.0)
    ext = np.abs(sites_um).max() + 3 * pitch_um
    e = int(ext / dx)
    region = stray[c - e : c + e + 1, c - e : c + e + 1]
    return {"max_stray_peak": float(region.max()), "zero_order_peak": float(intensity[c, c]) if mask[c, c] else float("nan")}
