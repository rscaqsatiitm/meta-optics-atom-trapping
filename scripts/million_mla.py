"""Metasurface lenslet-array (MLA) tweezer architecture for 10 to 10^6 traps.

Usage: python scripts/million_mla.py configs/million.yaml

A metasurface lenslet array in the front focal plane of a tube lens is relayed onto the atoms
by the meta-objective (demagnification M = f_T / f_obj). Each illuminated lenslet makes one
trap, so a flat-top beam of adjustable size sets the number of traps at constant depth per
trap (N proportional to the illuminated area and to the power; ref. Pause et al., Optica 11,
222 (2024), with refractive MLAs).

(a) Geometry. Atom pitch a, objective NA, demagnification M: lenslet pitch p = M a, lenslet NA
NA_L = NA / M (fills the objective pupil), focal length f_L = p / (2 NA_L) = M^2 a / (2 NA),
Fresnel number N_F = (p/2)^2 / (lambda f_L) = a NA / (2 lambda), independent of M (DERIVED).

(b) Per-lenslet aberration pre-compensation. A lenslet's aperture field is the angular
spectrum of its focal spot, so its aperture maps linearly onto the objective pupil
(p/2 <-> R1). A lenslet phase equal to minus the objective's wavefront error at that lenslet's
image position therefore corrects the field-dependent aberration of the objective, as the
exact per-site kernels do for the hologram. The residual wavefront W(x, y) at the stop is
obtained from the ray trace for each field point (common atom plane); the corrected trace adds
-W to the launch phase. We report Strehl ratios before and after, PV and rms of W, and the extra
local deflection the lenslet must produce, |grad W| x (2 R1 / p).

Output: results/million/mla.json.
"""

import json
import sys

import numpy as np
from scipy.interpolate import RegularGridInterpolator

from metatrap.metalens import doublet as D
from metatrap.metalens.design import _quality
from metatrap.metalens.raytrace import point_source_rays
from metatrap.propagation.layers import LayerStack
from metatrap.runio import REPO_ROOT, load_config, write_json


def geometry_table(mc):
    rows = []
    for a in mc["pitch_um"]:
        for M in mc["demag"]:
            na, lam = mc["na"], mc["wavelength_um"]
            p = M * a
            naL = na / M
            fL = p / (2 * naL)
            NF = (p / 2) ** 2 / (lam * fL)
            side_1e6 = 1000 * p
            rows.append({"a_um": a, "M": M, "lenslet_pitch_um": p, "lenslet_NA": naL, "lenslet_f_um": fL, "fresnel_number": NF,
                         "mla_side_mm_1e6": side_1e6 / 1000, "lenslets_per_cm2": 1e8 / p**2})
    return rows


def stop_grid(R, n):
    c = (np.arange(n) + 0.5) / n * 2 * R - R
    X, Y = np.meshgrid(c, c, indexing="ij")
    return c, X, Y


def corrected(s, R1, lam, th, dz, n=64):
    """Strehl before/after pre-compensation of the stop-plane wavefront error W(x, y)."""
    c, X, Y = stop_grid(R1, n)
    m = np.hypot(X, Y) <= R1
    px, py, o = point_source_rays(X, Y, th, None)
    xf, yf, pxf, pyf, opl, alive = s.trace(X.ravel(), Y.ravel(), px.ravel(), py.ravel(), lam, stop={"index": 0, "radius_um": R1}, opl0=o.ravel())
    q0 = _quality(xf, yf, pxf, pyf, opl, alive, lam, dz)
    W = np.full(X.size, np.nan)
    W[alive] = q0["residual"]
    W = W.reshape(X.shape)
    Wf = np.where(np.isfinite(W), W, 0.0)
    interp = RegularGridInterpolator((c, c), Wf, bounds_error=False, fill_value=0.0)
    # Corrected launch: subtract W (the lenslet supplies -W in pupil coordinates).
    o2 = o.ravel() - interp(np.column_stack([X.ravel(), Y.ravel()]))
    r2 = s.trace(X.ravel(), Y.ravel(), px.ravel(), py.ravel(), lam, stop={"index": 0, "radius_um": R1}, opl0=o2)
    q1 = _quality(*r2, lam, dz)
    gx, gy = np.gradient(Wf, c, c)
    # Interior points only (both neighbours valid along each axis): the zero fill outside the
    # pupil would otherwise create a spurious gradient at the rim.
    from scipy.ndimage import binary_erosion

    valid = binary_erosion(m & np.isfinite(W), iterations=2)
    grad = np.hypot(gx, gy)[valid]
    Wm = W[m & np.isfinite(W)]
    return {"strehl_before": q0["strehl"], "strehl_after": q1["strehl"], "W_pv_waves": float(np.ptp(Wm) / lam),
            "W_rms_waves": float(np.std(Wm) / lam), "max_grad_W": float(np.max(grad)), "x_um": q0["x_um"]}


def main(cfg_path):
    cfg = load_config(cfg_path)
    mc = cfg["mla"]
    path = REPO_ROOT / "results" / "million" / "mla.json"
    res = json.load(open(path)) if path.exists() else {}
    res["geometry"] = geometry_table(mc)
    for r in res["geometry"]:
        print("geometry", {k: (round(v, 3) if isinstance(v, float) else v) for k, v in r.items()}, flush=True)
    lam, na = mc["wavelength_um"], mc["na"]
    for des in mc["objectives"]:
        cell = [[des["air_gap_um"], "air"], [4000, "fused_silica"], [9000, "air"]]
        f0 = LayerStack.from_list(cell).aperture_for_na(na, lam) / na
        d = D.make(D.DoubletSpec(na, cell, [[f0, "air"]], n_pts=2001, ext=min(1.6, 0.95 / na)), lam)
        R1, f = d.info["R1_um"], d.info["f_um"]
        s = D.system(d, lam)
        dz, _ = D.common_plane(s, R1, lam, list(np.linspace(0, max(mc["field_um"]), 6)), f, n_rays=28)
        for rr in mc["field_um"]:
            key = f"{des['name']}:{rr:.0f}"
            if key in res:
                out = res[key]
            else:
                th = D.theta_for_height(s, R1, lam, rr, f)
                out = corrected(s, R1, lam, th, dz)
            # Extra lenslet deflection for each demagnification and site pitch (lenslet pitch M a);
            # earlier versions keyed it by M only, evaluated at the last pitch of the list.
            out.pop("extra_lenslet_deflection_note", None)
            for M in mc["demag"]:
                out.pop(f"extra_lenslet_deflection_M{M}", None)
                for a in mc["pitch_um"]:
                    out[f"extra_lenslet_deflection_M{M}_a{a:g}"] = out["max_grad_W"] * 2 * R1 / (M * a)
            out |= {"f_um": f, "R1_um": R1, "field_um": rr, "dz_common_um": dz}
            res[key] = out
            print(key, {k: (round(v, 4) if isinstance(v, float) else v) for k, v in out.items()}, flush=True)
            write_json(path, res)
    write_json(path, res)


if __name__ == "__main__":
    main(sys.argv[1])
