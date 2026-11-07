"""Imaging with per-atom aberration correction in image space (reciprocity with the lenslet relay).

Usage: python scripts/million_imaging_corrected.py configs/million.yaml

A corrector array behind an intermediate image (one element per atom, where the atoms' beams
are separated) can remove each atom's field-dependent aberration before the camera. By
reciprocity this is the time-reverse of per-lenslet pre-compensation in the trapping relay, and
the corrected point-spread function equals the one obtained by applying -W at the stop. For
the compact objectives (16 and 25 mm, NA 0.5, 852.3 nm) we compute, in wave optics
(superposition method), the ensquared energy in the 4 um site box and the largest
neighbour-box fraction (both normalised to the total PSF power, wave_offaxis.total_focal_power),
with and without correction, at the 10^6-site field radii.
Output: results/million/imaging_corrected.json.
"""

import json
import sys

import numpy as np
from scipy.interpolate import RegularGridInterpolator

from metatrap.metalens import doublet as D
from metatrap.metalens import wave_offaxis as WO
from metatrap.metalens.design import _quality
from metatrap.metalens.raytrace import point_source_rays
from metatrap.propagation.layers import LayerStack
from metatrap.runio import REPO_ROOT, load_config, write_json


def ee(fs, pf, amp, x0, dz, a, na, n=61):
    """Ensquared energy and largest neighbour-box fraction, normalised to the total PSF power."""
    g = np.linspace(-1.5 * a, 1.5 * a, n)
    GX, GY = np.meshgrid(x0 + g, g, indexing="ij")
    I = np.abs(fs.field(pf, amp, np.column_stack([GX.ravel(), GY.ravel()]), dz)).reshape(GX.shape) ** 2
    cx = np.sum(I * GX) / I.sum() - x0
    cy = np.sum(I * GY) / I.sum()
    ptot = WO.total_focal_power(fs, pf, amp, x0 + cx, cy, dz, na) / (g[1] - g[0]) ** 2
    box = (np.abs(GX - x0 - cx) <= a / 2) & (np.abs(GY - cy) <= a / 2)
    nb = []
    for ix, iy in ((1, 0), (-1, 0), (0, 1), (0, -1)):
        m = (np.abs(GX - x0 - cx - ix * a) <= a / 2) & (np.abs(GY - cy - iy * a) <= a / 2)
        nb.append(I[m].sum() / ptot)
    return float(I[box].sum() / ptot), float(max(nb))


def main(cfg_path):
    cfg = load_config(cfg_path)
    ic = cfg["imaging"]
    lam, na, a = ic["wavelength_um"], ic["na"], ic["pitch_um"]
    path = REPO_ROOT / "results" / "million" / "imaging_corrected.json"
    res = json.load(open(path)) if path.exists() else {}
    for des in ic["designs"][:2]:
        cell = [[des["air_gap_um"], "air"], [4000, "fused_silica"], [9000, "air"]]
        f0 = LayerStack.from_list(cell).aperture_for_na(na, lam) / na
        d = D.make(D.DoubletSpec(na, cell, [[f0, "air"]], n_pts=2001, ext=min(1.6, 0.95 / na)), lam)
        R1, f = d.info["R1_um"], d.info["f_um"]
        s = D.system(d, lam)
        dz, _ = D.common_plane(s, R1, lam, list(np.linspace(0, 2828, 6)), f, n_rays=28)
        fs = WO.FocalSum(d.spec.cell, lam, d.info["R2_um"] * 1.6 + 2900)
        for rr in ic.get("corrected_field_um", [1770, 2121, 2828]):
            key = f"{des['name']}:{rr}"
            if key in res:
                continue
            th = D.theta_for_height(s, R1, lam, rr, f)
            n = 96
            c = (np.arange(n) + 0.5) / n * 2 * R1 - R1
            X, Y = np.meshgrid(c, c, indexing="ij")
            px, py, o = point_source_rays(X, Y, th, None)
            r = s.trace(X.ravel(), Y.ravel(), px.ravel(), py.ravel(), lam, stop={"index": 0, "radius_um": R1}, opl0=o.ravel())
            q = _quality(*r, lam, dz)
            W = np.zeros(X.size)
            W[r[5]] = q["residual"]
            interp = RegularGridInterpolator((c, c), W.reshape(X.shape), bounds_error=False, fill_value=0.0)
            out = {"ray_strehl_uncorrected": q["strehl"]}
            for tag, extra in (("uncorrected", None), ("corrected", lambda x, y: -interp(np.column_stack([np.ravel(x), np.ravel(y)])).reshape(np.shape(x)))):
                pf = WO.plate_field(d, lam, th, 220, extra_opl=extra)
                e, nbm = ee(fs, pf, pf["amp"], q["x_um"], dz, a, na)
                out[tag] = {"ensquared_energy": e, "max_neighbour": nbm}
            res[key] = out
            print(key, out, flush=True)
            write_json(path, res)


if __name__ == "__main__":
    main(sys.argv[1])
