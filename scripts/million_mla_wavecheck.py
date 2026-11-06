"""Wave-optics check of per-lenslet pre-compensation (superposition method).

Usage: python scripts/million_mla_wavecheck.py configs/million.yaml
For the objectives of cfg['mla'] at the corner field (2.83 mm), the stop-plane wavefront error
W(x, y) from the ray trace (one common atom plane) is subtracted at the stop, and the focal spot
is computed by the superposition method (wave_offaxis), without and with the correction.
Output: results/million/mla_wavecheck.json.
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


def main(cfg_path):
    cfg = load_config(cfg_path)
    mc = cfg["mla"]
    lam, na = mc["wavelength_um"], mc["na"]
    out = {}
    for des in mc["objectives"]:
        cell = [[des["air_gap_um"], "air"], [4000, "fused_silica"], [9000, "air"]]
        f0 = LayerStack.from_list(cell).aperture_for_na(na, lam) / na
        d = D.make(D.DoubletSpec(na, cell, [[f0, "air"]], n_pts=2001, ext=min(1.6, 0.95 / na)), lam)
        R1, f = d.info["R1_um"], d.info["f_um"]
        s = D.system(d, lam)
        dz, _ = D.common_plane(s, R1, lam, list(np.linspace(0, max(mc["field_um"]), 6)), f, n_rays=28)
        for rr in mc.get("wavecheck_field_um", [max(mc["field_um"])]):
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
            res = {}
            for tag, extra in (("uncorrected", None), ("corrected", lambda x, y: -interp(np.column_stack([np.ravel(x), np.ravel(y)])).reshape(np.shape(x)))):
                pf = WO.plate_field(d, lam, th, 260, extra_opl=extra)
                # search the best focus near the common plane (the correction includes the defocus to it)
                w = WO.strehl_from_field(d, lam, pf, q["x_um"], dz, span_um=0.6, dz_span_um=6.0 if extra is None else 1.0)
                res[tag] = w
            res["ray_strehl_uncorrected"] = q["strehl"]
            out[f"{des['name']}:{rr:.0f}"] = res
            print(des["name"], rr, {k: (v if not isinstance(v, dict) else {a: round(b, 4) for a, b in v.items()}) for k, v in res.items()}, flush=True)
    write_json(REPO_ROOT / "results" / "million" / "mla_wavecheck.json", out)


if __name__ == "__main__":
    main(sys.argv[1])
