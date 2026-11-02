"""One shared metasurface doublet for several wavelengths (ray optics).

Usage: python scripts/metalens_multiwave.py configs/multiwave.yaml

For each shared design wavelength lambda_d and pre-corrector separation L: the air-spaced
aplanatic doublet is built at lambda_d. Each operating wavelength gets its two-surface
aplanatic pre-corrector (multiwave.aplanatic_precorrector; A a distance L before C, C relayed
onto the stop at the corrector metasurface). It is evaluated over the field in one common
focal plane per wavelength, with the stop at C. Output: results/<run_id>/optics.json.
"""

import json
import sys

import numpy as np

from metatrap.metalens import doublet as D
from metatrap.metalens import multiwave as M
from metatrap.metalens.design import _quality
from metatrap.metalens.raytrace import point_source_rays, pupil_grid
from metatrap.runio import load_config, run_dir, write_json, write_params


def qual(s, lam, th, rstop, RA, L, dz=None, n=48):
    x, y = pupil_grid(RA + abs(L * np.tan(th)) + 50, n)
    px, py, o = point_source_rays(x, y, th, None)
    q = _quality(*s.trace(x, y, px, py, lam, stop={"index": 2, "radius_um": rstop}, opl0=o), lam, dz)
    q.pop("residual", None)
    return q


def common(s, lam, rstop, RA, L, f, heights):
    ths = [np.arctan(h / f) for h in heights]
    free = [qual(s, lam, t, rstop, RA, L)["dz_um"] for t in ths]
    c = float(np.mean(free))
    span = 4 + 0.5 * np.ptp(free)
    grid = np.arange(c - span, c + span, 0.05)
    v = [min(qual(s, lam, t, rstop, RA, L, z, 32)["strehl"] for t in ths) for z in grid]
    z = float(grid[int(np.argmax(v))])
    return z, [qual(s, lam, t, rstop, RA, L, z) for t in ths], free


def main(cfg_path):
    cfg = load_config(cfg_path)
    out = run_dir(cfg["run_id"])
    write_params(out, cfg)
    path = out / "optics.json"
    res = json.load(open(path)) if path.exists() else {}
    spec = D.DoubletSpec(cfg["na"], cfg["cell"], [list(l) for l in cfg["spacer"]])
    for lam_d in cfg["design_um"]:
        d = D.make(spec, lam_d)
        R1, f = d.info["R1_um"], d.info["f_um"]
        for L in cfg["precorrector_L_um"]:
            for lam in cfg["operating_um"]:
                key = f"{lam_d}:{L}:{lam}"
                if key in res:
                    continue
                A, C, info = M.aplanatic_precorrector(d, lam, L)
                s = M.precorrected_system(d, lam, A, C, L)
                z, q, free = common(s, lam, info["h1_max_um"], R1 * 1.15, L, f, cfg["field_heights_um"])
                res[key] = {"design_um": lam_d, "L_um": L, "lambda_um": lam, "dz_common_um": z, "strehl": [x["strehl"] for x in q],
                            "image_x_um": [x["x_um"] for x in q], "free_focus_dz_um": free, "precorrector": info, "doublet": d.info}
                print(key, round(z, 2), [round(x["strehl"], 4) for x in q], "A/C pv waves", round(info["A_pv_waves"]), round(info["C_pv_waves"]), flush=True)
                write_json(path, res)


if __name__ == "__main__":
    main(sys.argv[1])
