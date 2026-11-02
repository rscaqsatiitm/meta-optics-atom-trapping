"""Off-axis wave check (superposition from the focusing plate) of the shared doublet with
per-wavelength aplanatic pre-correctors. Usage: python scripts/multiwave_wavecheck.py configs/multiwave.yaml
Output: results/<run_id>/wavecheck.json."""

import sys

import numpy as np

from metatrap.metalens import doublet as D
from metatrap.metalens import multiwave as M
from metatrap.metalens import wave_offaxis as WO
from metatrap.metalens.design import _quality
from metatrap.metalens.raytrace import Medium, point_source_rays, pupil_grid
from metatrap.runio import load_config, run_dir, write_json


def main(cfg_path):
    cfg = load_config(cfg_path)
    wc = cfg["wavecheck"]
    out = run_dir(cfg["run_id"])
    spec = D.DoubletSpec(cfg["na"], cfg["cell"], [list(l) for l in cfg["spacer"]])
    d = D.make(spec, wc["design_um"])
    R1, f, L = d.info["R1_um"], d.info["f_um"], wc["L_um"]
    res = {}
    for lam in wc["wavelengths_um"]:
        A, C, info = M.aplanatic_precorrector(d, lam, L)
        sysm = M.precorrected_system(d, lam, A, C, L)
        for h in wc["heights_um"]:
            th = np.arctan(h / f)
            RA = R1 * 1.15 + abs(L * np.tan(th)) + 50
            x, y = pupil_grid(RA, 48)
            px, py, o = point_source_rays(x, y, th, None)
            q = _quality(*sysm.trace(x, y, px, py, lam, stop={"index": 2, "radius_um": info["h1_max_um"]}, opl0=o), lam, None)
            pf = WO.plate_field_front(d, lam, th, [A, Medium(1.0, L), C], 2, info["h1_max_um"], RA, wc["n_pupil"])
            r = WO.strehl_from_field(d, lam, pf, q["x_um"], q["dz_um"])
            res[f"{lam}:{h}"] = {"ray_strehl": q["strehl"], "ray_x_um": q["x_um"], "ray_dz_um": q["dz_um"], "wave": r}
            print(lam, h, "ray", round(q["strehl"], 4), round(q["dz_um"], 3), "wave", {k: round(v, 4) for k, v in r.items()}, flush=True)
            write_json(out / "wavecheck.json", res)


if __name__ == "__main__":
    main(sys.argv[1])
