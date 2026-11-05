"""Imaging metric for very large arrays: ensquared energy per site and crosstalk.

Usage: python scripts/million_imaging.py configs/million.yaml

For site-resolved fluorescence imaging what matters is the fraction of each atom's light that
lands in its own pitch x pitch box in the image (ensquared energy, EE) and the fraction landing
in the neighbouring boxes (crosstalk), not the Strehl ratio. For each meta-objective design
(air-spaced aplanatic doublet, NA 0.5, 852.3 nm, construction extended to cover the field) the
point-spread function of a site at field radius r is computed by the superposition method of
docs/metalens_report.md Section 4b (rays to the focusing plate, exact point-to-point paths
through the cell), in one common focal plane chosen for the whole field, on a grid spanning
+/-1.5 pitch around the ray image. Output: results/million/imaging_ee.json.
"""

import json
import sys

import numpy as np

from metatrap.metalens import doublet as D
from metatrap.metalens import wave_offaxis as WO
from metatrap.runio import REPO_ROOT, load_config, write_json


def main(cfg_path):
    cfg = load_config(cfg_path)
    ic = cfg["imaging"]
    lam, na, a = ic["wavelength_um"], ic["na"], ic["pitch_um"]
    path = REPO_ROOT / "results" / "million" / "imaging_ee.json"
    res = json.load(open(path)) if path.exists() else {}
    for des in ic["designs"]:
        name = des["name"]
        cell = [[des["air_gap_um"], "air"], [4000, "fused_silica"], [9000, "air"]]
        from metatrap.propagation.layers import LayerStack

        f0 = LayerStack.from_list(cell).aperture_for_na(na, lam) / na
        d = D.make(D.DoubletSpec(na, cell, [[f0 * des.get("spacing_frac", 1.0), "air"]], n_pts=2001, ext=min(ic["ext"], 0.95 / na)), lam)
        R1, f = d.info["R1_um"], d.info["f_um"]
        s = D.system(d, lam)
        dz, _ = D.common_plane(s, R1, lam, list(np.linspace(0, max(ic["field_um"]), 6)), f, n_rays=28)
        for r in ic["field_um"]:
            key = f"{name}:{r:.0f}"
            if key in res:
                continue
            th = D.theta_for_height(s, R1, lam, r, f)
            q = D.quality(s, R1, lam, th, None, dz)
            pf = WO.plate_field(d, lam, th, ic["n_pupil"])
            fs = WO.FocalSum(d.spec.cell, lam, d.info["R2_um"] * 1.6 + r + 50)
            g = np.linspace(-1.5 * a, 1.5 * a, ic["n_grid"])
            GX, GY = np.meshgrid(q["x_um"] + g, g, indexing="ij")
            I = np.abs(fs.field(pf, pf["amp"], np.column_stack([GX.ravel(), GY.ravel()]), dz)).reshape(GX.shape) ** 2
            # Centre the box on the intensity centroid (distortion is calibrated in practice).
            cx = np.sum(I * GX) / I.sum() - q["x_um"]
            cy = np.sum(I * GY) / I.sum()
            inbox = (np.abs(GX - q["x_um"] - cx) <= a / 2) & (np.abs(GY - cy) <= a / 2)
            nb = {}
            for (ix, iy) in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                m = (np.abs(GX - q["x_um"] - cx - ix * a) <= a / 2) & (np.abs(GY - cy - iy * a) <= a / 2)
                nb[f"{ix},{iy}"] = float(I[m].sum() / I.sum())
            res[key] = {"design": name, "f_um": f, "field_um": r, "ray_strehl_common_plane": q["strehl"], "dz_common_um": dz,
                        "ensquared_energy": float(I[inbox].sum() / I.sum()), "neighbour_fraction": nb,
                        "max_neighbour": max(nb.values()), "peak_rel": float(I.max())}
            print(key, {k: (round(v, 4) if isinstance(v, float) else v) for k, v in res[key].items() if k != "neighbour_fraction"}, flush=True)
            write_json(path, res)


if __name__ == "__main__":
    main(sys.argv[1])
