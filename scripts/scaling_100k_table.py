"""Results table of docs/scaling_100k.md (working-distance study) from the stored results.

Usage: python scripts/scaling_100k_table.py configs/scaling_100k.yaml

From results/scaling_100k/scaling.json (exact-kernel site spots and meta-atom factors per WD,
imaging-lens field radius) and the v2 power budget (results/fullaperture_v2/metrics.json):
relative depth D(r) = I_rel(r) x eta(r) / eta(0); corner depths at 894 um (316 x 316 at 4 um)
and 670 um (3 um pitch); array power factor F = mean(1 / D_j) over the 316 x 316 sites (D
interpolated in r); power for 0.3 mK P = 0.3 x N x F / K, with K = depth x N / P of the v2
library-pillar design (NA 0.55); intensity P / (pi R^2); largest fully absorbing defect on fused
silica for a 50 K rise, a = 50 K x kappa / I (kappa 1.38 W/m/K); sites imageable N = 2 r^2 / a^2.
Output: results/scaling_100k/table.json (and the markdown rows on stdout).
"""

import json
import sys

import numpy as np

from metatrap.runio import REPO_ROOT, load_config, write_json


def main(cfg_path):
    cfg = load_config(cfg_path)
    res = json.load(open(REPO_ROOT / "results" / "scaling_100k" / "scaling.json"))
    v2 = json.load(open(REPO_ROOT / "results" / "fullaperture_v2" / "metrics.json"))
    K = v2["power_budget"]["library_lpa"]["depth_mK_at_50W"] * v2["n_sites"] / 50.0  # mK x site / W
    n_side, a4 = 316, 4.0
    c = (np.arange(n_side) - (n_side - 1) / 2) * a4
    R_sites = np.hypot(*np.meshgrid(c, c)).ravel()
    out = {"K_mK_site_per_W": K, "rows": {}}
    for key, g in res.items():
        if not key.startswith("WD"):
            continue
        WD, R = g["WD_um"], g["aperture_radius_um"]
        rs = np.array(sorted(float(k) for k in g["sites"]))
        s0 = g["sites"]["0"]
        D = np.array([g["sites"][f"{r:.0f}"]["peak"] / s0["peak"] * g["sites"][f"{r:.0f}"]["eta_site"] / s0["eta_site"] for r in rs])
        Dj = np.interp(R_sites, rs, D)
        F = float(np.mean(1 / Dj))
        P = 0.3 * n_side**2 * F / K
        I = P / (np.pi * (R * 1e-4) ** 2)  # W/cm^2
        a_def = 50 * 1.38 / (I * 1e4) * 1e6  # um
        r_img = g["lens"]["field_radius_strehl08_um"]
        corner = g["sites"]["894"] if "894" in g["sites"] else None
        row = {"WD_um": WD, "aperture_mm": 2 * R / 1000, "D_corner_4um": float(np.interp(894, rs, D)), "D_corner_3um": float(np.interp(670, rs, D)),
               "F": F, "P_0.3mK_W": P, "P_0.18mK_W": 0.6 * P, "intensity_W_cm2": I, "max_defect_um": a_def,
               "corner_w_radial_um": corner["w_radial_um"] if corner else None, "corner_w_tangential_um": corner["w_tangential_um"] if corner else None,
               "corner_s_max": corner["s_max"] if corner else None, "lens_field_radius_um": r_img,
               "sites_imageable_4um": 2 * r_img**2 / 16, "sites_imageable_3um": 2 * r_img**2 / 9}
        out["rows"][key] = row
        print(f"| {WD/1000:g} mm | {2*R/1000:.3g} mm | {row['D_corner_4um']:.2f} / {row['D_corner_3um']:.2f} | {F:.3f} | "
              f"{row['corner_w_radial_um']:.2f} / {row['corner_w_tangential_um']:.2f} um | {row['corner_s_max']:.2f} | {P:.0f} W ({0.6*P:.0f} W) | "
              f"{I:,.0f} W/cm^2 | {a_def:.3g} um | {r_img:,.0f} um | {row['sites_imageable_4um']:,.0f} / {row['sites_imageable_3um']:,.0f} |")
    print(f"K = {K:.1f} mK x site / W")
    write_json(REPO_ROOT / "results" / "scaling_100k" / "table.json", out)


if __name__ == "__main__":
    main(sys.argv[1])
