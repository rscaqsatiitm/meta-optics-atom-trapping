"""Absorbing defects on a thin plate, and the sensitivity of the global rise to the heat sinks.

Usage: python scripts/defect_plate.py configs/damage_v2.yaml

The defect table of docs/prefab_checks_report.md used the half-space result dT = I a / kappa for a
fully absorbing disc of radius a. That holds only for a << t (plate thickness). For a >~ t the heat
must spread laterally through the plate, dT ~ P / (4 pi kappa t) [1 + 2 ln(R / a)] with
P = I pi a^2, which grows as a^2 ln(R / a) instead of a. Here the defect is a top-hat absorber
of radius a on the film face of the axisymmetric finite-difference model (thermal/axisym_fd.py:
die radius 12.5 mm, rim held at ambient, h = 10 W/m^2/K on both faces), and the largest radius
with a rise <= 50 K is found by bisection, for fused silica and sapphire (0.5 mm) at the design
intensities. The rise is that of the defect alone; the global rise of the part adds to it.

Second part: the global peak rise of the reference case (fused silica 0.5 mm, Si-rich nominal,
100 W) against the face coefficient h and the rim condition.
Output: results/damage_v2/defect_plate.json.
"""

import sys

import numpy as np

from metatrap.propagation.layers import LayerStack
from metatrap.runio import REPO_ROOT, load_config, write_json
from metatrap.thermal import axisym_fd as A
from metatrap.thermal.damage import defect_rise


def defect_dT(sub, t, a, intensity, R, h=10.0):
    prob = A.ThermalProblem(radius_m=R, thickness_m=t, substrate=sub, beam="top_hat", beam_radius_m=a, h_faces=h, rim="fixed",
                            nr=300, nz=48, z_power=1.6)
    rc, zc, T = A.solve(prob, intensity * np.pi * a**2, 1.0)
    return float(T.max())


def largest(sub, t, intensity, R, dT_max):
    lo, hi = 1e-6, 0.99 * R
    for _ in range(40):
        mid = np.sqrt(lo * hi)
        if defect_dT(sub, t, mid, intensity, R) > dT_max:
            hi = mid
        else:
            lo = mid
        if hi / lo < 1.005:
            break
    return lo


def main(cfg_path):
    cfg = load_config(cfg_path)
    lam, na = cfg["wavelength_um"], cfg["na"]
    Rap = LayerStack.from_list(cfg["stack"]).aperture_for_na(na, lam) * 1e-6
    R = cfg["die_radius_mm"] * 1e-3
    t = 0.5e-3
    lim = cfg["limits"]["defect_dT_K"]
    res = {"aperture_radius_m": Rap, "die_radius_m": R, "thickness_m": t, "dT_limit_K": lim, "defects": {}}
    for P in (50, 100, 150, 200):
        I = P / (np.pi * Rap**2)
        row = {"intensity_W_cm2": I / 1e4}
        for sub in (A.FUSED_SILICA, A.SAPPHIRE):
            a = largest(sub, t, I, R, lim)
            a_hs = lim * sub.kappa / I
            row[sub.name] = {"a_max_plate_um": a * 1e6, "a_max_half_space_um": a_hs * 1e6,
                             "dT_at_half_space_radius_K": defect_dT(sub, t, min(a_hs, 0.99 * R), I, R),
                             "half_space_check_dT_K": defect_rise(I, a, sub.kappa)}
        res["defects"][str(P)] = row
        print(P, {k: (v if not isinstance(v, dict) else {a: round(b, 1) for a, b in v.items()}) for k, v in row.items()}, flush=True)

    # Global rise of the reference case against the heat-sink assumptions.
    sub, a_film, P = A.FUSED_SILICA, 5.9e-4, 100.0
    sens = {}
    for label, kw in (("h10_rim_fixed", dict(h_faces=10.0, rim="fixed")), ("h0_rim_fixed", dict(h_faces=0.0, rim="fixed")),
                      ("h5_rim_fixed", dict(h_faces=5.0, rim="fixed")), ("h20_rim_fixed", dict(h_faces=20.0, rim="fixed")),
                      ("h10_rim_convective_h10", dict(h_faces=10.0, rim="convective", h_rim=10.0))):
        prob = A.ThermalProblem(radius_m=R, thickness_m=t, substrate=sub, beam="top_hat", beam_radius_m=Rap, **kw)
        rc, zc, T = A.solve(prob, P, a_film)
        sens[label] = float(T.max())
        print("global", label, round(sens[label], 2), flush=True)
    res["global_rise_fs0.5_sirich_nominal_100W_K"] = sens
    write_json(REPO_ROOT / "results" / "damage_v2" / "defect_plate.json", res)


if __name__ == "__main__":
    main(sys.argv[1])
