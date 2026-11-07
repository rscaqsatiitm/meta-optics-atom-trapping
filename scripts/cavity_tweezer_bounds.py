"""S6: why a build-up cavity does not give 10x for the NA 0.5 tweezer array.
Usage:
    python scripts/cavity_tweezer_bounds.py configs/cavity_tweezer_bounds.yaml

(a) Loss bound: an element of per-pass efficiency eta inside the cavity limits G to
    1 / (1 - eta^n) for n passes per round trip (budget.intracavity_element_bound).
(b) Degenerate concentric cavity (meta-optics outside, only mirrors inside): two spherical
    mirrors share their centre C at the atom plane, so a focus at P is imaged to 2C - P and back.
    We ray-trace one reflection from a point at field radius r over an NA cone and compute the
    rms wavefront error after removing piston, tilt and focus (best image point), then the
    round-trip Strehl ratio exp(-(2 pi sigma_rt)^2) and the per-site round-trip phase spread,
    which must stay well below the cavity linewidth 2 pi / F for every site to build up.
Writes bounds.json.
"""

import sys

import numpy as np

from metatrap.cavity import budget as B
from metatrap.runio import load_config, run_dir, write_json, write_params


def concentric_wavefront_rms(R, r, na, n=81):
    """Rms wavefront error (same units as R) of one reflection from a sphere of radius R centred at
    the origin, for a point source at (r, 0, 0) and a cone of half-angle asin(na) around +z."""
    th_max = np.arcsin(na)
    u = np.linspace(-1, 1, n)
    U, Vv = np.meshgrid(u, u)
    sel = U**2 + Vv**2 <= 1
    sx, sy = (U[sel] * np.sin(th_max)), (Vv[sel] * np.sin(th_max))
    sz = np.sqrt(1 - sx**2 - sy**2)
    P = np.array([r, 0.0, 0.0])
    d = np.stack([sx, sy, sz], axis=1)
    # Intersection with the sphere |P + t d| = R.
    b = d @ P
    t = -b + np.sqrt(b**2 - (P @ P - R**2))
    Q = P + t[:, None] * d
    target = -P  # paraxial image 2C - P
    opl = t + np.linalg.norm(Q - target, axis=1)
    # Remove piston, tilt and focus: fit opl ~ c0 - (Q - target).delta / |Q - target| (linearised image shift).
    uq = (Q - target) / np.linalg.norm(Q - target, axis=1)[:, None]
    A = np.column_stack([np.ones(len(opl)), uq])
    coef, *_ = np.linalg.lstsq(A, opl, rcond=None)
    res = opl - A @ coef
    return float(np.std(res))


def main(cfg_path):
    cfg = load_config(cfg_path)
    out = run_dir(cfg["run_id"])
    write_params(out, cfg)
    lam = cfg["wavelength_um"]
    rep = {"loss_bound": {}, "concentric": []}
    for name, (eta, n) in cfg["intracavity_elements"].items():
        rep["loss_bound"][name] = {"eta": eta, "passes": n, "G_max": float(B.intracavity_element_bound(eta, n))}
    for R in cfg["mirror_radius_mm"]:
        for r in cfg["field_um"]:
            s = concentric_wavefront_rms(R * 1e3, r, cfg["na"])
            s_rt = 2 * s  # two reflections per round trip, same magnitude by symmetry
            rep["concentric"].append({"R_mm": R, "field_um": r, "rms_one_reflection_waves": s / lam,
                                      "rms_roundtrip_waves": s_rt / lam, "strehl_roundtrip": float(np.exp(-(2 * np.pi * s_rt / lam) ** 2)),
                                      "max_finesse_for_rms_below_linewidth_over_10": float(1 / (10 * s_rt / lam)) if s_rt > 0 else None})
    write_json(out / "bounds.json", rep)
    print(rep["loss_bound"])
    for c in rep["concentric"]:
        print(c)


if __name__ == "__main__":
    main(sys.argv[1])
