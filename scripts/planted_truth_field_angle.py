"""Mechanism of the field-angle part of the planted-truth error: obliquity of the off-axis foci.

Usage: python scripts/planted_truth_field_angle.py [--break]

scripts/planted_truth_pitch16_causes.py found that, at matched pupil sampling, the error of the
full-pixel field at the nominal sites grows with the field angle of the sites (0.62 % for the
8 um array on axis, 1.63 % for 16 um on axis, 2.39 % for the 8 um array centred at (16, 16) um;
no feedback). Hypothesis H, written before this script was run:

  The pupil model gives every site the same peak intensity per unit mode power: a single-site
  design has E_j = 1 wherever x_j lies. The physical field at the focus of an exact converging
  spherical wave to x_j is, by the Rayleigh-Sommerfeld integral with the phases cancelled,
      U(x_j) = (1 / 2 pi) sum_pix A(r) z (1 / R_j - i k) / R_j^2 dx^2,  R_j = |(r, 0) - (x_j, z)|,
  which falls with |x_j| (obliquity z / R_j and spreading 1 / R_j). The site-dependent gain
  G_j = |U(x_j)|^2 / mean_j |U(x_j)|^2 is ignored by the design and the evaluation compares
  against the plant directly, so the measured powers carry G_j.

Checks and thresholds (fixed before the first run):
  P0  primitive, no design: a single exact focus at x = (0,0), (16,0), (16,16), (24,24), (32,32)
      um propagated by the production angular-spectrum path; the nominal-point intensity ratio
      I(x) / I(0) agrees with the RS prediction |U(x)|^2 / |U(0)|^2 within 0.2 % at every point,
      and departs from 1 (the pupil-model value) by more than 1 % at (32, 32) um.
  P1  explanation: with the expectation plant_j G_j in place of plant_j, the max point-power
      error of each of the three cases (8 um on axis, 16 um on axis, 8 um at (16, 16)) is
      <= 1.0 %, and the three errors lie within 0.5 points of each other (now 0.62 to 2.39 %).
  P2  fix: designing with the target plant_j / G_j and comparing the full-pixel powers with the
      plant itself meets the same two conditions.
  Consistency: with G = 1 the three errors reproduce causes.json to 0.02 points.
  Declared break (--break): the compensation with the wrong sign (target plant_j G_j) must fail P2.
H is supported if P0, P1 and P2 pass; a failed check leaves H unsupported, not excluded.
Same design path, pupil grids (period 4.4 x field extent) and measurement as
planted_truth_pitch16_causes.py. Output: results/planted_truth_pitch16/field_angle.json
(field_angle_break.json with --break).
"""

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tests"))
import planted_truth_diagnostics as PD  # noqa: E402
import test_planted_hologram as T  # noqa: E402
from metatrap.hologram import exact as X  # noqa: E402
from metatrap.hologram.render import fourier_upsample, metasurface_field  # noqa: E402
from metatrap.hologram.wgs import PupilGrid  # noqa: E402
from metatrap.propagation.angular_spectrum import propagate  # noqa: E402
from metatrap.propagation.layers import LayerStack  # noqa: E402
from metatrap.runio import REPO_ROOT, write_json  # noqa: E402

P0_TOL, P0_MIN_DEPARTURE = 0.002, 0.01
P_ERR, P_SPREAD = 0.010, 0.005
CONSIST = 0.0002
CAUSES = {"8_axis": 0.006161203892405931, "16_axis": 0.016255644779854084, "8_off_axis": 0.023921619735622457}
CASES = {"8_axis": (8.0, (0.0, 0.0)), "16_axis": (16.0, (0.0, 0.0)), "8_off_axis": (8.0, (16.0, 16.0))}

STACK = LayerStack.from_list([[T.Z, "air"]])
RADIUS = STACK.aperture_for_na(T.NA, T.LAM)
COORDS = (np.arange(T.N) - T.N // 2) * T.DX
XX, YY = np.meshgrid(COORDS, COORDS, indexing="ij")
INSIDE = np.hypot(XX, YY) <= RADIUS
K0 = 2 * np.pi / T.LAM


def rel(a, ref):
    return float(np.max(np.abs((a / a.mean()) / (ref / ref.mean()) - 1)))


def rs_focus(x):
    """|U(x)|^2 of an exact converging spherical wave to (x, z), top-hat amplitude (RS sum)."""
    R = np.sqrt((XX[INSIDE] - x[0]) ** 2 + (YY[INSIDE] - x[1]) ** 2 + T.Z**2)
    u = np.sum(T.Z * (1 / R - 1j * K0) / R**2) * T.DX**2 / (2 * np.pi)
    return float(abs(u) ** 2)


def point_intensity(foc, pts, margin):
    """Nominal-point intensities on the 4x Fourier-upsampled focal window (as in causes)."""
    h = T.DX / 4
    w = int(2 * np.ceil((np.abs(pts).max() + margin) / T.DX)) + 8
    c0 = T.N // 2 - w // 2
    I = np.abs(fourier_upsample(foc[c0 : c0 + w, c0 : c0 + w], 4)) ** 2
    cU = I.shape[0] // 2
    return np.array([I[cU + int(round(x / h)), cU + int(round(y / h))] for x, y in pts])


def p0():
    pts = np.array([[0, 0], [16, 0], [16, 16], [24, 24], [32, 32]], float)
    meas, pred = [], []
    for x in pts:
        f = np.zeros((T.N, T.N), np.complex64)
        R = np.sqrt((XX[INSIDE] - x[0]) ** 2 + (YY[INSIDE] - x[1]) ** 2 + T.Z**2)
        f[INSIDE] = np.exp(-1j * K0 * R).astype(np.complex64)
        foc = propagate(f, T.DX, T.LAM, STACK)
        meas.append(point_intensity(foc, x[None, :], 16.0)[0])
        pred.append(rs_focus(x))
    meas, pred = np.array(meas) / meas[0], np.array(pred) / pred[0]
    dev = np.abs(meas / pred - 1)
    ok = bool(dev.max() <= P0_TOL and abs(meas[-1] - 1) > P0_MIN_DEPARTURE)
    return {"points_um": pts.tolist(), "field_angle_deg": np.degrees(np.arctan(np.hypot(pts[:, 0], pts[:, 1]) / T.Z)).tolist(),
            "as_ratio": meas.tolist(), "rs_ratio": pred.tolist(), "max_abs_dev": float(dev.max()), "pass": ok}


def design_eval(pitch, centre, target_factor):
    sites0, plant = PD.plant_for(pitch)
    sites = sites0 + np.array(centre)
    extent = 2 * np.abs(sites).max() + pitch
    g = PupilGrid.for_sites(T.LAM, pitch, 4, T.NA, extent, guard=4.4)
    ps = X.PupilSampling.build(g, STACK, T.NA, T.ILL)
    tab = X.OPLTable(STACK, T.LAM, RADIUS + 2 * np.abs(sites).max() + 20)
    K = X.kernels(ps, sites, tab, T.LAM)
    d = X.wgs_exact(ps, K, n_iter=80, n_fix=40, target=np.sqrt(plant * target_factor))
    holo = X.render_grid(ps, STACK, sites, d.coeffs, T.LAM, T.NA, kind="exact", table=tab, upsample=8)
    f = metasurface_field(T.N, T.DX, STACK, T.LAM, T.NA, ps, d.psi, T.ILL, holo=holo)
    foc = propagate(f.copy(), T.DX, T.LAM, STACK)
    return sites, plant, point_intensity(foc, sites, pitch)


def gain(sites):
    u = np.array([rs_focus(x) for x in sites])
    return u / u.mean()


def verdict(errs):
    e = np.array(list(errs.values()))
    return bool(e.max() <= P_ERR and e.max() - e.min() <= P_SPREAD)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--break", dest="brk", action="store_true", help="compensate with the wrong sign (P2 must fail)")
    a = ap.parse_args()
    res = {"thresholds": {"P0_tol": P0_TOL, "P0_min_departure": P0_MIN_DEPARTURE, "err": P_ERR, "spread": P_SPREAD, "consistency": CONSIST}}
    if not a.brk:
        res["P0"] = p0()
        print("P0", res["P0"], flush=True)
    raw, expl, fix = {}, {}, {}
    for name, (pitch, centre) in CASES.items():
        sites, plant, pt = design_eval(pitch, centre, 1.0)
        G = gain(sites)
        raw[name], expl[name] = rel(pt, plant), rel(pt, plant * G)
        _, _, pt2 = design_eval(pitch, centre, G if a.brk else 1.0 / G)
        fix[name] = rel(pt2, plant)
        res[name] = {"sites_um": sites.tolist(), "gain": G.tolist(), "gain_spread": float(G.max() - G.min()),
                     "err_raw": raw[name], "err_vs_plant_times_gain": expl[name], "err_compensated_design": fix[name]}
        print(name, {k: (round(v, 5) if isinstance(v, float) else None) for k, v in res[name].items() if not isinstance(v, list)}, flush=True)
    res["consistency_pass"] = bool(all(abs(raw[k] - CAUSES[k]) <= CONSIST for k in CASES))
    res["P1_pass"] = verdict(expl)
    res["P2_pass"] = verdict(fix)
    if a.brk:
        res["note"] = "declared break: target plant * G (wrong sign); P2 must fail"
    else:
        res["H_supported"] = bool(res["P0"]["pass"] and res["P1_pass"] and res["P2_pass"])
    print({k: v for k, v in res.items() if k.endswith("pass") or k == "H_supported"})
    out = REPO_ROOT / "results" / "planted_truth_pitch16"
    out.mkdir(parents=True, exist_ok=True)
    write_json(out / ("field_angle_break.json" if a.brk else "field_angle.json"), res)


if __name__ == "__main__":
    main()
