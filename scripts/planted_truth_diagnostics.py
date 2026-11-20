"""Diagnostics for the planted-truth recovery failure of the hologram design path.

Usage: python scripts/planted_truth_diagnostics.py

tests/test_planted_hologram.py designs a 3 x 3 array with planted powers (0.58 to 1.50, 8 um
pitch, NA 0.3, free space 200 um) and recovers positions to 0.12 to 0.13 um (tolerance 0.05 um)
and normalised core powers to 6.6 % (tolerance 5 %), with two independent propagators agreeing,
while the pupil model of the design reports a power error of 1.4e-7. Three diagnostics separate
the candidate causes. All measurements use the production route (metasurface_field, angular
spectrum, fourier_upsample x 4, spot_metrics) unless stated.

Interpretation rules, written before the first run:
  D1  Field at the nominal site points. |E(x_j)|^2 on the upsampled grid (the sites lie on grid
      points), normalised, against the planted powers. If the max relative error is <= 0.01, the
      design delivers what the WGS optimises (the field at the site points) and the recovery
      failure is a matter of what is measured (peak position and core power), not of the design.
  D2  Twice the pitch (16 um), same planted pattern in units of the pitch, same design path.
      If position error <= 0.05 um and core-power error <= 0.05 at 16 um (while failing at
      8 um), interference between neighbouring spots is the dominant cause.
  D3  Pupil-phase feedback at 8 um (as scripts/fullaperture.py, feedback_mode pupil_phase,
      gain 1.0): coefficient amplitudes scaled by sqrt(target / measured peak), normalised,
      phase recomputed as arg(K c), three rounds. If after <= 3 rounds the peak-power error is
      <= 0.01 and the recovery passes the test tolerances (0.05 um, 0.05), the gap is the
      pupil-model residual that feedback removes; if the peak error falls but the core-power or
      position error does not, the remaining gap is the metric definition (D1, D2).
These rules only label causes; they change no threshold of the test.
Output: results/planted_truth_diagnostics/metrics.json.
"""

import sys
from pathlib import Path

import numpy as np

from metatrap.hologram import evaluate as E
from metatrap.hologram import exact as X
from metatrap.hologram.render import fourier_upsample, metasurface_field
from metatrap.hologram.wgs import PupilGrid
from metatrap.propagation.angular_spectrum import propagate
from metatrap.propagation.layers import LayerStack
from metatrap.runio import REPO_ROOT, write_json

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tests"))
import test_planted_hologram as T  # noqa: E402  (plant definition and constants of the test)


def plant_for(pitch):
    c = np.array([-pitch, 0.0, pitch])
    X_, Y_ = np.meshgrid(c, c, indexing="ij")
    sites = np.column_stack([X_.ravel(), Y_.ravel()])
    u, v = sites[:, 0] / pitch, sites[:, 1] / pitch
    return sites, 1.0 + 0.30 * u + 0.12 * v + 0.5 * ((abs(u) < 1e-9) & (abs(v) < 1e-9))


def build(pitch, sites, plant):
    stack = LayerStack.from_list([[T.Z, "air"]])
    g = PupilGrid.for_sites(T.LAM, pitch, 4, T.NA, 2 * np.abs(sites).max() + pitch, guard=1.6)
    ps = X.PupilSampling.build(g, stack, T.NA, T.ILL)
    radius = stack.aperture_for_na(T.NA, T.LAM)
    tab = X.OPLTable(stack, T.LAM, radius + 2 * np.abs(sites).max() + 20)
    K = X.kernels(ps, sites, tab, T.LAM)
    d = X.wgs_exact(ps, K, n_iter=80, n_fix=40, target=np.sqrt(plant))
    return stack, ps, tab, K, d


def measure(stack, ps, tab, sites, coeffs, psi, pitch, w):
    holo = X.render_grid(ps, stack, sites, coeffs, T.LAM, T.NA, kind="exact", table=tab, upsample=8)
    f = metasurface_field(T.N, T.DX, stack, T.LAM, T.NA, ps, psi, T.ILL, holo=holo)
    foc = propagate(f.copy(), T.DX, T.LAM, stack)
    c0 = T.N // 2 - w // 2
    I = np.abs(fourier_upsample(foc[c0 : c0 + w, c0 : c0 + w], 4)) ** 2
    per = E.spot_metrics(I, T.DX / 4, sites, pitch, T.LAM, T.NA, 1.0)["per_site"]
    cU, h = I.shape[0] // 2, T.DX / 4
    point = np.array([I[cU + int(round(x / h)), cU + int(round(y / h))] for x, y in sites])
    return per, point


def errors(per, point, sites, plant):
    q = plant / plant.mean()
    rel = lambda a: float(np.max(np.abs((a / a.mean()) / q - 1)))  # noqa: E731
    return {"position_um": float(np.max(np.hypot(per["x"] - sites[:, 0], per["y"] - sites[:, 1]))),
            "core_power": rel(per["p_core"]), "peak_power": rel(per["peak"]), "point_power": rel(point)}


def main():
    res = {}
    # D1 and the 8 um baseline.
    sites, plant = plant_for(8.0)
    stack, ps, tab, K, d = build(8.0, sites, plant)
    per, point = measure(stack, ps, tab, sites, d.coeffs, d.psi, 8.0, 96)
    res["pitch8_baseline"] = errors(per, point, sites, plant)
    res["D1_point_power_error"] = res["pitch8_baseline"]["point_power"]
    res["D1_rule_met"] = res["D1_point_power_error"] <= 0.01
    print("8 um baseline", res["pitch8_baseline"], flush=True)
    # D3: pupil-phase feedback on the measured peaks, gain 1.0, three rounds.
    rounds = []
    c, psi = d.coeffs.astype(np.complex64), d.psi
    for r in range(1, 4):
        pk = per["peak"] / per["peak"].mean()
        corr = (plant / plant.mean()) / pk
        c = (c * np.sqrt(corr)).astype(np.complex64)
        psi = np.angle(K @ c).astype(np.float64)
        per, point = measure(stack, ps, tab, sites, c, psi, 8.0, 96)
        e = errors(per, point, sites, plant)
        rounds.append(e)
        print("feedback round", r, e, flush=True)
    res["D3_feedback_rounds"] = rounds
    best = rounds[-1]
    res["D3_peak_rule_met"] = best["peak_power"] <= 0.01
    res["D3_recovery_after_feedback"] = best["position_um"] <= T.POS_TOL and best["core_power"] <= T.POW_TOL
    # D2: 16 um pitch.
    sites2, plant2 = plant_for(16.0)
    stack2, ps2, tab2, K2, d2 = build(16.0, sites2, plant2)
    per2, point2 = measure(stack2, ps2, tab2, sites2, d2.coeffs, d2.psi, 16.0, 128)
    res["D2_pitch16"] = errors(per2, point2, sites2, plant2)
    res["D2_rule_met"] = res["D2_pitch16"]["position_um"] <= 0.05 and res["D2_pitch16"]["core_power"] <= 0.05
    print("16 um", res["D2_pitch16"], flush=True)
    out = REPO_ROOT / "results" / "planted_truth_diagnostics"
    out.mkdir(parents=True, exist_ok=True)
    write_json(out / "metrics.json", res)
    print({k: v for k, v in res.items() if k.endswith("rule_met") or k.endswith("feedback")}, flush=True)


if __name__ == "__main__":
    main()
