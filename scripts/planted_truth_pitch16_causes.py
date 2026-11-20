"""The remainder of the 16 um pitch excess at matched pupil sampling: aperture rim or field angle?

Usage: python scripts/planted_truth_pitch16_causes.py

scripts/planted_truth_pitch16.py: at a focal period of 4.4 times the field extent, the full-pixel
field at the nominal sites misses the planted powers by 1.63 % at 16 um pitch against 0.61 % at
8 um (no feedback; same design path as tests/test_planted_hologram.py). Two candidates, each
tested at a matched period-to-field ratio (~4.4):
  (a) the hard aperture rim (partially covered rim cells of the pupil model against the sharp
      edge of the metasurface): 16 um pitch with a raised-cosine edge of 10 % of the aperture
      radius (illumination top_hat, edge_um);
  (b) the field angle of the outer sites (6.4 deg at the 16 um corners): the 8 um array displaced
      off axis to centre (16, 16) um, so that its outer sites reach the same or larger angles.
Interpretation rules, written before the first run:
  (a) is supported if the soft edge brings the 16 um error to <= 1.0 % (from 1.63 %);
  (b) is supported if the displaced 8 um array misses by >= 1.3 % (from 0.61 % on axis).
Neither, both, or one may hold; a rule not met leaves that candidate unsupported, not excluded.
The 16 um hard-edge and 8 um on-axis cases are rerun as references.
Output: results/planted_truth_pitch16/causes.json.
"""

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


def rel(a, plant):
    q = plant / plant.mean()
    return float(np.max(np.abs((a / a.mean()) / q - 1)))


def case(pitch, centre, ill, ratio=4.4):
    sites0, plant = PD.plant_for(pitch)
    sites = sites0 + np.array(centre)
    stack = LayerStack.from_list([[T.Z, "air"]])
    radius = stack.aperture_for_na(T.NA, T.LAM)
    extent = 2 * np.abs(sites).max() + pitch  # field extent seen by the pupil grid (from the axis)
    g = PupilGrid.for_sites(T.LAM, pitch, 4, T.NA, extent, guard=ratio)
    ill = dict(ill)
    if "edge_frac" in ill:
        ill["edge_um"] = ill.pop("edge_frac") * radius
    ps = X.PupilSampling.build(g, stack, T.NA, ill)
    tab = X.OPLTable(stack, T.LAM, radius + 2 * np.abs(sites).max() + 20)
    K = X.kernels(ps, sites, tab, T.LAM)
    d = X.wgs_exact(ps, K, n_iter=80, n_fix=40, target=np.sqrt(plant))
    holo = X.render_grid(ps, stack, sites, d.coeffs, T.LAM, T.NA, kind="exact", table=tab, upsample=8)
    f = metasurface_field(T.N, T.DX, stack, T.LAM, T.NA, ps, d.psi, ill, holo=holo)
    foc = propagate(f.copy(), T.DX, T.LAM, stack)
    h = T.DX / 4
    w = int(2 * np.ceil((np.abs(sites).max() + pitch) / T.DX)) + 8
    c0 = T.N // 2 - w // 2
    I = np.abs(fourier_upsample(foc[c0 : c0 + w, c0 : c0 + w], 4)) ** 2
    cU = I.shape[0] // 2
    pt = np.array([I[cU + int(round(x / h)), cU + int(round(y / h))] for x, y in sites])
    angle = float(np.degrees(np.arctan(np.hypot(sites[:, 0], sites[:, 1]).max() / T.Z)))
    return {"pitch_um": pitch, "centre_um": list(centre), "illumination": ill, "n": g.n, "period_over_extent": g.period_um / extent,
            "max_field_angle_deg": angle, "point_power_error": rel(pt, plant)}


def main():
    res = {
        "ref_16_hard": case(16.0, (0, 0), {"kind": "top_hat"}),
        "a_16_soft_edge": case(16.0, (0, 0), {"kind": "top_hat", "edge_frac": 0.10}),
        "ref_8_axis": case(8.0, (0, 0), {"kind": "top_hat"}),
        "b_8_off_axis": case(8.0, (16.0, 16.0), {"kind": "top_hat"}),
    }
    for k, v in res.items():
        print(k, {a: (round(b, 4) if isinstance(b, float) else b) for a, b in v.items()}, flush=True)
    res["a_supported"] = res["a_16_soft_edge"]["point_power_error"] <= 0.010
    res["b_supported"] = res["b_8_off_axis"]["point_power_error"] >= 0.013
    print("rim (a) supported:", res["a_supported"], "; field angle (b) supported:", res["b_supported"])
    out = REPO_ROOT / "results" / "planted_truth_pitch16"
    out.mkdir(parents=True, exist_ok=True)
    write_json(out / "causes.json", res)


if __name__ == "__main__":
    main()
