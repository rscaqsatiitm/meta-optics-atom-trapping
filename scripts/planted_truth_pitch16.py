"""Why the planted-truth trap depths are less accurate at 16 um than at 8 um pitch (no feedback).

Usage: python scripts/planted_truth_pitch16.py

scripts/planted_truth_diagnostics.py found the field at the nominal sites within 0.61 % of the
planted powers at 8 um pitch but 3.07 % at 16 um, with the pupil model and the exact pupil sum
both reproducing the plant to ~1e-7 at both pitches, and with the render upsampling (8 or 16),
the AS grid (512 or 1024) and the measurement window making <= 0.2 points of difference. Both
designs use the same pupil grid (the 64-point floor of PupilGrid.for_sites), whose focal period
is 4.3 times the field extent at 8 um but 2.1 times at 16 um. This script varies the guard
(focal period over field extent) at both pitches with the same design path and reports the
error of the full-pixel field at the nominal site points. A diagnostic, no pass/fail.
Output: results/planted_truth_pitch16/metrics.json.
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


def main():
    stack = LayerStack.from_list([[T.Z, "air"]])
    radius = stack.aperture_for_na(T.NA, T.LAM)
    rows = []
    for pitch, guard, w in ((16.0, 1.6, 128), (16.0, 3.2, 128), (16.0, 4.3, 128), (8.0, 1.6, 96), (8.0, 8.6, 96)):
        sites, plant = PD.plant_for(pitch)
        field = 2 * np.abs(sites).max() + pitch
        g = PupilGrid.for_sites(T.LAM, pitch, 4, T.NA, field, guard=guard)
        ps = X.PupilSampling.build(g, stack, T.NA, T.ILL)
        tab = X.OPLTable(stack, T.LAM, radius + 2 * np.abs(sites).max() + 20)
        K = X.kernels(ps, sites, tab, T.LAM)
        d = X.wgs_exact(ps, K, n_iter=80, n_fix=40, target=np.sqrt(plant))
        holo = X.render_grid(ps, stack, sites, d.coeffs, T.LAM, T.NA, kind="exact", table=tab, upsample=8)
        f = metasurface_field(T.N, T.DX, stack, T.LAM, T.NA, ps, d.psi, T.ILL, holo=holo)
        foc = propagate(f.copy(), T.DX, T.LAM, stack)
        c0 = T.N // 2 - w // 2
        I = np.abs(fourier_upsample(foc[c0 : c0 + w, c0 : c0 + w], 4)) ** 2
        cU, h = I.shape[0] // 2, T.DX / 4
        pt = np.array([I[cU + int(round(x / h)), cU + int(round(y / h))] for x, y in sites])
        rows.append({"pitch_um": pitch, "guard": guard, "n": g.n, "period_over_field": g.period_um / field,
                     "pupil_samples": len(ps.amp), "point_power_error": rel(pt, plant)})
        print(rows[-1], flush=True)
    out = REPO_ROOT / "results" / "planted_truth_pitch16"
    out.mkdir(parents=True, exist_ok=True)
    write_json(out / "metrics.json", {"rows": rows})


if __name__ == "__main__":
    main()
