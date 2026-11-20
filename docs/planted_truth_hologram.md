# Planted truth for the hologram design path: result (2026-11-19)

Test: tests/test_planted_hologram.py (tolerances declared in its docstring, committed before
the first run). Case: 3 x 3 sites, 8 um
pitch, planted powers 0.58 to 1.50 (asymmetric), free space 200 um, NA 0.3, 1075 nm, 0.5 um
pixels, production design path (exact kernels, wgs_exact with targets, 80 iterations,
render_grid x 8, metasurface_field) and two propagators: production angular spectrum (AS) with
production spot metrics, and a direct Rayleigh-Sommerfeld sum (RS) written in the test.

| Check | Tolerance (declared first) | AS | RS | Verdict |
|---|---|---|---|---|
| Spot position against planted site | 0.05 um | 0.129 um | 0.124 um | FAIL |
| Normalised core power against planted power | 0.05 | 0.066 | 0.066 | FAIL |
| AS against RS: peak intensity (absolute) | 0.01 | 4.6e-3 | | PASS |
| AS against RS: position | 0.01 um | 0.0138 um | | FAIL |
| AS against RS: core power | 0.01 | 7.2e-4 | | PASS |
| Negative controls (x <-> y transposed plant; RS input shifted one pixel; RS at z x 1.01) | must fail | all fail | all fail | PASS |

The same numbers were obtained in CI (GitHub Actions run 44585582210) and locally.

Findings (reported, not fixed; no reported number was changed):
1. The design pipeline does not recover the planted array within the declared tolerances, and
   both propagators agree on this (positions 0.12 to 0.13 um, core powers 6.6 %). The pupil model
   of the same design reports a power error of 1.4e-7, so the gap is between the pupil-model
   design and the full-pixel field, i.e. the design residual that the feedback loops of
   docs/scaling_100k.md exist to remove. Peak-based powers are within 1.2 % (diagnostic computed
   after the run, not a declared check). The position offsets are of the size of the reported
   position errors of the full-aperture designs (~0.09 um rms), consistent with the WGS
   controlling the field at the nominal site points rather than the location of the intensity
   maximum (Provisional: not separated from other causes). Before tightening anything, decide
   whether the design claim is the field at the site point or the peak location and power, and
   whether a planted-truth check of a no-feedback design should use the feedback loop.
2. AS and RS agree on the field (diagnostic after the run: intensity maps on a common
   0.125 um grid differ by at most 0.52 % of the peak; summed intensity ratio 0.99995), but the
   declared position-agreement check fails at 0.0138 um against 0.01 um. The two use different
   peak estimators (2D quadratic fit in spot_metrics, 1D fits in the test) on a 0.125 um grid;
   the check as declared therefore tests the estimators as much as the propagators.
3. Found while reading spot_metrics (hologram/evaluate.py, line 56): the core and cell discs are
   centred with the sub-pixel offsets swapped between the axes (`yy, xx = np.mgrid[...]`, then
   `hypot(xx - ox, yy - oy)` with ox the first-axis offset). The displacement is at most
   0.7 grid step; its effect on reported powers was not quantified. Not changed.
