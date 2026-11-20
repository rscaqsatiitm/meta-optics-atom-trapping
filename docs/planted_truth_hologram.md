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

## Follow-up (2026-11-20)

**Finding 3 fixed.** `evaluate.disc_distances` (now shared by `spot_metrics` and
`tiled._site_metrics`) pairs each sub-pixel offset with its own axis. The new unit test
(tests/test_spot_metrics.py, thresholds declared first) failed on the legacy formula (disc
centroid (-0.245, 0.348) px for a peak at (0.3, -0.2) px) and passes with the fix (centroid
error 0.048 px against 0.05 px; the margin is the pixelation of a 15-pixel disc). Impact on the
same intensity maps, legacy against fixed (scripts/disc_fix_impact.py,
results/disc_fix_impact/metrics.json; decision rule declared first: erratum if a reported
quantity changes by more than its quoted precision, 0.0005 for efficiencies):

| Case | Array efficiency | Core efficiency | Cell-power non-uniformity | Max per-site change, core / cell |
|---|---|---|---|---|
| 10^5-site design, round 26, unfolded | 0.80896 -> 0.80916 (+1.9e-4) | 0.66387 -> 0.66354 (-3.3e-4) | 4.472 -> 4.480 % | 1.8 % / 1.7 % |
| 16 x 16 render-check design | +1.5e-4 | +4.9e-5 | -0.004 points of % | 0.35 % / 0.70 % |
| Planted 3 x 3 | +8e-5 (relative) | +2e-6 (relative) | -0.013 points of % | 0.08 % / 0.08 % |

No reported number changes at its quoted precision, so no erratum. Per-site powers move by up
to 1.8 %, but they are not reported individually.

**Finding 2 resolved.** With one measurement chain for both fields (same window, upsampling and
estimator; thresholds unchanged), AS and RS agree: peak 4.3e-3, position 0.0069 um, core power
7.2e-5. The negative controls (RS input shifted by one pixel, RS at z x 1.01) still fail on the
new chain.

**Finding 1 diagnosed** (scripts/planted_truth_diagnostics.py, results/planted_truth_diagnostics;
interpretation rules declared first):

| Diagnostic | Result | Rule |
|---|---|---|
| D1 field at the nominal site points against the plant (8 um pitch) | 0.61 % | <= 1 %: met. The design delivers what the WGS optimises. |
| D3 three pupil-phase feedback rounds on the measured peaks | peak 1.17 -> 0.21 %; core power 6.6 -> 5.9 %; position 0.129 -> 0.130 um | peak rule met; recovery still fails, so the rest is the metric definition |
| D2 pitch 16 um, same plant | position 0.030 um, core power 4.1 % | <= 0.05 um and <= 5 %: met. Neighbour interference is the dominant cause. |

Hence the recovery failure is not a defect of the design or of the propagation: at 8 um pitch
(3.6 Airy radii at NA 0.3) the Airy tails of neighbouring spots shift the intensity maxima by
0.13 um and add to each core integral, while the WGS controls only the field at the nominal
points, and feedback on the peak values does not move the peaks. The full-aperture designs sit
at a similar pitch-to-spot ratio (4 um at NA 0.55, 3.4 Airy radii) and report ~0.09 um rms
position errors, the same effect (Provisional). Not explained by these rules: at 16 um the
point and peak power errors rise to 3.1 % (0.6 % and 1.2 % at 8 um); left open.

Decision needed before the recovery tests are revised: which quantities the design claims. If
the trap depth (peak intensity) is the claim, the recovered quantity is the peak power, which
the design meets after feedback (0.21 %). If trap positions matter at the 0.05 um level, the
design must control the location of the maxima (e.g. position terms in the feedback), which it
does not today; the recovery tolerances would then stay and the design would change.

## Revision for the trap-depth claim (2026-11-20)

The PI decided that the design claims trap depth. tests/test_planted_hologram.py now recovers the
normalised peak intensity after the production feedback (three pupil-phase rounds, gain 1.0, on
the AS-measured peaks) with DEPTH_TOL = 1 %, measured by the production AS route (in-sample) and
by the RS route with its own estimator (independent of the loop and of the production metrics).
Positions and core powers are reported, not asserted. Provenance of the threshold: chosen after
diagnostics D1 to D3 had shown 1.2 % (AS, no feedback) and 0.21 % (AS, three rounds); the RS value
after feedback had not been computed; the level is that of the claim.

| Quantity | No feedback | After 3 rounds | Tolerance |
|---|---|---|---|
| Trap-depth error, AS (in-sample) | 1.17 % | 0.21 % | 1 %: PASS |
| Trap-depth error, RS (independent) | - | 0.56 % | 1 %: PASS |
| Position error, AS / RS (reported only) | 0.129 um | 0.130 / 0.124 um | - |
| Core-power error, AS / RS (reported only) | 6.6 % | 5.9 % / 5.9 % | - |
| AS against RS, one measurement chain: peak / position / core power | - | 0.43 % / 0.007 um / 7e-5 | 1 % / 0.01 um / 1 %: PASS |

Negative controls (all fail as required): planted powers transposed x <-> y; one planted power
changed by 2 %; RS input shifted by one pixel; RS at z x 1.01.

The independent route confirms the trap-depth claim with a factor ~1.8 margin. Its own resolution
is set by the AS-RS difference (0.43 % in absolute peak, 0.35 points between the two recovery
values), so the independent route cannot confirm relative trap depths below ~0.5 % in this
case; the in-loop 0.21 % is not independently confirmed at that level. Trap positions are not
part of the claim; they are offset by up to 0.13 um here by neighbour interference and are not
controlled by the design.
