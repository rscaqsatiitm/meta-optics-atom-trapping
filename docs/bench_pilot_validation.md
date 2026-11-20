# Validation of the bench_gpu.py pilot projection

`scripts/bench_gpu.py --pilot CONFIG` projects the time and memory of a full-aperture run from
scaled-down pilots (CLAUDE.md, compute planning). Before trusting it for the 10^6-site run, it is
checked against a run whose cost was measured: the unfolded 10^5-site feedback run
(configs/fullaperture_100k_full.yaml, results/fullaperture_100k_full), CPU only, this container
(4 cores, 15 GB).

Measured (the run's log): one feedback round (one kernel product plus one render, fold, unfolded
propagation and tiled metrics) took 724 to 787 s; the run completed within the container's 15 GB
(peak not recorded); the complex64 fold buffer alone is 36,864^2 x 8 B = 10.9 GB.

Acceptance criteria, written before the pilot output was seen:
1. Projected time per feedback round (kernel product + render and fold + propagation and metrics)
   within +/-30 % of 755 s, i.e. 530 to 980 s.
2. Projected peak memory between 10.9 GB (the fold buffer) and 15 GB (the run fitted in the
   container).
3. With three pilot sizes, the largest relative fit residual of each model <= 10 %.
If a criterion fails, the projection is not used for the 10^6 run until the cause is found.

## Result of the first validation: FAIL

Pilots at scales 0.1 / 0.2 / 0.3 (results/bench_pilot, first version of the pilot mode):

| Criterion | Projected | Required | Verdict |
|---|---|---|---|
| 1. Time per feedback round | 393 s (kernel product 51 + render and fold 145 + propagation and metrics 196) | 530 to 980 s | FAIL |
| 2. Peak memory | 11.5 GB | 10.9 to 15 GB | PASS |
| 3. Fit residuals | render 183 %, propagation 6.5 %, memory 0.1 % | <= 10 % | FAIL (render) |

Diagnosis: (i) the render-and-fold model a + b n^2 used the padded grid size n; the 0.2 pilot was
padded to twice its fold window, while the render cost follows the pixels inside the aperture,
hence the 183 % residual; (ii) the pilot rendered ideal pixels, whereas the measured run renders
library pillars (lookup per pixel). The projection is therefore not used.

Changes: the render model becomes a + b N_ap, with N_ap = pi (R / dx)^2 the aperture pixel count,
and the pilot renders library pillars whenever the config uses them (feedback_pillars: library or
pillar cases). A re-check against the same 755 s run is no longer blind (the expected value can be
estimated by hand from the pilot numbers above), so the fixed projector is validated blind against
a different measured run, declared here before it is run:

Blind validation (declared before running): configs/fullaperture_100k_pplib_g08.yaml (10^5 sites,
2.048 mm fold window, complex128 fold buffer, library pillars). Measured: one evaluation (render,
fold, propagation, tiled metrics; the `seconds` of results/fullaperture_100k_pplib_g08,
library_lpa) took 473 s. Criterion: projected render-and-fold plus propagation-and-metrics time per
evaluation within +/-30 % of 473 s (331 to 615 s), with fit residuals <= 10 %. The rerun of the
100k_full case is reported as well, labelled post hoc.
