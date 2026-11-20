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

## Result

(Appended after the run.)
