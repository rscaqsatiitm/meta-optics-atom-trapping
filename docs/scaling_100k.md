# Scaling to 100,000 sites: does a shorter working distance help?

Technical note. Question (USER): would moving the metasurface closer to the atoms (in vacuum,
e.g. 3 mm) help to reach 100,000 traps? All numbers are SIM (review before submission) or
DERIVED. Code: scripts/scaling_100k.py, configs/scaling_100k.yaml; supercells: scripts/multiwave_supercells.py configs/scaling_sc_*.yaml;
results: results/scaling_100k/{scaling,table}.json (table.json from scripts/scaling_100k_table.py
holds every number of the Results table), results/scaling_100k_sc_si_rich.

## Result

No. A shorter working distance does not remove any of the constraints that bind at
100,000 sites, and it tightens two of them. The trapping hologram tolerates short distances
well, because exact-kernel design has no field aberration. At 3 mm the corner site of a
316 x 316 array at 4 um pitch keeps 0.89 of the axial depth, which costs 4 % more power. The
binding constraints are three. (i) Total power: ~210 W for 0.3 mK, ~126 W for 0.18 mK. This is
independent of the distance, but the intensity on the metasurface grows as 1 / WD^2, to
1.7 kW/cm^2 at 3 mm. (ii) The field of the imaging and addressing lens: at 3 mm it images only
~25,000 sites at 4 um pitch. (iii) The design computation (~30 to 50 GB of kernels), which is
also independent of the distance. A working distance of >= 10 mm, the present 15.5 mm cell
geometry included, supports 100,000 sites on both sides.

## Method

Trapping hologram. The metasurface lies at a distance WD from the atom plane, in vacuum, at
NA 0.55 for the axial site and 1075 nm. Its aperture radius is R = 0.659 WD. For a site at
radius r we compute the ideal exact-kernel spot (the full aperture phased to focus at that
site) by a Rayleigh-Sommerfeld sum over a 360 x 360 aperture grid, with the maximum searched
laterally and axially. We report: (i) the peak relative to the axial site at equal power;
(ii) the 1/e^2 half-widths; (iii) the local deflection s the meta-atoms must produce; (iv) the
meta-atom efficiency factor of the aperture, (sum w sqrt(eta) / sum w)^2 with w = z / L^2, for
eta(s) from supercell RCWA of one pillar family (Si-rich SiN p 0.50 um, h 1.5 um, the 1075 nm
selection; eta_LPA 0.963): library_analysis_1075 for s <= 0.54 and a new N = 3 supercell
(results/scaling_100k_sc_si_rich_p050), 0.70 at s = 0.72. (An earlier version mixed this family
with p 0.45 um, h 1.3 um supercells, which gave a non-monotonic curve, and weighted eta by w;
the site factors change by < 1.3 %.) The relative depth
is D(r) = I_rel x eta(r) / eta(0). Uniform depth over a square array needs power in
proportion to 1 / D_j, so the array power factor is F = mean(1 / D_j) (DERIVED). The absolute
power is scaled from the simulated v2 design: 0.369 mK for 20,164 sites at 50 W with library
pillars, i.e. 148.7 mK x site / W at NA 0.55, including its array efficiency of 0.82 (power
budget with the corrected ideal core fraction 0.783 and vector factor 0.917; previously 0.355 mK
and 143). Table: results/scaling_100k/table.json, scripts/scaling_100k_table.py.

Imaging lens. For each WD the air-spaced aplanatic doublet (NA 0.5, 852 nm) is designed, with
the spacing scanned around f. The field radius with a free-focus Strehl ratio >= 0.8 then
gives the number of sites inside the inscribed square, N = 2 r^2 / a^2 (DERIVED).

## Results

| WD (vacuum) | Aperture | Corner depth D, 4 um / 3 um pitch | Power factor F | Corner spot half-widths (radial / tangential) | Max deflection, corner | Power for 0.3 mK (0.18 mK) | Intensity on metasurface | Max absorbing defect (FS, 50 K) | Imaging-lens field radius | Sites imageable, 4 / 3 um pitch |
|---|---|---|---|---|---|---|---|---|---|---|
| 3 mm | 3.95 mm | 0.89 / 0.94 | 1.043 | 0.84 / 0.80 um | 0.69 | 210 W (126 W) | 1,710 W/cm^2 | 4.0 um | 450 um | 25,000 / 45,000 |
| 6 mm | 7.9 mm | 0.97 / 0.98 | 1.011 | 0.80 / 0.80 um | 0.63 | 204 W (122 W) | 415 W/cm^2 | 17 um | 760 um | 73,000 / 129,000 |
| 10 mm | 13.2 mm | 0.99 / 0.99 | 1.004 | 0.80 / 0.78 um | 0.60 | 202 W (121 W) | 148 W/cm^2 | 47 um | 1,160 um | 167,000 / 296,000 |
| 15.5 mm | 20.4 mm | 1.00 / 1.00 | 1.002 | 0.78 / 0.78 um | 0.58 | 202 W (121 W) | 62 W/cm^2 | 112 um | 1,610 um | 324,000 / 577,000 |
| 25 mm | 32.9 mm | 1.00 / 1.00 | 1.001 | 0.78 / 0.78 um | 0.57 | 202 W (121 W) | 24 W/cm^2 | 292 um | 2,020 um | 510,000 / 907,000 |

The axial spot is 0.78 um (Airy 1/e^2 radius 0.41 lambda / NA = 0.80 um). The 100,000-site
array corner lies at 891 um (4 um pitch) or 669 um (3 um pitch). The defect column uses the half-space formula
a = 50 K x kappa / I; for radii approaching the 0.5 mm plate thickness the lateral spreading
lowers it (147 um instead of 184 um at 37 W/cm^2; docs/prefab_checks_report.md, Section 5), so
the 15.5 and 25 mm entries are upper bounds by ~15 to 30 %.

## Full 100,000-site design at the cell geometry

Correction (review 3, 2026-11-10). Every depth non-uniformity, efficiency and axis-light value
in the tables of this section up to "Feedback on the unfolded grid" was measured in a 2.048 mm
fold window of the folded propagation (propagation/folded.py). Folding makes the focal field
periodic with that period, so all light landing outside the window, the WGS design noise
(the coarse pupil grid spans ~2.2 mm here) and the stray light of the library pillars (zero
order, ghost foci, scatter), wraps back onto the sites; a smooth wide field adds coherently over
the (18.4 mm / W)^2 copies. The feedback loops measured in the same window and corrected the
alias pattern, so their reported convergence (0.015 % with ideal pillars, 0.45 % with library
pillars, 3.5 % with 5 nm width noise, and the "library floor") is in-sample and wrong. The same
designs re-measured with larger windows (scripts/fold_window_study.py,
configs/fold_window_100k.yaml, results/fold_window_100k):

| Design (pillars) | 2.048 mm window (as reported) | 4.6 / 6.1 / 9.2 mm windows | Array efficiency, 2.048 / 9.2 mm | Axis light / mean peak, 2.048 / 9.2 mm |
|---|---|---|---|---|
| Round 2 (ideal) | 7.14 % | 7.83 / 7.75 / 7.74 % | 0.857 / 0.840 | 0.102 / 0.060 |
| Round 8, pupil-phase feedback (ideal) | 0.015 % | 6.33 / 6.18 / 6.22 % | 0.857 / 0.841 | 0.108 / 0.060 |
| Round 17, library pillars in the loop (library) | 0.45 % | 6.81 / 6.60 / 6.72 % | 0.828 / 0.809 | 0.106 / 0.055 |
| Round 17, sigma 5 nm width noise | 3.47 % | 6.93 / 6.67 / 6.75 % | 0.802 / 0.778 | 0.120 / 0.053 |

The three larger windows agree to 0.2 percentage points, and the per-site errors of the round-8
and round-17 designs in the 2.048 and 9.2 mm windows are uncorrelated (r = -0.02 and 0.03), so
the in-loop values carry no information about the true errors. On the unfolded grid (36,864^2,
18.4 mm period, so nothing within the aperture footprint is aliased) the round-17 design gives
6.71 %, the 9.2 mm value. The array efficiency was 0.02 high in the 2.048 mm window, as was the
axis light (0.10 to 0.12 against 0.05 to 0.06). The v2 design (20,164 sites, 0.57 mm array) is
not affected: 5.33 / 5.48 / 5.76 % (ideal / library / sigma 5 nm) at 2.048 mm against 5.24 /
5.25 / 5.26 % at 9.2 mm, efficiency within 0.003 (configs/fold_window_v2.yaml,
results/fold_window_v2). The earlier folded validation (results/folded_validation_5mm) used ideal
pixels and an array filling 0.21 of the window, so it could not show the effect. The feedback
was therefore redone on the unfolded grid ("Feedback on the unfolded grid", below); the tables
that follow until then are kept as the record of what was run.

Run: 316 x 316 = 99,856 sites at 4.0 um pitch (array corner 891 um from the axis), NA 0.55,
1075 nm, 18.4 mm metasurface (1.36 x 10^9 pixels), present cell stack, array 50 um beyond the
lens focus. Settings as v2 (configs/fullaperture_100k.yaml; results/fullaperture_100k): 60
D4-symmetric WGS iterations with exact kernels, two feedback rounds of 30 iterations, then
library pillars (Si-rich SiN estimate, p 0.5 um, h 1.6 um) without and with 5 nm rms width
noise. The kernel (517,262 wedge pupil samples x 12,561 fundamental sites, ~50 GB dense) is
evaluated on the fly (src/metatrap/hologram/lazy.py; numba, vectorised sincos, relative
agreement with the dense kernel 1e-5): ~150 s per WGS iteration on 4 cores, 0.5 GB of memory.
The run survived two container restarts through the new checkpoints (every 10 WGS
iterations and before each feedback round). One deviation from v2: the render grid is
upsampled 5 x instead of 8 x, to fit in memory. The pupil-model row and round 0 are taken from
the run log (results/fullaperture_100k/run_log.txt); metrics.json holds rounds 1 and 2 and the
pillar cases.

| Case | Depth non-uniformity (rms) | Min / mean | Array efficiency | Mean waist | Position error (rms) | Interstitial p99 / max | Depth at 50 W | Power for 0.3 mK |
|---|---|---|---|---|---|---|---|---|
| Pupil model (exact kernels) | 0.24 % | - | mode efficiency 0.854 | - | - | - | - | - |
| Ideal phase, round 0 / 1 / 2 | 8.7 / 8.5 / 7.0 % | 0.66 / 0.69 / 0.76 | 0.852 / 0.855 / 0.857 | - | - | - | 0.076 mK | 197 W |
| Library pillars (LPA) | 7.5 % | 0.75 | 0.828 | 0.866 um | 0.095 um | 0.10 / 0.26 | 0.073 mK | 205 W |
| + 5 nm rms width noise | 8.5 % | 0.70 | 0.803 | 0.865 um | 0.098 um | 0.10 / 0.28 | 0.071 mK | 212 W |
| v2 (20,164 sites), library pillars, for comparison | 5.4 % | 0.84 | 0.826 | 0.869 um | 0.083 um | 0.09 / 0.12 | 0.370 mK | - |

What the full design shows (values in the 2.048 mm window; see the correction above for what
changes). (i) The physics scales as predicted above. The array efficiency
(0.83 with library pillars) and the spot size and position are those of the 20,164-site
design. The power for 0.3 mK, 205 W, matches the 202 W estimate. The depth shows no field
dependence: the radially binned mean (100 um bins) stays within 0.993 to 1.002 from the axis to
the corner, and
the outermost ring of sites averages 0.999 (v2: 0.987). (ii) The remaining non-uniformity is
site-to-site scatter, not a field effect. A smooth radial fit (quartic in r) explains 0.01 % of
the variance, and neighbouring sites are uncorrelated (correlation coefficient 0.008)
(results/fullaperture_100k/diagnostics.json, scripts/fullaperture_diagnostics.py). The feedback
converges more slowly than for v2 (8.7 to 7.0 % in two rounds, against 6.3 to 5.3 %). The number of
pupil samples per site is the same in both designs (8 x 103,459 / 20,164 = 41 for v2;
8 x 517,262 / 99,856 = 41 here), because the pupil grid is refined in proportion to the field.
We first attributed the slower convergence, provisionally, to the coarser render grid
(5 x instead of 8 x upsampling of a pupil field that now has 2.2 times finer structure per
axis). A later test rules this out. We re-rendered the final design with a render that
evaluates the hologram field at each metasurface pixel by a windowed-sinc interpolant of the
pupil grid, with no fine grid in memory (hologram/tiled.py). On 16 x 16-site tests (unfolded
grid, ideal pillars) this render matches an exact per-pixel kernel sum to 1.2e-4 rms in per-site
depth (array 50 um beyond focus) and 0.9e-4 (array in focus), against 0.4 % and 1.5 % for the
FFT-upsampled grid (results/tiled_render_check*/metrics.json, scripts/tiled_render_check.py;
tests/test_tiled.py checks the interpolant itself). The agreement requires the coarse pupil grid
to cover the array plus the defocus blur (+-33 um at 50 um and NA 0.55): with the default
guard of 1.25 times the 16 x 16 array diagonal it does not, and the exact render differs by
~6 % rms, because the design controls only the coarse samples. For the 10^5-site array the guard
covers 2.24 mm against 1.79 mm plus 0.07 mm. With it, the round-2 design gives 7.14 % rms (ideal
pillars; results/fullaperture_100k_tiled), against 7.05 % with the 5 x render grid, and 7.6 % against
7.5 % with library pillars (efficiency 0.828 in both). The per-site
depths of the two renders differ by 1.7 % rms (correlation 0.97). Hence the render adds an
uncorrelated ~1.7 % component, which does not change the total in quadrature, and the ~7 %
residual belongs to the design itself. Its origin, given that the pupil model predicts
0.24 %, is unresolved. Candidates are the pupil-model sampling of the aperture edge and too few
feedback rounds. More feedback rounds (~1.5 h each on 4 cores, minutes on a GPU) are the next
test. Further feedback rounds (results/fullaperture_100k_fb, configs/fullaperture_100k_fb.yaml;
ideal pillars, gain 0.5, 30 WGS iterations per round, tiled render in the loop, continuing
from the round-2 state). All of the residual is orbit-to-orbit: the scatter within the
8-site D4 orbits is <= 0.02 % rms, so the orbit targets can in principle remove it.

| Round | Depth non-uniformity rms | Min / mean | Array efficiency | Axis intensity / mean peak |
|---|---|---|---|---|
| 0 (fullaperture_100k) | 8.71 % | 0.66 | 0.852 | 0.020 |
| 1 | 8.52 % | 0.69 | 0.855 | 0.110 |
| 2 (FFT render / tiled render) | 7.05 / 7.14 % | 0.76 | 0.857 | 0.109 / 0.102 |
| 3 | 5.96 % | 0.75 | 0.858 | 0.147 |
| 4 | 5.47 % | 0.77 | 0.859 | 0.162 |
| 5 | 5.07 % | 0.82 | 0.860 | 0.172 |
| 6 | 4.78 % | 0.80 | 0.861 | 0.072 |
| Library pillars after round 6 | 5.48 % | 0.80 | 0.832 | 0.069 |

Convergence diagnostic (per-site depths saved every round). The per-orbit depth errors of
rounds 2 and 3 are nearly uncorrelated (correlation -0.13, regression slope -0.11;
results/fullaperture_100k_fb/diagnostics.json). A linear,
stable response at gain 0.5 would give a slope of +0.5. Hence each feedback round removes the
previous error pattern (with a slight overshoot), but the re-run WGS creates a new pattern of
almost the same size (~5.9 % rms uncorrelated with round 2). Round 4 behaves alike (5.47 %,
slope round 3 to 4: -0.04), and the axis light changes erratically from round to round
(0.15, 0.16, 0.17, then 0.07 in round 6). The mismatch between the pupil
model (0.24 %) and the full simulation is therefore not a fixed per-orbit bias that feedback
can learn; it changes with the design (Provisional, from one pair of rounds). The axis
intensity also grows with the rounds (0.02 to 0.15 of a mean site peak). More rounds of the
same scheme are expected to converge slowly. More effective options are a lower gain with the
WGS restarted closer to the previous solution (fewer iterations), or a feedback that acts on
the pupil phase directly rather than on the WGS targets (not tested).

Pupil-phase feedback (feedback_mode: pupil_phase in scripts/fullaperture.py;
configs/fullaperture_100k_pp.yaml). The WGS-target feedback above re-runs 30 WGS iterations
with free site phases and reset weights, which can carry the design to a different solution
each round. The alternative keeps the current design: the site phases and coefficients c_j are
retained, their amplitudes are scaled by the measured correction (field amplitude ~
(mean / orbit depth)^(g/2)), and the pupil phase is recomputed as arg(K c) with one kernel
product. The design then moves by a small step, so the measured error map stays valid, and a
round costs one simulation plus one product.

Test on the 16 x 16-site check case (configs/tiled_check.yaml, short stack, array 50 um beyond
focus, ideal pillars; depth non-uniformity rms, rounds 0 to 5; configs/tiled_check_fb_*.yaml,
results/tiled_check_fb_*):

| Feedback | Round 0 | 1 | 2 | 3 | 4 | 5 |
|---|---|---|---|---|---|---|
| WGS targets, gain 0.5 | 6.95 % | 2.36 % | 0.45 % | 0.12 % | 0.04 % | 0.01 % |
| Pupil phase, gain 0.5 | 6.95 % | 4.54 % | 3.16 % | 2.32 % | 1.77 % | 1.38 % |
| Pupil phase, gain 1.0 | 6.95 % | 2.57 % | 1.46 % | 0.92 % | 0.62 % | 0.42 % |
| Pupil phase, gain 1.4 | 6.95 % | 2.36 % | 1.06 % | 0.54 % | 0.31 % | 0.19 % |

At this size the WGS-target feedback converges fastest, so its slow convergence at 10^5 sites
is a large-array effect. The pupil-phase feedback is stable and keeps the error map
correlated (round-to-round orbit-error slopes 0.64 to 0.75 at gain 0.5); the phase-only
projection passes ~0.7 of an amplitude correction, hence the higher gain.

100k design, pupil-phase feedback at gain 1.2 from the same round-2 state (A/B comparison
with the WGS-target rounds above):

| Round | Depth non-uniformity rms | Min / mean | Array efficiency | Axis intensity / mean peak | Orbit-error slope |
|---|---|---|---|---|---|
| 2 (start) | 7.14 % | 0.76 | 0.857 | 0.102 | - |
| 3 | 1.77 % | 0.92 | 0.857 | 0.103 | -0.02 |
| 4 | 0.53 % | 0.98 | 0.857 | 0.106 | -0.04 |
| 5 | 0.18 % | 0.99 | 0.857 | 0.107 | -0.04 |
| 6 | 0.068 % | 0.99 | 0.857 | 0.108 | -0.04 |
| 7 | 0.029 % | 0.99 | 0.857 | 0.108 | -0.03 |
| 8 | 0.015 % | 0.999 | 0.857 | 0.108 | - |
| Library pillars after round 8 | 2.69 % | 0.91 | 0.828 | 0.122 | - |

One pupil-phase round removes three quarters of the residual (7.14 to 1.77 % rms, min / mean
0.76 to 0.92) with no loss of efficiency (0.857) and no increase of the axis light (0.103),
where one WGS-target round from the same state gave 5.96 % and raised the axis light to 0.147.
The orbit-error slope of -0.02 means the gain of 1.2 corrected the round-2 pattern almost
exactly; the remaining 1.77 % is uncorrelated with it and is the error introduced by the step
itself, which later rounds act on. Round 4 brings it to 0.53 % rms (min / mean 0.98), i.e.
the pupil-model level is approached in two rounds. Round 5 reaches 0.18 % rms (min / mean
0.99), below the 0.24 % of the pupil model itself. Rounds 6 to 8 continue at a factor ~2 to 3 per
round to 0.015 % rms (min / mean 0.999) with orbit-error slopes of -0.03 to -0.04, unchanged
efficiency 0.857 and axis light 0.108. One orbit (radius 768 um) lagged at 0.9945 in round 7 and
reached 0.9986 in round 8. Hence, with ideal pillars, the design residual of the 10^5-site
hologram is removed by pupil-phase feedback in five to six rounds (~14 min each on 4 cores);
with library pillars the same design gives 2.69 % rms (min / mean 0.91, efficiency 0.828),
against 7.5 % for the round-2 design. The remaining 2.7 % is the pillar library's own phase and
amplitude error under the local-periodic approximation. It is deterministic, so feedback with
the library pillars in the loop (the same scheme, simulated with library_lpa) should remove most
of it; this has not been run. Random width errors (sigma 5 nm) are not correctable this way.

Feedback with library pillars in the loop (feedback_pillars: library;
configs/fullaperture_100k_pplib.yaml). The same pupil-phase feedback, gain 1.2, but each round
simulates the metasurface with the library pillars (local-periodic model) instead of ideal
phase pixels, so the deterministic library error is corrected together with the design
residual. It starts from the round-8 pupil-phase design. On the 16 x 16 check case it converges
from 6.9 to 0.50 % rms in four rounds (orbit-error slopes 0.14 to 0.62), somewhat more slowly
than with ideal pillars. The final design is also evaluated with random pillar width errors
(sigma 5 nm, uncorrelated per pixel), which feedback cannot remove.

| Round | Depth non-uniformity rms (library pillars) | Min / mean | Array efficiency | Axis intensity / mean peak | Orbit-error slope |
|---|---|---|---|---|---|
| 8 (start; ideal-pillar design) | 2.69 % | 0.91 | 0.828 | 0.122 | - |
| 9 | 1.47 % | 0.94 | 0.827 | 0.119 | -0.10 |
| 10 | 1.07 % | 0.95 | 0.828 | 0.113 | -0.26 |
| 11 | 0.91 % | 0.97 | 0.827 | 0.107 | -0.37 |
| 12 | 0.84 % | 0.97 | 0.827 | 0.109 | -0.45 |
| 13 | 0.81 % | 0.97 | 0.827 | 0.113 | - |
| Final design, sigma 5 nm width errors | 3.53 % | 0.85 | 0.802 | 0.127 | - |

With library pillars in the loop the non-uniformity falls from 2.69 to 0.81 % rms in five rounds
(min / mean 0.97, efficiency 0.827, 0.073 mK at 50 W). It levels off at ~0.8 % while the
orbit-error slope grows from -0.10 to -0.45: the response of the library-pillar simulation to a
correction strengthens as the errors shrink, so the fixed gain of 1.2 over-corrects increasingly.
The lower-gain follow-up below separates the two. With independent pillar width errors of
sigma 5 nm the same design gives 3.53 % (efficiency 0.802), against 8.46 % for the original
round-2 design. Subtracting the 0.81 % in quadrature, the random width errors contribute ~3.4 %
rms; they are now the dominant term and set the fabrication requirement on width control.

Lower-gain follow-up (configs/fullaperture_100k_pplib_g08.yaml): library pillars in the loop,
gain 0.8, from the round-13 design.

| Round | Depth non-uniformity rms (library pillars) | Min / mean | Array efficiency | Axis intensity / mean peak | Orbit-error slope |
|---|---|---|---|---|---|
| 13 (start, gain-1.2 design) | 0.81 % | 0.97 | 0.827 | 0.113 | - |
| 14 | 0.56 % | 0.97 | 0.827 | 0.109 | -0.06 |
| 15 | 0.49 % | 0.98 | 0.827 | 0.107 | -0.03 |
| 16 | 0.46 % | 0.98 | 0.827 | 0.108 | -0.06 |
| 17 | 0.45 % | 0.98 | 0.827 | 0.106 | - |
| Final design, sigma 5 nm width errors | 3.47 % | 0.86 | 0.802 | 0.120 | - |

At gain 0.8 the orbit-error slopes stay at -0.03 to -0.06 and the non-uniformity falls from 0.81
to 0.45 % rms in four rounds (min / mean 0.98, efficiency 0.827), levelling off with a
per-round factor of 0.94 at the end. Hence about half of the 0.8 % plateau at gain 1.2 was
overshoot; the remaining ~0.45 % is the floor of this scheme with library pillars, set by the
new error that each correction introduces through the mapping of the phase onto discrete library
pillars (Provisional). With sigma 5 nm width errors the final design gives 3.47 %, i.e. a random
contribution of 3.44 % in quadrature, the same as for the gain-1.2 design.

Feedback on the unfolded grid (configs/fullaperture_100k_full.yaml;
results/fullaperture_100k_full). Pupil-phase feedback with library pillars in the loop, gain
1.0, from the round-17 design, with every round measured on the unfolded 36,864^2 grid
(fold_dtype complex64, the field propagated in place, ~11 GB; ~12.5 min per round on 4 cores).

| Round | Depth non-uniformity rms (library pillars) | Min / mean | Array efficiency | Axis intensity / mean peak | Orbit-error slope |
|---|---|---|---|---|---|
| 18 (round-17 design, unfolded) | 6.71 % | 0.75 | 0.809 | 0.054 | - |
| 19 | 1.46 % | 0.92 | 0.809 | 0.063 | 0.18 |
| 20 | 0.46 % | 0.97 | 0.809 | 0.064 | 0.29 |
| 21 | 0.17 % | 0.99 | 0.809 | 0.064 | 0.34 |
| 22 | 0.077 % | 0.995 | 0.809 | 0.064 | 0.37 |
| 23 | 0.040 % | 0.997 | 0.809 | 0.064 | 0.38 |
| 24 | 0.024 % | 0.997 | 0.809 | 0.064 | 0.33 |
| 25 | 0.017 % | 0.999 | 0.809 | 0.064 | 0.21 |
| 26 | 0.014 % | 0.997 | 0.809 | 0.064 | - |

Final design with fabrication errors (same grid; error models as in configs/checks_v2.yaml):

| Case | Depth non-uniformity rms | Min / mean | Array efficiency | Axis intensity / mean peak | Interstitial p99 / max |
|---|---|---|---|---|---|
| Library pillars (LPA) | 0.014 % | 0.997 | 0.809 | 0.064 | 0.093 / 0.22 |
| Uncorrelated width noise, sigma 5 nm (seeds 1 / 2 / 3) | 0.37 / 0.37 / 0.37 % | 0.982 / 0.984 / 0.985 | 0.778 | 0.063 to 0.065 | 0.093 / 0.22 |
| Correlated width error, 5 nm rms over 50 um | 5.5 % | 0.76 | 0.790 | 0.105 | 0.099 / 0.27 |
| Combined (100 um stitching fields, 20 nm stitch, 2 nm field bias, 3 nm over 5 um, 3 nm random) | [pending: run in progress] | | | | |

With an alias-free measurement the feedback converges as it did on the 16 x 16 test: the error
falls by a factor 3 to 4 per round at first and ~1.5 at the end, with positive orbit-error
slopes of 0.18 to 0.38 (each correction removes most of the previous pattern without overshoot),
and no loss of efficiency (0.809). The design residual and the deterministic library error are
both removed (0.014 % rms, in-sample but without aliasing); the "floor set by the discrete
library pillars" of the folded runs was an alias. The remaining limits are therefore (i) the
local-periodic model of the pillars, which these runs cannot test (a full-wave check of a
metasurface section is needed), and (ii) fabrication. Uncorrelated pixel-level width noise of
5 nm costs only 0.37 % rms and 0.031 in efficiency (identical for three draws), against 3.4 %
claimed from the folded runs. Smooth width errors of the same rms over 50 um (dose drift,
proximity effect) cost 5.5 % rms and drop the weakest site to 0.76 of the mean, because they
act like a slowly varying phase error that moves light between neighbouring regions of the
array. The fabrication specification is therefore on correlated width errors: to stay below
~1 % rms they must be held to ~1 nm over tens of um (DERIVED, assuming the rms scales linearly
with the error amplitude, untested), or measured after fabrication and corrected in a second
design iteration.

(iii) The brightest interstitial point is 0.22 of the mean site depth on the unfolded grid
(final design, at 568 um from the axis; 0.26 to 0.28 in the 2.048 mm window), against 0.12 for
v2, with the 99th percentile at 0.093. At 0.3 mK this is a ~65 uK spurious maximum at isolated
points. Interstitial suppression in the feedback (dark targets) would address it.

Conclusion. A single 18.4 mm metasurface at the present cell geometry generates 99,856 traps
with an array efficiency of 0.81 (library pillars, unfolded grid; 0.83 in the earlier folded
runs, which overstated it by 0.02). The folded runs with a 2.048 mm window aliased stray light
and design noise onto the array, and all their feedback results were in-sample artefacts: the
true non-uniformity of those designs is 6.2 to 7.7 %. Pupil-phase feedback measured on the
unfolded grid, with library pillars in the loop, removes the design residual and the
deterministic library error in eight rounds (6.71 to 0.014 % rms, min / mean 0.997, efficiency
unchanged). Uncorrelated 5 nm width noise then adds only 0.37 %, but smooth width errors of 5 nm
over 50 um give 5.5 % (min / mean 0.76): the uniformity is set by correlated fabrication errors
and by the validity of the local-periodic pillar model, which a full-wave check must confirm.
0.3 mK traps need ~210 W [pending: power budget of the unfolded run] of 1075 nm (Si-rich SiN on fused silica is then outside its thermal
limits; use LPCVD SiN or sapphire). A fold window of 9.2 mm reproduces the unfolded result for
this design (6.72 against 6.71 %), so the 10^6-site run (configs/fullaperture_1M.yaml) uses it.

## Discussion

Why the hologram does not care much. In the exact-kernel design every site receives its own
spherical wave through the stack, so there is no field aberration to grow with angle; this
is what made the 20,000-site design possible. At short distances the corner sites lose only
through obliquity, which shrinks the projected aperture: the peak falls to 0.90 at 891 um
for WD = 3 mm, the spot elongates radially by 5 %, and the steeper deflections
(s up to 0.69, 13 % of the aperture above 0.6) cost 2 % in meta-atom efficiency. Hence
the trapping side alone would tolerate 3 mm.

Why closer still does not help. (i) Power. The depth per site is set by NA and power, not by
the distance. 100,000 sites at 0.3 mK need ~200 W of 1075 nm at NA 0.55 and array efficiency
0.83, whatever the distance. A higher NA is the lever, since the depth per watt scales as
NA^2, so NA 0.7 would need ~125 W. However, at NA 0.7 the metasurface edge deflection rises
to 0.7 to 0.8, where eta ~ 0.66. (ii) Damage. The same power on a smaller aperture raises the
intensity as 1 / WD^2. At 3 mm a single absorbing particle larger than 4 um radius heats by
more than 50 K on fused silica. The global temperature rise of a thin plate depends little
on the aperture, so the film choice is the same at every distance: at ~200 W only LPCVD SiN,
or any film on sapphire, stays within the limits of docs/prefab_checks_report.md (Si-rich SiN
on fused silica at the upper absorption estimate is limited to 150 to 180 W). (iii) Imaging and
addressing. A lens images a field whose radius in angle is fixed by its design. At a
short distance that angle covers a small area, so a 3 mm doublet images only a quarter of a
100,000-site array, and a 6 mm doublet three quarters at 4 um pitch. (iv) Computation. The
kernel size scales with (field radius / wavelength)^2 x N and not with the distance: about
5 x 10^5 x 12,600 complex entries (~50 GB) at 4 um pitch, or ~30 GB at 3 um pitch (DERIVED from
the v2 kernel, 103,459 x 2,556). On this 15 GB machine this needs on-the-fly kernels, at an
estimated 10 to 20 h for the WGS loop. A workstation with 64 GB or a GPU would make it routine.

What would help. (i) Keep WD >= 10 mm; the present cell geometry (15.5 mm stack) is adequate on
both sides. (ii) Provide ~200 W of 1075 nm on LPCVD SiN or a sapphire substrate, with the
cleanliness specification of the damage model (absorbing particles < 100 um on fused silica
at 15.5 mm). (iii) Use a 3 um pitch if the trap-overlap budget allows. This shrinks the field
to 669 um and reduces the corner falloff, but the neighbouring Airy rings then overlap more,
and the interstitial background (Section 1 of docs/prefab_checks_report.md) should be
rechecked. (iv) A higher NA, at ~1.6 x less power per depth, needs steep-deflection meta-atoms
or a refractive-metasurface hybrid.

Limitations. Single-site ideal kernels: the multi-site WGS efficiency and uniformity are taken
from the 20,164-site v2 design, not recomputed for 100,000 sites. The meta-atom efficiency
versus deflection mixes two Si-rich families (p 0.50 / h 1.5 um below s = 0.54, p 0.45 / h 1.3 um
above) and assumes the local-grating picture. The lens field radius uses free focus per
field point; a common plane gives slightly smaller radii (Section 3 of docs/metalens_report.md).
Vacuum geometries without a window are used throughout; the 15.5 mm row approximates the
present cell. The full 100,000-site design is reported above.

## Imaging the 100,000-site array, and the capacity of a 20 mm chip

Imaging lens (852.3 nm, NA 0.5, ray trace, one common focal plane; results/scaling_100k/
imaging_field.json):

| Field radius covered | Dedicated 852 nm doublet (16 mm, Section 3 of docs/metalens_report.md): min Strehl | Shared 990 nm lens with 852 nm pre-corrector: min Strehl |
|---|---|---|
| 400 um (v2 corner) | 0.999 | 0.985 (wave optics, results/multiwave/wavecheck.json) |
| 670 um (100,000 sites, 3 um pitch) | 0.993 | - |
| 891 um (100,000 sites, 4 um pitch) | 0.977 | 0.686 (0.913 at 894.6 nm) |
| 1,265 um | 0.889 | 0.192 |
| 1,600 um | 0.677 | - |

Hence the 100,000-site array can be imaged with the dedicated 852 nm doublet (Strehl >= 0.977 in
one plane, focus spread 0.5 um), but not well through the shared multi-wavelength lens, whose
pre-corrected field is narrower. The field of the dedicated lens with Strehl >= 0.8 has a
radius of ~1.4 mm (interpolated), i.e. an inscribed square of 2.0 mm. Correction (docs/million_sites.md,
Section 2): these values used an aplanatic construction extended only to 1.1 R1, so the outer
focusing-plate zones reached by fields beyond ~0.8 mm carried a linearly extrapolated phase.
With the construction extended to 1.6 R1 and the plate spacing equal to f, the same NA 0.5
objective reaches a Strehl ratio >= 0.8 out to 1.62 mm (0.33 million sites at 4 um pitch).

Wave-optics check at the 100k corner (superposition method, Section 4b of
docs/metalens_report.md; 16 mm doublet with construction to 1.6 R1 and plate spacing f, NA 0.5,
852.3 nm; one common plane chosen over the 100k field; results/scaling_100k/imaging_ee_100k.json,
configs/imaging_100k.yaml). Entries: ensquared energy in the pitch box / largest neighbour-box
fraction, both relative to the total PSF power (ray Strehl ratio in the common plane). (A first
version normalised to the power inside the +/-1.5 pitch evaluation window, which holds only 0.968
(4 um) and 0.956 (3 um) of it; its values were 0.03 to 0.04 too high.)

| Pitch | 0 | ~0.5 corner | ~0.7 corner | Corner |
|---|---|---|---|---|
| 4 um (corner 891 um) | 0.900 / 0.013 (0.9995) | 0.900 / 0.014 at 445 um | 0.900 / 0.014 at 630 um | 0.897 / 0.017 (0.993) |
| 3 um (corner 670 um) | 0.876 / 0.015 (0.9998) | 0.876 / 0.015 at 335 um | 0.876 / 0.016 at 474 um | 0.875 / 0.017 (0.9985) |

The on-axis values are the diffraction limit of this pupil for each box (below the 0.920 of a
uniformly filled Airy pupil in the 4 um box, because the aplanatic pupil carries more power in
its outer zones; the 3 um box is smaller, so its EE is lower). Hence the present 16 mm meta-doublet
images the whole 100,000-site array, at either pitch, in one focal plane with at most 0.003 loss
of ensquared energy and at most 0.4 percentage points more neighbour light than on axis. The
poorer 16 mm values of docs/million_sites.md (0.76 on axis) come from a plane chosen over a
2.83 mm field, not from the lens.

Camera and photons (Provisional). Ref. [Manetsch et al., arXiv:2403.12021] images 6,100 atoms
at 7.2 um spacing on a qCMOS camera with magnification 5.1, 80 ms exposure and a 7 x 7 pixel box
per site (SOURCE, read in session). For 4 um pitch a magnification of ~5.75 gives ~5 pixels per
site and a 7.3 mm square image `[VALUE?: qCMOS pixel 4.6 um, sensor 4096 x 2304 (~18.8 x 10.6 mm);
recalled]`. The collected fraction at NA 0.5 is 0.067 of 4 pi (DERIVED) against 0.082 at NA 0.55.
With a lens efficiency of 0.80 (theory; ~0.6 if the fabrication penalty of the only measured
metalens tweezer applies), the detected photon number per atom is ~0.5 to 0.8 of that of a
refractive NA 0.55 objective, `[VALUE?: objective transmission]`, so the exposure time would
grow by ~1.3 to 2 x for the same fidelity.

Hologram field at the present geometry (single-site ideal spots, 15.5 mm vacuum model,
results/scaling_100k/hologram_large_field.json, scripts/hologram_large_field.py): relative depth
0.986 at 1.6 mm, 0.971 at 2.3 mm, 0.945 at 3.2 mm, 0.893 at 4.5 mm and 0.798 at 6.4 mm, with waists
<= 0.90 um and largest meta-atom deflection 0.61 to 0.73 (one-family efficiency curve; the first
version, from an uncommitted run with the mixed curve, gave 0.983 to 0.790).

Capacity of one 18.4 mm chip (DERIVED from the above; 4 um / 3 um pitch):

| Limit | Value | Sites |
|---|---|---|
| Hologram field (depth >= 0.89 of axis, radius 4.5 mm, inscribed square) | 6.4 mm square | 2.5 million / 4.5 million |
| Metasurface space-bandwidth | 1.36 x 10^9 pixels against ~1 um^2 per resolvable spot | not limiting |
| Imaging field of the dedicated 852 nm doublet (Strehl >= 0.8, radius 1.4 mm) | 2.0 mm square | 0.25 million / 0.44 million |
| Laser power, 0.2 mK (0.3 mK), N = 148.7 mK x site/W x P / U | 200 W: 149,000 (99,000); 1 kW: 744,000 (496,000) | |
| Design computation, ~ N^2 (100,000 sites: ~5 h on 4 cores) | 500,000 sites: ~5 days on 4 cores | GPU or cluster |

The chip is not the constraint: its field and space-bandwidth support millions of sites. The
binding limits are the imaging field of a single objective (~0.25 to 0.45 million sites) and the
trap power, ~1 kW for ~0.5 to 0.74 million sites at 0.3 to 0.2 mK. At 1 kW the intensity on an 18.4 mm
chip is 375 W/cm^2; this requires LPCVD SiN (global heating 0.4 K, SIM scaling) or a sapphire
substrate, and absorbing particles below ~18 um radius on fused silica (DERIVED, 50 K). Several
lasers can be added only at the same wavelength, because the array scales laterally with lambda:
1 nm of detuning moves a site at 891 um by 0.8 um, comparable to the waist. Two
polarisation-combined amplifiers at the same wavelength, as in ref. [Manetsch], double the power
through the polarisation-independent metasurface. Atom loading, rearrangement time and vacuum
lifetime at 10^5 to 10^6 atoms are outside this optical study.
