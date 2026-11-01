# Campaign 1b: full-aperture simulation, pillar errors, vector focusing, film estimates, and v2 (20,000 sites)

Technical report, Phase A modelling for ANRF/ARG/2026/009937, following
docs/campaign1_report.md. All numbers are from simulations executed in this repository
(type SIM, config and results directory given) unless marked otherwise; every SIM value
should be reviewed before it is quoted in a submission. This report supersedes the
pupil-design results of campaign 1 that were affected by the pupil-grid error described
in its erratum.

## Summary of results

(i) Correction. A pupil-grid error in the campaign 1 design code clipped every design and
produced the reported on-axis zero order; with it fixed, a plain hard-edged top hat
performs as well as a soft edge, and the 5 mm bench reaches 4.7 % rms after two feedback
rounds (Section 1).

(ii) Films (literature estimates). Si-rich SiN with n(852 nm) = 2.30 is estimated at
n(1075 nm) = 2.277 and k ~ 6.7e-5 (band 3e-5 to 3e-4); the selected pillar becomes 1.6 um tall,
absorbs 5.9e-4 of the light, and moves the focus by 0.11 um at 100 W on 0.5 mm fused silica
(0.50 um at the upper band), inside the 1 um budget. LPCVD SiN absorbs < 1e-5 (Section 2).

(iii) Vector optics at NA 0.55 lowers the peak by 7.7 %, makes the spots 7 % elliptical and
creates a fictitious magnetic field on Cs of up to 0.58 G per mK of depth; the bias field
should be along the tweezer polarisation (Section 3).

(iv) v1 at full aperture (48 x 48, 1.36 x 10^9 pixels): 8.0 % rms after two feedback rounds;
with library pillars and 5 or 10 nm rms width noise the uniformity is unchanged and the
efficiency falls from 0.84 to 0.80 or 0.75; 2,304 sites at 50 W reach 2.9 to 3.3 mK
(Section 4).

(v) v2, 20,164 sites at 4.0 um pitch from the same 18.4 mm metasurface: 5.3 to 5.7 % rms after
two feedback rounds, with library pillars and 5 nm noise; 0.34 to 0.37 mK at 50 W and 0.68 to
0.74 mK at 100 W. The exact-kernel design is mandatory (the standard design keeps a Strehl
ratio of 0.02 at the corners). Designing the array 50 um beyond the lens focus removes the
undiffracted light from the array plane and raises the weakest site from 0.78 to 0.84 of the
mean (Section 5).

## 1. Correction to campaign 1 (pupil-grid clipping)

The pupil-domain designs of campaign 1 sampled direction cosines only up to
sigma_max = lambda / (2 dx) = 0.478 (dx = pitch / 4 = 1.125 um), below the NA. The pupil was
clipped to a square, and at rendering the metasurface pixels outside it kept the plain lens
phase and focused onto the axis. The grid now covers at least 1.08 NA
(hologram/wgs.py, PupilGrid.for_sites), and pupil_amplitude refuses a grid that does not
cover the NA. Rerun results (configs unchanged except pupil_guard 2.0):

| Case (5 mm bench, 48 x 48, NA 0.5) | Before fix | After fix |
|---|---|---|
| Zero order on axis, hard edge (2 mm test, 24 x 24 off axis) | 0.463 of a site peak | 4e-4 |
| Array efficiency, same test | 0.858 | 0.878 |
| Exact design, hard edge: non-uniformity before / after 2 feedback rounds | 11.3 % / (not run) | 6.1 % / 4.7 % |
| Exact design, soft edge (10 %): before / after feedback | 7.8 % / 4.6 % | 6.9 % / 4.0 % |
| Debye design, hard edge | 23.6 % | 24.4 % |
| Folded vs full-grid propagation (non-uniformity) | 7.79 % vs 7.77 % | 7.71 % vs 7.67 % |

In-cell pupil model (18.4 mm, NA 0.55, top hat; results/hologram_b1_cell20mm): mode
efficiency 0.933 (exact) and 0.911 (Debye), non-uniformity 0.05 % and 5.0 %. The ideal spot
of the unclipped pupil has the peak intensity per watt of a uniformly filled NA 0.574
(previously 0.538 from the clipped pupil). The "edge-generated zero order" and the need for
a soft illumination edge are withdrawn: a hard-edged top hat performs as well as a soft
edge (4.7 % against 4.0 % rms after feedback, efficiency 0.879 against 0.889).

## 2. Film n and k: literature-based estimates

The films cannot be measured in this project phase, so n(lambda) and k are estimated from
published data (materials.py; results/library_analysis_est_1075).

Dispersion of Si-rich SiN. A single-oscillator (Wemple-DiDomenico) model,
n^2 - 1 = E0 Ed / (E0^2 - E^2), fitted to five measured SiN films (Luke et al., Opt. Lett.
40, 4823 (2015); Vogt, PhD thesis (2015), three PECVD films; Beliaev et al., Thin Solid Films
763, 139568 (2022); data via refractiveindex.info) gives an oscillator energy E0 that
decreases with index, from 8.0 to 9.3 eV near n(852) = 2.0 to 6.8 eV at 2.10 (one film,
Vogt #1 with n(852) = 1.89, lies below the trend). A linear fit of E0 against n(852),
extrapolated to the
n(852 nm) = 2.30 film quoted from ref. 14, E0 = 5.83 eV, Ed = 23.5 eV, and
n = 2.294 (895 nm), 2.278 (1064.4 nm), 2.277 (1075 nm), i.e. 0.023 below the constant value
used in campaign 1. The extrapolation is Provisional: the five films span n(852) = 1.89 to
2.10, and the recipe-dependent index of a CNNP Si-rich film is uncertain by about
+/-0.05, the band already simulated in the campaign 1 tolerances.

Absorption. Photothermal measurements at 1064 nm give an absorption coefficient of
(1.5 to 2.9) x 10^-2 cm^-1 for stoichiometric Si3N4 (k ~ 1.8 x 10^-7) and
7.94 +/- 0.50 cm^-1 for Si-rich low-stress SiNx (k ~ 6.7 x 10^-5) (Kumar et al., ACS Photonics
12, 6381 (2025), DOI 10.1021/acsphotonics.5c01886), consistent with 6.9 +/- 0.7 cm^-1 for
similar membranes (Steinlechner et al., Phys. Rev. D 96, 022007 (2017)). Those films are
LPCVD; PECVD and ICP-PECVD films contain more hydrogen, so we carry a band of
k = 3 x 10^-5 to 3 x 10^-4 for our Si-rich film (Provisional).

Consequences (RCWA with complex index, configs/library_est_1075.yaml):

| Quantity | Si-rich SiN (estimate) | LPCVD SiN |
|---|---|---|
| Selected family (p, h) | 0.50 um, 1.6 um | 0.60 um, 1.6 um |
| eta_16, +/-10 nm averaged | 0.959 | 0.960 |
| Pillar absorption, phase-averaged (max) | 5.9e-4 (1.1e-3) | < 1e-5 (numerical floor) |

The lower estimated index makes the selected Si-rich pillar 0.1 um taller than in
campaign 1. Thermal consequences at 100 W (configs/thermal_film_estimates.yaml):

| Absorption | 0.5 mm fused silica | 1 mm fused silica | 0.5 mm sapphire |
|---|---|---|---|
| LPCVD, 3.5e-6 | -0.001 um | -0.001 um | 0 |
| Si-rich, nominal 5.9e-4 | -0.11 um (dT = 6 K) | -0.14 um | -0.015 um |
| Si-rich, upper band 2.7e-3 | -0.50 um (dT = 27 K) | -0.65 um | -0.07 um |

All are inside the 1 um focal-shift budget. With these estimates, explaining the 19 %
deficit of ref. 4 thermally (a_film = 1.16 %, campaign 1) would need a film absorbing
~20 times more than the nominal Si-rich estimate, which strengthens the conclusion that
substrate heating is unlikely to be the cause.

## 3. Vector focusing at NA 0.55

Richards-Wolf (Debye) field of the in-cell soft-edge top-hat pupil with x polarisation,
two polarisation models for the metasurface output (meridional rotation as for an
aplanatic lens; projection of x onto the plane transverse to each ray), which agree to
~1 % (results/vector_check; scripts/vector_check.py):

| Quantity | Value |
|---|---|
| Peak intensity, vector / scalar | 0.923 |
| 1/e^2 waist along / across the polarisation (scalar 0.886 um) | 0.944 / 0.880 um |
| Peak \|E_z\|^2 relative to the central intensity | 0.04 |
| Circularity Im(E* x E)/\|E\|^2 within the half-maximum core | up to 0.48, along y |
| Cs vector-to-scalar ratio r_v at 1075 nm | -0.080 |
| Fictitious field in the half-maximum core | up to 0.58 G per mK of depth, along y |
| Its gradient at the centre | 1.07 G/um per mK |

Hence the scalar model overestimates the depth by 7.7 % (this factor now enters the power
budget), the spots are 7 % elliptical, and the longitudinal field in quadrature produces
local circular polarisation off-centre. For Cs F = 4 this acts as a fictitious magnetic
field along y with a 1 G/um-per-mK gradient. With the bias field along the polarisation axis
x, the fictitious field is perpendicular to it and shifts the levels only in second order,
by ~B_fict^2 / (2 B_0): for an atom at 10 % of the depth in a 1 mK trap (position spread
~0.07 um) and B_0 = 3 G, about 0.07^2 / 6 G ~ 1 mG (DERIVED). A bias along y would instead
give first-order shifts of ~0.07 G. The bias direction for the static array should therefore
be along the tweezer polarisation.

## 4. v1 at full aperture: 48 x 48 sites, 18.4 mm, with pillar errors

Folded propagation of the whole 18.4 mm metasurface (36,864^2 = 1.36 x 10^9 pixels at the
0.5 um pillar pitch) into a 2.048 mm focal window, after a D4-symmetric exact-kernel design
(results/fullaperture_v1; configs/fullaperture_v1.yaml). Library pillars use the Si-rich
estimate family (0.50 um, 1.6 um) in the local-periodic approximation per pixel; width
errors are a uniform bias or independent Gaussian noise per pillar.

| Case | Depth non-uniformity (rms) | Min / mean | Array efficiency | Metasurface transmission |
|---|---|---|---|---|
| Ideal phase, no feedback | 12.9 % | 0.70 | 0.862 | 1 |
| Ideal phase, 1 feedback round | 11.8 % | 0.62 | 0.868 | 1 |
| Ideal phase, 2 feedback rounds | 8.0 % | 0.83 | 0.868 | 1 |
| Library pillars (LPA) | 8.1 % | 0.82 | 0.836 | 0.971 |
| + 5 nm rms width noise | 8.1 % | 0.82 | 0.804 | 0.956 |
| + 10 nm rms width noise | 8.2 % | 0.80 | 0.750 | 0.956 |
| + 10 nm width bias | 8.1 % | 0.82 | 0.832 | 0.967 |

Array efficiency is the power within pitch/2 of each site divided by the power incident on
the metasurface. Three results follow. First, at full aperture the exact design meets
the 10 % target after two rounds of feedback (8.0 %), but starts higher than on the 5 mm
bench (12.9 % against 6.1 %) and converges more slowly. The residual is site-to-site
scatter of 7 to 8 % at every radius with no radial trend, and the largest outliers are
whole symmetry orbits (e.g. the four sites at (+/-92, +/-92) um at 1.42 of the mean), so
further rounds, or camera feedback in the experiment, continue to reduce it. Second,
pillar errors cost efficiency, not uniformity: 10 nm rms random width noise removes 10 % of
the light (into a diffuse background), and a uniform 10 nm bias only 0.5 %. Third, the
on-axis zero order is 3.9 to 5 % of a site peak, i.e. ~2e-5 of the power, harmless.

Power budget. Utilisation = core efficiency / 0.82 (core fraction of the ideal spot) x
supercell edge factor 0.879 x AR-coated substrate 0.995 `[VALUE?]` x AR-coated cell 0.99
`[CONFIRM?]` x top-hat shaper 0.95 `[VALUE?]` x vector peak factor 0.923, with NA_eff = 0.574:

| Case | Utilisation | Mean depth, 2,304 sites, 50 W | 100 W |
|---|---|---|---|
| Ideal phase | 0.673 | 3.31 mK | 6.62 mK |
| Library pillars | 0.648 | 3.19 mK | 6.37 mK |
| + 5 nm noise | 0.623 | 3.06 mK | 6.13 mK |
| + 10 nm noise | 0.581 | 2.86 mK | 5.71 mK |

## 5. v2: 20,000 sites

Target: 142 x 142 = 20,164 sites at 4.0 um pitch, from the same 18.4 mm, NA 0.55 metasurface
(results/fullaperture_v2; configs/fullaperture_v2.yaml). Our choice of 4.0 um is a
compromise between two considerations. First, a smaller pitch keeps the array corners
closer to the axis (399 um against 450 um at 4.5 um pitch), which limits the field the
hologram must correct and reduces the pupil sampling, and hence the design cost.
However, a smaller pitch increases the overlap between neighbouring traps (Airy first
zero 1.2 um, i.e. 0.3 of the pitch) and tightens the addressing tolerance of the 895 nm
beams; ref. 4 used 4.30 um.

Field of view. The Debye design is unusable at this size: its Strehl ratio is 0.82 at
200 um, 0.44 at 282 um (the array edge) and 0.02 at the corners (399 um), with up to 9 um of
field curvature (hologram/field_model.py). The exact-kernel design is therefore required,
not merely preferable. Without symmetry the exact problem needs a ~1.7 x 10^6 x 20,164
kernel (~270 GB); the D4 reduction (hologram/symmetric.py) needs 103,459 x 2,556 complex64
entries (2.1 GB, built in 7 min on 4 cores). The reduction is exact (tests/test_hologram.py).

Full-aperture results (1.36 x 10^9 pixels, folded into a 2.048 mm window):

| Case | Depth non-uniformity (rms) | Min / mean | Array efficiency | On-axis intensity / site peak |
|---|---|---|---|---|
| Pupil model (exact kernels) | 0.24 % | - | mode efficiency 0.869 | - |
| Ideal phase, no feedback | 7.0 % | 0.75 | 0.845 | 0.21 |
| Ideal phase, 2 feedback rounds | 5.3 % | 0.78 | 0.854 | 0.38 |
| Library pillars (LPA) | 5.5 % | 0.79 | 0.823 | 0.40 |
| + 5 nm rms width noise | 5.7 % | 0.77 | 0.793 | 0.40 |

Power budget (same factors as Section 4, NA_eff = 0.574):

| Case | Utilisation | Mean depth, 20,164 sites, 50 W | 100 W |
|---|---|---|---|
| Ideal phase | 0.655 | 0.37 mK | 0.74 mK |
| Library pillars | 0.631 | 0.35 mK | 0.71 mK |
| + 5 nm noise | 0.607 | 0.34 mK | 0.68 mK |

Hence 20,000 Cs traps at 0.34 to 0.37 mK need 50 W at 1075 nm, and 100 W gives 0.7 mK, inside
the 0.3 to 1 mK range the proposal considers for loading and imaging. The array is uniform
to 5 to 6 % rms after two feedback rounds, better than v1 at the same stage (8.0 %); we have
not established why.

Zero order. The on-axis spot, which lies between the four central sites, grows to 0.4 of a
site peak in v2 (0.05 in v1), since its power is fixed while each site's power scales as
1/N. It is a weak extra trap at ~40 % of the site depth, 2.8 um from four sites, and must be
removed. We attribute it to the residual difference between the sampled pupil model and
the continuous metasurface (Provisional): a 24 x 24 array on the 5 mm bench designed with
coarse pupil sampling (guard 1.25) shows 0.14 of a site peak on the axis, whereas a 24 x 24
array on a 2 mm aperture designed with fine sampling (guard 3, array off axis) shows 4e-4.
The remedy we adopt
is to design the array in a plane displaced from the lens focus: the exact kernels focus
the sites at the shifted plane, while undiffracted light still focuses at the original
plane, ~45 Rayleigh ranges away for a 50 um shift. On the 5 mm bench (24 x 24 sites) this
reduces the on-axis intensity from 0.10 to 0.010 of a site peak and, because the central
sites no longer interfere with it, the non-uniformity after one feedback round from 9.8 % to
4.4 %.

v2 with the array 50 um beyond the lens focus (results/fullaperture_v2_dz50;
configs/fullaperture_v2_dz50.yaml):

| Case | Depth non-uniformity (rms) | Min / mean | Array efficiency | On-axis intensity / site peak |
|---|---|---|---|---|
| Ideal phase, no feedback | 6.3 % | 0.68 | 0.846 | 0.026 |
| Ideal phase, 1 round | 6.0 % | 0.80 | 0.853 | 0.004 |
| Ideal phase, 2 rounds | 5.3 % | 0.84 | 0.858 | 0.11 |
| Library pillars (LPA) | 5.4 % | 0.84 | 0.826 | 0.11 |
| + 5 nm rms width noise | 5.7 % | 0.81 | 0.796 | 0.11 |

| Case | Utilisation | Mean depth, 20,164 sites, 50 W | 100 W |
|---|---|---|---|
| Ideal phase | 0.658 | 0.37 mK | 0.74 mK |
| Library pillars | 0.633 | 0.36 mK | 0.71 mK |
| + 5 nm noise | 0.609 | 0.34 mK | 0.68 mK |

The array defocus raises the weakest site from 0.78 to 0.84 of the mean at equal rms
non-uniformity and costs no efficiency. The on-axis intensity drops from 0.21 to 0.026 before
feedback and from 0.38 to 0.11 after two rounds. A zero order focused 50 um before the array
plane would be weaker there by ~(z_R / dz)^2 ~ 5e-4 (DERIVED, z_R ~ 1.1 um), so we expect the
residual 0.11 to be interstitial background rather than a focused spot. We have not checked
this: the interstitial intensity statistics were not saved in these runs. This is the first
check for the next campaign, together with an explicit dark-point constraint at the axis if
the residual proves to be a coherent spot.

## 6. Limitations

(i) Scalar propagation in all full-aperture runs; the vector check (Section 3) supplies a
7.7 % depth correction and the spot ellipticity for an on-axis spot only. Off-axis vector
effects at 400 um field were not computed.

(ii) Pillars enter through the local-periodic approximation per pixel. The coupling loss at
large deflection (supercell factor 0.879) is applied as an aperture-averaged multiplier,
not per pixel, so its effect on uniformity and stray light is not simulated. Width errors
are uncorrelated per pillar; correlated errors (proximity effect, stitching fields of
electron-beam writing) are not modelled and would act differently: a slowly varying width
error behaves like a phase aberration and would degrade spots locally rather than add a
diffuse background.

(iii) Feedback uses the simulated focal intensity as the measurement, with two rounds only.
In both v1 and v2 the non-uniformity was still falling after the second round.

(iv) The film n(lambda) and k are literature estimates (Section 2); the Si-rich index
extrapolation rests on five films that do not reach n(852) = 2.3.

(v) The residual on-axis intensity of 0.11 of a site peak in v2 is not yet identified as
background or as a coherent spot (Section 5).

(vi) The folded propagation aliases stray light from more than 1.0 mm off axis into the
window; on the 5 mm bench this changed the per-site depths by up to 2.5 %.

## Provenance ledger and open placeholders

All values below are SIM (review before submission) unless marked.

| Value | Section | Pointer |
|---|---|---|
| Corrected bench, folded-validation and in-cell numbers | 1 | results/fullpixel_bench5mm*, folded_validation_5mm, hologram_b1_cell20mm, fullpixel_zeroorder_2mm_* |
| Si-rich SiN n(lambda) (E0 = 5.83 eV, Ed = 23.5 eV) | 2 | DERIVED from SOURCE data (Luke 2015; Vogt 2015; Beliaev 2022); materials.py |
| k = 1.8e-7 (LPCVD), 6.7e-5 (Si-rich), band 3e-5 to 3e-4 | 2 | SOURCE: Kumar et al. 2025 (DOI 10.1021/acsphotonics.5c01886); Steinlechner et al. 2017; band is an assumption |
| Selected families, absorption 5.9e-4 / <1e-5 | 2 | results/library_analysis_est_1075 |
| Thermal focal shifts with estimated absorption | 2 | results/thermal_film_estimates |
| Vector factors (0.923, waists, circularity, fictitious field) | 3 | results/vector_check |
| v1 full-aperture metrics and power budget | 4 | results/fullaperture_v1, fullaperture_v1_dz50 |
| v2 full-aperture metrics and power budget | 5 | results/fullaperture_v2, fullaperture_v2_dz50 |
| Debye-design aberration at v2 radii | 5 | hologram/field_model.py (command in session log; reproduce with spot_aberrations) |

Open placeholders: cell geometry `[CONFIRM?]`; AR coatings (0.995, 0.99) and shaper
throughput (0.95) `[VALUE?]`; v2 pitch 4.0 um `[CONFIRM?: compatible with the 895 nm
addressing and Rydberg zone layout]`; CNNP etch limits `[CONFIRM?]`; Si-rich film index at
852 nm (2.30, from ref. 14 as quoted in the proposal) `[VALUE?]`.
