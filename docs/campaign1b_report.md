# Campaign 1b: full-aperture simulation, pillar errors, vector focusing, film estimates, and v2 (20,000 sites)

Technical report, Phase A modelling for ANRF/ARG/2026/009937, following
docs/campaign1_report.md. All numbers are from simulations executed in this repository
(type SIM, config and results directory given) unless marked otherwise; every SIM value
should be reviewed before it is quoted in a submission. This report supersedes the
pupil-design results of campaign 1 that were affected by the pupil-grid error described
in its erratum.

## Summary of results

<!-- SUMMARY -->

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

<!-- V2 -->

## 6. Limitations

<!-- LIMITATIONS -->

## Provenance ledger and open placeholders

<!-- LEDGER -->
