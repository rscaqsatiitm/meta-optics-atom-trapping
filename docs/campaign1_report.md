# Campaign 1: static meta-hologram (B1), meta-atom library and high-power thermal model

Technical report, Phase A modelling for ANRF/ARG/2026/009937. All numbers below are
from simulations executed in this repository (type SIM, config and results directory
given) unless marked otherwise; every SIM value should be reviewed before it is quoted
in a submission.

## Summary of results

<!-- SUMMARY -->

## 1. Scope and method

We model three things: (i) the meta-atom library for square SiN pillars on fused silica
at 1075, 1064.4, 895 and 852 nm (RCWA, fmmax, Jones-direct-Fourier formulation); (ii)
the direct-projection hologram that generates the static Cs array at 1075 nm through the
glass-cell wall (pupil-domain weighted Gerchberg-Saxton design, validated by exact
scalar angular-spectrum propagation of the full 10^8-pixel metasurface at 5 mm
aperture); (iii) steady-state heating of the metasurface substrate at 1 to 100 W and its
effect on the tweezers (axisymmetric finite volume plus a pupil wavefront analysis).
Equations, assumptions and the constants ledger are in `docs/model_notes.md`.

Assumed geometry (all `[CONFIRM?]`): 2.5 mm air, 4 mm fused-silica cell wall, 9 mm
vacuum to the atoms (modelled on ref. 4); metasurface on a 0.5 mm fused-silica wafer
diced to a 25 mm die, held at its rim; 48 x 48 square array at 4.5 um pitch; NA 0.55,
which needs an 18.4 mm aperture in this geometry (19.9 mm for NA 0.58, consistent with
the 19.8 mm of ref. 4).

## 2. Atomic model and power budget

The two-line polarisability model with Steck's matrix elements and measured static
polarisability gives alpha(1064 nm) = 1171 a.u. and alpha(1075 nm) = 1125 a.u., against
the proposal's 1,160 and 1,120 a.u. from ref. 36 (0.9 % and 0.4 % higher). Applied to
ref. 4 (Rb, 850 nm, 33 W, 18,225 sites, 37 % utilisation) with NA_eff = 0.467 derived
from the measured 1.11 um first dark radius, it gives 0.306 mK against the measured
0.305(17) mK. With the proposal's parameters (37 % utilisation, NA_eff 0.47) 50 W at
1075 nm supports 2,809 sites at 1 mK and 9,362 at 0.30 mK. The proposal's "~9,200 at
0.3 mK" corresponds to 0.305 mK; at 0.300 mK the figure is ~9,360. A transport tweezer
at 1.2 mK and 1064.4 nm needs 6.58 mW at NA 0.5 and 7.45 mW at NA 0.47 (proposal: 6.6
and 7.5 mW). (`tests/test_atoms.py`)

## 3. Meta-atom library

Numerics (results/rcwa_convergence, SIM). With the Jones-direct-Fourier formulation the
relative transmission phase converges to 0.003 rad at 197 Fourier terms. The
permittivity grid adds 0.01 to 0.03 rad between 256 and 1024 points, largest for
gaps below 60 nm, which the fabrication limits exclude. The plain-FFT formulation
converges about ten times more slowly (0.07 rad error at 317 terms), and grcwa, which
uses it, agrees with fmmax-FFT to 1e-4 rad at 437 terms. This is an independent check of
the solver wrapper and its phase convention. The library (data/library_main.h5,
configs/library_main.yaml) covers 5 pitches (0.40 to 0.60 um), 13 heights (0.8 to
2.0 um), 40 widths and 4 wavelengths for LPCVD SiN (Luke et al. dispersion) and Si-rich
SiN (n = 2.3, `[VALUE?]`), at 197 terms on a 256 grid.

Selection at 1075 nm (results/library_analysis_1075). Geometries are ranked by the
first-order efficiency of a 16-level blazed grating in the local-periodic approximation,
including the Fresnel loss at the exit face, averaged over a uniform width bias of -10,
0 and +10 nm. Fabrication limits are a minimum feature and gap of 100 nm and h/w <= 10
(`[CONFIRM?: CNNP etch capability]`). About 65 geometries per film are feasible, and the
best dozen lie within 0.5 % of each other, limited by the aspect ratio. Among those within
0.3 % of the best we take the shortest pillar:

| Film | Pitch (um) | Height (um) | Widths (nm) | eta_16 | eta_16, +/-10 nm averaged | Phase coverage (rad) | Min. T0 |
|---|---|---|---|---|---|---|---|
| LPCVD SiN | 0.60 | 1.6 | 168 to 492 | 0.967 | 0.960 | 6.13 | 0.93 |
| Si-rich SiN | 0.50 | 1.5 | 150 to 390 | 0.960 | 0.957 | 6.60 | 0.90 |

The Si-rich height agrees with the proposal's estimate (1.15 um at 852 nm from ref. 14,
~26 % taller at 1075 nm, i.e. ~1.45 um). An earlier ranking with cubic interpolation
returned eta > 1 for resonant geometries (minimum T0 of 0.02 to 0.3); this was spline
overshoot across under-sampled resonances, and those geometries also lose most under a
width bias. Linear interpolation and the bias-averaged metric remove them.

Tolerances of the selected families (8 phase levels, direct RCWA; the common phase
offset is removed, since a uniform offset does not affect the hologram):

| Perturbation | LPCVD: differential phase rms (rad) | LPCVD: efficiency ratio | Si-rich: differential phase rms (rad) | Si-rich: efficiency ratio |
|---|---|---|---|---|
| width +10 nm | 0.050 | 0.998 | 0.079 | 0.998 |
| width -10 nm | 0.071 | 0.997 | 0.064 | 0.991 |
| index +0.05 | 0.126 | 0.979 | 0.117 | 0.989 |
| index -0.05 | 0.146 | 0.984 | 0.109 | 0.986 |
| height +50 nm | 0.045 | 0.998 | 0.065 | 1.000 |

A ±0.05 index error costs 1 to 2 % of efficiency, several times more than a 10 nm width
bias or a 50 nm height error. Hence the film index at 1075 nm should be measured by
ellipsometry on each deposition run and fed into the lookup (materials.from_csv), rather
than taken from the literature.

<!-- SUPERCELL -->

## 4. Hologram design: field aberration of the standard (Debye) design

### 4.1 The Debye design is not adequate beyond ~100 um from the axis

The usual design (weighted Gerchberg-Saxton with an FFT between pupil and focal plane,
as for SLM arrays) assumes that the metasurface phase which focuses to site x_j is the
on-axis lens phase plus a tilt in direction-cosine space. For direct projection
through a 15 mm path this neglects a path difference of order |x_j|^2 / (2 WD) ~ 0.8 um
at 150 um from the axis (DERIVED). Ray-traced through the actual layer stack
(hologram/field_model.py), the residual after removing tilt and defocus gives, for
NA 0.55 and a top-hat beam (results/hologram_b1_cell20mm, aberration_*.png):

| Distance from axis (um) | 50 | 100 | 150 | 200 | 300 |
|---|---|---|---|---|---|
| Strehl ratio, Debye design | 0.999 | 0.988 | 0.939 | 0.819 | 0.346 |
| Axial focus shift (um) | -0.14 | -0.57 | -1.28 | -2.27 | -5.10 |

A 48 x 48 array at 4.5 um pitch has its corners at 150 um, so the Debye design costs ~6 %
of depth at the corners and has 1.3 um of field curvature, above the 1 um focal
tolerance of the error budget. A 100 x 100 array (corners at 318 um) would have
corner Strehl ratios near 0.3. We do not know how refs. 4 and 14 handled this (the
retrieved text of ref. 4 mentions block-wise FFT propagation), but ref. 4's measured
first dark radius of 1.11 um against the 0.9 um design is consistent in sign with
field-dependent aberration (Speculative).

### 4.2 Exact-kernel design

Replacing the tilt kernel by the exact path difference through the stack,
K_ij = exp(-i k0 [OPL(r_i -> x_j) - OPL0(r_i)]), removes the field dependence at no
extra cost for 2,304 sites (K is a 160,000 x 2,304 complex64 matrix, 16 s for 60 WGS
iterations on 4 CPU cores). For the in-cell geometry, evaluated with the exact kernels:

| Design | Mode efficiency, sum over sites | Depth non-uniformity (rms) | Corner / centre |
|---|---|---|---|
| Debye | 0.908 | 4.0 % | 0.82 |
| Exact | 0.928 | 0.41 % | 0.99 |

### 4.3 Full-pixel validation and two effects the pupil model misses

We rendered the designs onto a 5 mm, NA 0.5 metasurface (10,240^2 pixels at 0.5 um,
behind a 1 mm fused-silica plate) and propagated the full aperture with the exact scalar
angular-spectrum method (scripts/validate_fullpixel.py). The plain lens reproduces the
paraxial Airy peak intensity to 0.6 %. This geometry has a short working distance
(4.7 mm), so its field aberration is stronger than in the cell (Strehl 0.62 at 150 um
for the Debye design, against 0.94 in the cell). It is therefore a strict test of the
exact kernels.

The validation exposed two effects, both now handled in the code:

1. Rendering. Bilinear interpolation of the hologram phase from the design pupil grid
   acts as a focal-plane envelope sinc^2(x ds / lambda), which weakened corner sites
   by ~2x. The back-propagated pupil field is now evaluated slightly beyond the NA and
   Fourier-upsampled 8x before rendering (exact.render_grid). Direct per-pixel
   evaluation of the exact site superposition gives the same result to 0.1 %.
2. An edge-generated zero order. A hard-edged, phase-only aperture produces a focus
   on the optical axis carrying ~1e-3 of the power, even though the design places only
   2e-6 there. Resolved by radius, the discrepancy comes entirely from the outermost
   ~5 % of the aperture radius, where the zero-order integrand oscillates at full
   amplitude up to a discontinuity. Neither a finer pupil grid nor edge-coverage
   weighting removes it. A raised-cosine edge on the illumination over 10 % of the
   radius suppresses it ~30x (to 1.5 % of a site peak for a 24 x 24 array). Since the
   zero-order power is fixed while the site power scales as 1/N, a hard-edged aperture
   gives a zero-order spot brighter than a site for N ~ 1,000 or more; a soft-edged
   top-hat shaper is therefore a requirement, not an option.

Final full-pixel results, 48 x 48 sites, soft edge (results/fullpixel_bench5mm_softedge;
configs/fullpixel_bench5mm_softedge.yaml):

| Design | Depth non-uniformity (rms) | Min / mean depth | Mean Strehl | Array efficiency |
|---|---|---|---|---|
| Debye | 23.6 % | 0.28 | 0.91 | 0.877 |
| Exact | 7.8 % | 0.65 | 0.98 | 0.880 |
| Exact + 2 rounds of model feedback | 4.6 % | 0.83 | 0.98 | 0.887 |

Array efficiency is the power within pitch/2 of each site divided by the incident
power, for an ideal phase-only metasurface (no meta-atom losses). The feedback rounds
rescale the WGS targets by (mean / measured peak)^0.5 using the full-pixel result as
the measurement; experimentally this is camera or atom feedback. The residual 4.6 %
is within the 10 % rms target of Table 1, before meta-atom and fabrication errors.

<!-- INCELL -->

## 5. High-power thermal model (hypothesis H2a)

The film absorption a_film (fraction of the incident power absorbed in the SiN layer)
is unknown for our films, so it is the free parameter; all results are linear in
a_film x P. Results for the in-cell geometry at NA 0.55 (results/thermal_h2a,
config configs/thermal_h2a.yaml):

| Substrate | Beam | P (W) | a_film | dT_max (K) | Focal shift (um) | Strehl, refocused | Strehl, fixed plane |
|---|---|---|---|---|---|---|---|
| fused silica 0.5 mm | top hat | 50 | 1e-3 | 5.1 | -0.09 | 1.0000 | 0.9993 |
| fused silica 0.5 mm | top hat | 100 | 1e-3 | 10.2 | -0.18 | 1.0000 | 0.9974 |
| fused silica 1.0 mm | top hat | 100 | 1e-3 | 6.7 | -0.24 | 1.0000 | 0.9955 |
| sapphire 0.5 mm | top hat | 100 | 1e-3 | 0.75 | -0.025 | 1.0000 | 1.0000 |
| fused silica 0.5 mm | Gaussian, w = 2R | 100 | 1e-3 | 11.2 | -0.22 | 1.0000 | 0.9963 |

Three results follow. First, with the die heat-sunk at its rim the substrate thermal
lens is almost a pure defocus over the NA 0.55 pupil: the Strehl ratio after refocusing
stays above 0.9999 in every case run, so substrate heating does not degrade the spots,
it moves them. Second, the error-budget limit of 1 um focal shift at 100 W is reached
at a_film = 5.4e-3 (0.5 mm fused silica, top hat), 4.1e-3 (1 mm fused silica) and
4.0e-2 (0.5 mm sapphire). The sapphire margin is ~7x rather than the 24x ratio of
thermal conductivities for two reasons: face cooling already limits the fused-silica
temperature, so the temperature ratio is 13.6x (10.2 K against 0.75 K), and the
thermo-optic path coefficient chi = dn/dT + (n - 1)(1 + nu) alpha is ~1.8x larger for
sapphire (with its dn/dT still `[VALUE?]`).
Third, the proposal's footnote 1 scaling, dT ~ P_abs / (kappa w), holds for a
half-space; a 0.5 to 1 mm wafer under a 20 mm beam is in the opposite (thin-plate)
limit, where dT ~ P_abs / (4 pi kappa t) [1 + 2 ln(R/a)], nearly independent of beam
size at fixed absorbed power. Both limits are reproduced by the solver
(`tests/test_thermal.py`). The qualitative argument of Section 3 of the proposal (heating
grows with beam size at fixed intensity) survives, and is stronger in the plate limit,
since absorbed power grows as w^2 at fixed intensity.

Implication for the ref. 4 shortfall (12.2 W measured against 15.0 W predicted). If
the 19 % deficit were a Strehl loss of refocused spots, no film absorption up to 30 %
reproduces it in this model. As a fixed-plane loss it requires a_film = 1.65 % (dT =
32 K, focal shift -1.15 um), under our assumed 1 mm fused-silica substrate (ref. 4 does
not state its substrate). For a 1.15 um Si-rich film at 850 nm this corresponds to an
extinction coefficient of order 2e-3 [DERIVED: k ~ a lambda / (4 pi h f), with fill
fraction f ~ 0.5], well above typical values for SiN `[VALUE?: measured k of CNNP
films]`. With the same a_film, the ref. 14 case (2.5 W on 5 mm) gives a fixed-plane
Strehl of 0.997, consistent with its unchanged efficiency. Hence substrate heating
explains the ref. 4 deficit only if the film absorbs at the per-cent level; at the
0.01 to 0.1 % level it is negligible. We therefore regard the thermal explanation in the
proposal as Speculative, and the model points to two cheaper discriminating
measurements for H2a: the film absorption itself (photothermal or calorimetric, at
1075 nm), and the focal-plane position against power, which this model predicts to
move by -0.18 um per 100 W per 0.1 % absorption on 0.5 mm fused silica.
Thermo-optic changes inside the pillars are negligible by comparison: for LPCVD SiN
(dn/dT ~ 2.5e-5 /K `[VALUE?]`) at dT = 10 K, dn = 2.5e-4 changes the pillar phase by
~2 mrad (DERIVED, dphi/dn ~ k0 h ~ 9 rad); even the Si-rich value of 1.65e-4 /K
(Nejadriahi et al., Opt. Express 2020, DOI 10.1364/OE.396969) gives ~15 mrad.

## 6. Limitations

<!-- LIMITATIONS -->

## 7. Next steps

<!-- NEXT -->

## Provenance ledger and open placeholders

<!-- LEDGER -->
