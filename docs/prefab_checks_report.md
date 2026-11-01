# Pre-fabrication checks for the 18.4 mm, 20,000-site metasurface, and damage model

Technical report. Design under test: v2, 142 x 142 = 20,164 Cs sites at 4.0 um pitch, NA 0.55,
18.4 mm aperture, 1075 nm, array designed 50 um beyond the lens focus, two feedback rounds
(results/fullaperture_v2_dz50). All numbers are SIM (review before submission) unless marked.

## Summary

Result. Before the 18.4 mm, 20,164-site design (v2, array 50 um beyond the lens focus) is
committed to fabrication, the four checks and the damage model give the following.
(i) Light on the axis: with the defocused array the axial spot (0.11 of the mean site peak)
is no brighter than the interstitial background maxima (~0.1 to 0.15 everywhere). The
spurious-trap level of the array is therefore set by the background, not by a zero order.
(ii) Cell wall: thickness errors of +/-50 um add < 1 % rms. The relative tilt between the
metasurface and the wall must be <= 0.1 deg (7.4 % rms, against 5.3 % at zero tilt);
0.25 deg gives 14.8 %.
(iii) Fabrication: a realistic combination of stitching (20 nm), per-field bias (2 nm),
smooth (3 nm, 5 um) and random (3 nm) width errors gives 6.2 % rms non-uniformity, minimum
0.78, array efficiency 0.80. The most damaging error is long-range width variation (5 nm rms
over 50 um: 7.5 % rms).
(iv) Vector focusing: uniform from the axis to the array corner to 1e-4 in peak and 0.1 % in
waist (peak 0.917 of scalar, fictitious field 0.63 G per mK).
(v) Damage: thermo-elastic stress is never binding (<= 1.5 MPa at 100 W). The binding
limits are the thermal focal shift and the temperature, for Si-rich SiN on fused silica at the
upper absorption estimate (150 to 180 W). Sapphire or LPCVD SiN raises every limit above
1 kW. Absorbing particles larger than ~100 um radius must be excluded on fused silica above
150 W.

Main limitation. All non-uniformities are open-loop (no camera feedback after
fabrication), each error case is a single random realisation, and the meta-atoms are treated
in the local periodic approximation.

## 1. Light on the axis: coherent spot or background?

Intensity at every interstitial point (centre of four sites, 19,881 points) relative to the
mean site peak, compared with the on-axis value (results/checks_v2, axis_light):

| Design | Interstitial median | 99th percentile | Maximum (location) | On-axis |
|---|---|---|---|---|
| Array in the lens focus | 0.013 | 0.087 | 0.38 (axis) | 0.38 |
| Array 50 um beyond the focus | 0.013 | 0.086 | 0.12 (122 um off axis) | 0.11 |

Without the defocus the axis spot is a coherent zero order: it is the brightest interstitial
point, 4.4 times the 99th percentile. With the defocus it falls to 0.11, inside the tail of
the background (the brightest interstitial point, 0.12, lies elsewhere). Its radial profile
still falls to 0.01 within 0.75 um, i.e. it remains a diffraction-limited spot at 28 % of its
in-focus strength (results/checks_v2/axis_patch_*.npy). We interpret this as the part of the
zero order that originates in the hologram itself (the residual between the sampled design
and the continuous metasurface), which focuses in the array plane like any hologram
component, whereas the undiffracted lens light, which the defocus removes, focuses 50 um
away (Provisional).

Consequence. The spurious-trap level of the array is set by the interstitial maxima, about
0.1 of the site depth everywhere, not by the axis alone. At 50 W (site depth 0.36 mK) this is
~40 uK, comparable to the temperature of laser-cooled Cs, so loading into interstitial
maxima is possible in principle. Two remedies are available: feedback that includes the
interstitial points as dark targets (camera-based weighted GS with noise suppression), and a
larger pitch (the median background, 0.013, is dominated by the tails of the four
neighbouring Airy spots at 2.8 um).

## 2. Cell-wall tolerances

The design assumes a 4 mm fused-silica wall parallel to the metasurface. Errors in the wall
change only the propagation, so one folded render of the metasurface is re-propagated with the
perturbed transfer function: thickness errors through the layered stack, tilts with the exact
plane-wave transfer function of a tilted slab (propagation/folded.py,
transfer_function_tilted; zero tilt reproduces the layered result to 5e-8 rad). Each case is
evaluated at its best focus (scan) with the array re-centred and re-scaled by a least-squares
fit to the found peaks, and without repeating the feedback (results/checks_v2, cell:*).

| Wall error | Depth non-uniformity (rms) | Min / mean | Array efficiency | Array displacement |
|---|---|---|---|---|
| none | 5.3 % | 0.84 | 0.858 | - |
| thickness +50 um | 6.0 % | 0.81 | 0.855 | focus -33 um |
| thickness -50 um | 6.2 % | 0.80 | 0.855 | focus +33 um |
| thickness +100 um | 7.6 % | 0.75 | 0.850 | focus -66 um |
| thickness +200 um | 12.9 % | 0.58 | 0.824 | focus -132 um |
| tilt 0.1 deg | 7.4 % | 0.76 | 0.855 | 2.5 um lateral |
| tilt 0.25 deg | 14.8 % | 0.59 | 0.840 | 6.2 um lateral |
| tilt 0.5 deg | 24.7 % | 0.36 | 0.790 | 12.1 um lateral |
| tilt 1 deg | 31.7 % | 0.20 | 0.734 | array fails (peaks displaced ~1 um rms after fitting) |

The array follows the focus: a wall thicker by dt moves it by -dt cos(t)/sqrt(n^2 - sin^2 t)
~ -0.66 dt (rms ray), which the atoms follow, and a tilted wall shifts it laterally by
d sin(b)[1 - cos(b)/sqrt(n^2 - sin^2 b)]. The cost in uniformity comes from the
field-dependent part of the path through the wall, which the exact design pre-compensates
for the nominal wall only. Thickness errors of +/-50 um add < 1 % rms and +/-100 um ~2.3 %, so a
standard cell-wall tolerance is acceptable (`[VALUE?: thickness tolerance of the purchased
cell]`). Tilt is the critical tolerance: the relative tilt between the metasurface and the wall
must be kept to <= 0.1 deg (1.7 mrad) to stay within ~2 % of the nominal non-uniformity. This is
an alignment specification rather than a fabrication one: the metasurface mount needs tip/tilt
adjustment and an autocollimation (back-reflection) reference to the wall, which resolves
far below 0.1 deg. Alternatively the measured wall thickness and tilt can be included in the
design kernels, and camera feedback after installation recovers part of the uniformity (not
simulated here). The interstitial maxima rise with tilt (0.32 at 0.25 deg), so tilt also raises
the spurious-trap level.

## 3. Correlated fabrication errors

Method. Each case renders the full v2 metasurface with library pillars (LPA transmission of
the selected Si-rich family at 1075 nm, transmission and phase per pillar), with width errors
added before the lookup (hologram/render.py, pillar_errors). The result is propagated through
the cell by the folded transfer function and evaluated in the array plane without repeating
the feedback (scripts/fullaperture_checks.py; results/checks_v2, fab:*). The error models
are as follows:
(i) write-field stitching: each 100 um write field displaced rigidly by a Gaussian offset
of the stated rms (`[VALUE?: write-field size of the e-beam tool used]`);
(ii) per-field width bias, a dose or focus offset per write field;
(iii) smooth width noise, a Gaussian random field with the stated rms and correlation length;
(iv) a combination of all three with uncorrelated 3 nm noise.
Every case is one random realisation (one seed), so differences below ~0.3 % rms are not
significant `[VALUE?: seed-to-seed scatter not computed]`.

| Case | Depth non-uniformity (rms) | Min / mean | Array efficiency | Interstitial max | On-axis |
|---|---|---|---|---|---|
| ideal phases (reference) | 5.3 % | 0.84 | 0.858 | 0.12 | 0.11 |
| library pillars, no errors | 5.4 % | 0.84 | 0.826 | - | 0.11 |
| library, uncorrelated 5 nm (campaign 1b) | 5.7 % | 0.81 | 0.796 | - | 0.11 |
| stitching 20 nm rms | 5.6 % | 0.81 | 0.825 | 0.13 | 0.12 |
| stitching 50 nm rms | 6.7 % | 0.78 | 0.822 | 0.14 | 0.11 |
| per-field width bias 3 nm rms | 6.2 % | 0.78 | 0.812 | 0.14 | 0.11 |
| smooth 5 nm rms, correlation 5 um | 5.8 % | 0.81 | 0.796 | 0.13 | 0.11 |
| smooth 5 nm rms, correlation 50 um | 7.5 % | 0.75 | 0.806 | 0.15 | 0.12 |
| combined: stitch 20 nm, bias 2 nm, smooth 3 nm / 5 um, random 3 nm | 6.2 % | 0.78 | 0.799 | 0.14 | 0.11 |

The position error (0.083 to 0.086 um rms) and the mean waist (0.869 um) are unchanged by
every error model, so fabrication errors act on the depths only.

Interpretation. The correlation length decides where the error light goes. With a 5 um
correlation length the phase errors scatter into angles ~ lambda / L ~ 0.2. This carries the
light out of the array field: the efficiency drops (0.796, as for uncorrelated noise) and the
uniformity hardly changes (5.8 %). With a 50 um correlation length the scattering angles are
~ 0.02. This redistributes light within the array, over a few hundred um at the 15 mm working
distance: the efficiency barely changes (0.806), but the depths become non-uniform (7.5 %,
minimum 0.75) (Provisional; scaling argument consistent with both cases). Stitching and
per-field bias are periodic at the 100 um write field. They behave like weak gratings of
period 100 um, whose orders (sin theta ~ 0.011, i.e. ~160 um in the array plane) add faint
displaced copies of the array that beat with it. Hence 20 nm stitching costs 0.2 % rms,
50 nm costs 1.3 %, and a 3 nm per-field bias costs 0.8 % and 1.4 % in efficiency.

Specifications that follow (to stay within ~1 % rms of the library design, before camera
feedback): stitching <= 20 to 30 nm rms; per-field width (dose) uniformity <= 2 to 3 nm rms;
long-range width variation (tens of um and beyond) <= 3 to 5 nm rms. Long-range variation
also includes film thickness and etch-depth gradients across the 18.4 mm part, which change
the pillar phase like a width error `[VALUE?: SiN thickness uniformity across a 25 mm
substrate]`. The realistic combination gives 6.2 % rms and an efficiency of 0.80, against
5.4 % and 0.83 for perfect library pillars. These are within the 10 % rms target before any
feedback. Camera feedback (weighted GS on measured depths) acts on exactly these slowly
varying errors, since they are static and deterministic once the part is made; it is not
simulated here.

## 4. Vector focusing off axis

Richards-Wolf field for sites from the axis to the array corner, with every ray traced
exactly to the off-axis site through the cell stack (results/vector_offaxis;
scripts/vector_offaxis.py; top-hat pupil, projection polarisation model):

| Site (um) | Peak, vector / scalar | Waist along / across polarisation (um) | Peak \|E_z\|^2 | Max circularity in core | Fictitious field (G per mK) |
|---|---|---|---|---|---|
| (0, 0) | 0.9172 | 0.909 / 0.844 | 0.044 | 0.527 | 0.629 |
| (282, 0) | 0.9172 | 0.910 / 0.844 | 0.044 | 0.526 | 0.628 |
| (0, 282) | 0.9172 | 0.909 / 0.844 | 0.044 | 0.526 | 0.628 |
| (282, 282), corner | 0.9172 | 0.910 / 0.844 | 0.044 | 0.526 | 0.628 |

The vector corrections are uniform across the field to 1e-4 in peak and 0.1 % in waist: the
ray cones of the corner sites are tilted by less than 3 deg, too little to change the
polarisation structure. The vector-induced shift of the intensity maximum is below the
0.04 um grid step at every site. Hence the on-axis vector factor (0.917 for the top hat,
0.923 for the soft edge) can be applied to the whole array, and vector effects do not add
non-uniformity.

## 5. Damage model

Mechanism. At CW intensities of tens of W/cm^2 damage of dielectric optics is thermal
(vendor CW ratings of dielectric mirrors near 1.07 um are ~2 x 10^4 W/cm^2 for mm beams,
[Thorlabs](https://www.thorlabs.com/catalogpages/obsolete/2019/SSM1020A-K13.pdf)); the
intensity on the metasurface is P / 2.67 cm^2, i.e. 37 W/cm^2 at 100 W and 75 W/cm^2 at
200 W. The model (thermal/damage.py; configs/damage_v2.yaml; results/damage_v2) evaluates:
substrate thermo-elastic stress (free-edged disc, membrane plus through-thickness parts)
against strength / 7; peak temperature rise against 50 K (adhesive and film stability,
assumption); focal shift against 1 um; film thermal-mismatch stress; local heating at fully
absorbing defects; and a two-photon-absorption bound. Material data: fused silica E 72.7 GPa,
nu 0.16, abraded modulus of rupture 52.4 MPa (Corning 7980 sheet); sapphire E 335 GPa,
nu 0.25 (Alkor sheet), strength `[VALUE?: 400 MPa]`; Si-rich PECVD biaxial modulus 125 GPa
and CTE 2.6 ppm/K; LPCVD E 260 GPa; strength scatter of glass: Weibull modulus 4 to 5
([SPIE 7302](https://spiedigitallibrary.org/conference-proceedings-of-spie/7302/730210/Flexural-strength-of-fused-silica-Weibull-statistical-analysis/10.1117/12.817714.full)).

| Case (at 100 W) | Peak dT | Max tensile stress | Film mismatch stress | Focal shift | Power limit (binding criterion) |
|---|---|---|---|---|---|
| FS 0.5 mm, Si-rich nominal (a = 5.9e-4) | 6.0 K | 0.11 MPa | -1.6 MPa | -0.11 um | 832 W (temperature) |
| FS 0.5 mm, Si-rich upper (a = 2.7e-3) | 27.5 K | 0.48 MPa | -7.2 MPa | -0.50 um | 182 W (temperature) |
| FS 1.0 mm, Si-rich upper | 18.0 K | 0.31 MPa | -4.7 MPa | -0.65 um | 153 W (focal shift) |
| FS 0.5 mm, LPCVD (a = 3.5e-6) | 0.04 K | 0.001 MPa | -0.04 MPa | -0.001 um | > 10^5 W |
| Sapphire 0.5 mm, Si-rich nominal | 0.44 K | 0.34 MPa | +0.13 MPa | -0.015 um | 6.8 kW (focal shift) |
| Sapphire 0.5 mm, Si-rich upper | 2.0 K | 1.5 MPa | +0.6 MPa | -0.07 um | 1.5 kW (focal shift) |

Bulk absorption of the substrates (10 ppm/cm for fused silica, 30 ppm/cm for sapphire,
`[VALUE?]`) changes none of these numbers measurably, since the film absorbs 10^2 to 10^3 times
more. Stress is never binding: the largest value, 1.5 MPa, is below the conservative 3.4 MPa
guideline for fused-silica vacuum windows
([SPIE 2966](https://astronomicaltelescopes.spiedigitallibrary.org/conference-proceedings-of-spie/2966/0000/Laser-induced-damage-and-fracture-in-fused-silica-vacuum-windows/10.1117/12.274227.full))
and far below 7.5 MPa (fused silica) and 57 MPa (sapphire) at a safety factor of 7.

Local hot spots. A fully absorbing defect or contamination particle of radius a heats by
dT = I a / kappa. For dT <= 50 K:

| Power | Intensity | Largest tolerable absorbing defect, fused silica | Sapphire |
|---|---|---|---|
| 50 W | 19 W/cm^2 | 370 um radius | 8.8 mm |
| 100 W | 37 W/cm^2 | 184 um | 4.4 mm |
| 150 W | 56 W/cm^2 | 123 um | 2.9 mm |
| 200 W | 75 W/cm^2 | 92 um | 2.2 mm |

Nonlinear absorption. With the bulk-silicon two-photon coefficient (1.5 cm/GW at 1.06 um, an
upper bound for the wider-gap nitrides) and a local field enhancement of 20 inside the pillars,
two-photon absorption is 10^-7 of the linear absorption of Si-rich SiN at 300 W. Field
enhancement inside the pillars therefore matters only through linear absorption, which the
RCWA absorption values already include.

Design guidance from the damage model. (i) Global damage does not limit the design up to
several hundred watts for any film/substrate combination considered; the binding criteria are
the thermal focal shift and the temperature of the mount, and only for Si-rich films on fused
silica at the upper absorption estimate (150 to 180 W). (ii) Sapphire raises every limit above
1 kW and makes the part insensitive to contamination; with fused silica, absorbing particles
or residues (e.g. from the Cr hard mask) larger than ~100 um radius must be excluded above
150 W, which sets a cleanliness specification. (iii) LPCVD stoichiometric SiN removes thermal
limits altogether. (iv) Limitations: steady state only (no turn-on transients), rim held at
ambient (a poorly heat-sunk mount, which changes both the temperature and its profile, is
not modelled), no damage of the cell window coatings modelled.

## 6. Limitations

(i) One random realisation per fabrication case; seed-to-seed scatter not quantified. (ii)
LPA meta-atom response, with the aperture-averaged supercell edge factor accounted for in
campaign 1b but not in these renders. (iii) No camera feedback after fabrication or
installation. (iv) Thermal model steady-state only, rim at ambient. (v) Film absorption from
literature-based estimates (campaign 1b), not measured. Next logical step: a seed ensemble
for the two most damaging cases (smooth 50 um, stitching 50 um), then the 18.4 mm design
release with the tolerance list above.
