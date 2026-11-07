# Side-B metasurface objective: design study (ray optics with an axisymmetric wave check)

Technical report. Target: a flat (metasurface) objective outside the cell, NA 0.5 at the atoms,
imaging or addressing the v2 array (142 x 142 sites at 4.0 um, corner at 399 um from the axis)
through the cell stack 2.5 mm air | 4 mm fused silica | 9 mm vacuum `[CONFIRM?: side-B air gap,
wall, distance to the atoms]`. All numbers are SIM (review before submission) unless marked;
code: src/metatrap/metalens/{raytrace,aplanat,doublet,wave_axial}.py,
scripts/metalens_doublet.py, configs/metalens_doublet.yaml, results/metalens_doublet/.

## Summary

Result. An air-spaced aplanatic metasurface doublet works over the whole v2 field: corrector
metasurface (aperture stop, radius R1 = 8.106 mm), 16.0 mm of air, focusing metasurface, then
the cell. Its ray-traced Strehl ratio is >= 0.9992 at every field point out to 400 um in one
flat focal plane, at each of 852.3, 894.6 and 1064.4 nm, each with its own phase maps. The
field headroom is large: Strehl 0.96 (852 nm) and 0.975 (1064 nm) at 1.0 mm. An independent
axisymmetric wave-optics propagation (quasi-discrete Hankel transform through both plates and
the cell) gives on-axis Strehl 0.9989 / 0.9988 / 0.9980 at the three wavelengths, with the
focus on the nominal plane to < 1 nm and a focal-spot 1/e^2 radius equal to the Airy value
(0.71 um at 852 nm, 0.88 um at 1064 nm).

Main limitation. A metasurface phase map works at one wavelength. With the 852 nm maps,
the focus moves by -18.8 um per nm of detuning. If the input beam is re-focused, the Strehl
stays >= 0.8 only to about +/-1.8 nm. At 894.6 nm and 1064.4 nm the Strehl is 0.02 and 0.004
even with re-focusing. Co-focal operation at two or three of these wavelengths therefore needs
independent phase maps per wavelength (multi-wavelength meta-atoms covering the full 2 pi x 2 pi
phase torus) or separate optical paths.

Mechanical tolerances (852 nm, Strehl >= 0.92 over the field after refocusing):
(i) decentre between the plates <= 5 um (2 um gives 0.987);
(ii) plate spacing within +/-50 um (0.992);
(iii) cell-wall thickness within +/-50 um (0.988);
(iv) relative plate tilt <= 0.1 deg, after its coma is nulled by decentring the focusing plate;
uncorrected, a 0.1 deg tilt is equivalent to a 28 um decentre and fails;
(v) tilt of the cell relative to the objective <= 0.1 deg (0.93 at the corner).

Efficiency (Section 6): blazed-supercell RCWA of the single-pillar libraries gives an estimated
doublet focusing efficiency of 0.80 to 0.83 at 852 and 1064 nm (LPCVD or Si-rich SiN). The loss
sits in the outer zones of the focusing plate, where the deflection reaches 0.5 (ring
efficiency 0.73 to 0.81).

## 1. Why this architecture

| Architecture (852 nm, NA 0.5) | Strehl at image height 0 / 75 / 150 / 282 / 400 um (best focus per field) |
|---|---|
| Hyperbolic singlet, stop at the lens | 1.000 / 0.015 / 0.001 / 0.002 / 0.007 |
| Aplanatic doublet on a 2 mm silica spacer | 1.000 / 0.997 / 0.949 / 0.511 / 0.055 |
| Air-spaced aplanatic doublet, 16 mm | 1.000 / 1.000 / 1.000 / 0.9995 / 0.9992 |

The singlet fails because a single flat phase profile can make the axial focus stigmatic but
cannot also satisfy the Abbe sine condition: off-axis coma appears at once (Strehl 0.015 at
75 um). We also tried polynomial phase profiles optimised by ray tracing (singlet, singlet with
a displaced stop, doublets on 1 and 3 mm substrates; results/metalens_study). None reached a
useful field: the merit landscape at NA 0.5 and f = 16 mm, with up to ~2,400 waves of phase,
defeated the optimiser. We therefore abandoned that route for an analytic construction.

Construction (src/metatrap/metalens/aplanat.py). Take a ray entering the corrector at height
h1. We require that it leave the focusing plate with transverse direction cosine
s_f = h1 / f (sine condition, f = R1 / NA) and that it reach the axial focus through the cell
stack. This fixes its height at the focusing plate, h2 = sum_i d_i s_f / sqrt(n_i^2 - s_f^2).
The ray between the plates then has the transverse direction cosine p that satisfies
h2 - h1 = sum_j t_j p / sqrt(n_j^2 - p^2) over the spacer layers. The generalised Snell law
gives the two phase gradients, dphi1/dr (h1) = k0 p and dphi2/dr (h2) = -k0 (s_f + p).
Here, k0 = 2 pi / lambda, and d_i, n_i (t_j, n_j) are the thicknesses and indices of the cell
(spacer) layers. The result is stigmatic on axis with zero linear coma by construction. Only
astigmatism and field curvature remain, quadratic in the field.

The 2 mm silica-spaced version of this construction is good to only ~150 um. Spacing the
plates by about one focal length in air makes the focusing plate sit at the back focal plane
of the stop. The system becomes telecentric (chief-ray direction cosine at the atoms 3.6e-4
at 400 um), and every field point then uses nearly the same zone structure of the focusing
plate, shifted laterally. This is the flat-optics analogue of placing the stop at the front
focal plane, as in wide-field metalenses with a front aperture (`[REF?: Kalvach and Szabo,
J. Opt. Soc. Am. A 33, A66 (2016); unverified]`). Metasurface doublets corrected for
monochromatic aberrations by ray-traced optimisation were demonstrated by Arbabi et al.
(800 um aperture, f/0.9, > 60 deg field at 850 nm, ~70 % efficiency;
[Nat. Commun. 7, 13682 (2016)](https://doi.org/10.1038/ncomms13682), arXiv:1604.06160).
They did not treat a 16 mm aperture, a focus behind a 4 mm window, or a flat-field
requirement across a fixed atom plane, which is what this design addresses.

Spacing scan (one common focal plane over 0 to 400 um; minimum Strehl):

| Plate spacing | 852.3 nm | 894.6 nm | 1064.4 nm |
|---|---|---|---|
| 12 mm | 0.778 | 0.797 | 0.852 |
| 14 mm | 0.933 | 0.940 | 0.957 |
| 15 mm | 0.979 | 0.981 | 0.987 |
| 16 mm | 0.9992 | 0.9993 | 0.9995 |
| 17 mm | 0.991 | 0.992 | 0.995 |
| 18 mm | 0.957 | 0.961 | 0.972 |

The optimum, 16 mm, lies just below f = 16.21 mm, as expected from telecentricity.

## 2. Design data (16 mm spacing)

| Quantity | 852.3 nm | 894.6 nm | 1064.4 nm |
|---|---|---|---|
| Stop radius R1 (common) | 8.106 mm | 8.106 mm | 8.106 mm |
| Marginal ray height at the focusing plate R2 | 8.106 mm | 8.107 mm | 8.109 mm |
| Focal length f = R1 / NA | 16.21 mm | 16.21 mm | 16.21 mm |
| Corrector: max deflection / total phase | 0.025 / 155 waves | 0.025 / - | 0.025 / 125 waves |
| Focusing plate: deflection at R2 / at R2 + 0.4 mm | 0.500 / 0.525 | 0.500 / 0.525 | 0.500 / 0.525 |
| Focusing plate: local zone period at R2 + 0.4 mm | 1.62 um | - | 2.03 um |
| Focusing plate: total phase | 2,376 waves | - | 1,903 waves |

The corrector is weak (zone period >= 34 um at 852 nm). All the difficulty sits in the outer
zones of the focusing plate. To serve the 400 um field, the focusing plate needs a radius of
R2 + 0.4 mm ~ 8.5 mm, and its deflection there reaches 0.525. Section 6 gives the meta-atom
efficiency at this deflection. The plate substrates are neglected in this study. A plane
substrate before the corrector changes nothing for collimated input. A substrate between the
plates is handled by the layered-spacer option of the construction, but has not been
evaluated here.

## 3. Field performance

Free focus per field (852.3 nm; 1064.4 nm in parentheses), with the field curvature dz:

| Image height | 0 | 400 um | 600 um | 800 um | 1.0 mm | 1.3 mm | 1.6 mm |
|---|---|---|---|---|---|---|---|
| Strehl | 1.000 | 0.9992 (0.9995) | 0.996 (0.997) | 0.986 (0.991) | 0.961 (0.975) | 0.874 (0.917) | 0.670 (0.766) |
| Focus dz (um) | 0 | -0.11 | -0.24 | -0.42 | -0.66 | -1.11 | -1.94 |

In one flat plane (dz = -0.10 um) the Strehl ratio is >= 0.9992 everywhere from 0 to 400 um
at all three wavelengths. The residual is astigmatism and field curvature, quadratic in
field. The field curvature over the v2 array, 0.11 um, is 0.03 of the depth of focus
(lambda / NA^2 = 3.4 um at 852 nm). The field over which the Strehl ratio stays >= 0.8 has a
radius of ~1.4 mm, room for arrays about 3.5 times wider than v2.

## 4. Wave-optics check (on axis)

Method (src/metatrap/metalens/wave_axial.py): the radially symmetric scalar field is
propagated with the quasi-discrete Hankel transform (window 9 mm, spectral cut-off 0.70 / lambda,
N = 11,838 to 14,784 samples). The steps are: a top-hat beam with a 5 um tanh edge at the
stop, times exp(i phi1); 16 mm of air; times exp(i phi2); the cell layers, each by its exact
transfer function exp(i k0 d sqrt(n^2 - (lambda nu)^2)), without Fresnel interface losses. The
Strehl ratio is the on-axis intensity relative to that of the same angular spectrum with all
phases aligned, maximised over the focal plane. The transform was checked on a Gaussian beam
(spectrum error 4e-16, inverse transform error 3e-15).

| Wavelength | Strehl (wave) | Best focus relative to nominal plane | 1/e^2 radius | Spectral power beyond 1.1 NA |
|---|---|---|---|---|
| 852.3 nm | 0.9989 | +3e-5 um | 0.71 um | 1e-6 |
| 894.6 nm | 0.9988 | +4e-6 um | 0.74 um | 1e-6 |
| 1064.4 nm | 0.9980 | +2e-4 um | 0.88 um | 4e-6 |

The waists equal the Airy 1/e^2 radius 0.41 lambda / NA (0.70, 0.73 and 0.87 um) to within
the 0.01 um grid. The 0.1 to 0.2 % shortfall from the ray value (1.000) is consistent with
edge diffraction of the stop over the 16 mm spacer. This test is independent of the ray
tracer: it uses only the phase profiles and exact propagation. It confirms the phase
construction, the focal position through the 4 mm wall, and the stationary-phase (ray)
assumption at this Fresnel number (R1^2 / (lambda L) ~ 4,800). It is on axis only; Section 4b
covers the field.

## 4b. Off-axis check by superposition from the focusing plate

Method (src/metatrap/metalens/wave_offaxis.py; doublet.json, wave_offaxis). Rays from a
uniformly filled stop (300 x 300 grid, 70,688 samples) are traced through the corrector and
the spacer, where geometrical optics holds (Fresnel number ~4,800). Each sample then
radiates from the focusing plate with amplitude sqrt(J) dA and its accumulated phase. Here, J
is the local area magnification from the stop to the plate. The focal field is the sum of
these contributions through the cell, each with the exact point-to-point optical path of the
layer stack, the same kernel as the exact hologram model, which was validated against
full-pixel angular-spectrum propagation in campaign 1. The image is located as an intensity
maximum by a lateral and axial search, not by fitting the eikonal. Hence this tests the
ray-traced Strehl ratios and image positions of Section 3 independently in the second half of
the system. In the "apodised" cases the meta-atom power efficiency versus local deflection
(Section 6) weights both plates.

| Wavelength | Image height | Strehl, superposition / ray trace | Focus, superposition / ray (um) | 1/e^2 half-width, radial / tangential |
|---|---|---|---|---|
| 852.3 nm | 0 | 1.000 / 1.000 | 0.000 / 0.000 | 0.70 / 0.70 um |
| | 400 um | 0.9992 / 0.9992 | -0.105 / -0.105 | 0.70 / 0.70 um |
| | 1.0 mm | 0.960 / 0.961 | -0.657 / -0.657 | 0.71 / 0.70 um |
| | 1.3 mm | 0.872 / 0.874 | -1.07 / -1.11 | 0.73 / 0.70 um |
| 1064.4 nm | 0 | 1.000 / 1.000 | 0.000 / 0.000 | 0.86 / 0.86 um |
| | 400 um | 0.9995 / 0.9995 | -0.105 / -0.105 | 0.86 / 0.86 um |
| | 1.0 mm | 0.975 / 0.975 | -0.657 / -0.657 | 0.87 / 0.86 um |
| | 1.3 mm | 0.915 / 0.917 | -1.07 / -1.11 | 0.89 / 0.88 um |

The two methods agree to <= 0.003 in Strehl ratio and to 1e-4 um in image position at every
field point. At 1.3 mm the intensity maximum lies 0.04 um from the least-squares plane of the
ray fit, as expected once astigmatism is present: the diffraction maximum falls between the
line foci. The lateral sampling of the profiles is 0.02 um.

Apodisation by the meta-atom efficiency (Section 6 curves):

| Wavelength, film | Transmitted fraction, axis / 400 um | Peak relative to lossless uniform pupil, axis / 400 um | Half-width, axis |
|---|---|---|---|
| 852.3 nm, LPCVD | 0.798 / 0.797 | 0.792 / 0.791 | 0.70 um (unchanged) |
| 852.3 nm, Si-rich | 0.820 / 0.820 | 0.815 / 0.814 | 0.70 um (unchanged) |
| 1064.4 nm, LPCVD | 0.807 / 0.807 | 0.804 / 0.803 | 0.88 um (+0.02) |
| 1064.4 nm, Si-rich | 0.828 / 0.828 | 0.825 / 0.824 | 0.88 um (+0.02) |

Hence the weaker outer zones cost the peak intensity mainly through the transmitted power: the
ratio of peak to transmitted fraction is 0.993 to 0.996, so the apodisation itself costs at
most 0.7 % of the peak. Within the 0.02 um sampling of the profiles, they broaden the spot by at
most one step at 1064 nm and not at all at 852 nm. A site at the array corner receives the same
peak as one on the axis. The transmitted fractions equal the area averages of Section 6 (0.798
at 852 nm, LPCVD). (An earlier version of this table listed fractions 0.005 to 0.009 lower and
attributed the difference to the ray-path weighting; it was a J-weighted mean, corrected in
wave_offaxis.transmitted_fraction and recomputed by scripts/metalens_transmission.py.)

## 5. Chromatic behaviour and the multi-wavelength requirement

A metasurface keeps its phase in radians as the wavelength changes (fixed 2 pi zones). Its
deflection therefore scales as lambda / lambda_d, as for a diffractive lens. We evaluated the
852.3 nm design at other wavelengths in two ways: (a) collimated input and the nominal atom
plane; (b) the input vergence re-adjusted (point source at distance L) to bring the best
axial focus back to the atom plane.

| Detuning | Focus shift, collimated input | Strehl (a), nominal plane | Strehl (b), refocused input: on axis / 400 um | Residual rms (b) |
|---|---|---|---|---|
| 0.01 nm | -0.19 um | 0.997 | 1.000 / 0.998 | 0.0004 waves |
| 0.1 nm | -1.9 um | 0.74 | 0.999 / 0.998 | 0.004 waves |
| 0.5 nm | -9.4 um | 0.05 | 0.983 / 0.982 | 0.021 waves |
| 1 nm | -18.8 um | < 0.01 | 0.933 / 0.934 | 0.042 waves |
| 2 nm | -37.6 um | < 0.01 | 0.755 / 0.755 | 0.084 waves |
| 5 nm | -94 um | < 0.01 | 0.19 / 0.20 | 0.23 waves |
| 42.3 nm (894.6 nm) | -762 um | < 0.01 | 0.02 / 0.006 | 2.0 waves |
| 212.1 nm (1064.4 nm) | -3,250 um | < 0.01 | 0.004 / 0.001 | 16.7 waves |

Consequences. (i) For the 852 nm imaging light alone the dispersion is harmless: the D2
fluorescence is 5 MHz wide, and laser wavelengths are stable to far below 0.01 nm. (ii)
Re-focusing removes defocus but not the residual spherical aberration, which grows linearly
with the detuning. A single phase map therefore serves only a band of about +/-1.8 nm
(Strehl >= 0.8; computed for positive detunings, interpolated as exp(-(2 pi W_rms)^2)). (iii) 852/895/1064 nm co-focal operation through one objective requires
independent phase maps per wavelength. With independent maps (common stop radius), each
wavelength reproduces the single-wavelength results above (Sections 3 and 4), since the
construction is exact per wavelength.

What the meta-atoms must supply. The two focusing-plate maps differ by 473 waves (852/1064)
and 112 waves (852/895) at the edge. The required phase pairs (phi_852 mod 2 pi,
phi_1064 mod 2 pi), weighted by area, fill the torus nearly uniformly: every bin of a 36 x 36
grid holds 0.84 to 1.18 times the mean (0.83 to 1.12 for 852/895). The library must therefore
provide every phase pair with high transmission; a partial region of the torus does not
suffice. Two-wavelength meta-molecules that achieve independent full coverage were
demonstrated at 1550 and 915 nm (NA 0.46, 65 % measured efficiency; Arbabi et al.,
[Optica 3, 628 (2016)](https://authors.library.caltech.edu/records/32je3-2pw19),
arXiv:1601.05847). For 852 and 1064 nm (ratio 1.25) this is a scaled version of that problem
(Provisional). For 852 and 895 nm (ratio 1.05) it is much harder. A propagation-phase
pillar of height h gives phases in nearly the ratio 1.05 at the two wavelengths, so the phase
difference available across a library is only ~ 0.05 x (2 pi x 1.3) ~ 0.4 rad, against the
2 pi required (Speculative; estimate, not simulated). Broadband dispersion engineering is not an
option at this size: the delay-bandwidth bound limits achromatic metalenses to apertures and
NAs orders of magnitude smaller than these
([Presutti and Monticone, arXiv:2001.10899](https://arxiv.org/abs/2001.10899); venue not
verified in session). Discrete-wavelength operation (wrapped phases) is not bound by that
limit.

Recommendation. Use the objective at one design wavelength, or at 852 + 1064 nm with a
two-wavelength meta-atom library at roughly half the single-wavelength efficiency (Section 5b). Bring 894.6 nm light through another
path (for example the hologram side, or a dichroic-combined beam focused by its own optic).
`[CONFIRM?: which side-B functions need 852, 895 and 1064 nm at a diffraction-limited focus]`.

## 5b. Two-wavelength (852 / 1064 nm) meta-atoms: RCWA screen

Method (src/metatrap/rcwa/molecules.py; scripts/run_dual_library.py,
analyse_dual_library.py; configs/dual_library.yaml; results/dual_library). Two C4v-symmetric
(polarisation-independent) families on a 0.55 um square lattice: square rings (outer width,
square hole) and two-pillar cells (one square at the cell centre and one at the corners). The
layer heights are 1.6, 2.0 and 2.4 um, for Si-rich and LPCVD SiN on fused silica, solved at
852 and 1064.4 nm by RCWA. The numerics are reduced for a screen: 121 Fourier terms and a
128 grid. A spot check against 197 terms and a 256 grid on 12 atoms gives |dphi| <= 0.08 rad,
except one near-resonant ring at 0.17 rad. Geometries are on a 0.02 um grid with features
>= 0.10 um. We report two aspect-ratio limits (height / smallest solid feature): 10, the
working assumption for the etch, and 15. The figure of merit is the focusing efficiency at
each wavelength when every required phase pair must be realised. The atom for a target pair
maximises Re[t_852 e^{-i theta_852}] + Re[t_1064 e^{-i theta_1064}], and the efficiency is
|<t e^{-i theta}>|^2 (times n_out / n_in) over the targets. We use two target sets: (a) a
uniform torus, and (b) the actual phase pairs of each doublet plate, with the global phase
offset of each wavelength optimised.

| Film, height, aspect limit | Atoms | Uniform torus: eta 852 / 1064 | Same atoms, one wavelength only | Doublet (corrector x focusing plate, LPA): 852 / 1064 |
|---|---|---|---|---|
| Si-rich, 1.6 um, AR <= 10 | 58 | 0.67 / 0.56 | 0.94 / 0.95 | 0.45 / 0.33 |
| Si-rich, 1.6 um, AR <= 15 | 137 | 0.73 / 0.71 | 0.97 / 0.97 | 0.55 / 0.52 |
| Si-rich, 2.4 um, AR <= 15 | 63 | 0.75 / 0.71 | 0.96 / 0.95 | 0.57 / 0.51 |
| LPCVD, 1.6 um, AR <= 10 | 58 | 0.53 / 0.55 | 0.93 / 0.96 | 0.29 / 0.31 |
| LPCVD, 2.4 um, AR <= 15 | 63 | 0.70 / 0.63 | 0.94 / 0.95 | 0.49 / 0.41 |

(Union of both families; results/dual_library/dual_efficiency.json and plate_efficiency.json
list every family, height and limit.)

Observations. (i) Within AR <= 10, heights of 2.0 and 2.4 um admit no rings and no two-pillar
cells (walls and pillars would have to be >= 0.2 to 0.24 um), so only solid pillars remain,
and these give 0.45 to 0.55. Solid pillars alone are not hopeless, because their phase wraps
by different amounts at the two wavelengths (11.6 against 8.4 rad at 1.6 um, Si-rich), but
they fill the torus only along a few lines. (ii) The second geometric parameter is what fills
the torus. It needs thin walls or close pillar pairs, hence the strong dependence on the
aspect-ratio limit (Si-rich 1.6 um: 0.56 to 0.71 at 1064 nm). (iii) The corrector does not
benefit from its weaker phase profile: its required pairs, sampled along the radius, still
fill the torus nearly as uniformly as those of the focusing plate. (iv) Including the
deflection losses of the focusing plate (Section 6, factor ~0.87 relative to the LPA), we
estimate a two-wavelength doublet at ~0.45 to 0.50 (852 nm) and ~0.45 (1064 nm) with AR 15,
and ~0.3 to 0.4 with AR 10, against 0.80 to 0.83 for single-wavelength plates (Provisional:
screen numerics, LPA, coarse geometry grid).

Implication. A two-wavelength SiN objective is possible but costs roughly half of the light
at each wavelength, or more if the aspect ratio stays at 10. For the 852 nm imaging path this
halves the collection efficiency. For the 1064 nm tweezer path it doubles the required
power. The alternatives, in order of our preference, are as follows. (i) Separate the
functions: a single-wavelength 852 nm objective on side B (0.80 to 0.83), with 1064 nm
delivered by another path. (ii) Other films (Section 5c: TiO2, GaN and crystalline Si
were simulated; none closes the gap). (iii) A denser screen (finer geometry grid, more families
such as rings with a central pillar), which may add several per cent but is unlikely to close
the gap.

## 5c. Other films (TiO2, GaN, crystalline Si) and a finer geometry grid

Films (literature data, retrieved in session from the refractiveindex.info database):
(i) ALD TiO2, n = 2.332 / 2.308 at 852 / 1064 nm, k ~ 0
([Siefke et al., Adv. Opt. Mater. 4, 1780 (2016)](https://doi.org/10.1002/adom.201600250));
(ii) wurtzite GaN, ordinary index 2.348 / 2.333
([Barker and Ilegems, Phys. Rev. B 7, 743 (1973)](https://doi.org/10.1103/PhysRevB.7.743));
birefringence neglected);
(iii) crystalline Si, n = 3.634 / 3.555, k = 3.4e-3 / 8.2e-5, i.e. alpha = 500 / 10 cm^-1
([Schinke et al., AIP Adv. 5, 067168 (2015)](https://doi.org/10.1063/1.4923379)).
TiO2 is on fused silica at 0.55 um pitch. GaN and Si are on sapphire at 0.45 um pitch,
below lambda / n_sapphire = 0.484 um at 852 nm. Numerics as in Section 5b, except crystalline
Si, for which 121 Fourier terms were not converged (errors up to 0.25 rad and 50 % in
amplitude). Si was therefore solved with 197 terms and a 256 grid; against 289 terms the
errors are <= 0.036 rad for two-pillar cells and <= 0.125 rad / 10 % for resonant rings. The
TiO2 and GaN screens agree with 197-term solutions to <= 0.081 rad and 1.6 %
(results/dual_library_hi/convergence_tio2_gan.json).

The geometry grid matters as much as the film. Refining it from 0.02 to 0.01 um raises
every library by 0.05 to 0.10 (more atoms fill the torus), so the table compares films on the
fine grid. Doublet efficiency = corrector x focusing plate, for the actual phase pairs, in the
LPA, before the deflection factor (~0.87) of Section 6:

| Film, height | Aspect <= 10: 852 / 1064 nm | Aspect <= 15: 852 / 1064 nm | One wavelength only | Absorbed at 852 nm (median atom) |
|---|---|---|---|---|
| Si-rich SiN, 1.6 um | 0.53 / 0.40 | 0.63 / 0.58 | 0.97 to 0.98 | ~0 |
| ALD TiO2, 1.6 um | 0.49 / 0.44 | 0.62 / 0.61 | 0.96 to 0.98 | ~0 |
| LPCVD SiN, 1.6 um | 0.31 / 0.36 | 0.43 / 0.43 | 0.97 | ~0 |
| crystalline Si on sapphire, 1.0 um | 0.46 / 0.51 | 0.46 / 0.51 | 0.91 to 0.92 | 6 to 9 % |
| crystalline Si on sapphire, 1.2 um | 0.30 / 0.54 | 0.48 / 0.52 | 0.92 to 0.95 | 7 to 11 % |
| GaN on sapphire, 1.6 um (0.02 um grid only) | 0.33 / 0.26 | 0.45 / 0.32 | 0.88 to 0.93 | ~0 |

(results/dual_library_fine, dual_library_si_fine, dual_library_hi; plate_efficiency.json.)

Findings. (i) TiO2 and GaN offer no index advantage over Si-rich SiN: their indices at
852 to 1064 nm are only 0.03 to 0.05 higher. TiO2 performs like Si-rich SiN. GaN is worse
because sapphire forces the smaller pitch. Their practical value would lie in the process
(templated ALD and GaN etching reach high aspect ratios), i.e. in moving from the left to the
middle column `[VALUE?: aspect ratios demonstrated for TiO2 and GaN metasurfaces at
~1.6 um height]`. (ii) Crystalline Si needs only 1.0 um, so the aspect-ratio rule never
binds (h / 10 = 0.10 um = the minimum feature). It is the best film within aspect ratio 10
at 1064 nm (0.51) and comparable to Si-rich SiN at 852 nm (0.46 against 0.53). However, it
absorbs 6 to 11 % at 852 nm in its resonant atoms, and its single-wavelength ceiling is lower
(0.91 to 0.95). (iii) No film closes the gap. After the deflection factor, the best
two-wavelength doublets reach ~0.54 / 0.51 (Si-rich SiN or TiO2, aspect ratio 15) or
~0.40 to 0.46 (aspect ratio 10), against 0.80 to 0.83 for a single-wavelength doublet. The
limit is set by how densely a two-parameter, fabricable library covers the phase torus, not
by the index (Provisional: LPA, screen numerics, two families, 0.01 um grid).

Updated recommendation. For a two-wavelength objective, Si-rich SiN or ALD TiO2 at 1.6 um with
the highest aspect ratio the process allows, on the 0.01 um grid or finer (~0.5 per
wavelength at aspect ratio 15). If aspect ratio 10 is a hard limit, crystalline Si on
sapphire is the alternative at 1064 nm, at the price of 852 nm absorption. The single-
wavelength 852 nm objective (0.80 to 0.83) remains our preferred option unless 1064 nm must
pass through the same optic.

## 6. Meta-atom efficiency at the required deflection

Method (scripts/metalens_efficiency.py; results/library_analysis_852, _1064;
results/metalens_doublet/efficiency.json). For each film we use the pillar family chosen from
the single-pillar library at that wavelength (square SiN pillars on fused silica, fabrication
limits 0.10 um feature and gap, aspect ratio <= 10; analyse_library.py). We compute the
first-order efficiency of blazed supercells of N = 3, 4, 5, 6 and 8 pillars by RCWA (fmmax,
80 Fourier terms per pillar, 48 x 48 grid per pillar), for x and y polarisation, and take the
mean of the two. The efficiency is interpolated linearly in the local deflection
s = |dphi/dr| / k0, with the LPA value at s = 0. It is then averaged over each plate, weighted
by a uniform pupil followed along the design rays.

| Wavelength, film | Family (pitch, height) | eta at s ~ 0.5 (x / y) | Corrector | Focusing plate (rings 0-0.5 / 0.9-1.0 R1) | Doublet |
|---|---|---|---|---|---|
| 852 nm, LPCVD (n ~ 2.0) | 0.55 um, 1.7 um | 0.76 / 0.60 (s = 0.52) | 0.967 | 0.825 (0.90 / 0.73) | 0.80 |
| 852 nm, Si-rich (n ~ 2.3) | 0.45 um, 1.3 um | 0.69 / 0.76 (s = 0.47) | 0.963 | 0.851 (0.91 / 0.74) | 0.82 |
| 1064 nm, LPCVD | 0.60 um, 1.6 um | 0.82 / 0.78 (s = 0.44) | 0.960 | 0.841 (0.90 / 0.79) | 0.81 |
| 1064 nm, Si-rich | 0.50 um, 1.6 um | 0.82 / 0.70 (s = 0.53) | 0.961 | 0.862 (0.92 / 0.81) | 0.83 |

Hence a single-wavelength doublet is expected to focus 0.80 to 0.83 of the transmitted light
(Provisional). The losses are reflection (4 to 28 % in the steepest supercells) and light
in the zero and -1 orders. Three limitations apply. (i) The supercells sample s coarsely
(N = 3 and 4 are the only points above s = 0.35), and the N = 3 values scatter between
polarisations and films by +/-0.08, so the outer-ring values are uncertain at that level.
(ii) The averaging ignores the correlation between the zones of the two plates and the
polarisation structure of a linearly polarised beam. (iii) Substrate Fresnel losses
(~3.5 % per uncoated silica face) are excluded; anti-reflection coatings are assumed. The
outer-ring deficit (0.73 to 0.81 against 0.90 to 0.92 near the axis) apodises the pupil
slightly. This lowers the peak intensity relative to the transmitted power and broadens the
spot by an amount not computed here. These numbers are for one wavelength per plate. A
two-wavelength meta-molecule library (Section 5) will have lower efficiency, by an amount not
known before that library is simulated; Arbabi et al. measured 65 % at NA 0.46.

## 7. Mechanical tolerances

852.3 nm; image heights -400, -200, 0, 200, 400 um along the perturbation; one common focal
plane re-optimised per case (refocusing allowed); same input angles as the nominal design.
Tilts are modelled as exact rotations of the downstream optics (raytrace.Frame).

| Perturbation | Strehl at -400 / -200 / 0 / 200 / 400 um | Best common plane dz | Note |
|---|---|---|---|
| none | 0.9992 / 0.9995 / 0.9992 / 0.9995 / 0.9992 | -0.10 um | |
| spacing +50 um | 0.992 / 0.997 / 0.998 / 0.997 / 0.992 | -0.03 um | magnification change 1.4e-4 |
| spacing +/-200 um | 0.90 to 0.91 / 0.96 / 0.98 / 0.96 / 0.90 to 0.91 | +0.17 / -0.39 um | magnification change 5.5e-4 |
| wall +50 um | 0.988 / 0.989 / 0.989 / 0.989 / 0.988 | -33.2 um | atoms follow the focus |
| wall +200 um | 0.85 (all fields) | -132.7 um | uncorrected spherical aberration |
| decentre 2 um | 0.986 / 0.987 / 0.987 / 0.987 / 0.988 | | image shifts 0.96 x decentre |
| decentre 5 um | 0.921 / 0.924 / 0.925 / 0.925 / 0.928 | | |
| decentre 10 um | 0.73 (all fields) | | |
| decentre 25 um | 0.12 | | |
| plate tilt 0.1 deg | 0.05 to 0.07 | | equivalent to a 28 um decentre |
| plate tilt 0.1 deg, coma nulled by a 31 um decentre | 0.962 / 0.985 / 0.999 / 0.994 / 0.962 | | |
| plate tilt 0.25 deg, nulled (78 um) | 0.78 / 0.93 / 0.999 / 0.95 / 0.78 | | |
| plate tilt 0.5 deg, nulled (155 um) | 0.34 / 0.76 / 0.998 / 0.80 / 0.34 | | |
| cell tilt 0.1 deg | 0.933 / 0.967 / 0.971 / 0.957 / 0.933 | -0.14 um | |
| cell tilt 0.25 deg | 0.65 / 0.80 / 0.84 / 0.78 / 0.65 | | |
| cell tilt 0.5 deg | 0.20 / 0.40 / 0.48 / 0.38 / 0.20 | | |

Interpretation. The decentre sensitivity comes from the corrector. Its ~155 waves of
aspheric phase cancel the spherical aberration of the focusing plate, and a lateral offset
delta of a profile with spherical aberration W ~ a r^4 produces coma ~ 4 a delta r^3
(Provisional; consistent with the 2 to 10 um results). A tilt eps of the focusing plate
displaces it, relative to the incoming beam, by about the spacing times eps (16 mm x
1.75 mrad = 28 um). A decentre of 31 um nulls the axial coma of a 0.1 deg tilt, and the residual
(astigmatism, 0.96 at the corner) is acceptable. Hence the alignment procedure is: centre the
focusing plate on the corrector while observing the axial focus through the cell (coma
nulling), with tip/tilt kept within ~0.1 deg (autocollimation resolves far less). The
required centring, <= 2 to 5 um, is comparable to that of high-NA refractive objectives
`[REF?: typical centring tolerances of NA 0.5 objectives]`. It is set by mount
adjustment, not by the lithography, whose overlay between plates is irrelevant once the
plates are aligned optically. The cell tilt (<= 0.1 deg) is the same specification as for the
hologram (docs/prefab_checks_report.md, Section 2), so one alignment of the cell normal
serves both sides. The wall thickness must be known to +/-50 um, or the design recomputed
for the measured value: a +200 um error leaves 0.85, uncorrectable by refocusing.

## 8. Lens set for the experiment: two lenses, or one lens for all wavelengths

Requirement (USER): a 1075 nm lens (side A), and a lens for 852, 895 and 1064 nm (side B),
in the two-sided architecture; preferably one lens for all four wavelengths. All results are
SIM (review before submission), in ray optics unless stated. Cell stack 2.5 mm air | 4 mm
fused silica | 9 mm vacuum on both sides `[CONFIRM?]`.

### 8.1 Lens A: 1075 nm, NA 0.55 (configs/metalens_1075.yaml; results/metalens_1075)

This is the same air-spaced aplanatic doublet at NA 0.55 (as the side-A hologram): R1 =
9.21 mm, f = 16.75 mm. The spacing scan gives an optimum at 17.0 mm. In one common plane the
Strehl ratio is >= 0.9992 over 0 to 400 um; at free focus it is 0.968 at 1.0 mm and 0.890 at
1.3 mm. The axisymmetric wave check gives 0.9974. Tolerances (Strehl over -400 / 0 / 400 um):
decentre 2 um 0.984 to 0.988, 5 um 0.913 to 0.921; spacing +50 um 0.990; wall +50 um 0.983;
plate tilt 0.1 deg with coma nulled 0.959 / 0.998; cell tilt 0.1 deg 0.929 / 0.965. The
focusing plate needs a deflection of 0.55 at R2 (0.63 at the field-extended radius). The
estimated doublet efficiency from the 1075 nm supercells is 0.80 (LPCVD) / 0.82 (Si-rich).
The LPCVD supercells reach only s = 0.45, so the outer ring is extrapolated and this value
is optimistic. The side-A meta-hologram already combines this lens with the array hologram.
Lens A is needed only if the 1075 nm pattern is generated upstream (SLM or AOD) instead.

### 8.2 One metasurface doublet for several wavelengths

A metasurface phase map works at one wavelength (Section 5), and multi-wavelength meta-atoms
lose about half of the light (Sections 5b, 5c). We therefore keep single-wavelength meta-atoms
and move the wavelength dependence out of the shared lens. The shared doublet is designed at
a wavelength lambda_d between the operating wavelengths. Each operating wavelength then
receives, in its own beam path before the beams are combined, a two-surface aplanatic
pre-corrector. Surface A sits a distance L before surface C, and C is relayed 1:1 onto the stop
at the corrector metasurface (src/metatrap/metalens/multiwave.py, aplanatic_precorrector).
The pre-corrector is constructed exactly by back-tracing from the focus through the
diffractively dispersed plates. It restores both the stigmatic axial focus and the sine
condition at that wavelength, so the shared lens works at its NA over the field. The residual
is astigmatism and field curvature from the wavelength mismatch. A single stop corrector
without surface A restores only the axis: at 895 nm it gives Strehl 0.90 at 100 um and
0.16 at 400 um, and at 1064 nm it fails at once.

Minimum Strehl ratio over 0 to 400 um in one common plane per wavelength
(results/multiwave/optics.json; scripts/metalens_multiwave.py):

| Shared lambda_d, L | 852.3 nm | 894.6 nm | 1064.4 nm | 1075 nm |
|---|---|---|---|---|
| 852.3 nm (no correction needed at 852), 5 mm | 0.999 | 0.993 | 0.37 | - |
| 900 nm, 5 mm | 0.999 | 0.999 | 0.74 | 0.68 |
| 950 nm, 5 mm | 0.995 | 0.999 | 0.94 | 0.92 |
| 960 nm, 5 mm | 0.993 | 0.999 | 0.959 | 0.943 |
| 990 nm, 3 mm | 0.953 | 0.988 | 0.987 | 0.979 |
| 990 nm, 5 mm | 0.985 | 0.996 | 0.987 | 0.980 |
| 990 nm, 10 mm | 0.999 | 0.999 | 0.978 | 0.970 |

The pre-correctors are mostly power: 200 to 1,500 waves peak-to-valley, e.g. 430 / 530 waves
(A / C) at 1064 nm for lambda_d = 990 nm, L = 5 mm. Each is a weak refractive asphere pair or a
pair of single-wavelength phase plates in that wavelength's path. Without a relay (C placed
20 to 40 mm before the shared lens, in front of the dichroic combiner) the Strehl ratio drops
to 0.5 to 0.8 or worse, so the pupil relay is required. A 1:1 relay that is achromatic over
852 to 1075 nm is assumed ideal here `[not simulated: relay aberrations at 16 mm pupil
diameter and +/-1.4 deg field]`.

Meta-atom efficiency at the operating wavelengths (results/multiwave/atoms.json;
scripts/multiwave_atoms.py; square pillars on fused silica, aspect ratio <= 10, 197-term
RCWA). Atoms are chosen for the phase at lambda_d; at another wavelength the first diffraction
order carries the intended phase profile.

| Family, lambda_d | First order per plate: 852 / 895 / 1064 / 1075 nm | Zero order per plate |
|---|---|---|
| LPCVD, p 0.50, h 1.6 um, lambda_d 960 nm | 0.82 / 0.89 / 0.88 / 0.87 | 1 to 4 % |
| LPCVD, p 0.50, h 1.6 um, lambda_d 990 nm | 0.76 / 0.84 / 0.92 / 0.91 | 2 to 5 % |
| Si-rich, p 0.45, h 1.3 um, lambda_d 960 nm | 0.74 / 0.88 / 0.87 / 0.86 | 1 to 5 % |

With the deflection factor of the focusing plate and the corrector (~0.84 in total, Section 6),
the doublet efficiency is ~0.84 eta^2. For LPCVD, p 0.50, h 1.6 um this gives 0.57 / 0.66 /
0.66 / 0.64 (lambda_d = 960 nm) or 0.48 / 0.60 / 0.71 / 0.69 (990 nm), against 0.80 for a
single-wavelength doublet (Provisional: LPA, supercell factor transferred from 852 nm,
pre-corrector and relay losses excluded). The zero order of each plate (1 to 5 %) is only
weakly focused and spreads over millimetres in the atom plane. Its intensity there is ~1e-6
of a focal spot (area ratio ~1 um^2 / 1 mm^2, DERIVED), negligible for trapping.

### 8.3 Options

| Option | Lens(es) | Strehl over +/-400 um | Doublet efficiency | Notes |
|---|---|---|---|---|
| (a) Two lenses, as proposed | A: 1075 nm, NA 0.55 (or the hologram); B: shared doublet lambda_d = 960 to 990 nm with pre-correctors at 852 / 895 / 1064 | A >= 0.999; B >= 0.96 (960) or >= 0.985 (990) | A 0.80; B 0.57 / 0.66 / 0.66 (960) or 0.48 / 0.60 / 0.71 (990) | B needs three pre-corrector modules and a pupil relay |
| (b) One lens for all four | shared doublet lambda_d = 990 nm, L = 5 mm, NA 0.5 | >= 0.980 at all four | 0.48 / 0.60 / 0.71 / 0.69 | 1075 nm at NA 0.5 instead of 0.55; high-power 1075 nm through the shared lens |
| (c) Two-wavelength meta-atoms, for comparison | doublet with dual atoms (852 / 1064) | as single-wavelength | ~0.5 (AR 15), 0.3 to 0.45 (AR 10) | no 895 / 1075 without further correction |

Our choice is (a) with lambda_d = 990 nm for lens B, or (b) if one lens on one side is
preferred. The trade-off is as follows. First, (b) gives full flexibility, since every
wavelength reaches every site through the same optic, and keeps the cell accessible from
the other side. However, it lowers the 1075 nm NA from 0.55 to 0.5 and its efficiency from
~0.80 to ~0.69. This costs ~1.4 x in power per trap depth, since the depth scales as
eta NA^2, i.e. (0.80 x 0.55^2) / (0.69 x 0.50^2) = 1.40 (DERIVED). In addition, 50 to 100 W of
1075 nm then passes through the shared lens, so it must be LPCVD SiN on fused silica or on
sapphire (Section 5 of the damage model in docs/prefab_checks_report.md). Second, in both
options the 852 nm collection efficiency of lens B is ~0.5 instead of ~0.8, the price of
sharing the lens. A dedicated 852 nm doublet (0.80) remains the choice if imaging fidelity
dominates. Third, the pre-correctors and the relay are conventional optics but add alignment
degrees of freedom. Their centring tolerance against the shared lens is expected to be similar
to the plate decentre tolerance (a few um, Section 7) `[not simulated]`.

Next steps: (i) a meta-atom library at lambda_d = 990 nm with supercells over the
focusing-plate deflections at each operating wavelength; (ii) tolerances of the pre-corrector
modules and the relay; (iii) the wave-optics check of the corrected system off axis
(Section 4b method) at 852 and 1075 nm.

### 8.4 Checks of the shared lens (lambda_d = 990 nm, L = 5 mm)

Wave optics off axis (results/multiwave/wavecheck.json; scripts/multiwave_wavecheck.py). We
used the superposition method of Section 4b, with the pre-corrector surfaces included in the
ray part. At 400 um it gives Strehl ratios 0.985 / 0.996 / 0.988 / 0.982 at 852 / 895 / 1064 /
1075 nm, against 0.987 / 0.997 / 0.989 / 0.984 from the ray trace. On axis it gives 1.000 at
all four, and the foci agree to 1e-4 um.

Efficiency from supercells (results/multiwave/supercells.json, doublet_efficiency.json;
scripts/multiwave_supercells.py). The family is LPCVD SiN, p 0.50 um, h 1.6 um. Supercells of
N = 3, 4 and 6 pillars, chosen at 990 nm, were solved by RCWA at each operating wavelength,
x and y polarisation averaged. The first-order efficiency versus deflection (LPA value at
s = 0) is:

| Wavelength | s = 0 (LPA) | N = 6 | N = 4 | N = 3 | Corrector | Focusing plate | Doublet |
|---|---|---|---|---|---|---|---|
| 852.3 nm | 0.76 | 0.63 (s 0.28) | 0.50 (0.43) | 0.69 (0.57) | 0.75 | 0.61 | 0.46 |
| 894.6 nm | 0.84 | 0.69 (0.30) | 0.65 (0.45) | 0.73 (0.60) | 0.84 | 0.70 | 0.59 |
| 1064.4 nm | 0.92 | 0.79 (0.35) | 0.70 (0.53) | 0.67 (0.71) | 0.91 | 0.78 | 0.71 |
| 1075 nm | 0.91 | 0.78 (0.36) | 0.69 (0.54) | 0.64 (0.72) | 0.90 | 0.77 | 0.69 |

These confirm the estimate of Section 8.2 (0.48 / 0.60 / 0.71 / 0.69) to within 0.02. The
supercells reflect 6 to 20 % and send 2 to 19 % into the zero order. At 1064 and 1075 nm,
where the blaze is short of 2 pi, they send 4 to 13 % into the -1 order. The non-monotonic
852 nm values (N = 4 below N = 3) indicate a near-field resonance of this family at 852 nm.
A family chosen for the multi-wavelength figure of merit, rather than taken from the
single-wavelength ranking, should recover part of this (next step).

Pre-corrector module tolerances (results/multiwave/precorrector_tolerances.json; ideal
relay). The worst Strehl over -400 / 0 / 400 um after refocusing is as follows.

| Error | 852.3 nm | 1064.4 nm |
|---|---|---|
| none | 0.985 | 0.986 |
| module decentre 5 um | 0.966 | 0.983 |
| module decentre 20 um | 0.74 | 0.94 |
| module tilt 0.1 deg | 0.95 | 0.97 |
| module tilt 0.5 deg | 0.48 | 0.74 |
| A-C spacing +20 um | 0.984 (focus -4.0 um) | 0.986 (focus -1.6 um) |
| A-C spacing +100 um | 0.95 (focus -18 um) | 0.98 (focus -5.9 um) |

The 852 nm module carries the largest correction (~1,000 / 1,500 waves) and is the most
sensitive. Each module needs centring to ~5 um and tilt to <= 0.1 deg against the shared lens,
the same class as the plate tolerances of Section 7. The A-C spacing acts as a fine focus
adjustment (0.18 um of focus per um at 852 nm), which brings all wavelengths into the atom
plane. Relay aberrations remain unmodelled.

A comparison with infinity-corrected refractive objectives is given in
docs/objective_comparison.md.

## 9. Limitations and next steps

Limitations: (i) scalar optics throughout; off axis, the field between the plates is ray
optical (Section 4b), and no vector effects are included. (ii)
Ideal phase maps in Sections 1 to 5 and 7: the meta-atom efficiency enters only through the
aperture averages of Section 6, not through the imaging (no quantisation or local-periodic-
approximation phase error is propagated). (iii) Plate substrates and
their tolerances, thermal effects of the 1064 nm beam on the plates, and polarisation (vector
focusing: an 8 % peak reduction at NA 0.55 on side A, results/vector_check) are not included. (iv)
Collection efficiency for imaging (NA 0.5 collects 6.7 % of isotropic fluorescence, DERIVED:
(1 - cos asin 0.5) / 2) is set by the NA, not modelled here.

Next steps: (i) decide between a single-wavelength objective and a two-wavelength one at
roughly half the efficiency (Sections 5b, 5c), and obtain the achievable aspect ratio; (ii) denser supercell sampling above s = 0.35 (non-integer
blazes or larger supercells), to firm up the outer-ring efficiency of Section 6; (iii) a substrate-bearing
version of the design (layered spacer) and its tolerances.

## Provenance ledger

| Value | Location | Type | Pointer |
|---|---|---|---|
| Strehl, focus, deflections, tolerances, chromatic data | Sections 1 to 3, 5, 7 | SIM, review before submission | results/metalens_doublet/doublet.json |
| Plate and doublet efficiencies | Section 6, Summary | SIM, review before submission | results/metalens_doublet/efficiency.json; results/library_analysis_852, _1064/supercell.json |
| Off-axis superposition Strehl, foci, widths, apodisation | Section 4b | SIM, review before submission | doublet.json, wave_offaxis |
| Film indices TiO2, GaN, c-Si | Section 5c | SOURCE | Siefke 2016; Barker 1973; Schinke 2015 (refractiveindex.info database files, data/materials/) |
| Other-film and fine-grid efficiencies | Section 5c | SIM, review before submission | results/dual_library_fine, _si_fine, _hi |
| Lens A design, tolerances, efficiency | Section 8.1 | SIM, review before submission | results/metalens_1075 |
| Shared-lens Strehl ratios with pre-correctors | Section 8.2 | SIM, review before submission | results/multiwave/optics.json; 900 / 950 nm rows from an inline scan (not stored) |
| Shared-lens wave check, supercells, module tolerances | Section 8.4 | SIM, review before submission | results/multiwave/wavecheck.json, supercells.json, doublet_efficiency.json, precorrector_tolerances.json |
| Cross-wavelength atom efficiencies | Section 8.2 | SIM, review before submission | results/multiwave/atoms.json |
| Two-wavelength library efficiencies | Section 5b | SIM, review before submission | results/dual_library/dual_efficiency.json, plate_efficiency.json |
| Wave-optics Strehl, waists | Section 4 | SIM, review before submission | doublet.json, wave_axial |
| Chief-ray direction cosine 3.6e-4 at 400 um | Section 1 | SIM (inline evaluation, not stored) | raytrace of the 16 mm design |
| Polynomial-optimisation failures | Section 1 | SIM | results/metalens_study/metalens_study.json |
| Airy 1/e^2 radius 0.41 lambda / NA | Section 4 | DERIVED | (2 J1(v)/v)^2 = e^-2 at v = 2.58 |
| Depth of focus lambda / NA^2 = 3.4 um | Section 3 | DERIVED | 0.852 / 0.25 |
| Fresnel number 4,800 | Section 4 | DERIVED | 8106^2 / (0.852 x 16000) |
| Arbabi doublet: 800 um, f/0.9, > 60 deg, ~70 % | Section 1 | SOURCE | Nat. Commun. 7, 13682 (2016), abstract (search result) |
| Meta-molecules: 1550/915 nm, NA 0.46, 65 % | Section 5 | SOURCE | Optica 3, 628 (2016), Caltech repository record |
| Delay-bandwidth bound | Section 5 | SOURCE (arXiv) | arXiv:2001.10899 |
| Cs D2 natural linewidth 5.2 MHz | Section 5 | recalled `[VALUE?: verify, Steck]` | |

Open placeholders: `[CONFIRM?: side-B air gap, wall, distance to the atoms]`;
`[CONFIRM?: which side-B functions need 852, 895 and 1064 nm]`; `[REF?: Kalvach and Szabo]`;
`[REF?: centring tolerances of NA 0.5 objectives]`; `[VALUE?: Cs D2 linewidth]`;
`[VALUE?: aspect ratios demonstrated for TiO2 / GaN metasurfaces at ~1.6 um]`;
`[CONFIRM?: silicon-on-sapphire film thickness, ALD TiO2 height]`.
