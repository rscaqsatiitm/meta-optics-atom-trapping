# Scaling to 100,000 sites: does a shorter working distance help?

Technical note. Question (USER): would moving the metasurface closer to the atoms (in vacuum,
e.g. 3 mm) help to reach 100,000 traps? All numbers are SIM (review before submission) or
DERIVED. Code: scripts/scaling_100k.py, configs/scaling_100k.yaml, configs/scaling_sc_*.yaml;
results: results/scaling_100k/{scaling,summary}.json, results/scaling_100k_sc_si_rich.

## Result

No. A shorter working distance does not remove any of the constraints that bind at
100,000 sites, and it tightens two of them. The trapping hologram tolerates short distances
well, because exact-kernel design has no field aberration. At 3 mm the corner site of a
316 x 316 array at 4 um pitch keeps 0.88 of the axial depth, which costs 5 % more power. The
binding constraints are three. (i) Total power: ~210 W for 0.3 mK, ~125 W for 0.18 mK. This is
independent of the distance, but the intensity on the metasurface grows as 1 / WD^2, to
1.8 kW/cm^2 at 3 mm. (ii) The field of the imaging and addressing lens: at 3 mm it images only
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
meta-atom efficiency, averaged over the aperture with eta(s) from supercell RCWA (Si-rich SiN
at 1075 nm). The supercell data are library_analysis_1075 for s <= 0.54 and new N = 3, 4 and 5
supercells for p 0.45 um, h 1.3 um: 0.70 at s = 0.60 and 0.66 at s = 0.80. The relative depth
is D(r) = I_rel x eta(r) / eta(0). Uniform depth over a square array needs power in
proportion to 1 / D_j, so the array power factor is F = mean(1 / D_j) (DERIVED). The absolute
power is scaled from the simulated v2 design: 0.355 mK for 20,164 sites at 50 W, i.e.
143 mK x site / W at NA 0.55, including its array efficiency of 0.85.

Imaging lens. For each WD the air-spaced aplanatic doublet (NA 0.5, 852 nm) is designed, with
the spacing scanned around f. The field radius with a free-focus Strehl ratio >= 0.8 then
gives the number of sites inside the inscribed square, N = 2 r^2 / a^2 (DERIVED).

## Results

| WD (vacuum) | Aperture | Corner depth D, 4 um / 3 um pitch | Power factor F | Corner spot half-widths (radial / tangential) | Max deflection, corner | Power for 0.3 mK (0.18 mK) | Intensity on metasurface | Max absorbing defect (FS, 50 K) | Imaging-lens field radius | Sites imageable, 4 / 3 um pitch |
|---|---|---|---|---|---|---|---|---|---|---|
| 3 mm | 3.95 mm | 0.88 / 0.93 | 1.047 | 0.84 / 0.80 um | 0.69 | 219 W (131 W) | 1,790 W/cm^2 | 3.9 um | 450 um | 25,000 / 45,000 |
| 6 mm | 7.9 mm | 0.97 / 0.98 | 1.013 | 0.80 / 0.80 um | 0.63 | 212 W (127 W) | 430 W/cm^2 | 16 um | 760 um | 73,000 / 129,000 |
| 10 mm | 13.2 mm | 0.99 / 0.99 | 1.005 | 0.80 / 0.78 um | 0.60 | 211 W (127 W) | 155 W/cm^2 | 45 um | 1,160 um | 167,000 / 296,000 |
| 15.5 mm | 20.4 mm | 0.99 / 1.00 | 1.002 | 0.78 / 0.78 um | 0.58 | 210 W (126 W) | 64 W/cm^2 | 108 um | 1,610 um | 324,000 / 577,000 |
| 25 mm | 32.9 mm | 1.00 / 1.00 | 1.001 | 0.78 / 0.78 um | 0.57 | 210 W (126 W) | 25 W/cm^2 | 280 um | 2,020 um | 510,000 / 907,000 |

The axial spot is 0.78 um (Airy 1/e^2 radius 0.41 lambda / NA = 0.80 um). The 100,000-site
array corner lies at 891 um (4 um pitch) or 669 um (3 um pitch).

## Discussion

Why the hologram does not care much. In the exact-kernel design every site receives its own
spherical wave through the stack, so there is no field aberration to grow with angle; this
is what made the 20,000-site design possible. At short distances the corner sites lose only
through obliquity, which shrinks the projected aperture: the peak falls to 0.90 at 891 um
for WD = 3 mm, the spot elongates radially by 5 %, and the steeper deflections
(s up to 0.69, 13 % of the aperture above 0.6) cost 2 % in meta-atom efficiency. Hence
the trapping side alone would tolerate 3 mm.

Why closer still does not help. (i) Power. The depth per site is set by NA and power, not by
the distance. 100,000 sites at 0.3 mK need ~210 W of 1075 nm at NA 0.55 and array efficiency
0.85, whatever the distance. A higher NA is the lever, since the depth per watt scales as
NA^2, so NA 0.7 would need ~130 W. However, at NA 0.7 the metasurface edge deflection rises
to 0.7 to 0.8, where eta ~ 0.66. (ii) Damage. The same power on a smaller aperture raises the
intensity as 1 / WD^2. At 3 mm a single absorbing particle larger than 4 um radius heats by
more than 50 K on fused silica. The global temperature rise of a thin plate depends little
on the aperture, so the film choice is the same at every distance: at ~210 W only LPCVD SiN,
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
present cell. Next step: a full 100,000-site exact-kernel design at the cell geometry with
on-the-fly kernels, as a ~1 day job or on a larger machine.
