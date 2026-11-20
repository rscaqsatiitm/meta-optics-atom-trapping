# Metasurface lenses versus infinity-corrected refractive objectives for the Cs tweezer platform

Comparative analysis. Requirement set (from this project): NA 0.5 to 0.55 at the atoms; focus
behind the cell stack 2.5 mm air | 4 mm fused silica | 9 mm vacuum `[CONFIRM?]`; flat field
over at least +/-400 um (v2 array, 20,164 sites); wavelengths 852.3 nm (imaging), 894.6 nm
(addressing), 1064.4 nm and 1075 nm (trapping); 50 to 100 W of trap light. Numbers are labelled
SOURCE (read in this session), SEARCH (taken from a search-result summary, page not read),
SIM or DERIVED (this project), or `[VALUE?]` (not verified).

## Summary

A catalogue infinity-corrected objective cannot be used through this cell at the required NA.
An air-corrected objective focusing through 4 mm of fused silica carries 0.47 waves rms of
spherical aberration at NA 0.4 and 1.25 waves rms at NA 0.5 at 852 nm (DERIVED, exact plate
path minus best defocus). The Strehl ratio stays above 0.8 only below NA ~0.25. The realistic
refractive alternative is therefore a custom, window-corrected, multi-element objective. An
example is the Special Optics NA 0.65 objective of the 6,100-atom Caltech array, with a field
of view of 1.5 mm diameter (SOURCE). Against that benchmark, the metasurface lenses of this
project offer the following. (i) In ideal-phase models (unit transmission, angle-independent
meta-atom phase, scalar optics) a field quality at NA 0.5 to 0.55 comparable to that benchmark
(ray-traced Strehl >= 0.999 over +/-400 um, >= 0.8 to +/-1.4 mm, SIM). The angular response of
the meta-atoms, which a built objective includes, is not modelled, so this comparison is an upper
bound for the metasurface. (ii) A power-handling advantage
that targets the limitation reported for the refractive objective (thermal aberrations above
~35 W). (iii) Direct generation of 20,000 traps without an SLM. They pay for this with lower
efficiency (0.80 for one wavelength, 0.5 to 0.7 when shared between wavelengths, SIM, against
a measured 0.58 for the only metalens tweezer so far, SOURCE) and with intrinsic single-wavelength
operation. Our recommendation: a metasurface on the high-power 1075 nm side, where its
advantages are decisive. On the multi-wavelength imaging side, a window-corrected achromatic
objective is the lower-risk choice, unless size, in-vacuum mounting or the
shared-lens flexibility of Section 8 of docs/metalens_report.md is required.

## 1. The comparison set

| Property | Catalogue NIR plan-apochromat (Mitutoyo M Plan Apo NIR) | Custom window-corrected objective (Special Optics, Avantier) | Metasurface singlet, demonstrated (Hsu et al. 2022) | Metasurface doublet, this project |
|---|---|---|---|---|
| NA | 0.40 (20X), 0.42 (50X), 0.70 (100X HR) (SEARCH) | 0.65 at 1055/1061 nm, 0.55 at 852 nm (SOURCE, ref. 1); 0.7 (SEARCH, ref. 4) | 0.55 (SOURCE, ref. 2) | 0.50 (side B), 0.55 (side A) |
| Working distance | 20 / 17 / 10 mm (SEARCH) | 14 mm with 3 mm quartz window (SEARCH, ref. 4) | 3 mm (f = 3 mm) | focusing plate 2.5 mm from the cell; 15.5 mm stack to the atoms |
| Window correction | none (cover-glass versions to 0.7 mm only, SEARCH) | designed for the specific window | none (lens inside the cell) | exact for the specified stack; redesign for a measured wall is a new mask |
| Field (Strehl >= 0.8) | not specified at NA > 0.25 through 4 mm glass (DERIVED: unusable) | 1.5 mm diameter (SOURCE, ref. 1) | +/-11 um (+/-0.2 deg) (SOURCE) | +/-1.4 mm (852 nm, ray), >= 0.999 over +/-400 um (SIM) |
| Wavelengths | 480 to 1800 nm corrected (SEARCH) | NA and field differ by wavelength; 852 nm at reduced NA 0.55 (SOURCE) | one (852 nm design; 780 nm focus shifted by +300 um, SOURCE) | one per phase map; four with per-wavelength pre-correctors (SIM) |
| Focusing efficiency / transmission | `[VALUE?: ~0.7 to 0.9 at 1064 nm]` | measured, value not given in ref. 1 | 0.58 measured (0.78 theory; 14 % reflection, 6 % zero order, 2 % higher orders), zero order 13 % measured (SOURCE) | 0.80 to 0.82 single wavelength; 0.48 to 0.71 shared (SIM, LPA) |
| Power handling | `[VALUE?]` | ~35 W into the objective; "worsening optical aberrations at tweezer powers greater than that in the present study due to thermal heating of the objective. This is the main limitation on atom number" (SOURCE, ref. 1) | not tested at high power | damage model: not stress-limited; 150 W (Si-rich on FS, upper absorption) to > 1 kW (LPCVD or sapphire) (SIM) |
| Size, mass | ~ 34 mm diameter barrels `[VALUE?]` | larger than standard (SEARCH) | 4 x 4 mm, 0.66 um a-Si on 0.5 mm FS | two plates, 17 to 19 mm diameter, 0.5 to 1 mm thick, 16 to 17 mm apart |
| Vacuum compatibility | not for in-vacuum use `[VALUE?]` | in-vacuum five-singlet designs exist (SEARCH, ref. 5) | UHV compatible and bakeable (SOURCE) | SiN on FS or sapphire: same class of materials |
| Array generation | needs SLM / AOD upstream | SLM: 65 % diffraction efficiency at low power, ~45 % at full power, zero order filtered (SOURCE, ref. 1) | AOD (4 traps) | static hologram in the lens (side A, 20,164 sites, SIM); metasurface holograms with > 1000 atoms shown (SEARCH, ref. 3) |
| Alignment | factory-aligned barrel; tilt to window < 1 mrad (SEARCH, ref. 6) | same | single element | plate centring 2 to 5 um; cell tilt <= 0.1 deg (1.7 mrad) (SIM) |
| Availability | catalogue | custom, months `[VALUE?]` | research | research; 18 mm EBL write (stitching specs in docs/prefab_checks_report.md) |

## 2. Advantages of the metasurface lenses

Power handling. The refractive objective of ref. 1 limits the atom number through thermal
aberration at ~35 W, even after SLM correction. A multi-element objective absorbs in several
glasses, cements and coatings over centimetres of path, in a mount whose thermal design is not
the optical designer's. The metasurface puts 1.6 um of film on a 0.5 to 1 mm substrate. Our
damage model gives a thermal focal shift of 0.1 to 0.5 um at 100 W for Si-rich SiN on fused
silica, and 0.001 um for LPCVD SiN (SIM; docs/prefab_checks_report.md, Section 5). The
comparison is Provisional: we have not modelled the heating of a refractive objective, and
ref. 1 gives no absorption figures.

Field and window correction by construction. The aplanatic construction is exact for the
actual stack, so the design compensates the 4 mm wall at full NA without residual (Strehl
1.000 on axis, wave optics 0.998 to 0.999, SIM). A different wall thickness or a measured
value enters as a new phase map, i.e. a new mask, not a new glass design.

Array generation without an SLM. On side A the hologram and the lens are one element. This
removes the SLM diffraction efficiency (45 % at full power in ref. 1), the zero-order filter
and the LCOS heating. It is the largest single efficiency gain of the metasurface approach for
20,000 sites (0.80 to 0.86 array efficiency, SIM, against ~0.45 x objective transmission in
ref. 1).

Size and placement. Two sub-millimetre plates can sit 2.5 mm from the cell, or inside the
vacuum system (ref. 2), and they leave the other side of the cell free.

## 3. Disadvantages of the metasurface lenses

Efficiency. For one wavelength the theoretical doublet efficiency is 0.80 to 0.82 (SIM, LPA,
supercell edge factor). The only measured metalens tweezer reached 0.58 against 0.78 in
theory, which the authors attribute to pillar shape and sidewall errors (ref. 2). A similar
fabrication penalty would bring our 0.80 to ~0.6 `[estimate, scaling ref. 2]`. In imaging this
is a direct loss of collected photons. In trapping it costs power, which the high damage
threshold partly offsets.

Single-wavelength nature. Any phase map focuses one wavelength. Achromatic dispersion
engineering is excluded at a 16 mm aperture by the delay-bandwidth bound (Section 5 of the
metalens report). Two-wavelength meta-atoms lose about half of the light (Sections 5b, 5c), and
a shared lens needs per-wavelength pre-correctors with a pupil relay (Section 8). A refractive
apochromat covers 480 to 1800 nm with one element.

Zero order and scattering. 1 to 5 % per plate in our design and 13 % measured in ref. 2 is not
focused. At the traps its intensity is ~1e-6 (ours, DERIVED) to 1e-8 (measured, ref. 2) of the
focal intensity, so it does not perturb the traps, but it is lost for imaging.

Fabrication and maturity. An 18 mm metasurface is a large-area electron-beam write with
stitching, dose and film-uniformity requirements (docs/prefab_checks_report.md, Section 3).
There is no supplier, and a design error costs a fabrication cycle.

Alignment. The two plates must be centred to 2 to 5 um, with coma nulling for tilt (metalens
report, Section 7). A commercial objective is aligned in its barrel. The cell-tilt
specification (0.1 deg) is comparable to the < 1 mrad reported for a high-NA refractive
objective against its window (ref. 6, SEARCH).

## 4. Recommendation

Side A (1075 nm, 50 to 100 W, 20,000 sites): metasurface. Its advantages, power handling and
SLM-free array generation, address the limitations reported for the refractive state of the
art (ref. 1), and the single-wavelength penalty does not apply. Use LPCVD SiN, or Si-rich SiN
on sapphire, for the thermal margin.

Side B (852 / 895 / 1064 nm): a custom window-corrected achromatic objective at NA 0.5 to 0.65
is the lower-risk choice for imaging fidelity and wavelength flexibility. The metasurface
shared lens with pre-correctors (Strehl >= 0.98 at all four wavelengths, efficiency 0.48 to
0.71, SIM) is preferable if one or more of the following matter more than collection
efficiency: (i) compactness, (ii) in-vacuum mounting, (iii) a single lens for all four
wavelengths, (iv) cost of iterations once a process exists. A practical path is to start
side B with a refractive objective and qualify a metasurface version in parallel. The two are
interchangeable because both are infinity-conjugate at the stop.

## References

1. Manetsch et al., "A tweezer array with 6100 highly coherent atomic qubits", Nature 647, 60
   (2025); [arXiv:2403.12021](https://arxiv.org/abs/2403.12021) (read in session: objective NA,
   field, power, SLM efficiency, thermal limitation).
2. Hsu et al., "Single-atom trapping in a metasurface-lens optical tweezer", PRX Quantum 3,
   030316 (2022); [arXiv:2110.11559](https://arxiv.org/abs/2110.11559) (read in session).
3. Holman et al., "Trapping of single atoms in metasurface optical tweezer arrays",
   [arXiv:2411.05321](https://arxiv.org/abs/2411.05321) (search result only; page not read).
4. Avantier, "High NA and long working distance microscope objective lenses" (vendor case
   study, [link](https://www.vision-systems.com/directory/optics/blog/55294871/avantier-inc-high-na-and-long-working-distance-microscope-objective-lenses);
   search result only).
5. "A ship-in-a-bottle quantum gas microscope setup for magnetic mixtures",
   [arXiv:2306.05404](https://arxiv.org/abs/2306.05404) (search result only).
6. "A modular quantum gas platform", [arXiv:2501.08211](https://arxiv.org/abs/2501.08211)
   (search result only).
7. Mitutoyo M Plan Apo NIR specifications via distributor listings
   ([Edmund Optics](https://www.edmundoptics.com/f/mitutoyo-nir-nuv-and-uv-lcd-infinity-corrected-objectives/39905/);
   search result only; distributor pages returned HTTP 403).

## Verification notes

Provenance ledger:

| Value | Location | Type | Pointer |
|---|---|---|---|
| Window aberration 0.47 / 1.25 waves rms (NA 0.4 / 0.5, 852 nm) | Summary, Table | DERIVED | exact plate OPL minus best defocus, uniform pupil; command in session log |
| NA 0.65 / 0.55, FOV 1.5 mm, 35 W, 1.2 mW per tweezer, SLM 65 % / 45 %, thermal limitation | Table, Sections 2, 4 | SOURCE | ref. 1, main text and Methods |
| NA 0.55, f 3 mm, 58 % measured, 78 % theory, zero order 13 %, FOV +/-11 um, a-Si 660 nm | Table, Section 3 | SOURCE | ref. 2 |
| Mitutoyo NA and WD | Table | SEARCH | ref. 7 |
| Avantier NA 0.7, WD 14 mm, 3 mm window | Table | SEARCH | ref. 4 |
| Metasurface doublet field, efficiency, tolerances, damage limits | all | SIM, review before submission | docs/metalens_report.md; docs/prefab_checks_report.md |

Open placeholders: `[CONFIRM?: cell stack]`; `[VALUE?: transmission of NIR objectives at
1064 nm]`; `[VALUE?: barrel sizes, vacuum compatibility of catalogue objectives]`;
`[VALUE?: custom objective lead time]`; `[estimate: fabrication penalty scaled from ref. 2]`.
Items marked SEARCH should be read at the source before they are quoted externally.
