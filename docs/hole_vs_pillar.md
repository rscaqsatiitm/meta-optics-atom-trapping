# Hole (inverse, "Mie-void") meta-atoms versus pillars for the 1075 nm hologram

Technical note. Question: can circular air holes etched into a continuous Si-rich SiN film
replace the free-standing pillars of the B1 meta-hologram? The question was prompted by
reflective UV Mie-void metasurfaces in silicon (Schwab, Arslan et al., Advanced Photonics 2026,
doi:10.1117/1.AP.9.1.016002; read through the press release and a review, not the paper) and by
transmission-mode Mie-void metasurfaces in silicon-rich nitride (Goldberg, Mazurski and Levy,
arXiv:2604.10488, abstract only). The UV devices use the geometric phase of rotated elliptical
voids, which needs circular polarisation, and rely on the near-metallic reflectance of silicon at
325 nm; neither carries over to linearly polarised 1075 nm tweezers (bare silicon reflects ~32 %
at 1075 nm, DERIVED). What can carry over is the inverse geometry used for propagation phase in
transmission, which is what we test here.

## Result

Holes are not competitive. With the same solver settings and fabrication limits, the best hole
family reaches a 16-level efficiency (local-periodic approximation, averaged over a +/-10 nm width
bias) of 0.861, against 0.958 for pillars, needs a 2.8 um film instead of 1.6 um, does not reach
2 pi of phase coverage within the fabrication limits, and absorbs 5.7 times more power in the film.
The continuous film was expected to help with heat; it does the opposite, because the hole family
keeps two thirds of the film in the beam and is nearly twice as thick.

| Family (Si-rich SiN estimate, n = 2.277, on fused silica, 1075 nm) | Pitch / height | eta_16, +/-10 nm averaged | eta_16 | Phase coverage | Min. T0 | Smallest feature | Film in the beam | Absorbed fraction (k = 6.7e-5) |
|---|---|---|---|---|---|---|---|---|
| Pillars, aspect <= 10 | 0.50 / 1.6 um | 0.958 | 0.964 | 6.69 rad | 0.85 | 0.17 um pillar | 0.32 | 6.3e-4 |
| Holes, aspect <= 10 | 0.60 / 2.8 um | 0.861 | 0.903 | 5.47 rad | 0.78 | 0.28 um hole | 0.64 | 3.6e-3 (max 9.4e-3) |
| Pillars, aspect <= 15 | 0.50 / 1.6 um | 0.961 | 0.968 | 7.35 rad | 0.85 | 0.11 um pillar | 0.27 | 5.2e-4 |
| Holes, aspect <= 15 | 0.50 / 3.0 um | 0.875 | 0.884 | 5.99 rad | 0.77 | 0.21 um hole | 0.69 | 2.4e-3 |

(Film in the beam: material fill fraction averaged over the 16 selected meta-atoms. Absorbed
fraction: 1 - T - R of the 16 selected meta-atoms with the complex index, averaged; an unpatterned
film of 1.6 um absorbs 1.25e-3 and one of 2.8 um 2.19e-3. Selected family: the shortest within
0.3 % of the best, as for the pillar library.)

## Method

Square lattice, light incident from the fused-silica substrate, exit into air, normal incidence.
Holes: circular air holes through a film of thickness h (film continuous around them). Pillars:
square pillars of height h, as in the existing library. Pitches 0.40 to 0.60 um, heights 0.8 to
3.0 um, 40 widths or diameters from 0.06 um to p - 0.06 um (results/hole_vs_pillar/sweep.json;
configs/hole_vs_pillar.yaml; scripts/hole_vs_pillar.py). Both families are computed with fmmax in
the pol_fourier formulation, which shows no spurious cross-polarisation for C4v cells
(tests/test_rcwa.py): t_cross <= 8e-4 for the selected families and <= 9e-3 over the whole sweep
(largest at resonant holes, p 0.55, h 2.6 um), at 197 Fourier terms on a 256 grid. For holes
at p 0.50, h 1.6 um this agrees with 317 terms to <= 0.02 rad. Fabrication limits: minimum
feature and gap (wall between holes) 100 nm, aspect ratio (depth over smallest hole diameter, or
height over smallest pillar width) <= 10, or 15 as a relaxed case. Ranking:
scripts/hole_vs_pillar_rank.py (results/hole_vs_pillar/ranking.json); absorption:
scripts/hole_vs_pillar_absorption.py (absorption.json).

Consistency check. The best pillar family found here, p 0.50 um and h 1.6 um with eta_16 0.958
(+/-10 nm averaged) and 6.69 rad coverage, is the family selected earlier with the
jones_direct_fourier formulation (0.959, 6.74 rad; results/library_analysis_est_1075). The
spurious cross-polarisation of that formulation therefore did not affect the pillar selection.

## Why holes lose

The phase of a propagation-phase meta-atom spans k0 h (n_eff,max - n_eff,min). A pillar runs from
almost air (thin pillar) to almost film (wide pillar), a contrast of ~1.2 in effective index. A
hole cannot remove more than ~50 % of the film at p 0.5 um while keeping a 100 nm wall, and cannot
be narrower than h / 10, so its effective index stays between ~1.6 and ~2.2, a contrast of ~0.6
(Provisional, from the fill fractions of the selected families). Hence the film must be about
twice as thick for the same coverage. The thicker film raises the aspect ratio of the smallest
hole, which in turn removes the small holes that carried the high-index end of the range, so the
coverage stays below 2 pi within the fabrication limits. The selected hole families also have
lower and more strongly varying T0 (0.77 to 0.78 minimum), from Fabry-Perot and slab-guided
resonances of the thick continuous film, as Goldberg et al. note for their visible devices.

Thermal and damage consequences. Global heating of the metasurface scales with the absorbed
power (docs/campaign1_report.md, Section 5): at 200 W the hole family absorbs ~0.7 W against
~0.13 W for pillars (DERIVED, k = 6.7e-5 estimate). Lateral conduction in a 3 um nitride film is
negligible next to the substrate (conductance per unit width ~kappa_film x h ~ 3e-5 W/K for
kappa ~ 10 W/(m K), against ~1.4e-3 W/K for 1 mm of fused silica; DERIVED, kappa_film
`[VALUE?]`), so the continuous film does not offset the extra absorption. The largest per-atom
absorption (9.4e-3) also marks field concentration in the film at resonant hole sizes, relevant
for local damage.

## Limitations and remaining uses

This is the local-periodic approximation at normal incidence; supercell (edge) factors were not
computed for holes, and they penalise steep phase gradients in both families. Only circular holes
in the estimated Si-rich SiN were tested; square holes, LPCVD SiN, and multi-depth (grayscale)
voids, where the etch depth rather than the diameter sets the phase, were not. A multi-depth hole
family would remove the coverage limit but needs grayscale lithography with depth control at the
~10 nm level over an 18 mm aperture (Goldberg et al., arXiv:2603.27574, show multi-depth voids for
structural colour, without phase or efficiency data). The one advantage of holes that survives is
mechanical: a continuous film has no free-standing pillars to collapse or trap residues during
cleaning. That does not outweigh a ~10 % efficiency loss and a five-fold absorption increase for a
50 to 200 W element, so the pillar design remains the baseline.

## Verification notes

| Value | Location | Type | Pointer |
|---|---|---|---|
| eta_16, coverage, T0, fill fractions | Result table | SIM, review before submission | results/hole_vs_pillar/ranking.json |
| Absorbed fractions | Result table | SIM, review before submission | results/hole_vs_pillar/absorption.json |
| Convergence 197 vs 317 terms | Method | SIM | scratch check, values in configs/hole_vs_pillar.yaml comment |
| Pillar selection agreement 0.958 vs 0.959 | Method | SIM | results/library_analysis_est_1075/selection.json |
| Bare-silicon reflectance 32 % | Introduction | DERIVED | Fresnel, n = 3.56 (recalled, unverified) |
| UV Mie-void device parameters | Introduction | SOURCE (press release, review) | EurekAlert release; arXiv:2610.32494 |
| Absorbed power at 200 W, lateral conductance | Thermal | DERIVED | in text |

Open placeholders: `[VALUE?: nitride film thermal conductivity]`.
