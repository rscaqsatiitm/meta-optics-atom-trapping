# Lineage audit of the meta-atom solver checks

Question: which errors can the existing fmmax-against-grcwa cross-check rule out, and which can
it not? Written before the FDTD comparison was run (scripts/rcwa_fdtd_check.py, thresholds
declared in its docstring and committed first, commit 85ff125).

## The existing cross-check

`scripts/rcwa_convergence.py` (results/rcwa_convergence) compares fmmax (Fourier modal method in
JAX, the production solver: `rcwa/solver.py`, formulation `pol_fourier` for the recomputed
libraries) with grcwa (`solve_pillars_grcwa`). The reported agreement, 1e-4 rad at 437 terms,
is between grcwa and fmmax in its plain `fft` formulation.

## Shared layers

From the deepest up:

1. **Method.** Both are rigorous coupled-wave analysis (Fourier modal method): the
   permittivity of the patterned layer is expanded in plane waves, each layer is solved as an
   eigenproblem in that basis, and layers are joined by scattering matrices. Errors of the
   method itself (slow convergence of discontinuous permittivity, Gibbs ringing at pillar
   edges, the Fourier factorisation rule) can be common to both.
2. **Fourier factorisation.** grcwa uses the direct (Laurent) product rule for the permittivity.
   fmmax `fft` uses the same rule, which is why the two agree to 1e-4 rad: the agreement checks
   the wrapper and the phase convention, not the factorisation. The production formulation
   `pol_fourier` uses a different (vector) rule, and the `fft` formulation is itself still
   0.07 rad from its converged value at 317 terms (docs/campaign1_report.md, Section 3), so the
   production numbers are not covered by the grcwa agreement at all.
3. **Truncation.** Both truncate the plane-wave basis (fmmax circular, grcwa its own ordering)
   at a few hundred terms. A truncation error that converges slowly in both looks like
   agreement at equal term counts.
4. **Permittivity grid.** Both receive the same staircased permittivity map from
   `rcwa/cell.py` (`PillarCell.mask` / `permittivity`), sampled on 256 or 128 points per period.
   A wrong mask (orientation, centring, sub-pixel fill) is shared exactly. The grid error itself
   is 0.016 to 0.027 rad at fabricable widths (results/rcwa_grid_convergence.json).
5. **Material values.** Both use `materials.py` (the same n, k). An error there is shared
   exactly and no solver comparison can see it.
6. **Geometry and conventions.** Substrate incidence, exit into air, normal incidence, the
   zeroth-order extraction and the phase reference are set by the same wrapper code.
7. **Numerics.** numpy/JAX linear algebra (LAPACK eigensolvers). Errors here are unlikely to
   correlate with the physics and are not the feared failure.

Hence the fmmax-grcwa agreement certifies the wrapper (layer 6) against an independent RCWA
implementation, at the Laurent-rule level. It cannot detect a failure in layers 1 to 5.

## Checks that do not share the Fourier layer

- **Uniform-slab limit against the analytic thin-film result** (already present:
  tests/test_rcwa.py, `test_homogeneous_slab_is_airy_thin_film`, fmmax). A homogeneous layer has
  only the zeroth order, so this tests layer 6 (phase reference, normalisation, Fresnel
  factors) against an analytic formula, independent of the solver. It does not test layers 1
  to 4, because a uniform layer has no Fourier content.
- **Effective-medium limit** (new, check C4 of scripts/rcwa_fdtd_check.py). A lamellar grating
  with a period far below the wavelength (0.05 um at 1.075 um) must behave as a uniform
  anisotropic layer with the zeroth-order effective indices. This exercises the Fourier
  treatment of a discontinuous permittivity (layers 1 to 3) against an analytic limit, and it
  fixes which wrapper output is the field parallel to the ridges (a convention error of layer
  6 would otherwise be invisible). Shared with fmmax: the grid built in the check script and the
  index values.
- **FDTD** (new, checks C1 to C3 of scripts/rcwa_fdtd_check.py). The third-party `fdtd`
  package (version 0.3.5, pure numpy, Yee scheme with PML and periodic boundaries) solves
  Maxwell's equations in the time domain on a staggered grid. It shares no Fourier expansion,
  no eigenproblem and no scattering matrix with RCWA (layers 1 to 3 and 6 are not shared). It
  shares the problem definition (the index values 1.45 / 2.0 / 1.0 and the geometry, both
  written into the check script, not taken from `materials.py` or `rcwa/cell.py`) and numpy.
  Its own planted truth is the uniform layer against the analytic Airy result (C1), and its own
  sampling check is 10 nm against 5 nm (C2).
- Not added: a second FDTD or finite-element code (meep is not installable with pip in this
  environment), and a pillar (2D-periodic) FDTD case, which is a 3D simulation beyond the
  budget of this container; the check is on a 1D-periodic grating with the field parallel to
  the ridges.

## Results (results/rcwa_fdtd_check/metrics.json, metrics_break.json)

| Check | Threshold (declared first) | Result | Verdict |
|---|---|---|---|
| C1 FDTD uniform layer against Airy, 10 nm / 5 nm grid | amp 0.01, phase 0.03 rad | amp 3.2e-4 / 7.8e-5, phase 0.0051 / 0.0013 rad | PASS |
| C2 FDTD 10 nm against 5 nm, grating ratios (w 0.15 / 0.35 um) | amp 0.01, phase 0.02 rad | amp 0.0013 / 0.0023, phase 0.057 / 0.040 rad | FAIL |
| C3 fmmax (pol_fourier, 197 terms) against FDTD 5 nm (w 0.15 / 0.35 um) | amp 0.02, phase 0.05 rad | amp 6.6e-4 / 1.2e-3, phase 0.0027 / 0.0029 rad | within thresholds, but INCONCLUSIVE (C2 failed) |
| C4 effective-medium limit, p 0.05 um | phase 0.10 rad | TE 0.0074 rad, TM 0.016 rad; E parallel to the ridges is wrapper output t_y | PASS |

Verdict: OPEN. The FDTD route reproduces the analytic layer (C1), and fmmax agrees with the
5 nm FDTD grating results to 0.003 rad and 0.1 % in amplitude (C3). However, the FDTD grating
result itself changes by 0.04 to 0.06 rad between the 10 and 5 nm grids (C2), more than the
declared 0.02 rad, so the 5 nm FDTD is not shown to be converged to the level at which C3
agrees, and by the rule declared beforehand C3 does not certify the RCWA. The uniform layer
converges much faster than the grating (C1: 0.005 to 0.001 rad), which points at the
ridge-edge discretisation of the FDTD (field discontinuities at the vertical air-ridge
interfaces, not present in a uniform layer) as the slow part. Next step, not run: a 2.5 nm
FDTD grid (about 8 times the 5 nm cost, ~3 h on 4 cores, so a restructure is preferable, e.g.
sub-cell permittivity averaging at the ridge edges, or Richardson extrapolation from three
grids with its own threshold declared first).

Deliberate breaks (all four checks must fail; recorded in metrics_break.json): C1 with n = 2.05
in the analytic layer, phase error 0.28 rad, FAIL; C2 comparing different widths, FAIL; C3 with
the polarisation perpendicular to the ridges, 0.59 rad at w = 0.35 um, FAIL (at w = 0.15 um the
two polarisations differ by only 0.047 rad, inside the threshold, so the check relies on the
wider ridge to see a polarisation error); C4 with ridge index 2.1, 0.37 rad, FAIL. The first
version of the C3 and C4 breaks could not fail (C4 infers the orientation from the data, which
also undid the C3 swap); this was found by the break run and the breaks were replaced (see the
script docstring). C4 therefore tests the effective indices, and reports the orientation
(t_y), rather than testing an orientation fixed in advance.

## Superseding check with converged references (2026-11-20; scripts/rcwa_fdtd_convergence.py)

The results above are not evidence either way: the FDTD set-up of scripts/rcwa_fdtd_check.py had
a period of p - d (fdtd.PeriodicBoundary identifies the first and last cell) and, on the 10 nm
grid, ridges one cell too narrow. With the period exact, every width and layer on whole cells
(grids 25 / 12.5 / 6.25 nm), both polarisations, and RCWA at 437 terms on a 400-point grid
(thresholds declared in the script and committed before the run; results/rcwa_fdtd_convergence):

| Check | Threshold | Result | Verdict |
|---|---|---|---|
| C1 FDTD uniform layer against Airy, 25 / 12.5 / 6.25 nm | 0.03 rad, 0.01 amp | 0.0320 / 0.0079 / 0.0020 rad; 2.1e-3 / 4.9e-4 / 1.2e-4 amp | FAIL (25 nm grid, by 0.002 rad) |
| C2 FDTD sampling, Richardson over the three grids | order 0.7 to 3.5; correction 0.02 rad, 0.01 amp | TE converged (finest two grids within 1e-3); TM order 2.00 / 2.07, correction 5e-4 / 4e-4 rad | PASS |
| C5 RCWA digits: 317 vs 437 terms, grid 200 vs 400 | 0.005 rad, 0.002 amp | <= 4.8e-4 rad (TM, terms), <= 4e-5 rad (grid) | PASS |
| C3 fmmax (437 terms, grid 400) against extrapolated FDTD, TE and TM | 0.05 rad, 0.02 amp | <= 3e-4 rad, <= 2e-4 amp | PASS, formally INCONCLUSIVE (C1) |

Deliberate breaks (metrics_break.json): all four checks fail (C1 with n = 2.05: 0.26 rad; C2 with
the coarsest ratio of another width: fails on TM, while the TE entries are skipped as converged,
so the C2 break does not exercise the TE order test; C5 with 41 terms; C3 with RCWA TM against
FDTD TE).

Verdict: OPEN by the rule declared beforehand, because C1 fails on the coarsest grid. The failure
is understood: the uniform-layer error falls exactly as d^2 (ratios 4.05 and 3.95), the grid
dispersion of the Yee scheme, and the two finer grids pass with margin; the 25 nm grid enters
only the Richardson estimate, whose observed order (2.0) is the scheme's. Within that
qualification, fmmax (pol_fourier) and an FDTD route that shares no Fourier expansion,
eigenproblem or scattering matrix agree on the lamellar grating to 3e-4 rad and 2e-4 in
amplitude for both polarisations, with each side converged to better than 5e-4 rad. The earlier
C2 failure (0.04 to 0.06 rad between 10 and 5 nm) is explained by the two set-up defects: on
aligned grids the change between the two finer grids is <= 0.0015 rad. Not covered: 2D-periodic
pillars (a 3D FDTD case) and the material values, which both routes share.

## Finer grid sequence (2026-11-20; `scripts/rcwa_fdtd_convergence.py --fine`)

Declared after the run above and before this one: the same checks and thresholds on the three
finest grids, 12.5 / 6.25 / 3.125 nm (3.125 nm runs made four at a time in parallel, ~17 min
each; results/rcwa_fdtd_convergence/metrics_fine.json).

| Check | Threshold | Result | Verdict |
|---|---|---|---|
| C1 FDTD uniform layer against Airy, 12.5 / 6.25 / 3.125 nm | 0.03 rad, 0.01 amp | 0.0079 / 0.0020 / 0.0005 rad; 4.9e-4 / 1.2e-4 / 3e-5 amp | PASS |
| C2 FDTD sampling | finest two grids within 1e-3, else Richardson order and correction | every case within 1e-3 (max 3.9e-4, TM w 0.15 um); no extrapolation applied | PASS |
| C5 RCWA digits (as before) | 0.005 rad, 0.002 amp | <= 4.8e-4 rad | PASS |
| C3 fmmax (437 terms, grid 400) against FDTD at 3.125 nm, TE and TM | 0.05 rad, 0.02 amp | <= 2.1e-4 rad, <= 6e-5 amp | PASS, conclusive |

Deliberate breaks (metrics_fine_break.json): C1, C5 and C3 fail as required. The declared C2 break
(coarsest-grid ratio of another width) did not fire: in the converged branch C2 does not read the
coarsest grid, so that corruption cannot reach it. A control on the middle grid, the clause C2
actually used, was added after this result (`--break-c2-mid`, metrics_fine_break_c2_mid.json)
and makes C2 fail. Both the non-firing declared break and the added control are on record.

Verdict: CLOSED for the stated scope. fmmax with the pol_fourier formulation (the solver and
formulation of the libraries) agrees with an FDTD route that shares no Fourier expansion,
eigenproblem or scattering matrix to 2.1e-4 rad and 6e-5 in amplitude, for both polarisations,
on a lamellar grating (substrate 1.45, ridges 2.0, period 0.5 um, height 1.0 um, 1.075 um,
widths 0.15 to 0.35 um), with each route converged to better than 5e-4 rad. Not covered:
2D-periodic pillars (a 3D FDTD case), oblique incidence, other indices and heights, and the
material values, which both routes share.
