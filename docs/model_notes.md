# Model notes and provenance ledger

This file records the equations, assumptions and parameter provenance of the Phase A
simulation code. Provenance types: USER (supplied by the PI or the proposal text),
SOURCE (retrieved and read on 2026-10-31, with pointer), SIM (simulation executed in
this repository, with config), DERIVED (calculation shown here), and placeholders
`[VALUE?]`, `[CONFIRM?]`, `[REF?]` for items still to be verified.

## 1. Atomic model

Scalar ground-state polarisability of Cs (and Rb for the ref. 4 check):

    alpha(w) = (2/3) sum_{D1,D2} |d_k|^2 w_k / (w_k^2 - w^2) + alpha_rest,

with d_k the reduced matrix elements in Steck's convention and alpha_rest fixed so
that alpha(0) equals Steck's measured static polarisability. Trap depth
U/k_B = (2 pi a0^3 / c k_B) alpha I. Code: `src/metatrap/atoms/alkali.py`.

Array power budget (proposal Section 4 (ii)): I0 = pi NA_eff^2 P_site / lambda^2
(paraxial peak of a uniformly filled focus), P_site = utilisation x P_laser / N.
Code: `src/metatrap/atoms/power_budget.py`.

Check against ref. 4 (Rb, 850 nm, 18,225 sites, 33 W, 12.2 W effective): with
NA_eff = 0.61 x 0.85 um / 1.11 um = 0.467 (DERIVED from the measured first dark
radius) the model gives 0.306 mK against the measured 0.305(17) mK.

## 2. Meta-atom library

RCWA with fmmax 1.7.1, Jones-direct-Fourier formulation (vector FMM), square pillars on
fused silica, light incident from the substrate, exit into air, normal incidence.
Convergence (results/rcwa_convergence, SIM): relative phase converged to 0.003 rad at
197 Fourier terms; permittivity-grid dependence 0.01 to 0.03 rad between 256 and 1024
points (largest at gaps below 60 nm). Library numerics: 197 terms, 256 grid. The
plain-FFT (Laurent) formulation converges slowly (0.07 rad error at 317 terms); grcwa,
which uses it, agrees with fmmax-FFT to 1e-4 rad at 437 terms, an independent check of
the wrapper and phase conventions.

Numerical validity filter: entries with cross-polarised amplitude above 0.05 (forbidden
by symmetry; found only at narrow gaps and resonances) are excluded from the lookup.

Fabrication limits (default, `[CONFIRM?: CNNP etch capability]`): minimum feature
100 nm, minimum gap 100 nm, aspect ratio h/w <= 10.

## 3. Hologram design

Pupil-domain WGS with two kernels:

* Debye: K_ij = exp(-i k0 sigma_i . x_j) (FFT-equivalent, standard).
* Exact: K_ij = exp(-i k0 [OPL(r_i -> x_j) - OPL0(r_i)]) through the layer stack.

Here, sigma_i are pupil direction cosines, r_i their metasurface positions under the
ray mapping r(s) = sum_i d_i s / sqrt(n_i^2 - s^2), and OPL the optical path
sum_i n_i^2 d_i / sqrt(n_i^2 - s^2) along the actual ray. Mode efficiency of site j:
|sum_i A_i^2 e^{i psi_i} conj(K_ij)|^2 / (sum_i A_i^2)^2.

Validation: `scripts/validate_fullpixel.py` propagates the rendered 5 mm metasurface
(10^8 pixels at 0.5 um) with the exact scalar angular-spectrum method.

Limitations: scalar optics throughout (no vector focusing at NA 0.55; the depth of a
linearly polarised focus at this NA differs from the scalar value by a few per cent and
the spots become slightly elliptical [VALUE?: quantify with a vector Debye-Wolf model]);
local-periodic approximation for the metasurface (supercell check at the rim).

## 4. Thermal model

Axisymmetric steady-state finite volume, film absorption as a surface source,
optional bulk absorption, convective faces (h = 10 W/m^2/K), rim held at ambient.
Path change per unit thickness chi = dn/dT + (n - 1)(1 + nu) alpha (free plate,
photo-elastic terms neglected). Validated against the thin-plate disc result
P/(4 pi kappa t)[1 + 2 ln(R/a)] (2 %) and the half-space Gaussian result
P/(sqrt(2 pi) kappa w) (5 %).

## 5. Provenance ledger (constants)

| Value | Location | Type | Pointer |
|---|---|---|---|
| Cs D1 894.592 959 86 nm, D2 852.347 275 82 nm | constants.py | SOURCE | Steck, Cs D Line Data rev. 2.3.4, Tables 3, 4 |
| Cs d(D1) = 3.1869 e a0, d(D2) = 4.4837 e a0 | constants.py | SOURCE | Steck Cs Table 7 |
| Cs alpha0/h = 0.099 77 Hz/(V/cm)^2 (401.0 a.u.) | constants.py | SOURCE | Steck Cs Table 6 |
| Cs mass 132.905 451 931 u | constants.py | SOURCE | Steck Cs Table 1 |
| Rb-87 D1 794.978 851 156 nm, D2 780.241 209 686 nm, d = 2.9931, 4.227 52 e a0, alpha0/h = 0.0794 | constants.py | SOURCE | Steck, Rb-87 D Line Data rev. 2.3.4 |
| SiN (LPCVD) Sellmeier | materials.py | SOURCE | Luke et al., Opt. Lett. 40, 4823 (2015), DOI 10.1364/OL.40.004823, via refractiveindex.info DB |
| Fused silica Sellmeier | materials.py | SOURCE | Malitson, JOSA 55, 1205 (1965), DOI 10.1364/JOSA.55.001205 |
| Fused silica n(1064.00 nm) = 1.449633 (test) | test_materials.py | SOURCE | Corning 7980 product sheet (uqgoptics copy) |
| Sapphire (o) Sellmeier | materials.py | SOURCE | Malitson and Dodge (1972); Dodge, CRC Handbook IV (1986) |
| Si-rich SiN n = 2.3 | materials.py | USER | proposal Section 4 quoting ref. 14; `[VALUE?: verify, dispersion]` |
| Pyrex n = 1.47 | materials.py | `[VALUE?]` | used only for the ref. 4 geometry |
| Fused silica kappa 1.38 W/m/K | axisym_fd.py | SOURCE (secondary) | Corning HPFS figure quoted in search; sheet read gives 1.30 |
| Fused silica dn/dT 9.6 ppm/K (1064 nm), alpha 0.52 ppm/K | axisym_fd.py | SOURCE | Corning 7980 sheet |
| Sapphire kappa 33.0 W/m/K (perp. c), alpha 5.0 ppm/K | axisym_fd.py | SOURCE | Alkor sapphire data sheet |
| Sapphire dn/dT 13 ppm/K | axisym_fd.py | `[VALUE?]` | vendor value at 0.546 um, verify at 1075 nm |
| Poisson ratios 0.16 (FS), 0.25 (sapphire) | axisym_fd.py | `[VALUE?]` | from memory |
| Si-rich SiN dn/dT up to 1.65e-4 /K | report | SOURCE | Nejadriahi et al., Opt. Express (2020), DOI 10.1364/OE.396969 (abstract) |
| Ref. 4: 850 nm, 18,225 sites (135 x 135), 4.30(7) um pitch, 33 W in, 12.2 W effective, 15.0 W predicted, 67.5 % diffraction efficiency, 13.2 % rms power non-uniformity, dark radius 1.11(6) um, 0.305(17) mK, 3.5 mm Pyrex + 2.5 mm air + 9 mm vacuum, 19.8 mm, NA 0.58 | power_budget.py, configs | SOURCE | arXiv:2606.02715 (HTML), retrieved 2026-10-31 |
| Cell: 2.5 mm air, 4 mm fused silica, 9 mm vacuum | configs | `[CONFIRM?]` | default modelled on ref. 4 |
| Metasurface substrate thickness 0.5 mm, 25 mm die, rim heat-sunk | configs | `[CONFIRM?]` | proposal Section 4 (i) gives the 25 mm die |
