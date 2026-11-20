# metatrap: Phase A simulations of SiN meta-optics for cesium tweezer arrays

Simulation code for Phase A of the ANRF proposal ANRF/ARG/2026/009937
("integrated meta-optics and PIC control engine for neutral-atom qubits").
The first campaign targets element B1, the direct-projection SiN meta-hologram that
generates the static 1075 nm cesium array, and the high-power thermal behaviour of
SiN meta-optics (hypothesis H2a). Later work, on the same library and propagation code, added
the doublet meta-objective (B2), the pre-fabrication checks, full-aperture designs with 2,304,
20,164 and 99,856 sites, the 10^6-site scaling study (direct hologram and lenslet-array relay),
and a comparison of hole and pillar meta-atoms (see the docs rows of the layout table).

## Install

```bash
pip install -e ".[test]"
```

Everything is open-source Python and runs on a CPU: RCWA via
[fmmax](https://github.com/facebookresearch/fmmax) (JAX) with
[grcwa](https://github.com/weiliangjinca/grcwa) as an independent cross-check,
numpy/scipy for propagation and the thermal solver.

## Layout

| Path | Content |
|---|---|
| `src/metatrap/materials.py` | Dispersion of SiN (LPCVD, Si-rich), fused silica, sapphire; CSV loader for measured films |
| `src/metatrap/atoms/` | Cs/Rb polarisability (Steck data), trap depth, array power budget |
| `src/metatrap/rcwa/` | Pillar unit cell, RCWA solvers, meta-atom library sweep and lookup, supercell deflection efficiency |
| `src/metatrap/propagation/` | Angular-spectrum propagation, pupil mapping, glass-cell plate, illumination |
| `src/metatrap/hologram/` | Site targets, weighted Gerchberg-Saxton, rendering to pillar widths, spot metrics |
| `src/metatrap/thermal/` | Analytic and axisymmetric finite-difference heating, thermo-optic focal shift |
| `configs/` | One YAML per run |
| `scripts/` | Entry points; each writes `results/<run_id>/params.json` (config, git hash, versions) |
| `docs/model_notes.md` | Equations, assumptions and the provenance ledger |
| `docs/campaign1_report.md` | Results of the first campaign |
| `docs/campaign1b_report.md`, `docs/prefab_checks_report.md` | Literature film estimates and their consequences; pre-fabrication checks and damage model |
| `docs/metalens_report.md`, `docs/objective_comparison.md` | B2 doublet meta-objective; comparison with refractive objectives |
| `docs/scaling_100k.md`, `docs/million_sites.md`, `docs/gpu_run.md` | 10^5-site design, 10^6-site routes, running the 10^6 design on a GPU |
| `docs/hole_vs_pillar.md` | Hole (inverse) against pillar meta-atoms at 1075 nm |
| `src/metatrap/cavity/`, `docs/cavity_concept.md` | Cavity-enhanced large-period lattice with graded-phase flat-top modes for 10^5 to 10^6 sites: power budget, Fox-Li model, tolerances, defects |
| `docs/report/` | LaTeX editions of the cavity report: `cavity_report.tex` (edited, readable) and `cavity_concept_faithful.tex` (generated from the Markdown); build with `bash docs/report/build.sh` (needs TeX Live) |

## Provenance rules

Every constant in the code carries a source in a comment or docstring. Values that
were not retrieved and read are marked `[VALUE?: ...]`. Every simulated number in
the report points to the config and results directory that produced it.

## Tests

```bash
python -m pytest            # all tests (about 1 min on 4 cores); tests/test_lazy_gpu.py runs only with a CUDA GPU
```

No test is currently marked `slow` (the marker is defined in pyproject.toml for future use). The
RCWA convergence and solver cross-checks are scripts, not tests: `scripts/rcwa_grid_convergence.py`,
`scripts/compare_library_formulations.py` and the convergence configs.

```bash
python -m pytest tests/test_rcwa.py   # RCWA conventions against thin-film results
```
