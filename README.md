# metatrap: Phase A simulations of SiN meta-optics for cesium tweezer arrays

Simulation code for Phase A of the ANRF proposal ANRF/ARG/2026/009937
("integrated meta-optics and PIC control engine for neutral-atom qubits").
The first campaign targets element B1, the direct-projection SiN meta-hologram that
generates the static 1075 nm cesium array, and the high-power thermal behaviour of
SiN meta-optics (hypothesis H2a). The doublet meta-objective (B2) follows in a second
campaign and reuses the same library and propagation code.

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
| `src/metatrap/cavity/`, `docs/cavity_concept.md` | Cavity-enhanced large-period lattice with graded-phase flat-top modes for 10^5 to 10^6 sites: power budget, Fox-Li model, tolerances, defects |
| `docs/report/` | LaTeX editions of the cavity report: `cavity_report.tex` (edited, readable) and `cavity_concept_faithful.tex` (generated from the Markdown); build with `bash docs/report/build.sh` (needs TeX Live) |

## Provenance rules

Every constant in the code carries a source in a comment or docstring. Values that
were not retrieved and read are marked `[VALUE?: ...]`. Every simulated number in
the report points to the config and results directory that produced it.

## Tests

```bash
python -m pytest            # fast anchors
python -m pytest -m slow    # RCWA cross-checks and convergence
```
