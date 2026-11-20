# CLAUDE.md

Rules for working in this repository. The first two sections restate rules the repository
already follows (sources in brackets); the last two were added with the vendored skills.

## Provenance [README.md, "Provenance rules"; docs/model_notes.md]

- Every constant in the code carries its source in a comment or docstring. A value that was not
  retrieved and read is marked `[VALUE?: ...]`.
- Every number in a report has a provenance type: USER (supplied by the PI or the proposal
  text), SOURCE (retrieved and read, with a pointer), SIM (simulation executed in this
  repository, with its config) or DERIVED (calculation shown). Open items are marked
  `[VALUE?]`, `[CONFIRM?]` or `[REF?]`. Reports end with a provenance ledger
  (value | location | type | pointer).
- Every simulated number points to the config and the results directory that produced it.
  One YAML per run in `configs/`; each entry-point script writes
  `results/<run_id>/params.json` (config, git hash, package versions) through
  `metatrap.runio`.

## Checks [scripts/check_report_numbers.py; README.md, "Tests"]

- Every number in a LaTeX edition of a report must occur, as a value, in its Markdown source,
  or be listed with a reason in `scripts/report_numbers_allowlist.txt`.
  `python scripts/check_report_numbers.py` exits non-zero on an unmatched number or a stale
  allowlist entry.
- `python -m pytest` must pass. CI (`.github/workflows/ci.yml`) runs the tests and the
  report-number check on every push.

## Verification skills [.claude/skills/, vendored, see .claude/skills/UPSTREAM.md]

`planted-truth`, `acceptance-gate`, `independence-bookkeeping` and `timing-discipline` apply to
the simulations here with these meanings:

- **Digits** (acceptance-gate: count agreeing digits, rerun at raised precision): a result
  counts only to the extent it is stable when the discretisation is refined. For RCWA, more
  Fourier orders (e.g. 197 against 317 or 437 terms) and a finer permittivity grid; for FDTD,
  propagation and hologram rendering, a halved sampling step (grid, pixel, render or focal
  grid). Quote the change, not an adjective.
- **Reserved points** (acceptance-gate, independence-bookkeeping): held-out cases that entered
  no design, feedback round, fit or threshold choice: widths, wavelengths, sites, designs or
  fold windows kept back for the check. A metric measured by the same model that the feedback
  loop optimised is in-sample and is labelled so.
- **Independent route** (independence-bookkeeping): a second method, not a second wrapper of
  the same one, with the shared layers named in writing. Example: fmmax and grcwa are both
  RCWA and share the Fourier treatment of the permittivity, the staircased grid from
  `rcwa/cell.py` and the material values; FDTD and analytic limits do not share the Fourier
  layer (docs/rcwa_lineage.md). The angular-spectrum propagator and a direct
  Rayleigh-Sommerfeld sum share only the sampled field and the pixel convention
  (tests/test_planted_hologram.py).
- **Planted truth** (planted-truth): a design pipeline is first run on a target whose answer
  is known and recovered with a propagator that shares no code with the design loop.
- Thresholds are fixed before the check they apply to is run and are not loosened after a
  result. Every new check is broken once on purpose and seen to fail, and that is recorded. If
  a control fails on the existing pipeline, stop and report it as a finding; do not change a
  reported number silently.

## Compute planning [timing-discipline]

- Before any full-aperture run, run `scripts/bench_gpu.py` on a small aperture and extrapolate
  the time and the memory of the full run; write the projection, with its basis, in the run's
  config header or log before launching. Use `python scripts/bench_gpu.py --pilot CONFIG [--gpu]`
  (three scaled pilots, each in its own process; time and peak-memory projection with its basis,
  results/bench_pilot/). Its validation is open (docs/bench_pilot_validation.md): measured
  evaluation times were reproduced to +7 to +24 % on two 10^5-site configurations, but the
  fit-quality criterion was not met, so quote its projections with a +/-30 % allowance.
