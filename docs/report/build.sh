#!/usr/bin/env bash
# Build both LaTeX editions of the cavity-lattice report (needs a TeX Live with latexmk).
#   bash docs/report/build.sh
set -euo pipefail
cd "$(dirname "$0")/../.."
python scripts/cavity_report_tables.py latex docs/report/tables   # tables from the results JSON
python scripts/report_faithful.py                                  # faithful body from the Markdown
cd docs/report
for f in cavity_report cavity_concept_faithful; do
  latexmk -pdf -interaction=nonstopmode -halt-on-error "$f.tex"
done
latexmk -c
