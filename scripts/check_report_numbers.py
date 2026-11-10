"""Check that every number in the edited LaTeX report of the cavity study also occurs in its Markdown source.
Usage:
    python scripts/check_report_numbers.py

Scope: docs/report/cavity_report.tex (the hand-edited cavity-lattice report) against
docs/cavity_concept.md (the working report it was written from). The faithful LaTeX edition is
generated from the Markdown by scripts/report_faithful.py and needs no check; the other reports of the
repository are not covered.

Numbers are compared as signed values (1.4e-8, 1.4\\times10^{-8}, \\SI{1.4e-8}{...} and $1.4 x 10^-8$ are
equal; a dropped minus sign is a mismatch). Layout parameters (column widths, figure widths, font
sizes) are stripped first. A number in the LaTeX that is absent from the Markdown is either a
deliberate restatement, listed in ALLOWED with its reason, or a transcription error to fix.
Exit status: 0 if every flagged number is in ALLOWED, 1 otherwise.
"""

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
NUM = re.compile(r"(?<![\w.])-?\d+(?:\.\d+)?(?:(?:e|E)[+-]?\d+)?(?![\w.]*[A-Za-z]{2,}\d)")

# Deliberate restatements in cavity_report.tex: value (as keyed by values()) -> reason.
ALLOWED = {
    "-1": "exponent of (a, b)^{-1/2} in Section 6.1; the Markdown writes (a, b)^-0.5",
    "25.4": "the 1 inch mirror diameter of the thermal model, written in millimetres",
}


def normalise(text):
    text = re.sub(r"(?<=\d),(?=\d{3}\b)", "", text)  # 1,177 -> 1177
    text = re.sub(r"(\d(?:\.\d+)?)\s*\\times\s*10\^\{?(-?\d+)\}?", r"\1e\2", text)
    text = re.sub(r"(?<![\d.])10\^\{?(-?\d+)\}?", r"1e\1", text)
    return re.sub(r"(\d(?:\.\d+)?) x 10\^(-?\d+)", r"\1e\2", text)


def values(text):
    """Map value key (signed, 10 significant digits) -> (token, context of its first occurrence)."""
    text = normalise(text)
    out = {}
    for m in NUM.finditer(text):
        try:
            v = float(m.group(0))
        except ValueError:
            continue
        key = f"{v:.10g}"
        if key not in out:
            ctx = " ".join(text[max(0, m.start() - 50): m.end() + 30].split())
            out[key] = (m.group(0), ctx)
    return out


def strip_layout(tex):
    tex = re.sub(r"(?m)^\s*%.*$", "", tex)
    tex = re.sub(r"\\begin\{tabularx?\}(\{[^}]*\})?\{.*$", "", tex, flags=re.M)
    tex = re.sub(r"\\includegraphics\[[^\]]*\]", "", tex)
    tex = re.sub(r"\\(label|ref|eqref|cite|input|bibliography\w*|documentclass)\{[^}]*\}", "", tex)
    tex = re.sub(r"\\section|\\subsection|\[\d+ex\]|\\\\\[\d+ex\]", "", tex)
    tex = re.sub(r"(?:\d+(?:\.\d+)?)\\(?:linewidth|columnwidth|tabcolsep)", "", tex)
    tex = re.sub(r"\\real\{[\d.]+\}", "", tex)
    return tex


def main(tex_path=None):
    md = values((ROOT / "docs" / "cavity_concept.md").read_text())
    path = Path(tex_path) if tex_path else ROOT / "docs" / "report" / "cavity_report.tex"
    missing = {k: v for k, v in values(strip_layout(path.read_text())).items() if k not in md}
    unexplained = {k: v for k, v in missing.items() if k not in ALLOWED}
    print(f"{path.name}: {len(missing)} numbers not in docs/cavity_concept.md, {len(unexplained)} unexplained")
    for k, (tok, ctx) in sorted(missing.items(), key=lambda kv: float(kv[0])):
        tag = f"allowed: {ALLOWED[k]}" if k in ALLOWED else "UNEXPLAINED"
        print(f"  {tok:>10}  [{tag}]  ...{ctx}...")
    return 0 if not unexplained else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else None))
