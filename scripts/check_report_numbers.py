"""List numbers in the LaTeX editions that do not occur in docs/cavity_concept.md.
Usage:
    python scripts/check_report_numbers.py

Numbers are compared as values (1.4e-8, 1.4\\times10^{-8}, \\SI{1.4e-8}{...} and $1.4 x 10^-8$ are
equal). Layout parameters (column widths, figure widths, font sizes) are stripped first. Every
reported number must be a deliberate restatement (e.g. a ratio derived in the text) or a
transcription error to fix.
"""

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
NUM = re.compile(r"(?<![\w.])-?\d+(?:\.\d+)?(?:(?:e|E)[+-]?\d+)?(?![\w.]*[A-Za-z]{2,}\d)")


def values(text):
    text = re.sub(r"(?<=\d),(?=\d{3}\b)", "", text)  # 1,177 -> 1177
    text = re.sub(r"(\d(?:\.\d+)?)\s*\\times\s*10\^\{?(-?\d+)\}?", r"\1e\2", text)
    text = re.sub(r"(?<![\d.])10\^\{?(-?\d+)\}?", r"1e\1", text)
    text = re.sub(r"(\d(?:\.\d+)?) x 10\^(-?\d+)", r"\1e\2", text)
    out = {}
    for m in NUM.finditer(text):
        try:
            v = float(m.group(0))
        except ValueError:
            continue
        out.setdefault(f"{v:.10g}", m.group(0))  # signed: a dropped minus sign is a mismatch
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


def main():
    md = values((ROOT / "docs" / "cavity_concept.md").read_text())
    bad = 0
    for name in ("cavity_report.tex",):
        tex = (ROOT / "docs" / "report" / name).read_text()
        missing = {k: v for k, v in values(strip_layout(tex)).items() if k not in md}
        print(f"{name}: {len(missing)} numbers not in docs/cavity_concept.md")
        for k, v in sorted(missing.items(), key=lambda kv: float(kv[0])):
            i = tex.find(v)
            ctx = " ".join(tex[max(0, i - 50): i + 30].split())
            print(f"  {v:>10}  ...{ctx}...")
        bad += len(missing)
    return 0 if bad == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
