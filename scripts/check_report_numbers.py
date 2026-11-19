"""List numbers in the LaTeX editions that do not occur in docs/cavity_concept.md.
Usage:
    python scripts/check_report_numbers.py

Exit status 1 if any number of the LaTeX edition is neither in the Markdown source nor in the
explicit allowlist (scripts/report_numbers_allowlist.txt: one entry per line,
"<file> | <number as printed> | <reason>"), and also if an allowlist entry is no longer
needed (so the list cannot grow stale). tests/test_report_numbers.py checks that a planted
mismatch fails. Run in CI (.github/workflows/ci.yml).

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
    text = re.sub(r"(\d(?:\.\d+)?) x 10\^(-?\d+)", r"\1e\2", text)  # before the bare 10^ rule, which would split it
    text = re.sub(r"(?<![\d.])10\^\{?(-?\d+)\}?", r"1e\1", text)
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


def load_allowlist(path):
    """{(file, normalised value): reason} from '<file> | <number> | <reason>' lines ('#' comments)."""
    out = {}
    if not Path(path).exists():
        return out
    for line in Path(path).read_text().splitlines():
        line = line.split("#", 1)[0].strip() if not line.lstrip().startswith("#") else ""
        if not line:
            continue
        name, num, reason = (x.strip() for x in line.split("|", 2))
        if not reason:
            raise ValueError(f"allowlist entry without a reason: {line}")
        out[(name, f"{float(num):.10g}")] = reason
    return out


def check(md_path, tex_paths, allow_path):
    """Return the number of failures (unmatched numbers plus stale allowlist entries)."""
    md = values(Path(md_path).read_text())
    allow = load_allowlist(allow_path)
    used, bad = set(), 0
    for tex_path in tex_paths:
        name = Path(tex_path).name
        tex = Path(tex_path).read_text()
        missing = {k: v for k, v in values(strip_layout(tex)).items() if k not in md}
        allowed = {k: v for k, v in missing.items() if (name, k) in allow}
        unmatched = {k: v for k, v in missing.items() if (name, k) not in allow}
        used |= {(name, k) for k in allowed}
        print(f"{name}: {len(unmatched)} numbers not in {Path(md_path).name} and not allowlisted "
              f"({len(allowed)} allowlisted)")
        for k, v in sorted(allowed.items(), key=lambda kv: float(kv[0])):
            print(f"  allowlisted {v:>10}  ({allow[(name, k)]})")
        for k, v in sorted(unmatched.items(), key=lambda kv: float(kv[0])):
            i = tex.find(v)
            ctx = " ".join(tex[max(0, i - 50): i + 30].split())
            print(f"  UNMATCHED {v:>10}  ...{ctx}...")
        bad += len(unmatched)
    stale = [key for key in allow if key not in used and key[0] in {Path(t).name for t in tex_paths}]
    for name, k in stale:
        print(f"  STALE allowlist entry: {name} | {k} ({allow[(name, k)]}) is no longer needed; remove it")
    bad += len(stale)
    print("report numbers: " + ("OK" if bad == 0 else f"FAIL ({bad})"))
    return bad


def main():
    bad = check(ROOT / "docs" / "cavity_concept.md", [ROOT / "docs" / "report" / "cavity_report.tex"],
                ROOT / "scripts" / "report_numbers_allowlist.txt")
    return 0 if bad == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
