"""Faithful LaTeX edition of docs/cavity_concept.md (docs/report/cavity_concept_faithful_body.tex).
Usage:
    python scripts/report_faithful.py

The Markdown is converted by pandoc after a notation pass that changes only typography, never
content: file paths, result keys and URLs become code; the report's ASCII notation (lambda, theta,
urad, um, rtHz, 10^6, 1.4e-8, +-, <=, 101 x 101, w_z, ...) becomes LaTeX math. Placeholders such
as [VALUE?: ...] are kept verbatim and highlighted. docs/report/cavity_concept_faithful.tex wraps the
body with the shared preamble.
"""

import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "docs" / "cavity_concept.md"
OUT = ROOT / "docs" / "report" / "cavity_concept_faithful_body.tex"

GREEK = {"lambda": r"\lambda", "theta": r"\theta", "kappa": r"\kappa", "phi": r"\phi", "psi": r"\psi", "alpha": r"\alpha",
         "sigma": r"\sigma", "zeta": r"\zeta", "xi": r"\xi", "eta": r"\eta", "nu": r"\nu", "pi": r"\pi", "delta": r"\delta",
         "dphi": r"\mathrm{d}\phi", "dtheta": r"\mathrm{d}\theta", "dnu": r"\mathrm{d}\nu", "eps": r"\epsilon"}
SUBSCRIPTED = {"w_z": "w_z", "theta_z": r"\theta_z", "a_z": "a_z", "G_max": r"G_\mathrm{max}", "P_circ": r"P_\mathrm{circ}",
               "L_rt": r"L_\mathrm{rt}", "S_eps": r"S_\epsilon", "S_x": "S_x", "S_nu": r"S_\nu", "eta_return": r"\eta_\mathrm{return}",
               "delta_0": r"\delta_0", "dnu_r": r"\mathrm{d}\nu_r", "phi_G": r"\phi_G", "I_A": "I_A", "I_B": "I_B", "W_eq": r"W_\mathrm{eq}",
               "H_eq": r"H_\mathrm{eq}", "f_z": "f_z", "T_I": "T_I", "T_x": "T_x", "T_in": r"T_\mathrm{in}", "T_out": r"T_\mathrm{out}",
               "K_in": r"K_\mathrm{in}", "K_out": r"K_\mathrm{out}", "E_A": "E_A", "E_B": "E_B", "E_in": r"E_\mathrm{in}",
               "Ex_A": "E_{x,A}", "Ex_B": "E_{x,B}", "Ez_A": "E_{z,A}", "Ez_B": "E_{z,B}", "zeta_A": r"\zeta_A", "zeta_B": r"\zeta_B",
               "xi_A": r"\xi_A", "xi_B": r"\xi_B", "W_x": "W_x", "f1z": "f_{1z}", "f2z": "f_{2z}"}
UNITS = {"urad": r"\upmu\mathrm{rad}", "um": r"\upmu\mathrm{m}", "uK": r"\upmu\mathrm{K}", "rtHz": r"\sqrt{\mathrm{Hz}}",
         "nrad": r"\mathrm{nrad}"}
PROTECT = re.compile(
    r"`[^`]*`"                                         # existing code spans
    r"|\[(?:VALUE|REF|CONFIRM|BUDGET)\?[^\]]*\]"       # placeholders
    r"|!?\[[^\]]*\]\([^)]*\)"                          # links and images
    r"|https?://\S+|www\.\S+|\S+\.(?:gov|org|edu)/\S*"  # URLs
    r"|10\.\d{4,}/[^\s,;)]+"                           # DOIs
    r"|arXiv:\S+?(?=[\s,;)]|$)")
PATHLIKE = re.compile(r"^[(\[]*((?:[\w.*-]+/)+[\w.*-]*|[\w.*-]*_[\w.*-]*)([).,;:\]]*)$")


def codeify(token):
    """Wrap file paths and identifiers with underscores (not in SUBSCRIPTED) in backticks."""
    m = PATHLIKE.match(token)
    if not m:
        return token
    core = m.group(1)
    if core in SUBSCRIPTED or re.fullmatch(r"[\w]+/[\d.]+", core) or not re.search(r"[a-z]", core):
        return token
    if "/" in core and not re.search(r"\.(py|json|yaml|md|h5|png|npz|pdf|tex)$|^(results|scripts|configs|src|docs|data|tests)/|/$", core):
        return token  # a ratio such as W/cm^2 or lambda/1000, not a path
    lead = token[: token.index(core)]
    return f"{lead}`{core}`{m.group(2)}"


def mathify(text):
    """Typeset the ASCII notation of one unprotected text segment (may span lines)."""
    text = "".join(codeify(w) if not w.isspace() else w for w in re.split(r"(\s+)", text))
    parts = re.split(r"(`[^`]*`)", text)
    out = []
    for i, seg in enumerate(parts):
        if i % 2:
            out.append(seg)
            continue
        for k, v in SUBSCRIPTED.items():
            seg = re.sub(rf"(?<![\w$\\]){re.escape(k)}(?![\w])", lambda m, v=v: f"${v}$", seg)
        seg = re.sub(r"(\d(?:\.\d+)?)e([+-]?)0*(\d+)\b", lambda m: f"${m.group(1)}\\times10^{{{'-' if m.group(2) == '-' else ''}{m.group(3)}}}$", seg)
        seg = re.sub(r"\b10\^(-?\d+)", r"$10^{\1}$", seg)
        seg = re.sub(r"(?<=[\w)|$])\^(-?\d+)", r"$^{\1}$", seg)
        seg = re.sub(r"(?<=\d) x (?=\d)", r" $\\times$ ", seg)
        seg = re.sub(r"(?<=\d)x\b", r"$\\times$", seg)
        seg = re.sub(r"\+-(\d[\d.]*)", r"$\\pm\1$", seg)
        seg = re.sub(r"~(\d[\d.]*)", r"$\\sim\1$", seg)
        seg = seg.replace("+-", r"$\pm$").replace(">=", r"$\geq$").replace("<=", r"$\leq$").replace(" ~ ", r" $\sim$ ")
        for k, v in UNITS.items():
            seg = re.sub(rf"(?<![\w$\\]){k}\b", lambda m, v=v: f"${v}$", seg)
        for k, v in GREEK.items():
            seg = re.sub(rf"(?<![\w$\\]){k}(?![\w])", lambda m, v=v: f"${v}$", seg)
        seg = seg.replace("~$", r"$\sim$$").replace("$$", "")
        out.append(seg)
    return "".join(out)


def protect(tok):
    """Form of a protected token in the Markdown passed to pandoc."""
    if tok.startswith("[") and "?" in tok[:10]:
        return f"`{' '.join(tok.split())}`"
    m = re.fullmatch(r"(https?://\S+?|www\.\S+?|\S+\.(?:gov|org|edu)/\S*?)([.,;:)]*)", tok)
    if m:
        url = m.group(1) if m.group(1).startswith("http") else "https://" + m.group(1)
        return f"<{url}>{m.group(2)}"
    return tok


def convert(text):
    pieces, last = [], 0
    for m in PROTECT.finditer(text):
        pieces.append(mathify(text[last:m.start()]))
        pieces.append(protect(m.group(0)))
        last = m.end()
    pieces.append(mathify(text[last:]))
    return "".join(pieces)


def main():
    md = SRC.read_text().replace("](../results/", "](../../results/")
    title, body = md.split("\n", 1)
    assert title.startswith("# ")
    tmp = OUT.with_suffix(".md")
    tmp.write_text(convert(body))
    subprocess.run(["pandoc", str(tmp), "-f", "markdown-implicit_figures+pipe_tables+tex_math_dollars", "-t", "latex",
                    "--shift-heading-level-by=-1", "--wrap=preserve", "-o", str(OUT)], check=True)
    tmp.unlink()
    tex = OUT.read_text()
    tex = re.sub(r"\\texttt\{\{\[\}((?:VALUE|REF|CONFIRM|BUDGET)\?.*?)\{\]\}\}", r"\\placeholder{[\1]}", tex)
    def breakable(m):
        raw = m.group(1).replace("\\_", "_").replace("\\textasciitilde{}", "~")
        return f"\\nolinkurl{{{raw}}}" if re.fullmatch(r"[A-Za-z0-9_./*~+-]+", raw) else m.group(0)
    tex = re.sub(r"\\texttt\{([^{}]*)\}", breakable, tex)  # file paths and keys: allow line breaks
    OUT.write_text("% generated by scripts/report_faithful.py from docs/cavity_concept.md; do not edit\n" + tex)
    print("wrote", OUT.relative_to(ROOT))


if __name__ == "__main__":
    main()
