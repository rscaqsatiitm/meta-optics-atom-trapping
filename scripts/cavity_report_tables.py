"""Tables of docs/cavity_concept.md (Sections 4 and 4b) from results/cavity_foxli/foxli.json and
results/cavity_correction/correction.json, as Markdown or as LaTeX for docs/report/.
Usage:
    python scripts/cavity_report_tables.py                 # Markdown tolerance table (Section 4)
    python scripts/cavity_report_tables.py correction      # Markdown correction tables (Section 4b)
    python scripts/cavity_report_tables.py latex docs/report/tables   # both, as LaTeX tabulars
Entries: barrier-depth non-uniformity over the 10^6-site array (rms, median over realisations)
and, in parentheses, the lowest enhancement relative to the analytic value among those realisations.
"""

import json
import re
from collections import defaultdict

import numpy as np

from metatrap.runio import REPO_ROOT

VARIANTS = [("plane150", "Plane folds, 150 mm"), ("plane300", "Plane folds, 300 mm"), ("telK_M0.2", "Telescope 0.2x K"),
            ("telK_M0.1", "Telescope 0.1x K"), ("telG_M0.2", "Telescope 0.2x G")]


def cell(entries):
    if not entries:
        return "-"
    d = np.median([e["array"]["depth_rms"] for e in entries]) * 100
    g = min(e["g_ratio"] for e in entries)
    return f"{d:.2g} % ({g:.2f})"


def tolerance_rows():
    """Header and rows of the Section 4 tolerance table."""
    r = json.loads((REPO_ROOT / "results/cavity_foxli/foxli.json").read_text())
    rows = []

    def group(prefix, field, values, label, unit):
        for val in values:
            out = [f"{label}, {val:g} {unit}"]
            for v, _ in VARIANTS:
                out.append(cell([e for k, e in r.items() if k.startswith(prefix.format(v=v)) and abs(e.get(field, -1) - val) < 1e-12]))
            rows.append(out)

    group("figure_{v}_rms", "rms_waves", [0.001, 0.003, 0.01], "G figure, rms", "waves")
    group("figure_primary_{v}_rms", "rms_waves", [0.001, 0.003, 0.01], "P1 figure, rms", "waves")
    group("curvature_{v}_sag", "sag_waves_at_beam_halfwidth", [0.001, 0.003, 0.01], "G x^2 sag at beam half-width", "waves")
    group("tilt_{v}_", "tilt_urad", [0.3, 1.0, 3.0, 10.0], "G tilt", "urad")
    group("tilt_primary_{v}_", "tilt_urad", [0.3, 1.0, 3.0, 10.0], "P1 tilt", "urad")
    return ["Perturbation"] + [lab for _, lab in VARIANTS], rows, r


def main():
    header, rows, r = tolerance_rows()
    head = "| " + " | ".join(header) + " |"
    print(head)
    print("|" + "---|" * (len(VARIANTS) + 1))
    for row in rows:
        print("| " + " | ".join(row) + " |")
    # Single-variant studies, for the text.
    print()
    for k, e in r.items():
        if k.startswith(("input_", "meta_", "finesse_", "reference_")):
            print(f"{k}: depth {e['array']['depth_rms']*100:.3g} % G/G0 {e['g_ratio']:.3f} dy_affine_res {e['array']['dy_affine_residual_rms_um']*1e3:.1f} nm")
    worst = defaultdict(float)
    for k, e in r.items():
        if isinstance(e, dict) and "array" in e and isinstance(e.get("variant"), str) and e["variant"].startswith("tel") and e["array"]["depth_rms"] < 0.05:
            worst[e["variant"]] = max(worst[e["variant"]], e["array"]["dy_affine_residual_rms_um"] * 1e3)
    print("largest non-affine site displacement (nm rms) among telescope cases with depth rms < 5 %:", dict(worst))


CORRECTION_CASES = [("fig1000", "lambda/1000 figure, 2 mm"), ("fig1000_c0.5", "lambda/1000 figure, 0.5 mm"), ("fig300", "lambda/300 figure, 2 mm"),
                    ("tilt1", "tilts, 1 urad rms"), ("fig1000_tilt0.3", "lambda/1000 + 0.3 urad tilts")]
CORRECTION_COLS = [("uncorrected", "Uncorrected"), ("openloop_noise0.0nm", "Perfect metrology"), ("openloop_noise0.1nm", "0.1 nm metrology"),
                   ("openloop_noise0.3nm", "0.3 nm metrology"), ("openloop_noise1.0nm", "1 nm metrology"), ("closedloop_noise0.01", "Atom loop, 6 rounds"),
                   ("closedloop_camera_rand0.5pct", "Camera loop, 0.5 % noise"), ("closedloop_camera_rand0.5pct_prnu1pct", "+ 1 % pixel response"),
                   ("closedloop_camera_rand0.5pct_flat1pct", "+ 1 % flat field")]


def correction_rows(v, r=None):
    """Header and rows of the Section 4b correction table for ring variant v (columns without data dropped)."""
    r = r or json.loads((REPO_ROOT / "results/cavity_correction/correction.json").read_text())
    rows = []
    for cname, clab in CORRECTION_CASES:
        row = [clab]
        for ckey, _ in CORRECTION_COLS:
            vals = []
            for seed in (1, 2):
                e = r.get(f"{v}_{cname}_seed{seed}_{ckey}")
                if e is None:
                    continue
                d = e["final"]["depth_rms"] if "final" in e else e["depth_rms"]
                vals.append(f"{d * 100:.2g}")
            row.append(" / ".join(vals) + " %" if vals else "-")
        rows.append(row)
    keep = [0] + [j + 1 for j in range(len(CORRECTION_COLS)) if any(row[j + 1] != "-" for row in rows)]
    header = ["Errors on every mirror"] + [c for _, c in CORRECTION_COLS]
    return [header[j] for j in keep], [[row[j] for j in keep] for row in rows]


def _tex(cell):
    """LaTeX form of a table cell written in the report's ASCII notation."""
    out = cell.replace("%", r"\%").replace("x^2", r"$x^2$")
    out = re.sub(r"lambda/(\d+)", r"$\\lambda/\1$", out)
    out = re.sub(r"(\d(?:\.\d+)?)e\+0?(\d+)", r"$\1\\times10^{\2}$", out)
    out = out.replace(" urad", r"~\si{\micro\radian}").replace(" um", r"~\si{\micro\metre}")
    return out


def latex_tabular(header, rows, first_width="3.2cm"):
    n = len(header) - 1
    spec = rf"@{{}}>{{\raggedright\arraybackslash}}p{{{first_width}}}" + r">{\centering\arraybackslash}X" * n + "@{}"
    lines = [rf"\begin{{tabularx}}{{\linewidth}}{{{spec}}}", r"\toprule",
             " & ".join(r"\hspace{0pt}" + _tex(h) for h in header) + r" \\", r"\midrule"]
    lines += [" & ".join(_tex(c) for c in row) + r" \\" for row in rows]
    lines += [r"\bottomrule", r"\end{tabularx}", ""]
    return "\n".join(lines)


def write_latex(outdir):
    out = REPO_ROOT / outdir
    out.mkdir(parents=True, exist_ok=True)
    header, rows, _ = tolerance_rows()
    (out / "tolerance.tex").write_text("% generated by scripts/cavity_report_tables.py; do not edit\n" + latex_tabular(header, rows, "3.6cm"))
    for v, name in (("telK_M0.2", "correction_telescope.tex"), ("plane150", "correction_plain.tex")):
        h, rr = correction_rows(v)
        (out / name).write_text("% generated by scripts/cavity_report_tables.py; do not edit\n" + latex_tabular(h, rr, "2.1cm"))
    print("wrote", sorted(p.name for p in out.glob("*.tex")))


def correction_table():
    """Markdown table of results/cavity_correction/correction.json (Section 4b): depth rms per case,
    seeds separated by ' / '; closed loop after the last round (1 % per-site measurement noise)."""
    r = json.loads((REPO_ROOT / "results/cavity_correction/correction.json").read_text())
    cases = [("fig1000", "lambda/1000 figure, 2 mm"), ("fig1000_c0.5", "lambda/1000 figure, 0.5 mm"), ("fig300", "lambda/300 figure, 2 mm"),
             ("tilt1", "tilts, 1 urad rms"), ("fig1000_tilt0.3", "lambda/1000 + 0.3 urad tilts")]
    cols = [("uncorrected", "Uncorrected"), ("openloop_noise0.0nm", "Perfect metrology"), ("openloop_noise0.1nm", "0.1 nm metrology"),
            ("openloop_noise0.3nm", "0.3 nm metrology"), ("openloop_noise1.0nm", "1 nm metrology"), ("closedloop_noise0.01", "Atom loop, 6 rounds"),
            ("closedloop_camera_rand0.5pct", "Camera loop, 0.5 % noise"), ("closedloop_camera_rand0.5pct_prnu1pct", "+ 1 % pixel response"),
            ("closedloop_camera_rand0.5pct_flat1pct", "+ 1 % flat field")]
    for v in ("telK_M0.2", "plane150"):
        print(f"\n{v}\n")
        print("| Errors on every mirror | " + " | ".join(c for _, c in cols) + " |")
        print("|" + "---|" * (len(cols) + 1))
        for cname, clab in cases:
            row = [clab]
            for ckey, _ in cols:
                vals = []
                for seed in (1, 2):
                    e = r.get(f"{v}_{cname}_seed{seed}_{ckey}")
                    if e is None:
                        continue
                    d = e["final"]["depth_rms"] if "final" in e else e["depth_rms"]
                    vals.append(f"{d * 100:.2g}")
                row.append(" / ".join(vals) + " %" if vals else "-")
            print("| " + " | ".join(row) + " |")
    st = [(k, e["stroke_G1"]["rms_nm"], e["stroke_G1"]["pv_nm"]) for k, e in r.items() if "stroke_G1" in e and k.startswith("telK")]
    print("\nstroke rms nm range", min(s[1] for s in st), max(s[1] for s in st), "pv max", max(s[2] for s in st))


if __name__ == "__main__":
    import sys
    if len(sys.argv) == 1:
        main()
    elif sys.argv[1] == "correction":
        correction_table()
    elif sys.argv[1] == "latex":
        write_latex(sys.argv[2] if len(sys.argv) > 2 else "docs/report/tables")
