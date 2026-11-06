"""Markdown tolerance table of docs/cavity_concept.md (Section 4) from results/cavity_foxli/foxli.json.
Usage:
    python scripts/cavity_report_tables.py
Entries: barrier-depth non-uniformity over the 10^6-site array (rms, median over realisations)
and, in parentheses, the lowest enhancement relative to the analytic value among those realisations.
"""

import json
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


def main():
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
    head = "| Perturbation | " + " | ".join(lab for _, lab in VARIANTS) + " |"
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


if __name__ == "__main__":
    main()
