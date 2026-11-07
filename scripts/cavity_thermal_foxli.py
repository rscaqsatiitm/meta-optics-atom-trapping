"""Absorption-induced bulge of the graded mirrors G1, G2 in the Fox-Li ring (Section 5 check).
Usage:
    python scripts/cavity_thermal_foxli.py configs/cavity_thermal_foxli.yaml

The bulge profile across the strip is taken proportional to the circulating intensity on G (a
sharper profile than the thermal one, hence a conservative choice) and scaled to the centre-to-edge
bulge of results/cavity_thermal/thermal.json. Phase 4 pi h / lambda on both graded mirrors.
Writes thermal_foxli.json: depth rms and enhancement relative to the ideal design.
"""

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

from metatrap.runio import load_config, run_dir, write_json, write_params

ROOT = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location("cavity_foxli", ROOT / "scripts" / "cavity_foxli.py")
S = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(S)


def main(cfg_path):
    cfg = load_config(cfg_path)
    out = run_dir(cfg["run_id"])
    write_params(out, cfg)
    fcfg = load_config(ROOT / cfg["foxli_config"])
    s = S.build(fcfg, fcfg["variants"][cfg["variant"]])
    ring, ph = s["ring"], s["ph"]
    iw = np.abs(ring.to_mirror(s["F"])) ** 2
    shape = iw / iw.max()  # a constant (piston) part would be absorbed by the lock
    th = {d["name"]: d for d in json.load(open(ROOT / cfg["thermal_json"]))["deformation"]}
    k = 2 * np.pi / s["lam"]
    res = {}
    ev, _, _ = S.evaluate(fcfg, s, ring)
    res["ideal"] = {"depth_rms": ev["array"]["depth_rms"], "g_ratio": ev["g_ratio"]}
    print("ideal", res["ideal"])
    runs = [(n, 1.0) for n in cfg["cases"]] + [(cfg["cases"][0], f) for f in cfg["extra_scales"]]
    for name, f in runs:
        h_um = f * th[name]["bulge_centre_minus_edge_nm"] * 1e-3
        extra = 2 * k * h_um * shape
        ev, _, _ = S.evaluate(fcfg, s, S.with_phases(ring, ph + extra, ph + extra))
        key = name if f == 1.0 else f"{name} x{f:g}"
        res[key] = {"bulge_nm": h_um * 1e3, "depth_rms": ev["array"]["depth_rms"], "g_ratio": ev["g_ratio"],
                    "min_over_mean": ev["array"]["depth_min_over_mean"]}
        print(key, res[key])
    write_json(out / "thermal_foxli.json", res)


if __name__ == "__main__":
    main(sys.argv[1])
