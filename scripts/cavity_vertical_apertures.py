"""Round-trip design loss of the vertical-lattice ring's 1D modes versus mirror aperture (strip size).
Usage:
    python scripts/cavity_vertical_apertures.py configs/cavity_vertical_apertures.yaml

For each case the ideal design is rebuilt (scripts/cavity_vertical.py:build) with the given half-
apertures on P1 (ap1), P2 (ap2) and G (ap_g) and the round-trip loss of the design target,
1 - |<F, R F>|^2 / <F, F>^2, is reported, with the driven enhancement relative to the lossless value
when requested. Writes apertures.json.
"""

import sys

import numpy as np

from metatrap.runio import REPO_ROOT, load_config, run_dir, write_json, write_params

sys.path.insert(0, str(REPO_ROOT / "scripts"))
import cavity_vertical as V  # noqa: E402


def main(cfg_path):
    cfg = load_config(cfg_path)
    out = run_dir(cfg["run_id"])
    write_params(out, cfg)
    rep = {}
    for case in cfg["cases"]:
        vcfg = load_config(REPO_ROOT / case["vertical_config"])
        d = dict(vcfg[case["direction"]], **case.get("overrides", {}))
        s = V.build(vcfg, d)
        F = s["F"]
        lam = np.vdot(F, s["ring"].roundtrip(F)) / np.vdot(F, F)
        res = {"vertical_config": case["vertical_config"], "direction": case["direction"], **case.get("overrides", {}),
               "design_loss_per_roundtrip": float(1 - abs(lam) ** 2)}
        if case.get("solve"):
            res["g_ratio"] = float(V.solve(vcfg, s, s["ring"])[2])
        rep[case["name"]] = res
        print(case["name"], {k: v for k, v in res.items() if k not in ("vertical_config",)}, flush=True)
        write_json(out / "apertures.json", rep)


if __name__ == "__main__":
    main(sys.argv[1])
