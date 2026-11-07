"""Recompute the power budget of finished full-aperture runs from their stored metrics.

Usage: python scripts/recompute_power_budget.py

For every results/<run_id>/metrics.json written by scripts/fullaperture.py, the power-budget
inputs are read from configs/<run_id>.yaml (the current values) and the utilisation and trap
depths are recomputed from the stored core efficiencies, so a corrected factor does not require
re-running the simulations. The inputs used are stored as power_budget_inputs.
"""

import importlib.util
import json
from pathlib import Path

from metatrap.runio import REPO_ROOT, load_config, write_json

_spec = importlib.util.spec_from_file_location("fa", Path(__file__).with_name("fullaperture.py"))
FA = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(FA)


def main():
    for mp in sorted((REPO_ROOT / "results").glob("fullaperture*/metrics.json")):
        run_id = mp.parent.name
        cp = REPO_ROOT / "configs" / f"{run_id}.yaml"
        res = json.load(open(mp))
        if not cp.exists() or "runs" not in res:
            continue
        cfg = load_config(cp)
        if not cfg.get("power_budget") or not res["runs"]:
            continue
        old = res.get("power_budget", {})
        res["power_budget"], res["power_budget_inputs"] = FA.power_budget(cfg["power_budget"], res["runs"], int(res["n_sites"]), cfg["wavelength_um"])
        write_json(mp, res)
        for name, row in res["power_budget"].items():
            o = old.get(name, {}).get("depth_mK_at_50W", float("nan"))
            print(f"{run_id:30s} {name:20s} depth at 50 W: {o:.4f} -> {row['depth_mK_at_50W']:.4f} mK (utilisation {row['utilisation']:.4f})")


if __name__ == "__main__":
    main()
