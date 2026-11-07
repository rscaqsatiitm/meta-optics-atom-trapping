"""Recompute two stored diagnostics of metalens_doublet runs after code fixes.

Usage: python scripts/metalens_recheck.py configs/metalens_doublet.yaml
       python scripts/metalens_recheck.py configs/metalens_1075.yaml

(i) Transmitted fraction of the apodised wave checks (below). (ii) wave_axial.on_axis_focus
used 1.1 x 0.5 instead of 1.1 x the design NA for "spectrum_fraction_beyond_1.1NA"; it is
recomputed for every stored wave_axial entry (the other wave_axial values are unchanged).

wave_offaxis.strehl_offaxis returned a J-weighted mean of eta1 eta2 as "transmitted_weight"
(fixed in wave_offaxis.transmitted_fraction). The Strehl ratios and peaks were unaffected. This
recomputes the fraction for the stored apodised cases (same designs, field angles, pupil
sampling and efficiency curves as scripts/metalens_doublet.py) without the focal sums, and
updates results/metalens_doublet/doublet.json (key transmitted_weight; the previous value is
kept as transmitted_weight_jweighted).
"""

import json
import sys

import numpy as np

from metatrap.metalens import doublet as D
from metatrap.metalens import wave_offaxis as WO
from metatrap.runio import REPO_ROOT, load_config, write_json

sys.path.insert(0, str(REPO_ROOT / "scripts"))
from metalens_doublet import spec_for  # noqa: E402


def main(cfg_path):
    cfg = load_config(cfg_path)
    path = REPO_ROOT / "results" / cfg["run_id"] / "doublet.json"
    res = json.load(open(path))
    ep = REPO_ROOT / "results" / cfg["run_id"] / "efficiency.json"
    eff = json.load(open(ep)) if ep.exists() else {}
    lams = cfg["wavelengths_um"]
    spec = spec_for(cfg)
    d0 = D.make(spec, lams[0])
    R1, f = d0.info["R1_um"], d0.info["f_um"]
    from metatrap.metalens import wave_axial as W

    for lam in lams:
        key = str(lam)
        if key in res.get("wave_axial", {}) and "wave_axial" in cfg:
            dd = d0 if lam == lams[0] else D.make(spec, lam, R1_um=R1)
            r = W.on_axis_focus(dd, lam, **cfg["wave_axial"])
            old = res["wave_axial"][key].get("spectrum_fraction_beyond_1.1NA")
            res["wave_axial"][key]["spectrum_fraction_beyond_1.1NA"] = r["spectrum_fraction_beyond_1.1NA"]
            print("wave_axial", key, "beyond 1.1 NA:", old, "->", r["spectrum_fraction_beyond_1.1NA"], flush=True)
    wo = cfg.get("wave_offaxis", {"wavelengths_um": []})
    for lam in wo["wavelengths_um"]:
        dd = d0 if lam == lams[0] else D.make(spec, lam, R1_um=R1)
        sysm = D.system(dd, lam)
        for h in wo.get("apodised_heights_um", []):
            th = D.theta_for_height(sysm, R1, lam, h, f)
            pf = WO.plate_field(dd, lam, th, wo["n_pupil"])
            for mat, e in eff.get(str(lam), {}).items():
                key = f"{lam}:{h}:{mat}"
                if key not in res["wave_offaxis"]:
                    continue
                sg, eg = np.array(e["sin_theta"]), np.array(e["eta_xy_mean"])
                amp = pf["amp"] * np.sqrt(np.interp(pf["s1"], sg, eg) * np.interp(pf["s2"], sg, eg))
                r = res["wave_offaxis"][key]
                if "transmitted_weight_jweighted" in r:
                    continue  # already recomputed
                r["transmitted_weight_jweighted"] = r["transmitted_weight"]
                r["transmitted_weight"] = WO.transmitted_fraction(pf["amp"], amp)
                print(key, "transmitted", round(r["transmitted_weight_jweighted"], 4), "->", round(r["transmitted_weight"], 4),
                      "peak_rel_uniform_lossless", round(r["peak_rel_uniform_lossless"], 4), flush=True)
    write_json(path, res)


if __name__ == "__main__":
    main(sys.argv[1])
