"""Compare a meta-atom library computed with jones_direct_fourier against its pol_fourier recompute.

Usage: python scripts/compare_library_formulations.py data/library_main.h5 data/library_main_pf.h5 OUT_NAME

For every (material, pitch, height, wavelength) curve: the co-polarised transmission phase change
relative to the curve's mean offset (only relative phases matter for a hologram), the change of
|t|, and the cross-polarised amplitude, over all widths and over the fabricable widths (min feature
and gap 100 nm, aspect <= 10). Output: results/library_formulation_check/<OUT_NAME>.json.
"""

import sys

import numpy as np

from metatrap.rcwa import lookup as L
from metatrap.rcwa.library import Library
from metatrap.runio import REPO_ROOT, write_json


def main(old_path, new_path, name):
    a, b = Library.load(REPO_ROOT / old_path), Library.load(REPO_ROOT / new_path)
    lim = L.FabLimits(min_feature_um=0.10, min_gap_um=0.10, max_aspect=10.0)
    mats = [m.decode() if isinstance(m, bytes) else str(m) for m in a.materials]
    out = {"old": old_path, "new": new_path, "per_material_wavelength": {}}
    for im, mat in enumerate(mats):
        for il, lam in enumerate(a.wavelength_um):
            rows = []
            for ip, p in enumerate(a.pitch_um):
                for ih, h in enumerate(a.height_um):
                    ta_, tb_ = a.t[im, ip, ih, il], b.t[im, ip, ih, il]
                    if not (np.all(np.isfinite(ta_)) and np.all(np.isfinite(tb_)) and np.any(ta_ != 0) and np.any(tb_ != 0)):  # computed in both
                        continue
                    w = a.width_frac * p
                    ta, tb = a.t[im, ip, ih, il], b.t[im, ip, ih, il]
                    ok = L.feasible(w, p, h, lim)
                    d = np.angle(tb * np.conj(ta))
                    for tag, m in (("all", np.ones_like(ok)), ("fab", ok)):
                        if m.sum() < 2:
                            continue
                        dm = d[m] - np.angle(np.mean(np.exp(1j * d[m])))
                        rows.append({"tag": tag, "pitch": float(p), "height": float(h), "dphi_rel_max": float(np.max(np.abs(np.angle(np.exp(1j * dm))))),
                                     "dabs_max": float(np.max(np.abs(np.abs(tb[m]) - np.abs(ta[m])))),
                                     "xc_old_frac_above_0.01": float(np.mean(a.t_cross[im, ip, ih, il][m] > 0.01)),
                                     "xc_new_frac_above_0.01": float(np.mean(b.t_cross[im, ip, ih, il][m] > 0.01)),
                                     "xc_old_max": float(a.t_cross[im, ip, ih, il][m].max()), "xc_new_max": float(b.t_cross[im, ip, ih, il][m].max())})
            if not rows:
                continue
            s = {}
            for tag in ("all", "fab"):
                rr = [r for r in rows if r["tag"] == tag]
                dp = np.array([r["dphi_rel_max"] for r in rr])
                s[tag] = {"n_curves": len(rr), "dphi_rel_max_median": float(np.median(dp)), "dphi_rel_max_p90": float(np.percentile(dp, 90)),
                          "dphi_rel_max_max": float(dp.max()), "dabs_max_median": float(np.median([r["dabs_max"] for r in rr])),
                          "xc_old_frac_mean": float(np.mean([r["xc_old_frac_above_0.01"] for r in rr])),
                          "xc_new_frac_mean": float(np.mean([r["xc_new_frac_above_0.01"] for r in rr])),
                          "xc_old_max": float(max(r["xc_old_max"] for r in rr)), "xc_new_max": float(max(r["xc_new_max"] for r in rr))}
            out["per_material_wavelength"][f"{mat}:{lam:.4f}"] = s
            print(mat, lam, {k: {a_: round(b_, 4) if isinstance(b_, float) else b_ for a_, b_ in v.items()} for k, v in s.items()}, flush=True)
    (REPO_ROOT / "results" / "library_formulation_check").mkdir(parents=True, exist_ok=True)
    write_json(REPO_ROOT / "results" / "library_formulation_check" / f"{name}.json", out)


if __name__ == "__main__":
    main(*sys.argv[1:4])
