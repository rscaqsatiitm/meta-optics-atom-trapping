"""Per-site comparison of two evaluations of the same design (e.g. fold window against unfolded grid).

Usage: python scripts/compare_per_site.py A.npz KEY_A B.npz KEY_B OUT.json

Each file holds 'sites' and a per-site peak array under the given key ('peak' for
scripts/fullaperture.py outputs, 'peak_<n_fold>' for scripts/fold_window_study.py). The peaks
are normalised to their means; reported are the rms and maximum of the difference, the
correlation of the two error patterns, and each one's own rms. The rms difference is the
floor below which a feedback loop measured with model A cannot be trusted to improve model B.
"""

import json
import sys

import numpy as np


def main(a_path, ka, b_path, kb, out):
    A, B = np.load(a_path), np.load(b_path)
    if not np.allclose(A["sites"], B["sites"]):
        raise ValueError("site lists differ")
    ea, eb = A[ka] / A[ka].mean() - 1, B[kb] / B[kb].mean() - 1
    d = ea - eb
    res = {"a": f"{a_path}:{ka}", "b": f"{b_path}:{kb}", "n_sites": int(len(ea)), "rms_a": float(ea.std()), "rms_b": float(eb.std()),
           "rms_difference": float(d.std()), "max_abs_difference": float(np.abs(d).max()),
           "correlation": float(np.corrcoef(ea, eb)[0, 1]) if ea.std() > 0 and eb.std() > 0 else None,
           "mean_ratio_b_over_a": float(B[kb].mean() / A[ka].mean())}
    print(res)
    json.dump(res, open(out, "w"), indent=1)


if __name__ == "__main__":
    main(*sys.argv[1:6])
