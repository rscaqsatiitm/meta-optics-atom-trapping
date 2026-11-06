"""Vector focusing check at NA 0.55 for the in-cell pupil.

Usage: python scripts/vector_check.py configs/vector_check.yaml
"""

import sys

import numpy as np

from metatrap import vector_focus as V
from metatrap.atoms.alkali import alpha_au
from metatrap.hologram.wgs import PupilGrid, pupil_amplitude
from metatrap.propagation.layers import LayerStack
from metatrap.runio import load_config, run_dir, write_json, write_params


def main(cfg_path):
    cfg = load_config(cfg_path)
    out = run_dir(cfg["run_id"])
    write_params(out, cfg)
    lam, na = cfg["wavelength_um"], cfg["na"]
    stack = LayerStack.from_list(cfg["stack"])
    R = stack.aperture_for_na(na, lam)
    ill = dict(cfg["illumination"])
    if "edge_frac" in ill:
        ill["edge_um"] = ill.pop("edge_frac") * R
    g = PupilGrid(n=cfg["pupil_n"], ds=2.4 * na / cfg["pupil_n"], lam_um=lam)
    A = pupil_amplitude(g, stack, na, ill)
    s1 = g.s_coords()
    SX, SY = np.meshgrid(s1, s1, indexing="ij")
    m = A > 0
    sx, sy, amp = SX[m], SY[m], A[m]
    h = cfg["grid_step_um"]
    xs = np.arange(-cfg["grid_half_um"], cfg["grid_half_um"] + h / 2, h)
    X, Y = np.meshgrid(xs, xs, indexing="ij")
    pts = np.column_stack([X.ravel(), Y.ravel(), np.zeros(X.size)])
    c = len(xs) // 2
    res = {"r_v": V.vector_ratio(lam * 1e3), "alpha_s_au": float(alpha_au(lam * 1e3))}
    Es = V.focal_field(sx, sy, amp, lam, pts, None)[0].reshape(X.shape)
    Is = np.abs(Es) ** 2
    res["scalar"] = {"peak": float(Is[c, c]), "waist_x_um": V.curvature_waist(Is[:, c], h), "waist_y_um": V.curvature_waist(Is[c, :], h)}
    for model in ("rotation", "projection"):
        Ex, Ey, Ez = (e.reshape(X.shape) for e in V.focal_field(sx, sy, amp, lam, pts, model))
        I = np.abs(Ex) ** 2 + np.abs(Ey) ** 2 + np.abs(Ez) ** 2
        cross = np.stack([np.imag(np.conj(Ey) * Ez - np.conj(Ez) * Ey), np.imag(np.conj(Ez) * Ex - np.conj(Ex) * Ez), np.imag(np.conj(Ex) * Ey - np.conj(Ey) * Ex)])
        C = cross / np.maximum(I, 1e-30)  # circularity vector Im(E* x E)/|E|^2
        core = I >= 0.5 * I[c, c]
        entry = {
            "peak_over_scalar": float(I[c, c] / Is[c, c]),
            "waist_x_um": V.curvature_waist(I[:, c], h),
            "waist_y_um": V.curvature_waist(I[c, :], h),
            "Ez2_max_over_peak": float((np.abs(Ez) ** 2).max() / I[c, c]),
            "Ez2_fraction_power_in_window": float(np.sum(np.abs(Ez) ** 2) / np.sum(I)),
        }
        for axis, name in ((0, "x"), (1, "y"), (2, "z")):
            entry[f"C{name}_max_in_half_max_core"] = float(np.abs(C[axis][core]).max())
        # Fictitious field per mK of trap depth, along y (C_y is the dominant component for x pol).
        b_per_mK = abs(res["r_v"]) * V.KB_HZ_PER_K * 1e-3 / V.MU_B_HZ_PER_G  # G per unit circularity per mK
        cy_line = C[1][:, c]
        grad = np.gradient(cy_line, h)[c]
        entry["B_fict_y_max_core_G_per_mK"] = float(b_per_mK * np.abs(C[1][core]).max())
        entry["dBy_dx_centre_G_per_um_per_mK"] = float(b_per_mK * grad)
        res[model] = entry
        print(model, entry, flush=True)
    print("scalar", res["scalar"], "r_v", res["r_v"])
    write_json(out / "vector.json", res)


if __name__ == "__main__":
    main(sys.argv[1])
