"""Vector focusing of off-axis sites (up to the array corners) at NA 0.55.

Usage: python scripts/vector_offaxis.py configs/vector_offaxis.yaml

For a site at x_j the exact-kernel hologram focuses every ray from metasurface point r_i
onto x_j. The ray direction in the atom-plane medium is found exactly through the layer stack
(transverse distance |x_j - r_i| inverted with the ray mapping). The focal field near x_j is
the plane-wave superposition sum_i A_i p_i exp(i k0 [sigma'_i . d + cz'_i dz]), with A_i the
on-axis pupil amplitude (its change across the field is neglected) and p_i the polarisation
of an x-polarised field projected onto the plane transverse to each ray. The scalar field uses
p = 1 with the same directions. Reported per site: vector/scalar peak, waists along x and y,
longitudinal field, circularity in the half-maximum core, fictitious magnetic field on Cs,
and the vector-induced shift of the intensity maximum.
"""

import sys

import numpy as np

from metatrap import vector_focus as V
from metatrap.hologram.field_model import ExactOPL
from metatrap.hologram.wgs import PupilGrid, pupil_amplitude
from metatrap.propagation.layers import LayerStack
from metatrap.runio import load_config, run_dir, write_json, write_params


def site_field(sx, sy, amp, rx, ry, xj, ex, lam, pts, vector):
    rho = np.hypot(xj[0] - rx, xj[1] - ry)
    _, s = ex(rho)
    with np.errstate(invalid="ignore", divide="ignore"):
        ux = np.where(rho > 0, (xj[0] - rx) / rho, 0.0)
        uy = np.where(rho > 0, (xj[1] - ry) / rho, 0.0)
    sxp, syp = ux * s, uy * s
    cz = np.sqrt(1 - s**2)
    k0 = 2 * np.pi / lam
    ph = np.exp(1j * k0 * (np.outer(sxp, pts[:, 0]) + np.outer(syp, pts[:, 1]) + np.outer(cz, pts[:, 2])))
    if not vector:
        return (amp @ ph)[None, :]
    px, py, pz = V.polarisation(sxp, syp, "projection")
    W = np.stack([amp * px, amp * py, amp * pz])
    return W @ ph


def main(cfg_path):
    cfg = load_config(cfg_path)
    out = run_dir(cfg["run_id"])
    write_params(out, cfg)
    lam, na = cfg["wavelength_um"], cfg["na"]
    stack = LayerStack.from_list(cfg["stack"])
    g = PupilGrid(n=cfg["pupil_n"], ds=2.4 * na / cfg["pupil_n"], lam_um=lam)
    A = pupil_amplitude(g, stack, na, {"kind": "top_hat"})
    s1 = g.s_coords()
    SX, SY = np.meshgrid(s1, s1, indexing="ij")
    m = A > 0
    sx, sy, amp = SX[m], SY[m], A[m]
    s = np.minimum(np.hypot(sx, sy), na)
    rho = stack.ray_height(s, lam)
    with np.errstate(invalid="ignore", divide="ignore"):
        u = np.where(s > 0, 1 / np.where(s > 0, np.hypot(sx, sy), 1), 0.0)
    rx, ry = -sx * u * rho, -sy * u * rho
    ex = ExactOPL(stack, lam)
    h = cfg["grid_step_um"]
    xs = np.arange(-cfg["grid_half_um"], cfg["grid_half_um"] + h / 2, h)
    X, Y = np.meshgrid(xs, xs, indexing="ij")
    pts = np.column_stack([X.ravel(), Y.ravel(), np.zeros(X.size)])
    c = len(xs) // 2
    b_per_mK = abs(V.vector_ratio(lam * 1e3)) * V.KB_HZ_PER_K * 1e-3 / V.MU_B_HZ_PER_G
    res = {}
    for xj in cfg["sites_um"]:
        Es = site_field(sx, sy, amp, rx, ry, xj, ex, lam, pts, False)[0].reshape(X.shape)
        Is = np.abs(Es) ** 2
        Ex, Ey, Ez = (e.reshape(X.shape) for e in site_field(sx, sy, amp, rx, ry, xj, ex, lam, pts, True))
        I = np.abs(Ex) ** 2 + np.abs(Ey) ** 2 + np.abs(Ez) ** 2
        cross_y = np.imag(np.conj(Ez) * Ex - np.conj(Ex) * Ez) / np.maximum(I, 1e-30)
        cross_x = np.imag(np.conj(Ey) * Ez - np.conj(Ez) * Ey) / np.maximum(I, 1e-30)
        core = I >= 0.5 * I.max()
        iv, jv = np.unravel_index(np.argmax(I), I.shape)
        isc, jsc = np.unravel_index(np.argmax(Is), Is.shape)
        entry = {
            "peak_over_scalar": float(I.max() / Is.max()),
            "scalar_peak_over_onaxis_scalar": None,
            "waist_x_um": V.curvature_waist(I[:, jv], h), "waist_y_um": V.curvature_waist(I[iv, :], h),
            "Ez2_max_over_peak": float((np.abs(Ez) ** 2).max() / I.max()),
            "Cy_max_core": float(np.abs(cross_y[core]).max()), "Cx_max_core": float(np.abs(cross_x[core]).max()),
            "B_fict_max_core_G_per_mK": float(b_per_mK * np.hypot(cross_x, cross_y)[core].max()),
            "vector_minus_scalar_peak_shift_um": (float((iv - isc) * h), float((jv - jsc) * h)),
        }
        res[f"{xj[0]:.0f},{xj[1]:.0f}"] = entry | {"_scalar_peak": float(Is.max())}
        print(xj, {k: (round(v, 4) if isinstance(v, float) else v) for k, v in entry.items() if v is not None}, flush=True)
    p0 = res[next(iter(res))]["_scalar_peak"]
    for k in res:
        res[k]["scalar_peak_over_onaxis_scalar"] = res[k].pop("_scalar_peak") / p0
    write_json(out / "vector_offaxis.json", res)


if __name__ == "__main__":
    main(sys.argv[1])
