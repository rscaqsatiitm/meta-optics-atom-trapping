"""S4: geometry of the telescope ring and its vertical (Gaussian) mode; beam sizes on every mirror.
Usage:
    python scripts/cavity_layout.py configs/cavity_layout.yaml

Horizontal: the 1D Fox-Li design of S3 (flat-top width at X, compression by the telescope to M
at the secondary and the graded mirror). Vertical: a Gaussian mode with waist w_z at the
crossing, shaped by a vertical curvature on the primary (focal length f_z) and on the graded
mirror (radius R_z); the half trip is symmetric about G, so the vertical round trip is that of
a two-element system unfolded about G. For each f_z we solve R_z for the target waist and
report the stability, the round-trip Gouy phase and the vertical beam size on each mirror.
Ring layout in the atom plane: arms cross at X at +-theta; the telescope folds outward and the
graded mirrors sit on the bisector normal. Writes layout.json and layout.png.
"""

import sys

import numpy as np

from metatrap.cavity import budget as B
from metatrap.runio import load_config, run_dir, write_json, write_params


def abcd_prop(d):
    return np.array([[1.0, d], [0.0, 1.0]])


def abcd_lens(f):
    return np.array([[1.0, 0.0], [-1.0 / f, 1.0]])


def vertical_roundtrip(d0, D, f_z, R_z):
    """Half trip X -> d0 -> lens f_z -> D -> mirror R_z -> D -> lens f_z -> d0 -> X; round trip = half^2."""
    half = abcd_prop(d0) @ abcd_lens(f_z) @ abcd_prop(D) @ abcd_lens(R_z / 2) @ abcd_prop(D) @ abcd_lens(f_z) @ abcd_prop(d0)
    return half @ half, half


def waist_from_half(half, lam):
    """Eigen q at X of the symmetric half trip (waist at X by symmetry): returns w0 (m) or nan."""
    A, Bm, C, Dm = half.ravel()
    m = (A + Dm) / 2
    if abs(m) >= 1:
        return np.nan, m
    # q = i zR at X with zR = |B| / sqrt(1 - m^2) (A = D for a symmetric system)
    zr = abs(Bm) / np.sqrt(1 - m**2)
    return np.sqrt(lam * zr / np.pi), m


def beam_size(q, M):
    A, Bm, C, Dm = M.ravel()
    q2 = (A * q + Bm) / (C * q + Dm)
    return q2


def main(cfg_path):
    cfg = load_config(cfg_path)
    out = run_dir(cfg["run_id"])
    write_params(out, cfg)
    lam = cfg["wavelength_nm"] * 1e-9
    t = cfg["telescope"]
    f1, M = t["f1_m"], t["mag"]
    f2 = M * f1
    d0, L = t["d0_m"], t["L_m"]
    D = f1 + f2 + L  # primary -> graded mirror
    w0_target = cfg["w_z_m"]
    rows = []
    for f_z in cfg["f_z_scan_m"]:
        # Solve R_z (bisection on log R) for the target waist.
        best = None
        for R in np.geomspace(0.01, 100.0, 20001):
            _, half = vertical_roundtrip(d0, D, f_z, R)
            w, m = waist_from_half(half, lam)
            if np.isfinite(w):
                err = abs(w - w0_target) / w0_target
                if best is None or err < best[0]:
                    best = (err, R, w, m)
        if best is None or best[0] > 1e-3:
            rows.append({"f_z_m": f_z, "solution": None})
            continue
        _, R, w, m = best
        rt, half = vertical_roundtrip(d0, D, f_z, R)
        q0 = 1j * np.pi * w**2 / lam
        q_p1 = beam_size(q0, abcd_prop(d0))
        q_g = beam_size(q0, abcd_prop(D) @ abcd_lens(f_z) @ abcd_prop(d0))
        wz = lambda q: float(np.sqrt(-lam / (np.pi * np.imag(1 / q))))  # noqa: E731
        gouy_rt = float(2 * np.arccos(np.clip((rt[0, 0] + rt[1, 1]) / 2, -1, 1)))
        # Misalignment sensitivity: dR/R of 1e-3 -> relative waist change.
        _, half2 = vertical_roundtrip(d0, D, f_z, R * (1 + 1e-3))
        w2, _ = waist_from_half(half2, lam)
        rows.append({"f_z_m": f_z, "R_z_m": R, "w0_m": w, "half_trip_m": m, "gouy_roundtrip_rad": gouy_rt,
                     "w_z_at_primary_m": wz(q_p1), "w_z_at_graded_m": wz(q_g), "waist_change_per_1e-3_dR": float((w2 - w) / w)})
    # Ring layout in the atom plane (mm): X at origin, arms at +-theta.
    th = B.crossing_half_angle(cfg["period_um"], cfg["wavelength_nm"] * 1e-3)
    P1 = d0 * np.array([np.cos(th), np.sin(th)])
    P1b = d0 * np.array([-np.cos(th), np.sin(th)])
    turn = np.radians(t["primary_turn_deg"])  # outgoing direction from P1, measured from +x
    P2 = P1 + (f1 + f2) * np.array([np.cos(turn), np.sin(turn)])
    yG = P2[1] + np.sqrt(max(L**2 - P2[0] ** 2, 0.0))
    G = np.array([0.0, yG])
    P2b = np.array([-P2[0], P2[1]])

    def aoi(prev, here, nxt):
        a, b = prev - here, nxt - here
        return float(np.degrees(0.5 * np.arccos(np.dot(a, b) / np.linalg.norm(a) / np.linalg.norm(b))))

    layout = {"theta_deg": float(np.degrees(th)), "X": [0, 0], "P1": P1.tolist(), "P2": P2.tolist(), "G1": G.tolist(),
              "P2b": P2b.tolist(), "P1b": P1b.tolist(),
              "aoi_deg": {"P1": aoi(np.zeros(2), P1, P2), "P2": aoi(P1, P2, G), "G1": aoi(P2, G, P2b)},
              "half_trip_path_m": float(d0 + D), "round_trip_m": float(4 * (d0 + D)),
              "fsr_MHz": float(299792458.0 / (4 * (d0 + D)) / 1e6)}
    rep = {"vertical_designs": rows, "layout_m": layout}
    write_json(out / "layout.json", rep)
    for r in rows:
        print(r)
    print(layout)
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(6.5, 5.0))
        # Ring for the y-lattice, and its mirror image (lower half) for the second half trip (G2).
        pts_up = np.array([np.zeros(2), P1, P2, G, P2b, P1b, np.zeros(2)]) * 1e3
        pts_dn = pts_up * np.array([1, -1])
        for pts, c in ((pts_up, "C0"), (pts_dn, "C0")):
            ax.plot(pts[:, 0], pts[:, 1], "-", color=c, lw=1)
        rot = np.array([[0, -1], [1, 0]])
        for pts in (pts_up, pts_dn):
            q = pts @ rot.T
            ax.plot(q[:, 0], q[:, 1], "-", color="C1", lw=1)
        for name, p in (("P1", P1), ("P2", P2), ("G1", G)):
            ax.plot(*(p * 1e3), "ks", ms=4)
            ax.annotate(name, p * 1e3, textcoords="offset points", xytext=(4, 4), fontsize=8)
        ax.add_patch(plt.Rectangle((-15, -15), 30, 30, fill=False, ec="0.4", lw=1))
        ax.set_aspect("equal")
        ax.set_xlabel("x (mm)")
        ax.set_ylabel("y (mm)")
        ax.set_title("In-plane rings (y-lattice blue, x-lattice orange); cell outline 30 mm [CONFIRM?]", fontsize=8)
        fig.tight_layout()
        fig.savefig(out / "layout.png", dpi=150)
    except Exception as exc:  # plotting is optional
        print("plot failed", exc)


if __name__ == "__main__":
    main(sys.argv[1])
