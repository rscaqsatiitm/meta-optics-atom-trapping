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

import json
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
    # Ring layout in the atom plane: X at origin; each ring's arms cross at +-theta about its axis.
    th = B.crossing_half_angle(cfg["period_um"], cfg["wavelength_nm"] * 1e-3)
    rings = {name: ring_geometry(th, d0, f1, f2, L, np.radians(r["turn_deg"]), np.radians(r["rotation_deg"]), t["ap1_mm"], t["ap2_mm"])
             for name, r in cfg["rings"].items()}
    names = list(rings)
    layout = {"theta_deg": float(np.degrees(th)), "half_trip_path_m": float(d0 + D), "round_trip_m": float(4 * (d0 + D)),
              "fsr_MHz": float(299792458.0 / (4 * (d0 + D)) / 1e6), "rings": {}}
    for name, rg in rings.items():
        layout["rings"][name] = {"turn_deg": cfg["rings"][name]["turn_deg"], "rotation_deg": cfg["rings"][name]["rotation_deg"],
                                 "aoi_deg": rg["aoi_deg"], "mirrors_mm": {k: v.tolist() for k, v in rg["mirrors"].items()},
                                 "diagonal_clearance_mm": diagonal_clearance(rg, cfg["diagonal_min_radius_mm"])}
    layout["ring_clearance_mm"] = clearance(rings[names[0]], rings[names[1]], cfg["beam_margin_mm"])
    rep = {"vertical_designs": rows, "layout": layout}
    write_json(out / "layout.json", rep)
    for r in rows:
        print(r)
    print(json.dumps(layout, indent=1)[:1500])
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(6.5, 6.0))
        for (name, rg), col in zip(rings.items(), ("#2a78d6", "#eb6834")):
            for path in rg["paths"]:
                ax.plot(path[:, 0], path[:, 1], "-", color=col, lw=1, label=name if path is rg["paths"][0] else None)
            for k, (p_, hw) in rg["mirror_sizes"].items():
                ax.add_patch(plt.Circle(p_, hw, fill=False, ec=col, lw=0.8))
        ax.add_patch(plt.Rectangle((-15, -15), 30, 30, fill=False, ec="0.4", lw=1))
        ax.set_aspect("equal")
        ax.set_xlabel("x (mm)")
        ax.set_ylabel("y (mm)")
        ax.legend(fontsize=8, frameon=False)
        ax.set_title(f"In-plane rings, clearance {layout['ring_clearance_mm']:.1f} mm; cell outline 30 mm [CONFIRM?]", fontsize=8)
        fig.tight_layout()
        fig.savefig(out / "layout.png", dpi=150)
    except Exception as exc:  # plotting is optional
        print("plot failed", exc)


def ring_geometry(th, d0, f1, f2, L, turn, rot, ap1_mm, ap2_mm):
    """Mirror positions (mm) of a self-intersecting telescope ring: arms cross at X at +-theta; the
    primary P1 at d0 turns the beam by `turn` (direction of the outgoing ray from +x) towards the
    secondary P2 (f1 + f2 further); the graded mirror G lies on the ring's mirror axis at L from P2.
    The second half trip (G2) is the mirror image across the ring axis; `rot` rotates the ring."""
    P1 = d0 * np.array([np.cos(th), np.sin(th)])
    P2 = P1 + (f1 + f2) * np.array([np.cos(turn), np.sin(turn)])
    G = np.array([0.0, P2[1] + np.sqrt(max(L**2 - P2[0] ** 2, 0.0))])
    up = [np.zeros(2), P1, P2, G, P2 * [-1, 1], P1 * [-1, 1], np.zeros(2)]
    dn = [p * [1, -1] for p in up]
    R = np.array([[np.cos(rot), -np.sin(rot)], [np.sin(rot), np.cos(rot)]])

    def aoi(prev, here, nxt):
        a, b = prev - here, nxt - here
        return float(np.degrees(0.5 * np.arccos(np.dot(a, b) / np.linalg.norm(a) / np.linalg.norm(b))))

    names = ["P1", "P2", "G1", "P2'", "P1'"]
    mirrors, sizes = {}, {}
    for half, pts in (("", up), ("_lower", dn)):
        for i, n in enumerate(names):
            key = (n.replace("G1", "G2") if half else n) + half
            mirrors[key] = R @ pts[i + 1] * 1e3
            sizes[key] = (mirrors[key], ap1_mm if n.startswith("P1") else ap2_mm)
    paths = [np.array([R @ p * 1e3 for p in pts]) for pts in (up, dn)]
    return {"mirrors": mirrors, "mirror_sizes": sizes, "paths": paths,
            "aoi_deg": {"P1": aoi(np.zeros(2), P1, P2), "P2": aoi(P1, P2, G), "G": aoi(P2, G, P2 * [-1, 1])}}


def _segdist(p, a, b):
    t = np.clip(np.dot(p - a, b - a) / np.dot(b - a, b - a), 0, 1)
    return float(np.linalg.norm(a + t * (b - a) - p))


def clearance(r1, r2, margin):
    """Smallest gap (mm) between mirrors of two rings, or between a mirror and the other ring's beam
    path (beam half-width plus `margin` taken as 3 mm), using the mirror apertures as radii."""
    c = np.inf
    for p, h in r1["mirror_sizes"].values():
        for q, k in r2["mirror_sizes"].values():
            c = min(c, np.linalg.norm(p - q) - h - k)
    for ra, rb in ((r1, r2), (r2, r1)):
        for p, h in rb["mirror_sizes"].values():
            for path in ra["paths"]:
                for a, b in zip(path[:-1], path[1:]):
                    c = min(c, _segdist(p, a, b) - h - margin)
    return float(c)


def diagonal_clearance(rg, r_min):
    """Smallest distance (mm) of the ring's mirrors from the two array diagonals beyond radius r_min,
    the corridor reserved for the vertical-lattice ring."""
    c = np.inf
    for p, h in rg["mirror_sizes"].values():
        if np.linalg.norm(p) < r_min:
            continue
        for d in (np.array([1, 1]) / np.sqrt(2), np.array([1, -1]) / np.sqrt(2)):
            c = min(c, abs(d[0] * p[1] - d[1] * p[0]) - h)
    return float(c)


if __name__ == "__main__":
    main(sys.argv[1])
