"""Fold geometry of the vertical-lattice ring.
Usage:
    python scripts/cavity_vertical_layout.py configs/cavity_vertical_layout.yaml

The vertical ring lies in the vertical plane V through the atoms that contains the array diagonal
d = (1, 1, 0)/sqrt(2); in-plane coordinates (s along d, z up). Its two arms cross at X at +-theta_z.
A planar ring has no image or polarisation rotation: the wide (6.9 mm) flat-top axis and the
polarisation stay along u = (-1, 1, 0)/sqrt(2), perpendicular to V, i.e. s-polarised at every
mirror, and the cylindrical telescope mirrors curve along u, perpendicular to the plane of
incidence, where an (ideally parabolic) cylinder acts as an exact cylindrical lens of focal length
R / (2 cos AOI) (DERIVED).
Each half trip is symmetric under s -> -s (graded mirror G on the vertical axis s = 0):
  half 1: X -> P1 = d0 (cos t, +sin t) -> P2 -> G1 -> P2' -> P1' -> X  (exits arm A, returns on B)
  half 2: X -> P3 = d0 (cos t, -sin t) -> P4 -> G2 -> P4' -> P3' -> X  (exits arm B, returns on A)
Both are folded BELOW the atom plane, so nothing of the ring is above z = d0 sin t: the space
above the atoms stays free for the imaging objective. The fold directions alpha1, alpha2 (from
+s, degrees) and the distances d0 and L (P2 -> G) are chosen by a search that maximises the
smallest clearance between mirror faces, beams and foreign mirrors (in V), and against the
in-plane rings in 3D (results/cavity_layout/layout.json), subject to AOI <= max_aoi_deg and the
cell keep-out. Writes vertical_layout.json and vertical_layout.png.
"""

import itertools
import json
import sys

import numpy as np

from metatrap.cavity import budget as B
from metatrap.runio import REPO_ROOT, load_config, run_dir, write_json, write_params

D = np.array([1.0, 1.0, 0.0]) / np.sqrt(2)
U = np.array([-1.0, 1.0, 0.0]) / np.sqrt(2)
Z = np.array([0.0, 0.0, 1.0])


def to3(p):
    return p[0] * D + p[1] * Z


def unit(v):
    return v / np.linalg.norm(v)


def ring_v(theta, d0, f12, L, a1, a2, g1=-1, g2=-1, L2=None):
    """Mirror positions (mm, in V) and the two closed half-trip paths; None if G cannot reach s = 0.
    g = -1 puts G below P2, +1 above; L2 is the P4 -> G2 distance of half 2 (default L)."""
    c, s = np.cos(theta), np.sin(theta)
    out = {}
    for half, (sgn, a, gs, LL) in {1: (1.0, a1, g1, L), 2: (-1.0, a2, g2, L if L2 is None else L2)}.items():
        P1 = d0 * np.array([c, sgn * s])
        P2 = P1 + f12 * np.array([np.cos(np.radians(a)), np.sin(np.radians(a))])
        if abs(P2[0]) > LL:
            return None
        G = np.array([0.0, P2[1] + gs * np.sqrt(LL**2 - P2[0] ** 2)])
        names = ["P1", "P2", "G", "P2'", "P1'"] if half == 1 else ["P3", "P4", "G", "P4'", "P3'"]
        pts = [np.zeros(2), P1, P2, G, P2 * [-1, 1], P1 * [-1, 1], np.zeros(2)]
        out[half] = {"names": [n if n != "G" else f"G{half}" for n in names], "pts": pts}
    return out


def mirrors_of(r):
    """Per mirror: position, unit normal, AOI (deg)."""
    m = {}
    for half, h in r.items():
        pts = h["pts"]
        for i, n in enumerate(h["names"]):
            p_in, p, p_out = pts[i], pts[i + 1], pts[i + 2]
            a, b = unit(p - p_in), unit(p_out - p)
            nrm = unit(b - a)
            aoi = 0.5 * np.degrees(np.arccos(np.clip(np.dot(-a, b), -1, 1)))
            m[n] = {"p": p, "n": nrm, "aoi": float(aoi), "half": half}
    return m


def seg_seg(a, b, c, d, n=41):
    t = np.linspace(0, 1, n)
    P = a[None, :] + t[:, None] * (b - a)[None, :]
    Q = c[None, :] + t[:, None] * (d - c)[None, :]
    return float(np.min(np.linalg.norm(P[:, None, :] - Q[None, :, :], axis=2)))


def face(mi, h_r):
    t = np.array([-mi["n"][1], mi["n"][0]])
    return mi["p"] - h_r * t, mi["p"] + h_r * t


def evaluate(r, cfg, inplane):
    m = mirrors_of(r)
    lim = cfg["limits"]
    hr = {k: (cfg["mirror_half_height_mm"]["primary"] if k[:2] in ("P1", "P3") else cfg["mirror_half_height_mm"]["small"]) for k in m}
    wu = {k: (cfg["mirror_half_width_u_mm"]["primary"] if k[:2] in ("P1", "P3") else cfg["mirror_half_width_u_mm"]["small"]) for k in m}
    faces = {k: face(v, hr[k]) for k, v in m.items()}
    segs = []
    for half, h in r.items():
        for i in range(6):
            segs.append((h["names"][i - 1] if i > 0 else "X", h["names"][i] if i < 5 else "X", h["pts"][i], h["pts"][i + 1], i in (0, 5)))
    beam = cfg["beam_half_in_plane_mm"]
    c_mm, c_bm = np.inf, np.inf
    keys = list(m)
    for k1, k2 in itertools.combinations(keys, 2):
        c_mm = min(c_mm, seg_seg(*faces[k1], *faces[k2]))
    for e0, e1, a, b, arm in segs:
        for k in keys:
            if k in (e0, e1):
                continue
            c_bm = min(c_bm, seg_seg(a, b, *faces[k]) - beam)
    if min(c_mm, c_bm) < cfg.get("min_clearance_mm", 0.0):  # cheap reject before the 3D checks
        return {"ok": False, "score": float(min(c_mm, c_bm))}
    # Cell keep-out: folded segments (not the arms through X) stay away from X.
    c_cell = min(seg_seg(a, b, np.zeros(2), np.zeros(2) + 1e-9) for _, _, a, b, arm in segs if not arm)
    # Objective space: nothing above z_max_mm.
    z_top = max(max(v["p"][1] + hr[k] for k, v in m.items()), max(max(a[1], b[1]) for _, _, a, b, _ in segs))
    # In-plane rings (3D): their mirrors as spheres of the aperture radius; ours as rectangles h_r x w_u.
    c_ip = np.inf
    for ring in inplane["rings"].values():
        for kk, cxy in ring["mirrors_mm"].items():
            rad = 5.0 if kk.startswith("P1") else 2.5
            cen = np.array([cxy[0], cxy[1], 0.0])
            for k, v in m.items():
                t = np.array([-v["n"][1], v["n"][0]])
                t3, n3 = t[0] * D + t[1] * Z, v["n"][0] * D + v["n"][1] * Z
                rel = cen - to3(v["p"])
                q = np.array([np.clip(rel @ t3, -hr[k], hr[k]), np.clip(rel @ U, -wu[k], wu[k]), 0.0])
                c_ip = min(c_ip, float(np.linalg.norm(rel - (q[0] * t3 + q[1] * U))) - rad)
    # In-plane ring beams (z = 0 sheets, half-width 2.9 mm, half-height 0.6 mm) under/over our mirrors.
    c_ipb = np.inf
    for ring in inplane["rings"].values():
        mm_ = ring["mirrors_mm"]
        for half in ("", "_lower"):
            g = "G1" if half == "" else "G2"
            seq = [np.zeros(2), mm_["P1" + half], mm_["P2" + half], mm_[g + half], mm_["P2'" + half], mm_["P1'" + half], np.zeros(2)]
            for a2, b2 in zip(seq[:-1], seq[1:]):
                a3, b3 = np.array([*a2, 0.0]), np.array([*b2, 0.0])
                for k, v in m.items():
                    if np.linalg.norm(v["p"]) < 1e-9:
                        continue
                    t = np.array([-v["n"][1], v["n"][0]])
                    t3 = t[0] * D + t[1] * Z
                    for fa in np.linspace(-1, 1, 5):
                        for fb in np.linspace(-1, 1, 9):
                            q = to3(v["p"]) + fa * hr[k] * t3 + fb * wu[k] * U
                            tt = np.clip(np.dot(q - a3, b3 - a3) / np.dot(b3 - a3, b3 - a3), 0, 1)
                            dv = q - (a3 + tt * (b3 - a3))
                            horiz = np.hypot(dv[0], dv[1])
                            c_ipb = min(c_ipb, max(horiz - 2.9, abs(dv[2]) - 0.6))
    aoi_max = max(v["aoi"] for v in m.values())
    depth = min(v["p"][1] for v in m.values())
    ok = aoi_max <= lim["max_aoi_deg"] and c_cell >= lim["cell_keepout_mm"] and z_top <= lim["z_top_max_mm"] and depth >= -lim["max_depth_mm"]
    score = min(c_mm, c_bm, c_ip, c_ipb)
    return {"ok": bool(ok), "score": float(score), "mirror_mirror_mm": float(c_mm), "beam_mirror_mm": float(c_bm),
            "cell_keepout_mm": float(c_cell), "z_top_mm": float(z_top), "inplane_mirror_mm": float(c_ip),
            "inplane_beam_mm": float(c_ipb), "aoi_max_deg": float(aoi_max), "depth_mm": float(depth),
            "aoi_deg": {k: v["aoi"] for k, v in m.items()}}


def check_chosen(cfg, out, th, inplane):
    """Re-evaluate the chosen fold geometry with other mirror strip sizes and beam heights (no search)."""
    ch = json.loads((REPO_ROOT / cfg["check_chosen"]["layout_json"]).read_text())["chosen"]
    f12 = ch["f1_mm"] + ch["f2_mm"]
    h1, h2 = ch["half1"], ch["half2"]
    r = ring_v(th, ch["d0_mm"], f12, h1["L_mm"], h1["alpha_deg"], h2["alpha_deg"], h1["G_side"], h2["G_side"], L2=h2["L_mm"])
    rep = {}
    for case in cfg["check_chosen"]["cases"]:
        c = dict(cfg, mirror_half_height_mm=case["mirror_half_height_mm"], beam_half_in_plane_mm=case["beam_half_in_plane_mm"], min_clearance_mm=0.0)
        ev = evaluate(r, c, inplane)
        rep[case["name"]] = dict(case, **{k: v for k, v in ev.items() if k != "aoi_deg"})
        print(case["name"], {k: (round(v, 1) if isinstance(v, float) else v) for k, v in rep[case["name"]].items() if k not in ("name",)}, flush=True)
    write_json(out / "check.json", rep)


def main(cfg_path):
    cfg = load_config(cfg_path)
    out = run_dir(cfg["run_id"])
    write_params(out, cfg)
    if cfg.get("check_chosen"):
        th = B.crossing_half_angle(cfg["vertical_period_um"], cfg["wavelength_nm"] * 1e-3)
        inplane = json.loads((REPO_ROOT / cfg["inplane_layout_json"]).read_text())["layout"]
        check_chosen(cfg, out, th, inplane)
        return
    th = B.crossing_half_angle(cfg["vertical_period_um"], cfg["wavelength_nm"] * 1e-3)
    inplane = json.loads((REPO_ROOT / cfg["inplane_layout_json"]).read_text())["layout"]
    best = []
    sc = cfg["search"]
    counts = {}
    for f1 in sc.get("f1_mm", [cfg["f1_mm"]]):
        b_, c_ = search_f1(cfg, th, inplane, f1)
        best += b_
        counts[str(f1)] = c_
    # Selection: among layouts with at least min_clearance_mm, the longest graded-mirror distance L
    # (mode discrimination grows with L / M^2), then the largest clearance.
    good = [b for b in best if b[0] >= cfg["min_clearance_mm"]]
    good.sort(key=lambda b: (-b[2][1], -b[0]))
    best = good + sorted([b for b in best if b[0] < cfg["min_clearance_mm"]], key=lambda b: -b[0])
    write_report(cfg, out, th, inplane, best, counts)


def search_f1(cfg, th, inplane, f1):
    f12 = f1 * (1 + cfg["mag"])
    best = []
    sc = cfg["search"]
    # Each half separately first (cheap AOI / objective / cell filter), then pairs.
    halves = {1: [], 2: []}
    for d0 in sc["d0_mm"]:
        for L in sc["L_mm"]:
            for a in np.arange(*sc["alpha_deg"]):
                for gs in (-1, 1):
                    for h in (1, 2):
                        r = ring_v(th, d0, f12, L, a, a, gs, gs)
                        if r is None:
                            continue
                        m = {k: v for k, v in mirrors_of(r).items() if v["half"] == h}
                        if max(v["aoi"] for v in m.values()) > cfg["limits"]["max_aoi_deg"]:
                            continue
                        pts = r[h]["pts"]
                        if max(p[1] for p in pts) > cfg["limits"]["z_top_max_mm"] - 2.0:
                            continue
                        if min(seg_seg(pts[i], pts[i + 1], np.zeros(2), np.zeros(2) + 1e-9) for i in range(1, 5)) < cfg["limits"]["cell_keepout_mm"]:
                            continue
                        halves[h].append((d0, L, float(a), gs))
    print(f"f1 {f1}: feasible halves", len(halves[1]), len(halves[2]), flush=True)
    for h1 in halves[1]:
        for h2 in halves[2]:
            if h1[0] != h2[0] or (cfg.get("equal_L") and h1[1] != h2[1]):
                continue
            d0 = h1[0]
            r = ring_v(th, d0, f12, h1[1], h1[2], h2[2], h1[3], h2[3], L2=h2[1])
            if r is None:
                continue
            ev = evaluate(r, cfg, inplane)
            if ev["ok"]:
                best.append((ev["score"], d0, h1, h2, ev))
    return [(b[0], b[1], b[2], b[3], dict(b[4], f1_mm=f1)) for b in best], {k: len(v) for k, v in halves.items()}


def write_report(cfg, out, th, inplane, best, counts):
    rep = {"theta_z_deg": float(np.degrees(th)), "candidates_found": len(best), "feasible_halves_per_f1": counts,
           "best": [{"score_mm": b[0], "d0_mm": b[1], "half1": {"L_mm": b[2][1], "alpha_deg": b[2][2], "G_side": b[2][3]},
                     "half2": {"L_mm": b[3][1], "alpha_deg": b[3][2], "G_side": b[3][3]}, **b[4]} for b in best[:10]]}
    if best:
        _, d0, h1, h2, ev = best[0]
        f12 = ev["f1_mm"] * (1 + cfg["mag"])
        r = ring_v(th, d0, f12, h1[1], h1[2], h2[2], h1[3], h2[3], L2=h2[1])
        m = mirrors_of(r)
        Ls = (h1[1], h2[1])
        rep["chosen"] = {"f1_mm": ev["f1_mm"], "f2_mm": ev["f1_mm"] * cfg["mag"], "d0_mm": d0, "half1": {"L_mm": h1[1], "alpha_deg": h1[2], "G_side": h1[3]},
                         "half2": {"L_mm": h2[1], "alpha_deg": h2[2], "G_side": h2[3]}, **ev,
                         "half_trip_paths_mm": [d0 + f12 + x for x in Ls], "round_trip_mm": 2 * sum(d0 + f12 + x for x in Ls),
                         "fsr_MHz": 299792458.0 / (2 * sum(d0 + f12 + x for x in Ls) * 1e-3) / 1e6,
                         "mirrors_V_mm": {k: v["p"].tolist() for k, v in m.items()},
                         "mirrors_3d_mm": {k: to3(v["p"]).tolist() for k, v in m.items()}}
        print(json.dumps({k: v for k, v in rep["chosen"].items() if k not in ("mirrors_3d_mm",)}, indent=1, default=float))
    write_json(out / "vertical_layout.json", rep)
    for b in rep["best"][:6]:
        print({k: (round(v, 1) if isinstance(v, float) else v) for k, v in b.items() if k not in ("aoi_deg",)})


if __name__ == "__main__":
    main(sys.argv[1])
