"""Optical access around the cavity-lattice atoms: MOT beams, imaging, optical pumping / Raman beams,
and a grating-MOT alternative.
Usage:
    python scripts/cavity_access.py configs/cavity_access.yaml

Obstacles (mm, atoms at the origin, z up): every mirror of the two in-plane rings
(results/cavity_layout/layout.json) and of the vertical ring (results/cavity_vertical_layout/
vertical_layout.json), each as a sphere of its aperture radius plus a mount margin; the imaging
objective above the atoms as a vertical cylinder (radius, from the height of its first plate upwards);
the imaging light cone is not an obstacle (light may cross light). Beams through the atoms are cylinders of radius r_b along a line
(both directions, to a length R). The clearance of a direction is the smallest gap between the beam
and any obstacle (negative = blocked). Optical windows are assumed wherever needed (in-vacuum mirror
design, custom chamber; [CONFIRM?]). Ring beams are light only and may be crossed.
Outputs: an access map (polar angle from +z, azimuth) for the MOT and the pumping/Raman beam radii;
the best orthogonal six-beam MOT triad (search over orientations); the best pumping/Raman
directions; checks of a grating-MOT chip below the atoms with its input through the objective.
Writes access.json and access.png.
"""

import json
import sys

import numpy as np

from metatrap.runio import REPO_ROOT, load_config, run_dir, write_json, write_params

D = np.array([1.0, 1.0, 0.0]) / np.sqrt(2)
Z = np.array([0.0, 0.0, 1.0])


def obstacles(cfg):
    obs = []
    lay = json.loads((REPO_ROOT / cfg["inplane_layout_json"]).read_text())["layout"]
    for rn, ring in lay["rings"].items():
        for k, p in ring["mirrors_mm"].items():
            r = (5.0 if k.startswith("P1") else 2.5) + cfg["mount_margin_mm"]
            obs.append({"name": f"{rn}:{k}", "c": np.array([p[0], p[1], 0.0]), "r": r})
    vl = json.loads((REPO_ROOT / cfg["vertical_layout_json"]).read_text())["chosen"]
    vcfg = load_config(REPO_ROOT / cfg["vertical_layout_config"])
    for k, p in vl["mirrors_3d_mm"].items():
        prim = k[:2] in ("P1", "P3")
        hr = vcfg["mirror_half_height_mm"]["primary" if prim else "small"]
        wu = vcfg["mirror_half_width_u_mm"]["primary" if prim else "small"]
        obs.append({"name": f"vertical:{k}", "c": np.array(p), "r": float(np.hypot(hr, wu)) + cfg["mount_margin_mm"]})
    return obs, lay, vl


def line_clearance(e, rb, obs, ob, R):
    """Smallest gap between a beam (radius rb along +-e through the origin, length R each way) and obstacles."""
    c = np.inf
    for o in obs:
        t = np.clip(o["c"] @ e, -R, R)
        c = min(c, float(np.linalg.norm(o["c"] - t * e)) - o["r"] - rb)
    # Objective: cylinder radius ro above z0, and the cone of the field corners below z0.
    ts = np.linspace(-R, R, 801)
    pts = ts[:, None] * e[None, :]
    rho = np.hypot(pts[:, 0], pts[:, 1])
    z = pts[:, 2]
    above = z >= ob["z0_mm"]
    if np.any(above):
        c = min(c, float(np.min(rho[above] - ob["radius_mm"])) - rb)
    # The imaging light cone itself is not an obstacle: every beam through the atoms crosses it there.
    return c


def direction(theta, phi):
    return np.array([np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)])


def rot(a, b, g):
    """Rotation matrix from z-y-z Euler angles (rad)."""
    def rz(t):
        return np.array([[np.cos(t), -np.sin(t), 0], [np.sin(t), np.cos(t), 0], [0, 0, 1]])

    def ry(t):
        return np.array([[np.cos(t), 0, np.sin(t)], [0, 1, 0], [-np.sin(t), 0, np.cos(t)]])

    return rz(a) @ ry(b) @ rz(g)


def main(cfg_path):
    cfg = load_config(cfg_path)
    out = run_dir(cfg["run_id"])
    write_params(out, cfg)
    obs, lay, vl = obstacles(cfg)
    ob = cfg["objective"]
    R = cfg["beam_length_mm"]
    rep = {"assumptions": {"objective": ob, "mount_margin_mm": cfg["mount_margin_mm"], "beam_length_mm": R}}
    # Access maps.
    th = np.radians(np.arange(0.5, 90.0, 1.0))  # lines through the atoms: theta in [0, 90] covers all lines
    ph = np.radians(np.arange(0.0, 360.0, 1.0))
    maps = {}
    for label, rb in cfg["beam_radii_mm"].items():
        m = np.array([[line_clearance(direction(t, p), rb, obs, ob, R) for p in ph] for t in th])
        maps[label] = m
        free = m > 0
        rep[f"map_{label}"] = {"beam_radius_mm": rb, "free_fraction_of_lines": float(free.mean()),
                               "best": [{"theta_deg": float(np.degrees(th[i])), "phi_deg": float(np.degrees(ph[j])), "clearance_mm": float(m[i, j])}
                                        for i, j in zip(*np.unravel_index(np.argsort(-m, axis=None)[:5], m.shape))]}
    # Standard MOT (one vertical + two horizontal pairs) and the best orthogonal triad.
    rb = cfg["beam_radii_mm"]["mot"]
    std = [line_clearance(e, rb, obs, ob, R) for e in (np.array([1.0, 0, 0]), np.array([0, 1.0, 0]), Z)]
    std45 = [line_clearance(e, rb, obs, ob, R) for e in (D, np.array([-1.0, 1.0, 0]) / np.sqrt(2), Z)]
    rep["mot_standard_xyz"] = {"clearance_mm": std}
    rep["mot_standard_diagonals"] = {"clearance_mm": std45}
    best = None
    for a in np.radians(np.arange(0, 120, 3)):
        for b in np.radians(np.arange(0, 91, 3)):
            for g in np.radians(np.arange(0, 90, 3)):
                Rm = rot(a, b, g)
                cl = [line_clearance(Rm[:, i], rb, obs, ob, R) for i in range(3)]
                if best is None or min(cl) > best[0]:
                    best = (min(cl), Rm, cl)
    Rm = best[1]
    rep["mot_best_triad"] = {"min_clearance_mm": best[0], "clearance_mm": best[2],
                             "axes": [{"theta_deg": float(np.degrees(np.arccos(abs(Rm[2, i])))), "phi_deg": float(np.degrees(np.arctan2(Rm[1, i], Rm[0, i]) % (2 * np.pi)))} for i in range(3)]}
    # Extra pumping / Raman lines: greedy choice of the clearest lines at least min_sep_deg from the
    # MOT axes and from each other (lines, i.e. +-e).
    m = maps["pumping_raman"]
    taken = [Rm[:, i] for i in range(3)]
    extra = []
    order = np.argsort(-m, axis=None)
    for idx in order:
        i, j = np.unravel_index(idx, m.shape)
        if m[i, j] <= 0 or len(extra) >= cfg["n_extra_lines"]:
            break
        e = direction(th[i], ph[j])
        if all(np.degrees(np.arccos(min(1.0, abs(e @ t)))) >= cfg["min_sep_deg"] for t in taken):
            taken.append(e)
            extra.append({"theta_deg": float(np.degrees(th[i])), "phi_deg": float(np.degrees(ph[j])), "clearance_mm": float(m[i, j]),
                          "clearance_mot_size_mm": float(maps["mot"][i, j])})
    rep["extra_lines_pumping_raman"] = extra
    # Grating MOT: chip below the atoms, input beam along -z from above (through the objective).
    gm = cfg["grating_mot"]
    chip_z, half = gm["chip_z_mm"], gm["chip_half_mm"]

    def seg_to_chip(a, b):
        """Smallest distance from a segment (3D) to the chip square (z = chip_z, |x|,|y| <= half)."""
        t = np.linspace(0, 1, 401)[:, None]
        p = a[None, :] + t * (b - a)[None, :]
        q = np.stack([np.clip(p[:, 0], -half, half), np.clip(p[:, 1], -half, half), np.full(len(p), chip_z)], axis=1)
        return float(np.min(np.linalg.norm(p - q, axis=1)))

    dmin = {}
    for rn, ring in lay["rings"].items():
        mm = ring["mirrors_mm"]
        for hs in ("", "_lower"):
            g = "G1" if hs == "" else "G2"
            seq = [[0, 0], mm["P1" + hs], mm["P2" + hs], mm[g + hs], mm["P2'" + hs], mm["P1'" + hs], [0, 0]]
            seq = [np.array([p[0], p[1], 0.0]) for p in seq]
            dmin[f"{rn}{hs}"] = min(seg_to_chip(a, b) for a, b in zip(seq[:-1], seq[1:]))
    mv = {k: np.array(v) for k, v in vl["mirrors_3d_mm"].items()}
    for h, seq in {1: ["X", "P1", "P2", "G1", "P2'", "P1'", "X"], 2: ["X", "P3", "P4", "G2", "P4'", "P3'", "X"]}.items():
        pts = [np.zeros(3) if k == "X" else mv[k] for k in seq]
        dmin[f"vertical_half{h}"] = min(seg_to_chip(a, b) for a, b in zip(pts[:-1], pts[1:]))
    # Diffracted beams: from the chip, upwards at the diffraction angle, in n_beams azimuths.
    diff = []
    for k in range(gm["n_beams"]):
        phi = 2 * np.pi * k / gm["n_beams"] + np.radians(gm["azimuth0_deg"])
        e = np.array([np.sin(np.radians(gm["diffraction_deg"])) * np.cos(phi), np.sin(np.radians(gm["diffraction_deg"])) * np.sin(phi), np.cos(np.radians(gm["diffraction_deg"]))])
        start = np.array([-np.cos(phi), -np.sin(phi), 0.0]) * half * 0.5 + np.array([0, 0, chip_z])
        c = np.inf
        for o in obs:
            t = max(0.0, (o["c"] - start) @ e)
            c = min(c, float(np.linalg.norm(o["c"] - (start + t * e))) - o["r"] - gm["diffracted_half_width_mm"])
        ts = np.linspace(0, R, 801)
        p = start[None, :] + ts[:, None] * e[None, :]
        above = p[:, 2] >= ob["z0_mm"]
        if np.any(above):
            c = min(c, float(np.min(np.hypot(p[above, 0], p[above, 1]) - ob["radius_mm"])) - gm["diffracted_half_width_mm"])
        diff.append({"azimuth_deg": float(np.degrees(phi)), "clearance_mm": c})
    rep["grating_mot"] = {**gm, "chip_to_ring_beams_mm": dmin, "diffracted_beams": diff,
                          "retro_beam_through_chip_hole_to_first_vertical_mirror_mm": float(min(-v[2] for k, v in mv.items() if abs(v[0]) < 1e-6 and abs(v[1]) < 1e-6 and v[2] < 0) + chip_z) if any(abs(v[0]) < 1e-6 and abs(v[1]) < 1e-6 for v in mv.values()) else None}
    write_json(out / "access.json", rep)
    print(json.dumps({k: v for k, v in rep.items() if not k.startswith("map_")}, indent=1, default=float)[:3000])
    for k in maps:
        print(k, rep[f"map_{k}"]["free_fraction_of_lines"], rep[f"map_{k}"]["best"][:3])
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, axs = plt.subplots(1, len(maps), figsize=(5.5 * len(maps), 3.6), facecolor="#fcfcfb")
        for ax, (k, m) in zip(np.atleast_1d(axs), maps.items()):
            im = ax.imshow(np.clip(m, -30, 30), extent=[0, 360, 90, 0], aspect="auto", cmap="RdBu", vmin=-30, vmax=30)
            ax.contour(np.degrees(ph), np.degrees(th), m, levels=[0.0], colors="k", linewidths=0.8)
            for ax_ in rep["mot_best_triad"]["axes"]:
                for dphi in (0.0, 180.0):
                    ax.plot((ax_["phi_deg"] + dphi) % 360, ax_["theta_deg"], "k*", ms=9)
            if k == "pumping_raman":
                for ex in rep["extra_lines_pumping_raman"]:
                    ax.plot(ex["phi_deg"], ex["theta_deg"], "o", mfc="none", mec="k", ms=7)
            ax.set_xlabel("azimuth (deg; 0 = +x, 45 = vertical-ring plane)", fontsize=7)
            ax.set_ylabel("angle of the line from vertical (deg)", fontsize=7)
            ax.set_title(f"clearance of a {cfg['beam_radii_mm'][k]:.0f} mm-radius beam ({k}); black: 0 mm", fontsize=8)
            ax.tick_params(labelsize=7)
            fig.colorbar(im, ax=ax, label="clearance (mm)")
        fig.tight_layout()
        fig.savefig(out / "access.png", dpi=150)
    except Exception as exc:
        print("plot failed", exc)


if __name__ == "__main__":
    main(sys.argv[1])
