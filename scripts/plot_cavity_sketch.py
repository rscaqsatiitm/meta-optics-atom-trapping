"""Design sketches for docs/cavity_concept.md (Figures 3 and 4).
Usage:
    python scripts/plot_cavity_sketch.py

Figure 3 (results/cavity_layout/sketch_layout.png):
  (a) top view to scale, mirror positions from results/cavity_layout/layout.json;
  (b) the crossing to scale: the two sheets of one axis and the 10^6-site array, with the computed
      two-axis lattice intensity (plane-wave model) in an inset;
  (c) side view, schematic (angles and thicknesses exaggerated).
Figure 4 (results/cavity_layout/sketch_ring.png):
  (a) beam along one unfolded half trip of the adopted ring (Fox-Li, design telK_M0.2_d70):
      intensity normalised per plane, and the profiles at the crossing and on G;
  (b) optics and control chain (block diagram; powers from results/cavity_budget/budget.json).
"""

import json
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Polygon, Rectangle  # noqa: E402

from metatrap.cavity import budget as B  # noqa: E402
from metatrap.cavity import lattice as LA  # noqa: E402
from metatrap.runio import REPO_ROOT, load_config  # noqa: E402

sys.path.insert(0, str(REPO_ROOT / "scripts"))
import cavity_foxli as S  # noqa: E402

BLUE, ORANGE, AQUA, YELLOW, MAGENTA, GREEN, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6"
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
OUT = REPO_ROOT / "results" / "cavity_layout"


def style(ax, grid=True):
    ax.set_facecolor(SURFACE)
    if grid:
        ax.grid(True, color=GRID, lw=0.5)
    for s in ax.spines.values():
        s.set_color(INK2)
        s.set_linewidth(0.6)
    ax.tick_params(colors=INK2, labelsize=7)


def top_view(ax, lay):
    style(ax)
    for (name, rg), col in zip(lay["rings"].items(), (BLUE, ORANGE)):
        m = rg["mirrors_mm"]
        for half in ("", "_lower"):
            g = "G1" if half == "" else "G2"
            seq = [[0, 0], m["P1" + half], m["P2" + half], m[g + half], m["P2'" + half], m["P1'" + half], [0, 0]]
            seq = np.array(seq)
            ax.plot(seq[:, 0], seq[:, 1], "-", color=col, lw=1.0)
        for k, p in m.items():
            r = 5.0 if k.startswith("P1") else 2.5
            ax.add_patch(plt.Circle(p, r, fc=SURFACE, ec=col, lw=1.0, zorder=3))
        label = "y-lattice ring" if name.startswith("y") else "x-lattice ring"
        if name.startswith("y"):
            ax.annotate(f"{label}: G1 (top), G2 (bottom)\ngraded mirrors", m["G1"], xytext=(12, -14), textcoords="offset points",
                        fontsize=6.5, color=INK, arrowprops=dict(arrowstyle="-", color=INK2, lw=0.5))
        else:
            ax.annotate(f"{label}:\nG1, G2 at left, right", m["G1"], xytext=(4, 22), textcoords="offset points",
                        fontsize=6.5, color=INK, arrowprops=dict(arrowstyle="-", color=INK2, lw=0.5))
    # Labels on the y-lattice ring.
    m = lay["rings"]["y_lattice"]["mirrors_mm"]
    ax.annotate("P1 (f 100 mm)", m["P1"], xytext=(8, -12), textcoords="offset points", fontsize=6.5, color=INK)
    ax.annotate("P2 (f 20 mm)", m["P2"], xytext=(-20, 12), textcoords="offset points", fontsize=6.5, color=INK,
                arrowprops=dict(arrowstyle="-", color=INK2, lw=0.5))
    # Diagonal corridors for the vertical-lattice ring.
    for sgn in (1, -1):
        ax.plot([-150, 150], [-150 * sgn, 150 * sgn], ls=(0, (4, 3)), color=AQUA, lw=0.9)
    vl_path = OUT.parent / "cavity_vertical_layout" / "vertical_layout.json"
    if vl_path.exists():
        ch = json.loads(vl_path.read_text())["chosen"]
        d = np.array([1, 1]) / np.sqrt(2)
        sv = [p[0] for p in ch["mirrors_V_mm"].values()]
        for sv_ in sv:
            q = sv_ * d
            ax.plot([q[0] - 6 * d[1] * -1, q[0] + 6 * d[1] * -1], [q[1] - 6 * d[0], q[1] + 6 * d[0]], "-", color=AQUA, lw=2.0)
        ax.annotate("vertical-lattice ring: vertical plane on this\ndiagonal, folded below the atom plane (Figure 5)", (min(sv) * d), xytext=(-140, -142),
                    fontsize=6.5, color=INK, arrowprops=dict(arrowstyle="-", color=INK2, lw=0.5))
    else:
        ax.annotate("vertical-lattice ring\n(vertical plane, along a diagonal)", (118, -118), xytext=(30, -142), fontsize=6.5, color=INK,
                    arrowprops=dict(arrowstyle="-", color=INK2, lw=0.5))
    ax.add_patch(Rectangle((-15, -15), 30, 30, fill=False, ec=INK2, lw=0.8, ls="--"))
    ax.annotate("cell, 30 mm [CONFIRM?]", (15, -15), xytext=(22, -26), fontsize=6.5, color=INK2, arrowprops=dict(arrowstyle="-", color=INK2, lw=0.5))
    ax.add_patch(Rectangle((-2, -2), 4, 4, fc=VIOLET, ec="none"))
    # Input coupler and camera (transmission of a fold mirror), schematic arrows.
    p_in = np.array(m["P1'_lower"])
    ax.add_patch(FancyArrowPatch(p_in + [-35, -25], p_in + [-4, -2], arrowstyle="-|>", mutation_scale=8, color=GREEN, lw=1.0))
    ax.annotate("input: flat-top sheet\nfrom meta-optic converter", p_in + [-35, -25], xytext=(-20, -26), textcoords="offset points", fontsize=6.5, color=INK)
    p_c = np.array(m["P1'"])
    ax.add_patch(FancyArrowPatch(p_c + [-3, 2], p_c + [-30, 22], arrowstyle="-|>", mutation_scale=8, color=MAGENTA, lw=1.0))
    ax.annotate("leak -> camera\n(mode profile)", p_c + [-30, 22], xytext=(-30, 6), textcoords="offset points", fontsize=6.5, color=INK)
    ax.set_xlim(-150, 150)
    ax.set_ylim(-150, 150)
    ax.set_aspect("equal")
    ax.set_xlabel("x (mm)", fontsize=7, color=INK)
    ax.set_ylabel("y (mm)", fontsize=7, color=INK)
    ax.set_title("(a) Top view, to scale: two in-plane rings, mirrors in vacuum", fontsize=8, color=INK, loc="left")


def crossing(ax, fox):
    style(ax)
    lam = 1.075
    th = B.crossing_half_angle(4.0, lam)
    W = fox["equivalent_width_um"] * 1e-3
    L = 14.0
    for sgn, col in ((1, BLUE), (-1, BLUE)):
        d = np.array([np.cos(th), sgn * np.sin(th)])
        n = np.array([-d[1], d[0]])
        c = [(-L * d - W / 2 * n), (L * d - W / 2 * n), (L * d + W / 2 * n), (-L * d + W / 2 * n)]
        ax.add_patch(Polygon(np.array(c), closed=True, fc=col, alpha=0.18, ec=col, lw=0.8))
        ax.add_patch(FancyArrowPatch(-6 * d + 0 * n, 6 * d, arrowstyle="-|>", mutation_scale=7, color=col, lw=0.8))
    ax.add_patch(Rectangle((-2, -2), 4, 4, fc="none", ec=VIOLET, lw=1.2))
    ax.annotate("10^6 sites, 4 mm", (-2, 2), xytext=(-9, 7.2), fontsize=6.5, color=INK, arrowprops=dict(arrowstyle="-", color=INK2, lw=0.5))
    ax.annotate(f"arm A, arm B: flat-top sheets {W:.1f} mm wide (equivalent), crossing at 2 theta = {2 * np.degrees(th):.1f} deg", (-13.5, -9.2), fontsize=6.5, color=INK)
    ax.set_xlim(-14, 14)
    ax.set_ylim(-10, 10)
    ax.set_aspect("equal")
    ax.set_xlabel("x (mm)", fontsize=7, color=INK)
    ax.set_ylabel("y (mm)", fontsize=7, color=INK)
    ax.set_title("(b) The crossing of one axis, to scale", fontsize=8, color=INK, loc="left")
    ins = ax.inset_axes([0.70, 0.68, 0.27, 0.30])
    k = 2 * np.pi / lam
    c, s = np.cos(th), np.sin(th)
    fy = [(1.0, k * np.array([c, s, 0])), (1.0, k * np.array([c, -s, 0]))]
    fx = [(1.0, k * np.array([s, c, 0])), (1.0, k * np.array([-s, c, 0]))]
    g = np.linspace(-6, 6, 241)
    X, Y = np.meshgrid(g, g)
    pts = np.stack([X, Y, np.zeros_like(X)], axis=-1)
    inten = LA.plane_wave_intensity(pts, [fx, fy])
    ins.imshow(inten, extent=[-6, 6, -6, 6], origin="lower", cmap="Blues")
    ins.set_xticks([-4, 0, 4])
    ins.set_yticks([-4, 0, 4])
    ins.tick_params(labelsize=5.5, colors=INK2)
    ins.set_title("two axes, 4 um pitch (um)", fontsize=6, color=INK)


def side_view(ax):
    ax.set_facecolor(SURFACE)
    ax.set_xlim(-10, 10)
    ax.set_ylim(-6, 9)
    ax.axis("off")
    ax.add_patch(Rectangle((-6, -3.5), 12, 6.5, fill=False, ec=INK2, lw=0.8, ls="--"))
    ax.text(-5.8, -3.3, "glass cell / vacuum", fontsize=6.5, color=INK2)
    ax.plot([-9.5, 9.5], [0, 0], color=BLUE, lw=3.0, alpha=0.5)
    ax.text(-9.6, 1.6, "in-plane sheets (x and y rings),\n2 w_z = 150 um (thickness exaggerated)", fontsize=6.5, color=INK)
    for sgn in (1, -1):
        ax.plot([-9.5, 9.5], [-9.5 * np.tan(sgn * 0.12), 9.5 * np.tan(sgn * 0.12)], color=AQUA, lw=1.4)
    ax.text(1.0, -3.0, "vertical-lattice beams, +-3.85 deg (angle exaggerated);\n8 um planes, one plane loaded", fontsize=6.5, color=INK)
    ax.add_patch(Rectangle((-0.25, -0.12), 0.5, 0.24, fc=VIOLET, ec="none", zorder=4))
    # Objective (meta-doublet) above, NA 0.5 cone.
    ax.add_patch(Polygon(np.array([[0, 0.15], [-3.3, 5.4], [3.3, 5.4]]), closed=True, fc=YELLOW, alpha=0.2, ec="none"))
    ax.add_patch(Rectangle((-3.6, 5.4), 7.2, 0.35, fc=YELLOW, ec=INK2, lw=0.6))
    ax.add_patch(Rectangle((-3.6, 6.3), 7.2, 0.35, fc=YELLOW, ec=INK2, lw=0.6))
    ax.text(3.9, 5.8, "852 nm meta-objective, NA 0.5\n(imaging, addressing)", fontsize=6.5, color=INK)
    ax.plot([-3.0, 3.0], [7.8, 7.8], color=INK2, lw=0.8)
    ax.text(-9.6, 7.6, "dichroic: rearrangement\ntweezers (AOD) in, camera out", fontsize=6.5, color=INK)
    ax.set_title("(c) Side view, schematic", fontsize=8, color=INK, loc="left")


def caustic(ax_img, ax_prof, cfg):
    v = cfg["variants"]["telK_M0.2_d70"]
    s = S.build(cfg, v)
    ring, F, g = s["ring"], s["F"], s["g"]
    p = ring.p_um
    zs = np.linspace(0, 2 * p, 600)
    stops = ring.stop_positions()
    u = F.copy()
    img = []
    z_prev = 0.0
    for z in zs:
        if z_prev < p <= z:  # apply G
            u = ring._segment(u, z_prev, p, 1) * ring._ap_g * np.exp(1j * ring.phase_g1)
            u = ring._segment(u, p, z, 1)
        else:
            u = ring._segment(u, z_prev, z, 1)
        z_prev = z
        i = np.abs(u) ** 2
        img.append(i / i.max())
    img = np.array(img).T
    sel = np.abs(g.x) <= 4500
    ax_img.imshow(img[sel], extent=[0, 2 * p / 1e3, g.x[sel][0] / 1e3, g.x[sel][-1] / 1e3], aspect="auto", origin="lower", cmap="Blues")
    names = ["P1", "P2", "P2'", "P1'"]
    for z, n in zip(stops, names):
        ax_img.axvline(z / 1e3, color=INK2, lw=0.6, ls="--")
        ax_img.text(z / 1e3 + 2, 3.9, n, fontsize=6.5, color=INK)
    ax_img.axvline(p / 1e3, color=ORANGE, lw=0.9)
    ax_img.text(p / 1e3 + 2, 3.9, "G1", fontsize=6.5, color=INK)
    ax_img.text(2, -4.2, "X (arm A)", fontsize=6.5, color=INK)
    ax_img.text(2 * p / 1e3 - 2, -4.2, "X (arm B)", fontsize=6.5, color=INK, ha="right")
    ax_img.set_xlabel("path along the half trip (mm)", fontsize=7, color=INK)
    ax_img.set_ylabel("transverse, in plane (mm)", fontsize=7, color=INK)
    ax_img.tick_params(labelsize=7, colors=INK2)
    ax_img.set_title("(a) One half trip of the ring, unfolded (Fox-Li; intensity normalised per plane)", fontsize=8, color=INK, loc="left")
    style(ax_prof)
    iX = np.abs(F) ** 2
    iG = np.abs(ring.to_mirror(F, 1)) ** 2
    ax_prof.plot(g.x / 1e3, iX / iX.max(), color=BLUE, lw=1.5, label="at the crossing (flat-top)")
    ax_prof.plot(g.x / 1e3, iG / iG.max(), color=ORANGE, lw=1.5, label="on the graded mirror G")
    h = s["h"] / 1e3
    ax_prof.axvspan(-h, h, color=VIOLET, alpha=0.08, lw=0)
    ax_prof.text(0, 0.45, "array range", fontsize=6.5, color=INK, ha="center")
    ax_prof.set_xlim(-4.5, 4.5)
    ax_prof.set_xlabel("transverse, in plane (mm)", fontsize=7, color=INK)
    ax_prof.set_ylabel("intensity (norm.)", fontsize=7, color=INK)
    ax_prof.set_ylim(0, 1.35)
    ax_prof.legend(fontsize=6.5, frameon=False, loc="upper center", ncol=1)


def chain(ax, bud):
    ax.set_facecolor(SURFACE)
    ax.axis("off")
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 40)
    c = [r for r in bud["cases"] if r["n_sites"] == 1e6 and r["depth_mK"] == 0.2 and r["match"] == "equal_depth"][0]
    pa, pv = c["circulating"]["inplane_per_axis_W"], c["circulating"]["vertical_W"]
    g = 300.0
    pin_a, pin_v = pa / (g * 0.95), pv / (g * 0.95)

    def box(x, y, w, h, txt, col):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.4,rounding_size=1.2", fc=SURFACE, ec=col, lw=1.0))
        ax.text(x + w / 2, y + h / 2, txt, fontsize=6.3, color=INK, ha="center", va="center")

    def arrow(a, b, col=INK2):
        ax.add_patch(FancyArrowPatch(a, b, arrowstyle="-|>", mutation_scale=7, color=col, lw=0.8))

    box(1, 16, 14, 8, f"1075 nm laser\n+ fibre amplifier\n~{pin_v + 2 * pin_a:.1f} W at the cavities", INK2)
    box(20, 16, 14, 8, "AOMs: 3 frequencies,\n>= 100 MHz apart\n(axes incoherent)", INK2)
    box(39, 16, 16, 8, "meta-optic mode\nconverters: Gaussian\n-> flat-top sheet", GREEN)
    box(60, 27, 18, 9, f"x and y rings, G = {g:.0f}\n{pin_a:.2f} W in, {pa:.0f} W circ.\neach", BLUE)
    box(60, 4, 18, 9, f"vertical ring, G = {g:.0f}\n{pin_v:.1f} W in,\n{pv:.0f} W circulating", AQUA)
    box(83, 16, 16, 8, "atoms: 10^6 sites\n0.2 mK, 4 um pitch", VIOLET)
    arrow((15.6, 20), (19.4, 20))
    arrow((34.6, 20), (38.4, 20))
    arrow((55.6, 21), (59.4, 30))
    arrow((55.6, 19), (59.4, 9))
    arrow((78.6, 31), (82.4, 22))
    arrow((78.6, 8), (82.4, 18))
    ax.text(57, 38.5, "PDH lock of the laser to each ring", fontsize=6.3, color=INK2, ha="center")
    box(60, 15.5, 18, 8, "measure: cameras on leaked\nlight, atom depth maps;\nactuate: graded-mirror shape", MAGENTA)
    for a_, b_ in (((69, 26.4), (69, 24.3)), ((69, 13.6), (69, 15.1))):
        ax.add_patch(FancyArrowPatch(a_, b_, arrowstyle="<|-|>", mutation_scale=6, color=MAGENTA, lw=0.8))
    ax.text(1, 1.5, f"Powers for 10^6 sites at 0.2 mK (results/cavity_budget/budget.json), mode matching 0.95; in vacuum with G = 300.", fontsize=6.0, color=INK2)
    ax.set_title("(b) Optics and control chain", fontsize=8, color=INK, loc="left")


def vertical_ring(ax_full, ax_zoom, vl):
    ch = vl["chosen"]
    mv = {k: np.array(v) for k, v in ch["mirrors_V_mm"].items()}
    paths = {1: ["X", "P1", "P2", "G1", "P2'", "P1'", "X"], 2: ["X", "P3", "P4", "G2", "P4'", "P3'", "X"]}
    cols = {1: AQUA, 2: GREEN}

    def pt(k):
        return np.zeros(2) if k == "X" else mv[k]

    for ax in (ax_full, ax_zoom):
        style(ax)
        for h, seq in paths.items():
            arr = np.array([pt(k) for k in seq])
            ax.plot(arr[:, 0], arr[:, 1], "-", color=cols[h], lw=1.1, label=f"half trip {h}" if ax is ax_full else None)
            for i in range(1, 6):
                a_, p_, b_ = pt(seq[i - 1]), pt(seq[i]), pt(seq[i + 1])
                n = (b_ - p_) / np.linalg.norm(b_ - p_) - (p_ - a_) / np.linalg.norm(p_ - a_)
                n = n / np.linalg.norm(n)
                t = np.array([-n[1], n[0]]) * 2.0
                ax.plot([p_[0] - t[0], p_[0] + t[0]], [p_[1] - t[1], p_[1] + t[1]], "-", color=INK, lw=2.2, solid_capstyle="butt")
        ax.add_patch(Rectangle((-21.2, -15), 42.4, 30, fill=False, ec=INK2, lw=0.8, ls="--"))
        ax.axhline(0, color=BLUE, lw=0.6, alpha=0.5)
        ax.plot(0, 0, "s", color=VIOLET, ms=4)
        ax.set_aspect("equal")
        ax.set_xlabel("s, along the array diagonal (mm)", fontsize=7, color=INK)
        ax.set_ylabel("z (mm)", fontsize=7, color=INK)
    for k, p_ in mv.items():
        off = (5, 4) if p_[0] >= 0 else (-22, 4)
        ax_full.annotate(k, p_, xytext=off, textcoords="offset points", fontsize=6.5, color=INK)
    # Objective as in configs/cavity_access.yaml: 42 mm plates + mount (radius 26 mm) from 38 mm up [CONFIRM?];
    # NA 0.5 light cone from the field corners.
    ax_full.add_patch(Polygon(np.array([[-2.83, 0], [2.83, 0], [24.8, 38], [-24.8, 38]]), closed=True, fc=YELLOW, alpha=0.25, ec="none"))
    ax_full.add_patch(Rectangle((-26, 38), 52, 30, fc=YELLOW, alpha=0.45, ec="none"))
    ax_full.text(0, 52, "42 mm meta-objective + mount", fontsize=6.5, color=INK, ha="center")
    ax_full.text(-60, 2.5, "plane of the in-plane rings", fontsize=6.5, color=BLUE, ha="center")
    ax_full.text(23, -14, "cell, 30 mm [CONFIRM?]", fontsize=6.5, color=INK2)
    ax_full.set_xlim(-125, 125)
    ax_full.set_ylim(-max(110, -min(p[1] for p in mv.values()) + 12), 72)
    ax_full.legend(fontsize=6.5, frameon=False, loc="lower center")
    ax_full.set_title(f"(a) Vertical-lattice ring in its own vertical plane, to scale (f1 {ch['f1_mm']:.0f} mm, {ch['f2_mm'] / ch['f1_mm']:.2f}x, "
                      f"P1 at {ch['d0_mm']:.0f} mm, G at {ch['half1']['L_mm']:.0f} mm from P2)", fontsize=8, color=INK, loc="left")
    d0 = ch["d0_mm"]
    ax_zoom.set_xlim(d0 - 30, d0 + 12)
    ax_zoom.set_ylim(-16, 16)
    ax_zoom.set_title("(b) Exit and return mirrors near P1, P3 (to scale)", fontsize=8, color=INK, loc="left")
    for k in ("P1", "P3"):
        ax_zoom.annotate(k, mv[k], xytext=(6, 0), textcoords="offset points", fontsize=7, color=INK)


def main():
    lay = json.loads((OUT / "layout.json").read_text())["layout"]
    fox = json.loads((REPO_ROOT / "results/cavity_foxli/foxli.json").read_text())["design_telK_M0.2_d70"]
    bud = json.loads((REPO_ROOT / "results/cavity_budget/budget.json").read_text())
    cfg = load_config(REPO_ROOT / "configs/cavity_foxli.yaml")

    fig = plt.figure(figsize=(11.0, 7.6), facecolor=SURFACE)
    gs = fig.add_gridspec(2, 2, width_ratios=[1.05, 1.0], height_ratios=[1.0, 0.8])
    top_view(fig.add_subplot(gs[:, 0]), lay)
    crossing(fig.add_subplot(gs[0, 1]), fox)
    side_view(fig.add_subplot(gs[1, 1]))
    fig.tight_layout()
    fig.savefig(OUT / "sketch_layout.png", dpi=170, facecolor=SURFACE)

    vl_path = OUT.parent / "cavity_vertical_layout" / "vertical_layout.json"
    if vl_path.exists():
        vl = json.loads(vl_path.read_text())
        fig = plt.figure(figsize=(11.0, 5.0), facecolor=SURFACE)
        gs5 = fig.add_gridspec(1, 2, width_ratios=[2.6, 1.0])
        vertical_ring(fig.add_subplot(gs5[0, 0]), fig.add_subplot(gs5[0, 1]), vl)
        fig.tight_layout()
        fig.savefig(OUT / "sketch_vertical.png", dpi=170, facecolor=SURFACE)

    fig = plt.figure(figsize=(11.0, 6.4), facecolor=SURFACE)
    gs = fig.add_gridspec(2, 2, width_ratios=[2.2, 1.0], height_ratios=[1.0, 0.95])
    ax_img = fig.add_subplot(gs[0, 0])
    ax_prof = fig.add_subplot(gs[0, 1])
    caustic(ax_img, ax_prof, cfg)
    chain(fig.add_subplot(gs[1, :]), bud)
    fig.tight_layout()
    fig.savefig(OUT / "sketch_ring.png", dpi=170, facecolor=SURFACE)


if __name__ == "__main__":
    main()
