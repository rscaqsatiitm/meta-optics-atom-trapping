"""Rank meta-atom geometries and quantify sensitivities. Usage:
    python scripts/analyse_library.py configs/library_analysis.yaml

Outputs results/<run_id>/: ranking_<material>_<lam>.json, curves_<...>.png,
sensitivity.json, supercell.json, and selection.json (the chosen family per material).
"""

import sys

import numpy as np

from metatrap import materials
from metatrap.rcwa import lookup as L
from metatrap.rcwa.cell import PillarCell
from metatrap.rcwa.library import Library
from metatrap.rcwa.solver import solve_pillars
from metatrap.rcwa.supercell import deflection_efficiency
from metatrap.runio import REPO_ROOT, load_config, run_dir, write_json, write_params


def plot_curves(lib, material, lam, picks, path, n_in):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    for pk in picks:
        w, t, T0 = lib.curve(material, pk["pitch_um"], pk["height_um"], lam)
        lbl = f"p={pk['pitch_um']:.2f} h={pk['height_um']:.1f} (eta16={pk['eta']:.3f}, robust {pk['eta_robust']:.3f})"
        ax[0].plot(w * 1e3, np.unwrap(np.angle(t)) / (2 * np.pi), label=lbl)
        ax[1].plot(w * 1e3, T0, label=lbl)
    ax[0].set_xlabel("pillar width (nm)")
    ax[0].set_ylabel("transmission phase (waves)")
    ax[1].set_xlabel("pillar width (nm)")
    ax[1].set_ylabel("zeroth-order transmittance T0")
    ax[0].set_title(f"{material}, lambda = {lam*1e3:.1f} nm")
    ax[1].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)


def main(cfg_path):
    cfg = load_config(cfg_path)
    out = run_dir(cfg["run_id"])
    write_params(out, cfg)
    lib = Library.load(REPO_ROOT / cfg["library"])
    lim = L.FabLimits(**cfg["fab_limits"])
    lam = cfg["design_wavelength_um"]
    sub = materials.get(lib.attrs["substrate"])
    n_in = float(sub.n(lam))
    selection = {}
    for mat in lib.materials:
        rows = L.scan_geometries(lib, mat, lam, n_in, lim, n_levels=16)
        write_json(out / f"ranking_{mat}_{int(lam*1e3)}.json", rows)
        if not rows:
            print(mat, "no feasible geometry")
            continue
        top = rows[: cfg.get("n_plot", 5)]
        plot_curves(lib, mat, lam, top, out / f"curves_{mat}_{int(lam*1e3)}.png", n_in)
        # Tie-break: among geometries within ``tie_tol`` of the best robust efficiency,
        # take the shortest pillar (lowest etch depth and aspect ratio), then the best eta.
        tol = cfg.get("tie_tol", 0.003)
        near = [r for r in rows if r["eta_robust"] >= rows[0]["eta_robust"] - tol]
        hmin = min(r["height_um"] for r in near)
        best = max((r for r in near if r["height_um"] == hmin), key=lambda r: r["eta_robust"])
        selection[mat] = best
        print(mat, "best", best, flush=True)
    write_json(out / "selection.json", selection)

    # Absorption of the chosen family (A = 1 - T - R), averaged uniformly over phase
    # through the lookup table (each phase level is equally used by a hologram).
    absorb = {}
    for mat, best in selection.items():
        p, h = best["pitch_um"], best["height_um"]
        im, ip, ih, il = lib.index(mat, p, h, lam)
        w = lib.width_frac * p
        A = 1.0 - lib.T[im, ip, ih, il] - lib.R[im, ip, ih, il]
        lk = L.build(w, lib.t[im, ip, ih, il], p, h, n_in, 1.0, lim)
        a_tab = np.interp(lk.table_width, w, A)
        absorb[mat] = {"absorbed_fraction_phase_avg": float(np.mean(a_tab)), "absorbed_fraction_max": float(a_tab.max()),
                       "absorbed_fraction_min": float(a_tab.min())}
        print(mat, "absorption", absorb[mat], flush=True)
    write_json(out / "absorption.json", absorb)
    if not cfg.get("run_sensitivity", True):
        return

    # Sensitivity of the chosen family: +/- width bias and +/- index change (direct RCWA).
    sens = {}
    for mat, best in selection.items():
        p, h = best["pitch_um"], best["height_um"]
        w, t, _ = lib.curve(mat, p, h, lam)
        lk = L.build(w, t, p, h, n_in, 1.0, lim)
        targets = np.linspace(0, 2 * np.pi, 8, endpoint=False)
        wsel = lk.width_for_phase(targets)
        cell = PillarCell(p, h)
        n_p = complex(materials.get(mat).nk(lam))
        base = solve_pillars(cell, wsel, lam, n_p**2, n_in**2)
        entry = {"widths_um": wsel}
        for label, dw, dn, dh in [("width+10nm", 0.01, 0, 0), ("width-10nm", -0.01, 0, 0), ("n+0.05", 0, 0.05, 0), ("n-0.05", 0, -0.05, 0), ("height+50nm", 0, 0, 0.05)]:
            c2 = PillarCell(p, h + dh)
            r = solve_pillars(c2, wsel + dw, lam, (n_p + dn) ** 2, n_in**2)
            dphi = np.angle(r.t_x / base.t_x)
            dphi_rel = dphi - np.angle(np.mean(np.exp(1j * dphi)))  # remove common piston
            eta_ratio = np.abs(np.mean(r.t_x * np.exp(-1j * targets))) ** 2 / np.abs(np.mean(base.t_x * np.exp(-1j * targets))) ** 2
            entry[label] = {"dphi_rms_rad": float(np.sqrt(np.mean(dphi_rel**2))), "dphi_mean_rad": float(np.angle(np.mean(np.exp(1j * dphi)))), "eta8_ratio": float(eta_ratio)}
            print(mat, label, entry[label], flush=True)
        sens[mat] = entry
    write_json(out / "sensitivity.json", sens)

    # Edge efficiency: blazed supercells at the deflection angles of the outer aperture.
    sc = {}
    for mat, best in selection.items():
        p, h = best["pitch_um"], best["height_um"]
        w, t, _ = lib.curve(mat, p, h, lam)
        lk = L.build(w, t, p, h, n_in, 1.0, lim)
        n_p = complex(materials.get(mat).nk(lam))
        rows = []
        for nc in cfg["supercell_n"]:
            if lam / (nc * p) >= 1:
                continue
            for pol in ("x", "y"):
                r = deflection_efficiency(lk, nc, p, h, lam, n_p**2, n_in**2, n_grid=cfg.get("supercell_grid", 48), terms_per_cell=cfg.get("supercell_terms_per_cell", 100), polarisation=pol)
                r.pop("widths_um")
                rows.append(r)
                print(mat, r, flush=True)
        sc[mat] = rows
    write_json(out / "supercell.json", sc)


if __name__ == "__main__":
    main(sys.argv[1])
