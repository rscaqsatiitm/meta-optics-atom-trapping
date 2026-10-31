"""Design and evaluate the static B1 meta-hologram in the pupil domain.

Usage: python scripts/design_hologram.py configs/hologram_b1_cell20mm.yaml

For each illumination, designs the array with Debye kernels and with exact per-site
kernels, evaluates both with the exact kernels (mode efficiency per site), runs the
Debye-design aberration model (Strehl, axial and lateral shift per site), and converts
the result into trap depths and site counts with the Cs power budget, through an
explicit efficiency chain. Writes metrics.json, per-site npz and figures.
"""

import gc
import sys

import numpy as np

from metatrap.atoms.alkali import depth_K
from metatrap.hologram import exact as X
from metatrap.hologram import targets
from metatrap.hologram.field_model import spot_aberrations
from metatrap.hologram.wgs import PupilGrid
from metatrap.propagation.layers import LayerStack
from metatrap.runio import load_config, run_dir, write_json, write_params


def ideal_peak_per_watt(ps, lam):
    """Peak intensity (1/um^2) per unit power of the aberration-free focus of this pupil."""
    s = np.hypot(ps.sigma[:, 0], ps.sigma[:, 1])
    cz = np.sqrt(1 - s**2)
    return (ps.grid.ds / lam) ** 2 * np.sum(ps.amp) ** 2 / np.sum(ps.amp**2 * cz)


def main(cfg_path):
    cfg = load_config(cfg_path)
    out = run_dir(cfg["run_id"])
    write_params(out, cfg)
    lam, na, pitch = cfg["wavelength_um"], cfg["na"], cfg["site_pitch_um"]
    stack = LayerStack.from_list(cfg["stack"])
    sites = targets.make(**cfg["sites"])
    field_extent = np.ptp(sites[:, 0]) + pitch
    radius = stack.aperture_for_na(na, lam)
    chain = cfg["efficiency_chain"]
    chain_total = float(np.prod(list(chain.values())))
    res = {"aperture_diameter_mm": 2 * radius / 1000, "n_sites": len(sites), "efficiency_chain": chain, "chain_total": chain_total, "cases": {}}
    print(f"aperture {2*radius/1000:.2f} mm, {len(sites)} sites, chain {chain_total:.3f}", flush=True)
    for ill in cfg["illuminations"]:
        if ill["kind"] == "top_hat":
            ef = ill.get("edge_frac", 0.0)
            tag = "top_hat" if ef == 0 else f"top_hat_edge{ef:g}"
            ill_abs = {"kind": "top_hat", "edge_um": ef * radius}
        else:
            tag = f"gauss_w{ill['waist_over_radius']}"
            ill_abs = {"kind": "gaussian", "waist_um": ill["waist_over_radius"] * radius}
        clip = 0.0 if ill["kind"] == "top_hat" else float(np.exp(-2 / ill["waist_over_radius"] ** 2))
        g = PupilGrid.for_sites(lam, pitch, cfg.get("pupil_oversample", 4), na, field_extent, guard=cfg.get("pupil_guard", 1.6))
        ps = X.PupilSampling.build(g, stack, na, ill_abs)
        tab = X.OPLTable(stack, lam, radius + 2 * np.abs(sites).max() + 20)
        Kd = X.debye_kernels(ps, sites, lam)
        d_debye = X.wgs_exact(ps, Kd, n_iter=cfg["wgs_iter"], n_fix=cfg["wgs_fix"])
        del Kd
        gc.collect()
        K = X.kernels(ps, sites, tab, lam)
        d_exact = X.wgs_exact(ps, K, n_iter=cfg["wgs_iter"], n_fix=cfg["wgs_fix"])
        ab = spot_aberrations(stack, lam, na, sites, ill_abs, n_pupil=cfg.get("aberration_pupil", 160))
        i_per_w = ideal_peak_per_watt(ps, lam)
        case = {"clip_loss": clip, "pupil_samples": int(len(ps.amp)), "ideal_peak_per_W_per_um2": float(i_per_w),
                "debye_model": {"strehl_min": float(ab["strehl"].min()), "strehl_mean": float(ab["strehl"].mean()),
                                "dz_min_um": float(ab["dz_um"].min()), "dz_max_um": float(ab["dz_um"].max()),
                                "lateral_shift_max_um": float(np.hypot(ab["dx_um"], ab["dy_um"]).max())}}
        for name, d in (("debye", d_debye), ("exact", d_exact)):
            eta = np.abs(X.mode_amplitudes(ps, d.psi, K)) ** 2
            eta_site = eta * chain_total * (1 - clip)
            case[name] = {
                "mode_efficiency_sum": float(eta.sum()),
                "depth_nonuniformity_rms": float(eta.std() / eta.mean()),
                "depth_min_over_mean": float(eta.min() / eta.mean()),
                "delivered_fraction_chain": float(eta_site.sum()),
            }
            for p_w in cfg["laser_powers_W"]:
                i_peak = eta_site.mean() * p_w * i_per_w * 1e12  # W/m^2
                case[name][f"mean_depth_mK_at_{p_w}W"] = float(depth_K(i_peak, lam * 1e3) * 1e3)
            np.savez_compressed(out / f"per_site_{tag}_{name}.npz", sites=sites, eta=eta, psi=d.psi)
            print(tag, name, case[name], flush=True)
        # Equivalent parameters of the proposal's power model for direct comparison.
        eta_ex = case["exact"]["delivered_fraction_chain"]
        na_eff = float(np.sqrt(i_per_w * lam**2 / np.pi))
        case["proposal_model_equivalent"] = {"utilisation": eta_ex, "na_eff_ideal_spot": na_eff}
        res["cases"][tag] = case
        np.savez_compressed(out / f"aberration_{tag}.npz", sites=sites, **ab)
        plot(sites, ab, out / f"aberration_{tag}.png", tag)
        del K
        gc.collect()
    write_json(out / "metrics.json", res)


def plot(sites, ab, path, tag):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 2, figsize=(10, 4.2))
    sc = ax[0].scatter(sites[:, 0], sites[:, 1], c=ab["strehl"], s=6, cmap="viridis")
    fig.colorbar(sc, ax=ax[0], label="Strehl ratio (Debye design)")
    ax[0].set_aspect("equal")
    ax[0].set_xlabel("x (um)")
    ax[0].set_ylabel("y (um)")
    sc = ax[1].scatter(sites[:, 0], sites[:, 1], c=ab["dz_um"], s=6, cmap="coolwarm")
    fig.colorbar(sc, ax=ax[1], label="axial focus shift dz (um)")
    ax[1].set_aspect("equal")
    ax[1].set_xlabel("x (um)")
    ax[0].set_title(f"{tag}: Debye-design field aberration")
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)


if __name__ == "__main__":
    main(sys.argv[1])
