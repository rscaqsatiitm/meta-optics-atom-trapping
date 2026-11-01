"""Full-pixel angular-spectrum validation of the hologram design models.

Usage: python scripts/validate_fullpixel.py configs/fullpixel_bench5mm.yaml

Designs the array hologram in the pupil domain with (a) Debye kernels and (b) exact
per-site kernels, renders each onto the metasurface pixel grid, propagates the full
aperture through the layer stack with the exact scalar angular-spectrum method, and
measures every spot. A plain lens provides the Strehl reference. Predictions of the
pupil-domain models (mode efficiency per site) are stored next to the measurements.
"""

import gc
import sys
import time

import numpy as np

from metatrap.hologram import evaluate as E
from metatrap.hologram import exact as X
from metatrap.hologram import targets
from metatrap.hologram.field_model import spot_aberrations
from metatrap.hologram.render import fourier_upsample, metasurface_field
from metatrap.hologram.wgs import PupilGrid
from metatrap.propagation.angular_spectrum import propagate
from metatrap.propagation.layers import LayerStack
from metatrap.runio import load_config, run_dir, write_json, write_params


def focal_window(field_ms, dx, lam, stack, half_um, up):
    p_in = float(np.sum(np.abs(field_ms) ** 2) * dx * dx)
    out = propagate(field_ms, dx, lam, stack)
    n = out.shape[0]
    w = int(np.ceil(2 * half_um / dx / 2) * 2)
    c0 = n // 2 - w // 2
    win = np.array(out[c0 : c0 + w, c0 : c0 + w])
    del out
    gc.collect()
    f = fourier_upsample(win, up)
    return np.abs(f) ** 2, dx / up, p_in


def main(cfg_path):
    cfg = load_config(cfg_path)
    out = run_dir(cfg["run_id"])
    write_params(out, cfg)
    lam, na, pitch = cfg["wavelength_um"], cfg["na"], cfg["site_pitch_um"]
    stack = LayerStack.from_list(cfg["stack"])
    ill = cfg["illumination"]
    sites = targets.make(**cfg["sites"])
    half = np.abs(sites).max() + 4 * pitch
    radius = stack.aperture_for_na(na, lam)
    dx = cfg["pixel_um"]
    n = int(np.ceil(2 * radius / dx / 256 + 0.05) * 256)
    print(f"aperture {2*radius/1000:.3f} mm, grid {n}^2 at {dx} um", flush=True)

    # Pupil-domain designs.
    field_extent = 2 * np.abs(sites).max() + pitch  # covers arrays displaced from the axis
    g = PupilGrid.for_sites(lam, pitch, cfg.get("pupil_oversample", 4), na, field_extent, guard=cfg.get("pupil_guard", 1.6))
    ps = X.PupilSampling.build(g, stack, na, ill)
    tab = X.OPLTable(stack, lam, radius + 2 * np.abs(sites).max() + 20)
    designs = {}
    Kd = X.debye_kernels(ps, sites, lam)
    designs["debye"] = X.wgs_exact(ps, Kd, n_iter=cfg["wgs_iter"], n_fix=cfg["wgs_fix"])
    del Kd
    gc.collect()
    K = X.kernels(ps, sites, tab, lam)
    designs["exact"] = X.wgs_exact(ps, K, n_iter=cfg["wgs_iter"], n_fix=cfg["wgs_fix"])
    pred = {k: np.abs(X.mode_amplitudes(ps, d.psi, K)) ** 2 for k, d in designs.items()}
    up = cfg.get("render_upsample", 8)
    holos = {k: (X.render_grid(ps, stack, sites, d.coeffs, lam, na, kind=k, table=tab, upsample=up) if up > 1 else None) for k, d in designs.items()}
    del K
    gc.collect()
    ab = spot_aberrations(stack, lam, na, sites, ill, n_pupil=160)

    # Plain-lens reference.
    t0 = time.time()
    f = metasurface_field(n, dx, stack, lam, na, None, None, ill)
    I, dxu, p_in = focal_window(f, dx, lam, stack, 4 * pitch, cfg["upsample"])
    ref = E.single_spot_reference(I, dxu, lam, na, pitch)
    ref["p_in"] = p_in
    ref["core_fraction"] = ref["p_core"] / p_in
    print("reference", ref, f"{time.time()-t0:.0f}s", flush=True)
    del f, I
    gc.collect()

    results = {"reference": ref, "grid": {"n": n, "dx_um": dx, "pupil_n": g.n, "pupil_ds": g.ds, "pupil_samples": int(len(ps.amp))}}
    n_fb = cfg.get("feedback_rounds", 0)
    if n_fb:
        # Closed-loop correction with the full-pixel model as the 'measurement', as an
        # experiment would use camera or atom-based feedback (refs. 4, 14 use camera feedback
        # [REF?: verify]). Target weights are rescaled by (mean peak / peak)^gain and the WGS
        # is warm-started from the previous phase.
        K = X.kernels(ps, sites, tab, lam)
        tgt = np.ones(len(sites))
        d = designs["exact"]
        hist = []
        for rnd in range(n_fb + 1):
            holo = X.render_grid(ps, stack, sites, d.coeffs, lam, na, kind="exact", table=tab, upsample=up)
            f = metasurface_field(n, dx, stack, lam, na, ps, d.psi, ill, holo=holo)
            I, dxu, p_in = focal_window(f, dx, lam, stack, half, cfg["upsample"])
            del f
            m = E.spot_metrics(I, dxu, sites, pitch, lam, na, p_in, ref)
            pk = m["per_site"]["peak"]
            hist.append({"round": rnd, "depth_nonuniformity_rms": m["depth_nonuniformity_rms"], "depth_min_over_mean": m["depth_min_over_mean"], "array_efficiency": m["array_efficiency"], "strehl_mean": m.get("strehl_mean")})
            print("feedback", hist[-1], flush=True)
            if rnd == n_fb:
                break
            tgt = tgt * (pk.mean() / pk) ** cfg.get("feedback_gain", 0.5)
            d = X.wgs_exact(ps, K, n_iter=cfg.get("feedback_wgs_iter", 30), n_fix=None, psi0=d.psi, target=tgt)
            del I
            gc.collect()
        designs["exact_feedback"] = d
        pred["exact_feedback"] = np.abs(X.mode_amplitudes(ps, d.psi, K)) ** 2
        holos["exact_feedback"] = X.render_grid(ps, stack, sites, d.coeffs, lam, na, kind="exact", table=tab, upsample=up)
        del K
        gc.collect()
    for name, d in designs.items():
        t0 = time.time()
        f = metasurface_field(n, dx, stack, lam, na, ps, d.psi, ill, holo=holos[name])
        I, dxu, p_in = focal_window(f, dx, lam, stack, half, cfg["upsample"])
        del f
        m = E.spot_metrics(I, dxu, sites, pitch, lam, na, p_in, ref)
        per = m.pop("per_site")
        # Intensity on the optical axis relative to the mean site peak (zero-order check;
        # meaningful when no site lies within ~1 um of the axis).
        m["axis_intensity_over_mean_site_peak"] = float(I[I.shape[0] // 2, I.shape[1] // 2] / per["peak"].mean())
        m["predicted_mode_efficiency_sum"] = float(pred[name].sum())
        m["predicted_nonuniformity_rms"] = float(pred[name].std() / pred[name].mean())
        m["design_history_nonuniformity_final"] = d.history_nonuniformity[-1]
        results[name if up > 1 else f"{name}_coarse_render"] = m
        tag = name if up > 1 else f"{name}_nointerp"
        np.savez_compressed(out / f"per_site_{tag}.npz", sites=sites, pred_mode_eff=pred[name], strehl_debye_model=ab["strehl"], dz_debye_model=ab["dz_um"], **per)
        print(name, {k: (round(v, 4) if isinstance(v, float) else v) for k, v in m.items()}, f"{time.time()-t0:.0f}s", flush=True)
        if cfg.get("save_png", True):
            plot(I, dxu, sites, per, pred[name], ab, out / f"focal_{tag}.png", tag)
        del I
        gc.collect()
    if n_fb:
        results["feedback_history"] = hist
    write_json(out / "metrics.json", results)


def plot(I, dx, sites, per, pred, ab, path, title):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 3, figsize=(15, 4.5))
    n = I.shape[0]
    ext = np.array([-n // 2, n - n // 2]) * dx
    ax[0].imshow(I.T, origin="lower", extent=[*ext, *ext], cmap="magma")
    ax[0].set_title(f"{title}: focal intensity")
    ax[0].set_xlabel("x (um)")
    ax[0].set_ylabel("y (um)")
    r = np.hypot(*sites.T)
    pk = per["peak"] / per["peak"].mean()
    ax[1].plot(r, pk, ".", ms=3, label="full-pixel peak / mean")
    ax[1].plot(r, pred / pred.mean(), ".", ms=3, label="pupil model (exact kernels) / mean")
    ax[1].set_xlabel("site radius (um)")
    ax[1].set_ylabel("relative trap depth")
    ax[1].legend(fontsize=8)
    if "strehl" in per:
        ax[2].plot(r, per["strehl"], ".", ms=3, label="full-pixel Strehl")
    ax[2].plot(r, ab["strehl"], ".", ms=3, label="Debye-design aberration model")
    ax[2].set_xlabel("site radius (um)")
    ax[2].set_ylabel("Strehl ratio")
    ax[2].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)


if __name__ == "__main__":
    main(sys.argv[1])
