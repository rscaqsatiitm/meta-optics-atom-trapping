"""Fold-window convergence of the 100k full-aperture metrics (aliasing of scattered light).

Usage: python scripts/fold_window_study.py configs/fold_window_100k.yaml

The folded propagation (propagation/folded.py) makes the focal field periodic with period
W = n_fold dx, so light that lands more than W/2 from the axis wraps back onto the array.
Each case (a stored design state, ideal or library pillars, optional width noise) is rendered
once and folded onto several windows at the same time (fold_blocks_multi); every window is
then propagated and measured with the same tiled metrics as scripts/fullaperture_lazy.py
--tiled. If the per-site residual is an alias of light scattered far from the array, its
amplitude falls as 1/W and the rms obeys rms(W)^2 = r0^2 + (a / W)^2; r0 is the window-
independent (physical) value. A design whose feedback was run at one window is out of sample
at every other window.

With n_folds equal to [n] (configs/fold_window_100k_full.yaml) the metasurface field is not folded:
the full 36,864^2 grid (10.9 GB in complex64, 18.4 mm period) is propagated in place, the
reference with no aliasing of light that stays within the aperture footprint.

Also recorded per window: the mean background intensity in a band outside the array
(relative to the mean site peak), which for aliased stray light falls as 1/W^2.
Writes results/<run_id>/metrics.json and per-site peaks per case and window.
"""

import gc
import importlib.util
import sys
import time
from pathlib import Path

import numpy as np

from metatrap.hologram import symmetric as S
from metatrap.hologram import targets
from metatrap.hologram.exact import OPLTable
from metatrap.hologram.tiled import RolledView, focal_metrics_tiled, tiled_render_grid
from metatrap.hologram.wgs import PupilGrid
from metatrap.propagation.angular_spectrum import propagate
from metatrap.propagation.folded import fold_blocks_multi
from metatrap.propagation.layers import LayerStack
from metatrap.runio import REPO_ROOT, load_config, run_dir, write_json, write_params

_spec = importlib.util.spec_from_file_location("fal", Path(__file__).with_name("fullaperture_lazy.py"))
LZ = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(LZ)
FA = LZ.FA


def geometry(cfg):
    lam, na, pitch, dx = cfg["wavelength_um"], cfg["na"], cfg["site_pitch_um"], cfg["pixel_um"]
    stack = LayerStack.from_list(cfg["stack"])
    radius = stack.aperture_for_na(na, lam)
    layers = [list(l) for l in cfg["stack"]]
    layers[-1][0] = float(layers[-1][0]) + float(cfg.get("array_defocus_um", 0.0))
    stack_sites = LayerStack.from_list(layers)
    ill = dict(cfg["illumination"])
    sites = targets.make(**cfg["sites"])
    fs, fstab, inv = S.wedge_sites(sites)
    rmax = np.abs(sites).max() * np.sqrt(2)
    g = PupilGrid.for_sites(lam, pitch, 4, na, 2 * rmax, guard=cfg.get("pupil_guard", 1.25))
    tab = OPLTable(stack_sites, lam, radius + 2 * rmax + 20)
    n = int(np.ceil((2 * radius / dx + 2) / 4096) * 4096)
    return dict(lam=lam, na=na, dx=dx, stack=stack, stack_sites=stack_sites, ill=ill, sites=sites, fs=fs, fstab=fstab, g=g, tab=tab, n=n)


def band_background(foc, dx, x_in, x_out, rows=512):
    """Mean |E|^2 over the square band x_in < max(|x|, |y|) < x_out (um), axis at N//2.

    Evaluated in row blocks, so the unfolded 36,864^2 grid (a RolledView) is never copied whole."""
    c = foc.shape[0] // 2
    k = int(x_out / dx)
    ax = (np.arange(-k, k)) * dx
    total, count = 0.0, 0
    for r0 in range(0, 2 * k, rows):
        r1 = min(2 * k, r0 + rows)
        sub = foc[c - k + r0 : c - k + r1, c - k : c + k]
        m = np.maximum(np.abs(ax[r0:r1])[:, None], np.abs(ax)[None, :])
        sel = (m > x_in) & (m < x_out)
        total += float(np.sum(np.abs(sub[sel]) ** 2))
        count += int(sel.sum())
    return total / count


def main(cfg_path):
    cfg = load_config(cfg_path)
    out = run_dir(cfg["run_id"])
    write_params(out, cfg)
    G = geometry(cfg)
    KC = LZ._kernel_class(False)
    n_folds = cfg["n_folds"]
    tmp = Path(cfg.get("scratch_dir", out / "tmp"))
    tmp.mkdir(parents=True, exist_ok=True)
    mpath = out / "metrics.json"
    import json

    res = json.load(open(mpath)) if mpath.exists() else {}
    res.setdefault("cases", {})
    half_array = float(np.abs(G["sites"]).max())
    # Band outside the array, capped at the 2.048 mm window's band so all runs compare like with like.
    bg_band = (half_array + 40.0, min(min(n_folds) * G["dx"] / 2 - 20.0, 1004.0))
    res["background_band_um"] = bg_band
    res["n_metasurface"] = G["n"]
    for case in cfg["cases"]:
        name = case["name"]
        if name in res["cases"] and len(res["cases"][name].get("windows", {})) == len(n_folds):
            print("skip (done)", name, flush=True)
            continue
        t0 = time.time()
        st = np.load(REPO_ROOT / case["state"])
        coeffs = st["coeffs"]
        lookup = FA.build_lookup(cfg, G["lam"]) if case.get("pillars", "ideal") == "library" else None
        holo = tiled_render_grid(G["g"], G["stack"], G["na"], G["ill"], G["fs"], G["fstab"], coeffs, G["tab"], G["lam"],
                                 upsample=cfg.get("render_upsample", 5), kernel=lambda *a: KC(*a))
        blocks = FA.metasurface_blocks(G["n"], G["dx"], G["stack"], G["lam"], G["na"], None, None, G["ill"], holo=holo, lookup=lookup,
                                       pillar_errors=case.get("errors"), block=cfg.get("block_rows", 64), with_incident=True)
        Fs, p_pix, p_inc_pix = fold_blocks_multi(blocks, G["n"], n_folds)
        del holo
        in_memory = n_folds == [G["n"]]  # unfolded grid: keep in memory, propagate in place
        if not in_memory:
            for nf, F in zip(n_folds, Fs):
                np.save(tmp / f"F_{nf}.npy", F)
            del Fs, F
        gc.collect()
        print(f"{name}: rendered and folded in {time.time() - t0:.0f}s", flush=True)
        p_in = p_inc_pix * G["dx"] ** 2
        cres = {"state": case["state"], "pillars": case.get("pillars", "ideal"), "errors": case.get("errors"),
                "metasurface_transmission": float(p_pix / p_inc_pix), "windows": {}}
        peaks = {}
        for nf in n_folds:
            t1 = time.time()
            if in_memory:
                foc = RolledView(propagate(Fs.pop(), G["dx"], G["lam"], G["stack_sites"]))  # in place (complex64)
            else:
                F = np.load(tmp / f"F_{nf}.npy")
                foc = np.fft.fftshift(propagate(F, G["dx"], G["lam"], G["stack_sites"]))
                del F
            gc.collect()
            m, per, _ = focal_metrics_tiled(cfg, None, G["dx"], G["lam"], G["na"], G["sites"], p_in, foc=foc,
                                            tile_px=cfg.get("metrics_tile_px", 1024))
            bg = band_background(foc, G["dx"], *bg_band)
            del foc
            gc.collect()
            pk = per["peak"]
            m["background_band_over_mean_peak"] = bg / float(pk.mean())
            m["window_mm"] = nf * G["dx"] / 1000
            m["seconds"] = time.time() - t1
            cres["windows"][str(nf)] = m
            peaks[f"peak_{nf}"] = pk
            print(f"  W {m['window_mm']:.3f} mm: rms {m['depth_nonuniformity_rms']:.5f}, eta {m['array_efficiency']:.4f}, "
                  f"core {m['core_efficiency']:.4f}, bg {m['background_band_over_mean_peak']:.2e} ({m['seconds']:.0f}s)", flush=True)
        if len(n_folds) < 3:
            cres["seconds"] = time.time() - t0
            res["cases"][name] = cres
            np.savez_compressed(out / f"per_site_{name}.npz", sites=G["sites"], **peaks)
            write_json(mpath, res)
            continue
        # Fit rms^2 = r0^2 + (a / W)^2 over the windows (least squares in rms^2 versus 1/W^2).
        W = np.array([cres["windows"][str(nf)]["window_mm"] for nf in n_folds])
        r = np.array([cres["windows"][str(nf)]["depth_nonuniformity_rms"] for nf in n_folds])
        A = np.column_stack([np.ones_like(W), 1 / W**2])
        (c0, c1), *_ = np.linalg.lstsq(A, r**2, rcond=None)
        cres["fit_rms2_vs_inv_W2"] = {"r0": float(np.sqrt(max(c0, 0.0))), "r0_squared": float(c0), "a_mm": float(np.sqrt(max(c1, 0.0))),
                                      "max_abs_residual_rms": float(np.max(np.abs(np.sqrt(np.clip(A @ [c0, c1], 0, None)) - r)))}
        # Correlation of the per-site errors between the smallest and the largest window.
        e = {k: v / v.mean() - 1 for k, v in peaks.items()}
        ks = [f"peak_{nf}" for nf in n_folds]
        cres["error_correlation_smallest_largest"] = float(np.corrcoef(e[ks[0]], e[ks[-1]])[0, 1])
        cres["rms_of_difference_smallest_largest"] = float(np.std(e[ks[0]] - e[ks[-1]]))
        cres["seconds"] = time.time() - t0
        res["cases"][name] = cres
        np.savez_compressed(out / f"per_site_{name}.npz", sites=G["sites"], **peaks)
        write_json(mpath, res)
        print(name, "fit", cres["fit_rms2_vs_inv_W2"], "corr", cres["error_correlation_smallest_largest"], flush=True)
    for f in tmp.glob("F_*.npy"):
        f.unlink()


if __name__ == "__main__":
    main(sys.argv[1])
