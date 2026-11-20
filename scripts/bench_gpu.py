"""Benchmark and validate the GPU kernels on the fullaperture_100k geometry.

Usage: python scripts/bench_gpu.py [n_sites_per_side]      (default 316, the 100k design)

Builds the wedge pupil and fundamental sites as scripts/fullaperture.py does, then
(i) compares LazyKGGPU with the CPU LazyKG on 2,000 pupil rows and 500 sites (relative error),
(ii) times one full K @ x and one conj(v) @ K on the GPU and reports kernel evaluations per
second, and (iii) extrapolates the WGS time of the 10^6-site design (configs/fullaperture_1M.yaml:
5.2e12 evaluations per product, 120 iterations of two products).

Pilot mode (CLAUDE.md, compute planning):
    python scripts/bench_gpu.py --pilot configs/fullaperture_1M.yaml [--gpu] [--scales 0.1,0.2,0.3]

Runs two scaled-down versions of the full-aperture run described by the config, each in its own
process so that its peak memory (ru_maxrss) is measured separately. A pilot keeps the NA, the
site pitch, the pixel, the fold dtype and the render settings of the config, and scales the
stack (hence the aperture), the number of sites per side and the fold window by the factor s.
Each pilot times one kernel product (CPU LazyKG, or the GPU kernel with --gpu), the render and
fold of the whole metasurface, and the propagation plus tiled metrics, with random design
coefficients (the cost does not depend on their values). Projection to the config (printed with
its basis and written to results/bench_pilot/<config stem>.json):
  kernel product time  proportional to the number of kernel evaluations (8 x wedge pupil samples
                       x fundamental sites, counted exactly for the config), at the mean pilot rate;
                       the number of products and of full evaluations follows the config (fresh
                       WGS or initial_state, WGS-target or pupil-phase feedback, pillar cases)
  render and fold      a + b N_ap (pixels inside the aperture, pi (R / dx)^2), with library pillars
                       when the config renders them (an earlier version used the padded n^2 and
                       ideal pixels and failed its validation, docs/bench_pilot_validation.md)
  propagation, metrics a + b n_fold^2 log2(n_fold), fitted to the two pilots
  peak memory          a + b n_fold^2 x itemsize, fitted to the two pilots
These are extrapolations over the ratio of target to pilot size, stated in the output. With three
pilot sizes (the default) the largest relative fit residual is reported; pilots must be large
enough that the scaled terms dominate the fixed overheads (scales of 0.03 to 0.06 gave fits
dominated by the runtime's own memory and by padding, i.e. useless projections).
"""

import sys
import time
import types

import numpy as np

from metatrap.hologram import lazy as L
from metatrap.hologram import symmetric as S
from metatrap.hologram import targets
from metatrap.hologram.exact import OPLTable
from metatrap.hologram.wgs import PupilGrid
from metatrap.propagation.layers import LayerStack



def _pilot_geometry(cfg, s):
    from metatrap.hologram import targets as TG

    lam, na, pitch, dx = cfg["wavelength_um"], cfg["na"], cfg["site_pitch_um"], cfg["pixel_um"]
    layers = [[float(t) * s, m] for t, m in cfg["stack"]]
    stk = LayerStack.from_list(layers)
    dz = float(cfg.get("array_defocus_um", 0.0)) * s
    lsites = [list(l) for l in layers]
    lsites[-1][0] += dz
    stk_sites = LayerStack.from_list(lsites)
    R = stk.aperture_for_na(na, lam)
    nside = max(4, 2 * int(round(cfg["sites"]["nx"] * s / 2)))  # even: the D4 design has no sites on the axes
    sites_ = TG.make(kind="square", nx=nside, ny=nside, pitch_um=pitch)
    fs_, fstab_, _ = S.wedge_sites(sites_)
    rmax_ = np.abs(sites_).max() * np.sqrt(2)
    g_ = PupilGrid.for_sites(lam, pitch, 4, na, 2 * rmax_, guard=cfg.get("pupil_guard", 1.25))
    nf = max(256, int(round(cfg["n_fold"] * s / 64)) * 64)
    n = int(np.ceil((2 * R / dx + 2) / nf) * nf)
    return dict(lam=lam, na=na, pitch=pitch, dx=dx, stack=stk, stack_sites=stk_sites, R=R, sites=sites_, fs=fs_, fstab=fstab_,
                g=g_, n=n, nf=nf, rmax=rmax_)


def _pilot_worker(args):
    """One pilot in a fresh process; returns timings and peak RSS."""
    import resource

    cfg_path, s, gpu = args
    from metatrap.hologram import tiled as TL
    from metatrap.hologram.render import metasurface_blocks
    from metatrap.propagation.angular_spectrum import propagate as prop
    from metatrap.propagation.folded import fold_blocks
    from metatrap.runio import load_config

    cfg = load_config(cfg_path)
    G = _pilot_geometry(cfg, s)
    lam, na = G["lam"], G["na"]
    ill = dict(cfg["illumination"])
    tab_ = OPLTable(G["stack_sites"], lam, G["R"] + 2 * G["rmax"] + 20)
    wp_ = S.WedgePupil.build(G["g"], G["stack"], na, ill)
    rho = float(np.max(np.hypot(wp_.r[:, 0], wp_.r[:, 1]))) + float(np.max(np.abs(G["fs"]))) * 1.5 + 50.0
    if gpu:
        from metatrap.hologram.lazy_gpu import LazyKGGPU as KC
    else:
        KC = L.LazyKG
    K_ = KC(wp_, G["fs"], tab_, lam, rho)
    rng_ = np.random.default_rng(0)
    c = np.exp(2j * np.pi * rng_.random(len(G["fs"]))).astype(np.complex64)
    K_ @ c[:1].repeat(len(c))  # warm-up (compilation)
    reps, t0 = 0, time.time()
    while reps == 0 or time.time() - t0 < 1.0:  # at least 1 s of products, so the rate is not timer noise
        K_ @ c
        reps += 1
    t_k = (time.time() - t0) / reps
    pairs_ = 8.0 * len(G["fs"]) * len(wp_.opl0)
    lookup = None
    if cfg.get("library") and (cfg.get("feedback_pillars") == "library" or cfg.get("pillar_cases")):
        import importlib.util as _iu
        from pathlib import Path as _P

        _spec = _iu.spec_from_file_location("fa", _P(__file__).with_name("fullaperture.py"))
        _fa = _iu.module_from_spec(_spec)
        _spec.loader.exec_module(_fa)
        lookup = _fa.build_lookup(cfg, lam)  # the run renders library pillars: include the lookup cost
    t0 = time.time()
    holo = TL.tiled_render_grid(G["g"], G["stack"], na, ill, G["fs"], G["fstab"], c, tab_, lam, upsample=cfg.get("render_upsample", 5),
                                kernel=lambda *a: KC(*a))
    blocks = metasurface_blocks(G["n"], G["dx"], G["stack"], lam, na, None, None, ill, holo=holo, lookup=lookup,
                                block=cfg.get("block_rows", 64), with_incident=True)
    F, _, p_inc = fold_blocks(blocks, G["n"], G["nf"], dtype=np.dtype(cfg.get("fold_dtype", "complex128")))
    t_render = time.time() - t0
    del holo
    t0 = time.time()
    foc = TL.RolledView(prop(F.astype(np.complex64, copy=False), G["dx"], lam, G["stack_sites"]))
    TL.focal_metrics_tiled(cfg, None, G["dx"], lam, na, G["sites"], p_inc * G["dx"] ** 2, foc=foc, tile_px=cfg.get("metrics_tile_px", 1024))
    t_prop = time.time() - t0
    peak_mb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
    return {"scale": s, "n": G["n"], "n_aperture_px": float(np.pi * (G["R"] / G["dx"]) ** 2), "library_pillars": lookup is not None,
            "n_fold": G["nf"], "n_sites": len(G["sites"]), "pairs": pairs_, "t_kernel_s": t_k,
            "rate_per_s": pairs_ / t_k, "t_render_fold_s": t_render, "t_prop_metrics_s": t_prop, "peak_rss_mb": peak_mb,
            "itemsize": np.dtype(cfg.get("fold_dtype", "complex128")).itemsize}


def pilot(cfg_path, scales, gpu):
    import json as _json
    import multiprocessing as mp
    from pathlib import Path

    from metatrap.hologram import targets as TG
    from metatrap.runio import REPO_ROOT, load_config, write_json

    cfg = load_config(cfg_path)
    ctx = mp.get_context("spawn")
    rows = []
    for s_ in scales:
        with ctx.Pool(1) as pool:
            rows.append(pool.map(_pilot_worker, [(cfg_path, s_, gpu)])[0])
        print("pilot", {k: (round(v, 3) if isinstance(v, float) else v) for k, v in rows[-1].items()}, flush=True)
    # Target counts.
    lam, na, pitch, dx = cfg["wavelength_um"], cfg["na"], cfg["site_pitch_um"], cfg["pixel_um"]
    stk = LayerStack.from_list(cfg["stack"])
    R = stk.aperture_for_na(na, lam)
    sites_ = TG.make(**cfg["sites"])
    fs_, _, _ = S.wedge_sites(sites_)
    g_ = PupilGrid.for_sites(lam, pitch, 4, na, 2 * np.abs(sites_).max() * np.sqrt(2), guard=cfg.get("pupil_guard", 1.25))
    wp_ = S.WedgePupil.build(g_, stk, na, dict(cfg["illumination"]))
    pairs_t = 8.0 * len(fs_) * len(wp_.opl0)
    nf_t = int(cfg["n_fold"])
    n_t = int(np.ceil((2 * R / dx + 2) / nf_t) * nf_t)
    # Operation count of the run as configured (scripts/fullaperture.py): a fresh design runs the WGS
    # (two products per iteration); each feedback round is one full evaluation plus, for WGS-target
    # feedback, feedback_wgs_iter WGS iterations, or one product for pupil-phase feedback; the
    # pillar cases are one full evaluation each. A run resuming from initial_state starts at the
    # stored round (taken from the state file when it exists).
    start = 0
    if cfg.get("initial_state"):
        try:
            start = int(np.load(REPO_ROOT / cfg["initial_state"])["next_round"])
        except Exception:
            start = 0
    rounds = max(0, cfg.get("feedback_rounds", 0) - start)
    pp = cfg.get("feedback_mode", "wgs") == "pupil_phase"
    n_products = (0 if cfg.get("initial_state") else 2 * cfg.get("wgs_iter", 60)) + rounds * (1 if pp else 2 * cfg.get("feedback_wgs_iter", 30))
    n_evaluations = rounds + 1 + len(cfg.get("pillar_cases", []))
    rate = float(np.mean([r["rate_per_s"] for r in rows]))

    def fit(xs, ys):
        """Least squares a + b x; with three or more pilots also the largest relative residual."""
        xs, ys = np.asarray(xs, float), np.asarray(ys, float)
        A = np.column_stack([np.ones(len(xs)), xs])
        (a, b), *_ = np.linalg.lstsq(A, ys, rcond=None)
        resid = float(np.max(np.abs(A @ [a, b] - ys) / ys)) if len(xs) > 2 else float("nan")
        return float(a), float(b), resid

    xr = [r["n_aperture_px"] for r in rows]  # render cost follows the pixels inside the aperture, not the padded grid
    nap_t = float(np.pi * (R / dx) ** 2)
    xp = [r["n_fold"] ** 2 * np.log2(r["n_fold"]) for r in rows]
    xm = [r["n_fold"] ** 2 * r["itemsize"] / 2**20 for r in rows]
    ar, br, rr = fit(xr, [r["t_render_fold_s"] for r in rows])
    ap, bp, rp = fit(xp, [r["t_prop_metrics_s"] for r in rows])
    am, bm, rm = fit(xm, [r["peak_rss_mb"] for r in rows])
    itemsize = np.dtype(cfg.get("fold_dtype", "complex128")).itemsize
    t_eval = (ar + br * nap_t) + (ap + bp * nf_t**2 * np.log2(nf_t))
    proj = {
        "kernel_product_s": pairs_t / rate, "n_products": n_products, "products_total_h": n_products * pairs_t / rate / 3600,
        "render_fold_s_per_evaluation": ar + br * nap_t,
        "prop_metrics_s_per_evaluation": ap + bp * nf_t**2 * np.log2(nf_t),
        "n_evaluations": n_evaluations, "evaluations_total_h": n_evaluations * t_eval / 3600,
        "total_h": (n_products * pairs_t / rate + n_evaluations * t_eval) / 3600,
        "peak_memory_gb": (am + bm * nf_t**2 * itemsize / 2**20) / 1024,
        "fit_max_relative_residual": {"render": rr, "propagation": rp, "memory": rm},
        "basis": {"kernel": f"mean pilot rate {rate:.3g} evaluations/s ({'GPU' if gpu else 'CPU LazyKG'}), "
                            f"{pairs_t:.3g} evaluations per product (exact count), extrapolated x{pairs_t / max(r['pairs'] for r in rows):.0f}",
                  "render": f"a + b N_aperture fitted to {[round(x / 1e6, 1) for x in xr]} Mpx, target {nap_t / 1e6:.0f} Mpx "
                            f"(x{nap_t / max(xr):.0f}); library pillars: {rows[0]['library_pillars']}",
                  "propagation": f"a + b n_fold^2 log2 n_fold fitted to n_fold = {[r['n_fold'] for r in rows]}, target {nf_t}",
                  "memory": f"a + b n_fold^2 itemsize fitted to peak RSS {[round(r['peak_rss_mb']) for r in rows]} MB, target n_fold {nf_t}"},
    }
    print(_json.dumps(proj, indent=1), flush=True)
    out = REPO_ROOT / "results" / "bench_pilot"
    out.mkdir(parents=True, exist_ok=True)
    write_json(out / f"{Path(cfg_path).stem}.json", {"config": str(cfg_path), "gpu": gpu, "pilots": rows, "projection": proj})


def _default():
    from metatrap.hologram.lazy_gpu import LazyKGGPU

    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    nx = int(args[0]) if args else 316
    lam, na, pitch = 1.075, 0.55, 4.0
    stack = LayerStack.from_list([[2500, "air"], [4000, "fused_silica"], [9000, "air"]])
    sst = LayerStack.from_list([[2500, "air"], [4000, "fused_silica"], [9050, "air"]])
    R = stack.aperture_for_na(na, lam)
    sites = targets.make(kind="square", nx=nx, ny=nx, pitch_um=pitch)
    fs, fstab, _ = S.wedge_sites(sites)
    rmax = np.abs(sites).max() * np.sqrt(2)
    g = PupilGrid.for_sites(lam, pitch, 4, na, 2 * rmax, guard=1.25)
    wp = S.WedgePupil.build(g, stack, na, {"kind": "top_hat"})
    tab = OPLTable(sst, lam, R + 2 * rmax + 20)
    rho_max = float(np.max(np.hypot(wp.r[:, 0], wp.r[:, 1]))) + float(np.max(np.abs(fs))) * 1.5 + 50.0
    print(f"{len(sites)} sites, {len(fs)} fundamental, {len(wp.opl0)} wedge pupil samples", flush=True)

    # (i) accuracy on a subset
    rng = np.random.default_rng(0)
    rows = rng.choice(len(wp.opl0), min(2000, len(wp.opl0)), replace=False)
    cols = rng.choice(len(fs), min(500, len(fs)), replace=False)
    sub = types.SimpleNamespace(r=wp.r[rows], opl0=wp.opl0[rows])
    x = np.exp(2j * np.pi * rng.random(len(cols))).astype(np.complex64)
    v = np.exp(2j * np.pi * rng.random(len(rows))).astype(np.complex64)
    Kc, Kg = L.LazyKG(sub, fs[cols], tab, lam, rho_max), LazyKGGPU(sub, fs[cols], tab, lam, rho_max)
    for name, a, b in (("K @ x", Kc @ x, Kg @ x), ("conj(v) @ K", np.conj(v) @ Kc, np.conj(v) @ Kg)):
        print(f"accuracy {name}: max |gpu - cpu| / max |cpu| = {np.abs(a - b).max() / np.abs(a).max():.2e}", flush=True)

    # (ii) full-size timing
    K = LazyKGGPU(wp, fs, tab, lam, rho_max)
    x = np.exp(2j * np.pi * rng.random(len(fs))).astype(np.complex64)
    v = np.exp(2j * np.pi * rng.random(len(wp.opl0))).astype(np.complex64)
    pairs = 8.0 * len(fs) * len(wp.opl0)
    t0 = time.time()
    K @ x
    t1 = time.time()
    np.conj(v) @ K
    t2 = time.time()
    rate = pairs / (0.5 * (t2 - t0))
    print(f"K @ x {t1 - t0:.1f} s, conj(v) @ K {t2 - t1:.1f} s, {pairs:.2e} evaluations each: {rate:.2e} evaluations/s", flush=True)
    # (iii) extrapolation to 10^6 sites
    p1m = 5.22e12
    print(f"extrapolated 10^6-site design: {p1m / rate:.0f} s per product, {2 * p1m / rate / 60:.1f} min per WGS iteration, "
          f"{120 * 2 * p1m / rate / 3600:.1f} h for 120 iterations (render, propagation and metrics not included)")



if __name__ == "__main__":
    if "--pilot" in sys.argv:
        cfg_arg = sys.argv[sys.argv.index("--pilot") + 1]
        sc = [0.1, 0.2, 0.3]
        if "--scales" in sys.argv:
            sc = [float(x) for x in sys.argv[sys.argv.index("--scales") + 1].split(",")]
        pilot(cfg_arg, sc, "--gpu" in sys.argv)
    else:
        _default()
