"""Impact of the disc-centring fix (evaluate.disc_distances, 2026-11-20) on reported metrics.

Usage: python scripts/disc_fix_impact.py [small|full]

The cell and core discs of the spot metrics were centred with the sub-pixel peak offsets paired
with the wrong array axes (displacement up to 0.7 pixel). This script evaluates the same focal
intensity with the legacy and the fixed disc and reports the change of every disc-based metric
(array and core efficiency, cell-power non-uniformity, per-site core and cell power). Peaks,
positions and depth non-uniformity do not use the disc and cannot change.

Cases. small: the planted-truth design (tests/test_planted_hologram.py, untiled spot_metrics) and
the 16 x 16 render-check design (data/designs/tiled_render_check_dz50_g2_state.npz, tiled
metrics, unfolded 4096^2 grid). full: in addition the final 10^5-site design
(data/designs/fullaperture_100k_full_round26_state.npz, library pillars, unfolded 36,864^2 grid,
tiled metrics; ~13 min, ~12 GB).

Decision rule (written before the first run): if a reported quantity changes by more than the
precision it is quoted with (efficiencies: 0.0005 absolute; non-uniformities quoted in % to two
significant figures: 0.05 of the value), the affected runs are re-evaluated and the reports get
an erratum; otherwise the change is recorded here and in docs/planted_truth_hologram.md only.
Output: results/disc_fix_impact/metrics.json.
"""

import importlib.util
import sys
from pathlib import Path

import numpy as np

from metatrap.hologram import evaluate as E
from metatrap.hologram import tiled as TL
from metatrap.runio import REPO_ROOT, write_json

FIXED = E.disc_distances


def legacy(rr, ox, oy, dx):
    yy, xx = np.mgrid[-rr : rr + 1, -rr : rr + 1]
    return np.hypot((xx - ox) * dx, (yy - oy) * dx)


def use(fn):
    E.disc_distances = fn
    TL.disc_distances = fn


def diff(m_old, m_new, per_old, per_new):
    keys = ("array_efficiency", "core_efficiency", "cell_power_nonuniformity_rms")
    out = {k: {"legacy": float(m_old[k]), "fixed": float(m_new[k]), "change": float(m_new[k] - m_old[k])} for k in keys}
    for k in ("p_core", "p_cell"):
        out[f"{k}_per_site_max_rel_change"] = float(np.max(np.abs(per_new[k] / per_old[k] - 1)))
    return out


def planted_case():
    sys.path.insert(0, str(REPO_ROOT / "tests"))
    import test_planted_hologram as T

    # Rebuild the design exactly as the test fixture does, then measure twice.
    from metatrap.hologram import exact as X
    from metatrap.hologram.render import fourier_upsample, metasurface_field
    from metatrap.hologram.wgs import PupilGrid
    from metatrap.propagation.angular_spectrum import propagate
    from metatrap.propagation.layers import LayerStack

    sites, plant = T.sites_and_plant()
    stack = LayerStack.from_list([[T.Z, "air"]])
    g = PupilGrid.for_sites(T.LAM, T.PITCH, 4, T.NA, 2 * np.abs(sites).max() + T.PITCH, guard=1.6)
    ps = X.PupilSampling.build(g, stack, T.NA, T.ILL)
    radius = stack.aperture_for_na(T.NA, T.LAM)
    tab = X.OPLTable(stack, T.LAM, radius + 2 * np.abs(sites).max() + 20)
    K = X.kernels(ps, sites, tab, T.LAM)
    d = X.wgs_exact(ps, K, n_iter=80, n_fix=40, target=np.sqrt(plant))
    holo = X.render_grid(ps, stack, sites, d.coeffs, T.LAM, T.NA, kind="exact", table=tab, upsample=8)
    f = metasurface_field(T.N, T.DX, stack, T.LAM, T.NA, ps, d.psi, T.ILL, holo=holo)
    foc = propagate(f.copy(), T.DX, T.LAM, stack)
    w = 96
    c0 = T.N // 2 - w // 2
    I = np.abs(fourier_upsample(foc[c0 : c0 + w, c0 : c0 + w], 4)) ** 2
    use(legacy)
    mo = E.spot_metrics(I, T.DX / 4, sites, T.PITCH, T.LAM, T.NA, 1.0)
    use(FIXED)
    mn = E.spot_metrics(I, T.DX / 4, sites, T.PITCH, T.LAM, T.NA, 1.0)
    return diff(mo, mn, mo["per_site"], mn["per_site"])


def tiled_case(cfg_path, state, pillars):
    spec = importlib.util.spec_from_file_location("fws", Path(__file__).with_name("fold_window_study.py"))
    FWS = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(FWS)
    from metatrap.propagation.angular_spectrum import propagate
    from metatrap.propagation.folded import fold_blocks
    from metatrap.runio import load_config

    cfg = load_config(cfg_path)
    G = FWS.geometry(cfg)
    coeffs = np.load(REPO_ROOT / state)["coeffs"]
    lookup = FWS.FA.build_lookup(cfg, G["lam"]) if pillars == "library" else None
    KC = FWS.LZ._kernel_class(False)
    holo = TL.tiled_render_grid(G["g"], G["stack"], G["na"], G["ill"], G["fs"], G["fstab"], coeffs, G["tab"], G["lam"],
                                upsample=cfg.get("render_upsample", 5), kernel=lambda *a: KC(*a))
    blocks = FWS.FA.metasurface_blocks(G["n"], G["dx"], G["stack"], G["lam"], G["na"], None, None, G["ill"], holo=holo, lookup=lookup,
                                       block=cfg.get("block_rows", 64), with_incident=True)
    F, p_pix, p_inc = fold_blocks(blocks, G["n"], G["n"], dtype=np.complex64)
    del holo
    foc = TL.RolledView(propagate(F, G["dx"], G["lam"], G["stack_sites"]))
    p_in = p_inc * G["dx"] ** 2
    out = []
    for fn in (legacy, FIXED):
        use(fn)
        m, per, _ = TL.focal_metrics_tiled(cfg, None, G["dx"], G["lam"], G["na"], G["sites"], p_in, foc=foc,
                                           tile_px=cfg.get("metrics_tile_px", 1024))
        out.append((m, {k: v.copy() for k, v in per.items()}))
    return diff(out[0][0], out[1][0], out[0][1], out[1][1])


def main(mode):
    res = {"planted_3x3": planted_case()}
    print("planted", res["planted_3x3"], flush=True)
    res["render_check_16x16"] = tiled_case(REPO_ROOT / "configs" / "tiled_render_check.yaml",
                                           "data/designs/tiled_render_check_dz50_g2_state.npz", "ideal")
    print("16x16", res["render_check_16x16"], flush=True)
    if mode == "full":
        res["fullaperture_100k_round26"] = tiled_case(REPO_ROOT / "configs" / "fold_window_100k_full.yaml",
                                                      "data/designs/fullaperture_100k_full_round26_state.npz", "library")
        print("100k", res["fullaperture_100k_round26"], flush=True)
    out = REPO_ROOT / "results" / "disc_fix_impact"
    out.mkdir(parents=True, exist_ok=True)
    write_json(out / "metrics.json", res)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "small")
