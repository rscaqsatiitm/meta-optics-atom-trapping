"""Wave simulation of a patch of lenslet-array traps through the meta-objective.

Usage: python scripts/million_patch.py configs/million.yaml

A 7 x 7 patch of traps (pitch a) is centred on the axis or at the 10^6-site corner radius. Each
lenslet's light fills the objective pupil (circular stop, square lenslet circumscribing it) as a
plane wave at the field angle of its site. With pre-compensation it carries -W_j, the
objective's wavefront error for that site (ray trace, one common atom plane). All traps come
from one laser and add coherently with per-lenslet pistons psi_j (equal, random, or alternating
0 / pi). The objective is axisymmetric, so the plate field of a site at (rho cos phi,
rho sin phi) is that of (rho, 0) rotated by phi. The intensity is computed by the superposition
method on a grid over the patch.

Reported per case: depth (peak intensity) non-uniformity over the inner 5 x 5 traps (the outer
ring lacks neighbours on one side), the largest peak position error, and the brightest
interstitial point relative to the mean trap peak. Output: results/million/patch.json.
"""

import json
import sys

import numpy as np
from scipy.interpolate import RegularGridInterpolator

from metatrap.metalens import doublet as D
from metatrap.metalens import wave_offaxis as WO
from metatrap.metalens.design import _quality
from metatrap.metalens.raytrace import point_source_rays
from metatrap.propagation.layers import LayerStack
from metatrap.runio import REPO_ROOT, load_config, write_json


def site_field(d, s, lam, R1, f, rho, dz, correct, n_pupil):
    th = D.theta_for_height(s, R1, lam, rho, f) if rho > 0 else 0.0
    extra = None
    if correct:
        n = 96
        c = (np.arange(n) + 0.5) / n * 2 * R1 - R1
        X, Y = np.meshgrid(c, c, indexing="ij")
        px, py, o = point_source_rays(X, Y, th, None)
        r = s.trace(X.ravel(), Y.ravel(), px.ravel(), py.ravel(), lam, stop={"index": 0, "radius_um": R1}, opl0=o.ravel())
        q = _quality(*r, lam, dz)
        W = np.zeros(X.size)
        W[r[5]] = q["residual"]
        interp = RegularGridInterpolator((c, c), W.reshape(X.shape), bounds_error=False, fill_value=0.0)
        extra = lambda x, y: -interp(np.column_stack([np.ravel(x), np.ravel(y)])).reshape(np.shape(x))  # noqa: E731
    pf = WO.plate_field(d, lam, th, n_pupil, extra_opl=extra)
    x_img = D.quality(s, R1, lam, th, None, dz)["x_um"] if rho > 0 else 0.0
    return pf, x_img


def main(cfg_path):
    cfg = load_config(cfg_path)
    pc = cfg["patch"]
    lam, na, a = pc["wavelength_um"], pc["na"], pc["pitch_um"]
    cell = [[pc["air_gap_um"], "air"], [4000, "fused_silica"], [9000, "air"]]
    f0 = LayerStack.from_list(cell).aperture_for_na(na, lam) / na
    d = D.make(D.DoubletSpec(na, cell, [[f0, "air"]], n_pts=2001, ext=min(1.6, 0.95 / na)), lam)
    R1, f = d.info["R1_um"], d.info["f_um"]
    s = D.system(d, lam)
    dz, _ = D.common_plane(s, R1, lam, list(np.linspace(0, 2828, 6)), f, n_rays=28)
    fs = WO.FocalSum(d.spec.cell, lam, d.info["R2_um"] * 1.6 + 3000)
    path = REPO_ROOT / "results" / "million" / "patch.json"
    res = json.load(open(path)) if path.exists() else {}
    k = np.arange(-3, 4)
    rng = np.random.default_rng(5)
    for centre in pc["centres_um"]:
        for correct in (False, True):
            key0 = f"r{centre}:{'corrected' if correct else 'uncorrected'}"
            if all(f"{key0}:{p}" in res for p in pc["pistons"]):
                continue
            # Per-site plate fields (rotated copies of the radial computation).
            fields, nominal = [], []
            for i in k:
                for j in k:
                    x, y = centre + i * a, j * a
                    rho, phi = np.hypot(x, y), np.arctan2(y, x)
                    pf, x_img = site_field(d, s, lam, R1, f, rho, dz, correct, pc["n_pupil"])
                    cph, sph = np.cos(phi), np.sin(phi)
                    pf = dict(pf)
                    pf["x"], pf["y"] = cph * pf["x"] - sph * pf["y"], sph * pf["x"] + cph * pf["y"]
                    fields.append(pf)
                    nominal.append((x_img * cph, x_img * sph))
            nominal = np.array(nominal)
            # Distortion is calibrated by the lenslet positions: shift each site to its nominal grid
            # point (x, y) by referencing the evaluation grid to the imaged positions.
            g = np.arange(-3.5 * a, 3.5 * a + 1e-9, pc["grid_um"])
            cx, cy = nominal[:, 0].mean(), nominal[:, 1].mean()
            GX, GY = np.meshgrid(cx + g, cy + g, indexing="ij")
            P = np.column_stack([GX.ravel(), GY.ravel()])
            U = [fs.field(pf, pf["amp"], P, dz) for pf in fields]
            for pist in pc["pistons"]:
                key = f"{key0}:{pist}"
                if key in res:
                    continue
                if pist == "equal":
                    ps = np.zeros(len(U))
                elif pist == "random":
                    ps = rng.uniform(0, 2 * np.pi, len(U))
                else:  # alternating checkerboard 0 / pi
                    ps = np.array([np.pi * ((i + j) % 2) for i in k for j in k])
                tot = sum(u * np.exp(1j * p) for u, p in zip(U, ps))
                I = (np.abs(tot) ** 2).reshape(GX.shape)
                Isingle = np.array([np.abs(u).reshape(GX.shape) ** 2 for u in U])
                pk, pos_err = [], []
                for n_, (i, j) in enumerate([(i, j) for i in k for j in k]):
                    xs, ys = nominal[n_]
                    m = (np.abs(GX - xs) <= 0.6) & (np.abs(GY - ys) <= 0.6)
                    idx = np.argmax(np.where(m, I, -1))
                    pk.append(I.ravel()[idx])
                    pos_err.append(np.hypot(GX.ravel()[idx] - xs, GY.ravel()[idx] - ys))
                pk = np.array(pk).reshape(7, 7)
                inner = pk[1:-1, 1:-1]
                single_pk = np.array([Isingle[n_].max() for n_ in range(len(U))]).reshape(7, 7)[1:-1, 1:-1]
                # Interstitial points: centres of four inner traps.
                inter = []
                for i in range(1, 5):
                    for j in range(1, 5):
                        xs = 0.25 * (nominal[i * 7 + j, 0] + nominal[(i + 1) * 7 + j, 0] + nominal[i * 7 + j + 1, 0] + nominal[(i + 1) * 7 + j + 1, 0])
                        ys = 0.25 * (nominal[i * 7 + j, 1] + nominal[(i + 1) * 7 + j, 1] + nominal[i * 7 + j + 1, 1] + nominal[(i + 1) * 7 + j + 1, 1])
                        m = (np.abs(GX - xs) <= 0.8) & (np.abs(GY - ys) <= 0.8)
                        inter.append(I[m].max())
                res[key] = {"centre_um": centre, "corrected": correct, "pistons": pist,
                            "depth_nonuniformity_rms_inner": float(inner.std() / inner.mean()),
                            "min_over_mean_inner": float(inner.min() / inner.mean()),
                            "single_trap_nonuniformity_rms_inner": float(single_pk.std() / single_pk.mean()),
                            "interference_modulation_rms": float(np.std(inner / single_pk) / np.mean(inner / single_pk)),
                            "max_position_error_um": float(np.max(pos_err)),
                            "interstitial_max_over_mean_peak": float(np.max(inter) / inner.mean()),
                            "peak_rel_axis_single": float(single_pk.mean())}
                print(key, {a_: (round(b, 4) if isinstance(b, float) else b) for a_, b in res[key].items()}, flush=True)
                write_json(path, res)


if __name__ == "__main__":
    main(sys.argv[1])
