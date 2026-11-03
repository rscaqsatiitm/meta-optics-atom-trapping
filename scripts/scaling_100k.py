"""Scaling to 100,000 sites: does a shorter working distance help?

Usage: python scripts/scaling_100k.py configs/scaling_100k.yaml

Trapping hologram (side A). For each geometry (metasurface a distance WD from the atom plane
in vacuum, aperture radius R set by the NA for the axial site) and each site radius r, the
ideal exact-kernel spot of site j (all aperture light phased to focus at x_j) is computed by a
Rayleigh-Sommerfeld sum over the aperture,

    U(P) = sum_i (z / L_i^2) exp(i k (L_i(P) - L_i(x_j))) dA,  L_i(P) = |r_i - P|,

uniformly illuminated, searched for its maximum around x_j in a few planes. Reported: (i) the
peak intensity relative to the axial site at equal power, (ii) the 1/e^2 half-widths along and
across the radius, (iii) the largest local deflection the meta-atoms must produce for that
site, s = |r_i - x_j| / L_i, and (iv) the power-weighted meta-atom efficiency from the
supercell curve eta(s). The trap depth of site j relative to the axial site is
D_j = I_rel x eta_j / eta_0. For a square array of N sites at pitch a, uniform depth requires
power in proportion to 1 / D_j, so the power needed relative to an ideal flat field is
F = mean(1 / D_j) (DERIVED). The absolute power is scaled from the simulated v2 design.

Imaging/addressing lens (side B): for the same geometries the air-spaced aplanatic doublet is
designed (spacing scanned near f) and the field radius with a free-focus Strehl ratio >= 0.8
is found.

Output: results/<run_id>/scaling.json.
"""

import json
import sys
import time

import numpy as np

from metatrap.metalens import doublet as D
from metatrap.runio import REPO_ROOT, load_config, run_dir, write_json, write_params


def eta_curve(cfg):
    pts = {}
    for path in cfg["supercell_files"]:
        p = REPO_ROOT / path
        if not p.exists():
            continue
        d = json.load(open(p))
        rows = d.values() if isinstance(d, dict) and "sin_lpcvd" not in d else [r | {"eta_plus1": r["eta_rcwa"]} for r in d[cfg["supercell_family"]]]
        for r in rows:
            if "lambda_um" in r and abs(r["lambda_um"] - cfg["wavelength_um"]) > 1e-3:
                continue
            pts.setdefault(round(r["sin_theta"], 3), []).append(r.get("eta_plus1", r.get("eta_rcwa")))
    s = np.array(sorted(pts))
    e = np.array([np.mean(pts[k]) for k in s])
    return np.concatenate([[0.0], s]), np.concatenate([[cfg["eta_lpa"]], e])


def spot(WD, R, r_site, lam, n_ap, eta_s=None):
    k = 2 * np.pi / lam
    c = (np.arange(n_ap) + 0.5) / n_ap * 2 * R - R
    X, Y = np.meshgrid(c, c, indexing="ij")
    m = np.hypot(X, Y) <= R
    x, y = X[m], Y[m]
    dA = (c[1] - c[0]) ** 2
    Lj = np.sqrt((x - r_site) ** 2 + y**2 + WD**2)
    s_loc = np.hypot(x - r_site, y) / Lj

    def field(P, dz):
        z = WD + dz
        out = np.empty(len(P), complex)
        for a in range(0, len(P), 64):
            Pb = P[a : a + 64]
            L = np.sqrt((x[None, :] - Pb[:, 0:1]) ** 2 + (y[None, :] - Pb[:, 1:2]) ** 2 + z**2)
            out[a : a + 64] = np.sum((z / L**2) * np.exp(1j * k * (L - Lj[None, :])), axis=1) * dA
        return out

    best = None
    for dz in np.linspace(-1.0, 1.0, 9):
        g = np.linspace(-0.4, 0.4, 9)
        GX, GY = np.meshgrid(r_site + g, g, indexing="ij")
        P = np.column_stack([GX.ravel(), GY.ravel()])
        I = np.abs(field(P, dz)) ** 2
        j = int(np.argmax(I))
        if best is None or I[j] > best[0]:
            best = (I[j], dz, P[j])
    I0, dz0, p0 = best
    for _ in range(2):
        g = np.linspace(-0.06, 0.06, 7)
        cands = []
        for dz in dz0 + np.linspace(-0.12, 0.12, 5):
            GX, GY = np.meshgrid(p0[0] + g, p0[1] + g, indexing="ij")
            P = np.column_stack([GX.ravel(), GY.ravel()])
            I = np.abs(field(P, dz)) ** 2
            j = int(np.argmax(I))
            cands.append((I[j], dz, P[j]))
        I0, dz0, p0 = max(cands, key=lambda t: t[0])
    rr = np.linspace(-1.6, 1.6, 161)
    Ix = np.abs(field(np.column_stack([p0[0] + rr, np.full_like(rr, p0[1])]), dz0)) ** 2
    Iy = np.abs(field(np.column_stack([np.full_like(rr, p0[0]), p0[1] + rr]), dz0)) ** 2

    def hw(I):
        I = I / I.max()
        c0 = int(np.argmax(I))
        rgt = c0 + int(np.argmax(I[c0:] < np.exp(-2)))
        lft = c0 - int(np.argmax(I[: c0 + 1][::-1] < np.exp(-2)))
        return float(0.5 * (rr[rgt] - rr[lft]))

    w_ap = z_weight = (WD / Lj**2)
    out = {"peak": float(I0), "dz_um": float(dz0), "x_um": float(p0[0]), "w_radial_um": hw(Ix), "w_tangential_um": hw(Iy),
           "s_max": float(s_loc.max()), "s_mean": float(np.sum(w_ap * s_loc) / np.sum(w_ap)),
           "frac_s_above_0.6": float(np.sum(w_ap * (s_loc > 0.6)) / np.sum(w_ap))}
    if eta_s is not None:
        S, E = eta_s
        out["eta_site"] = float(np.sum(w_ap * np.interp(s_loc, S, E)) / np.sum(w_ap))
    return out


def lens_field_radius(WD, na, lam, f_guess):
    best = None
    for frac in (0.9, 0.95, 1.0, 1.05, 1.1):
        sp = frac * f_guess
        d = D.make(D.DoubletSpec(na, [[WD, "air"]], [[sp, "air"]], n_pts=2001), lam)
        R1, f = d.info["R1_um"], d.info["f_um"]
        s = D.system(d, lam)
        h_test = 0.12 * f
        q = D.quality(s, R1, lam, D.theta_for_height(s, R1, lam, h_test, f), n_rays=32)
        if best is None or q["strehl"] > best[0]:
            best = (q["strehl"], sp, d, s)
    _, sp, d, s = best
    R1, f = d.info["R1_um"], d.info["f_um"]
    hs = np.linspace(0, 0.4 * f, 41)
    st = []
    for h in hs:
        try:
            st.append(D.quality(s, R1, lam, D.theta_for_height(s, R1, lam, h, f), n_rays=32)["strehl"])
        except Exception:
            st.append(0.0)
    st = np.array(st)
    bad = np.where(st < 0.8)[0]
    r08 = float(hs[bad[0] - 1]) if len(bad) else float(hs[-1])
    return {"spacing_um": float(sp), "f_um": float(f), "R1_um": float(R1), "field_radius_strehl08_um": r08,
            "strehl_vs_height": list(zip(hs.tolist(), st.tolist()))}


def main(cfg_path):
    cfg = load_config(cfg_path)
    out = run_dir(cfg["run_id"])
    write_params(out, cfg)
    path = out / "scaling.json"
    res = json.load(open(path)) if path.exists() else {}
    lam, na = cfg["wavelength_um"], cfg["na"]
    eta_s = eta_curve(cfg)
    res["eta_curve"] = {"s": eta_s[0].tolist(), "eta": eta_s[1].tolist()}
    for WD in cfg["working_distances_um"]:
        R = WD * na / np.sqrt(1 - na**2)
        key = f"WD{WD:.0f}"
        g = res.get(key, {"WD_um": WD, "aperture_radius_um": R, "sites": {}})
        for r_site in cfg["site_radii_um"]:
            k = f"{r_site:.0f}"
            if k in g["sites"]:
                continue
            t0 = time.time()
            g["sites"][k] = spot(WD, R, r_site, lam, cfg["n_aperture"], eta_s)
            g["sites"][k]["seconds"] = time.time() - t0
            print(key, k, {a: (round(b, 4) if isinstance(b, float) else b) for a, b in g["sites"][k].items()}, flush=True)
            res[key] = g
            write_json(path, res)
        if "lens" not in g and cfg.get("lens", True):
            g["lens"] = lens_field_radius(WD, cfg["lens_na"], cfg["lens_wavelength_um"], WD / np.sqrt(1 - cfg["lens_na"] ** 2))
            print(key, "lens field radius (Strehl >= 0.8):", round(g["lens"]["field_radius_strehl08_um"]), "um", flush=True)
            res[key] = g
            write_json(path, res)


if __name__ == "__main__":
    main(sys.argv[1])
