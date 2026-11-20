"""RCWA against 3D FDTD for square pillars (the geometry of the meta-atom libraries).

Usage: python scripts/rcwa_fdtd_pillars.py [--break]

scripts/rcwa_fdtd_convergence.py closed the RCWA-FDTD comparison for a lamellar (1D-periodic)
grating. The libraries use square pillars on a square lattice, where the Fourier treatment of a
2D-discontinuous permittivity (the part of fmmax's pol_fourier formulation that differs most from
the 1D case) is exercised. This script repeats the comparison for square pillars with a 3D FDTD
(the `fdtd` package; periodic boundaries in x and y identify the first and last cell, so the
grid has p / d + 1 cells per period; PML above and below).

Case: substrate 1.45, pillars 2.0, air, square lattice 0.50 um, height 1.00 um, 1.075 um, normal
incidence from the substrate, E along x (square pillars are polarisation independent), pillar
widths 0.15 / 0.25 / 0.35 um, filling whole cells from a cell corner (a lateral shift leaves the
zeroth order unchanged). Compared quantity: t(w) / t(0.25 um). Grids 25, 16.67 and 12.5 nm, on
which every width, the period and every layer are whole numbers of cells (asserted). A finer grid
is out of reach here (pilot: 27 min per run at 12.5 nm, 16 times that at 6.25 nm).

Checks and thresholds (declared before the first run; never changed after a result):
  P1  FDTD uniform layer against Airy: on the finest grid |amp - 1| <= 0.01 and |dphi| <= 0.03 rad,
      and the phase error over the three grids must fall with an observed order within
      [1.5, 2.5] (the Yee scheme is second order). The 0.03 rad threshold is applied to the finest
      grid only, unlike rcwa_fdtd_convergence.py, where the coarsest grid's dispersion error
      (0.032 rad at 25 nm) failed it although the scheme was in its asymptotic regime; this choice
      is made knowing that result.
  P2  FDTD sampling per width (Richardson for unequal steps h1 > h2 > h3): observed order p from
      (q1 - q2) / (q2 - q3) = (h1^p - h2^p) / (h2^p - h3^p), within [0.7, 3.5], and correction
      |q* - q3| <= 0.01 in amplitude and <= 0.02 rad in phase with
      q* = q3 + (q3 - q2) h3^p / (h2^p - h3^p). If |q2 - q3| < 1e-3 the grids count as converged
      and q* = q3.
  P5  RCWA digits (fmmax pol_fourier): 317 against 437 terms (grid 200) and grid 200 against 400
      (437 terms), every width: |damp| <= 0.002, |dphi| <= 0.005 rad.
  P3  RCWA (437 terms, grid 400) against FDTD q*: |damp| <= 0.02, |dphi| <= 0.05 rad. Conclusive only
      if P1, P2 and P5 pass.
--break: P1 with n = 2.05 in the analytic layer; P2 with the middle-grid ratio of another width
(the clause used in the converged branch) and the coarsest-grid ratio of another width; P5 with
41 against 437 terms; P3 with RCWA evaluated at widths 0.025 um (one coarse cell) larger. Every
check must fail.
Output: results/rcwa_fdtd_pillars/{metrics,metrics_break}.json, FDTD phasors in fdtd_phasors.json.
"""

import json
import sys
import time

import numpy as np
from scipy.optimize import brentq

from metatrap.runio import REPO_ROOT, write_json

LAM, P, H = 1.075, 0.50, 1.00
N_SUB, N_PIL = 1.45, 2.0
WIDTHS = [0.15, 0.25, 0.35]
W_REF = 0.25
GRIDS = [0.025, 1 / 60, 0.0125]
THR = {"P1_amp": 0.01, "P1_phi": 0.03, "P1_order": [1.5, 2.5], "P2_order": [0.7, 3.5], "P2_amp": 0.01, "P2_phi": 0.02,
       "P2_converged": 1e-3, "P5_amp": 0.002, "P5_phi": 0.005, "P3_amp": 0.02, "P3_phi": 0.05}


def airy_t(n1, n2, n3, h, lam):
    t12, t23 = 2 * n1 / (n1 + n2), 2 * n2 / (n2 + n3)
    r12, r23 = (n1 - n2) / (n1 + n2), (n2 - n3) / (n2 + n3)
    d = 2 * np.pi / lam * n2 * h
    return t12 * t23 * np.exp(1j * d) / (1 + r12 * r23 * np.exp(2j * d))


def fdtd_t(width_um, d_um, layer):
    import fdtd

    fdtd.set_backend("numpy")
    lengths = {"pml": 0.8, "sub": 1.0, "h": H, "air": 1.0, "src": 0.3, "mon": 0.8, "p": P, "w": width_um}
    cells = {k: int(round(v / d_um)) for k, v in lengths.items()}
    assert all(abs(cells[k] * d_um - v) < 1e-9 for k, v in lengths.items()), "geometry not on the grid"
    n_per = int(round(LAM / (0.56 * d_um)))  # Courant <= 1 / sqrt(3)
    courant = LAM / (n_per * d_um)
    npml, nsub, nh, nair, ncell = cells["pml"], cells["sub"], cells["h"], cells["air"], cells["p"]
    nx, nz = ncell + 1, 2 * npml + nsub + nh + nair
    eps = np.ones((nx, nx, nz))
    eps[:, :, : npml + nsub] = N_SUB**2
    z0, z1 = npml + nsub, npml + nsub + nh
    if layer == "slab":
        eps[:, :, z0:z1] = N_PIL**2
    elif layer == "pillar":
        m = cells["w"]
        eps[1 : 1 + m, 1 : 1 + m, z0:z1] = N_PIL**2  # cells 0 and nx - 1 (identified) stay air
    assert np.allclose(eps[0], eps[-1]) and np.allclose(eps[:, 0], eps[:, -1])
    g = fdtd.Grid((nx, nx, nz), grid_spacing=d_um * 1e-6, permittivity=eps, courant_number=courant)
    g[0, :, :] = fdtd.PeriodicBoundary(name="xper")
    g[:, 0, :] = fdtd.PeriodicBoundary(name="yper")
    g[:, :, 0:npml] = fdtd.PML(name="pml_lo")
    g[:, :, -npml:] = fdtd.PML(name="pml_hi")
    z_src, z_mon = npml + cells["src"], z1 + cells["mon"]
    omega = 2 * np.pi / n_per
    n_ramp, n_win = 10 * n_per, 20 * n_per
    n_total = n_ramp + 30 * n_per + 2 * n_win
    acc = [0j, 0j]
    for q in range(n_total):
        g.step()
        ramp = np.sin(0.5 * np.pi * min(1.0, q / n_ramp)) ** 2
        g.E[:, :, z_src, 0] += ramp * np.sin(omega * q)
        k = (q - (n_total - 2 * n_win)) // n_win
        if 0 <= k < 2:
            acc[k] += g.E[1:, 1:, z_mon, 0].mean() * np.exp(1j * omega * q)  # unique cells, zeroth order
    a1, a2 = (2 / n_win) * acc[0], (2 / n_win) * acc[1]
    return a2, abs(a2 / a1 - 1)


def _run(job):
    key, w, d, layer = job
    v, steady = fdtd_t(w, d, layer)
    return key, [v.real, v.imag, steady]


def rcwa_ratios(n_terms, n_grid, widths):
    from metatrap.rcwa.solver import solve_eps

    xc = (np.arange(n_grid) + 0.5) / n_grid * P - P / 2
    eps = np.ones((len(widths), n_grid, n_grid))
    for i, w in enumerate(widths):
        inside = np.abs(xc) < w / 2 - 1e-9
        eps[i][np.ix_(inside, inside)] = N_PIL**2
    r = solve_eps(eps.astype(complex), P, H, LAM, N_SUB**2, 1.0, n_terms=n_terms, formulation="pol_fourier")
    return r.t_x


def cmp(a, b):
    return float(abs(abs(a) - abs(b))), float(abs(np.angle(a / b)))


def richardson_order(e1, e2, h):
    """Observed order p for unequal steps from the ratio e1 / e2 = (h1^p - h2^p) / (h2^p - h3^p)."""
    target = e1 / e2
    f = lambda p: (h[0] ** p - h[1] ** p) / (h[1] ** p - h[2] ** p) - target  # noqa: E731
    try:
        return float(brentq(f, 0.05, 8.0))
    except ValueError:
        return float("nan")


def main(brk):
    t0 = time.time()
    out = REPO_ROOT / "results" / "rcwa_fdtd_pillars"
    out.mkdir(parents=True, exist_ok=True)
    case = {"lam": LAM, "pitch": P, "height": H, "n_sub": N_SUB, "n_pillar": N_PIL, "widths": WIDTHS, "grids": GRIDS}
    cache_path = out / "fdtd_phasors.json"
    cache = json.load(open(cache_path)) if cache_path.exists() else {}
    if cache.get("key") != json.dumps(case, sort_keys=True):
        cache = {"key": json.dumps(case, sort_keys=True)}
    jobs = [(f"{d:.6f}:{name}", w, d, layer) for d in GRIDS
            for name, layer, w in [("slab", "slab", 0), ("air", "air", 0)] + [(f"pillar:{w}", "pillar", w) for w in WIDTHS]]
    jobs = [j for j in jobs if j[0] not in cache]
    jobs.sort(key=lambda j: -1 / j[2])  # finest (longest) first, for load balance
    if jobs:
        from multiprocessing import Pool

        with Pool(4) as pool:
            for key, val in pool.imap_unordered(_run, jobs):
                cache[key] = val
                json.dump(cache, open(cache_path, "w"), indent=1)
                print(f"FDTD {key} done ({time.time() - t0:.0f}s, steady-state change {val[2]:.1e})", flush=True)
    ph = {k: complex(v[0], v[1]) for k, v in cache.items() if k != "key"}
    gk = [f"{d:.6f}" for d in GRIDS]
    res = {"thresholds": THR, "case": case, "break": brk, "fdtd_max_steady_state_change": max(v[2] for k, v in cache.items() if k != "key")}
    # P1
    ana = airy_t(N_SUB, 2.05 if brk else N_PIL, 1.0, H, LAM) / airy_t(N_SUB, 1.0, 1.0, H, LAM)
    p1 = {}
    for g in gk:
        r = ph[f"{g}:slab"] / ph[f"{g}:air"]
        p1[g] = {"amp_err": float(abs(abs(r / ana) - 1)), "phi_err": float(abs(np.angle(r / ana)))}
    errs = [p1[g]["phi_err"] for g in gk]
    order1 = float(np.log(errs[0] / errs[2]) / np.log(GRIDS[0] / GRIDS[2])) if errs[2] > 0 else float("nan")
    p1_pass = (p1[gk[2]]["amp_err"] <= THR["P1_amp"] and p1[gk[2]]["phi_err"] <= THR["P1_phi"]
               and THR["P1_order"][0] <= order1 <= THR["P1_order"][1])
    # P2
    q = {g: {w: ph[f"{g}:pillar:{w}"] / ph[f"{g}:pillar:{W_REF}"] for w in WIDTHS} for g in gk}
    p2, qstar = {}, {W_REF: 1.0 + 0j}
    for i, w in enumerate(WIDTHS):
        if w == W_REF:
            continue
        other = WIDTHS[(i + 1) % 3] if WIDTHS[(i + 1) % 3] != W_REF else WIDTHS[(i + 2) % 3]
        q1 = q[gk[0]][other if brk else w]
        q2 = q[gk[1]][other if brk else w]
        q3 = q[gk[2]][w]
        e1, e2 = abs(q1 - q2), abs(q2 - q3)
        if e2 < THR["P2_converged"]:
            p, qs, ok = float("nan"), q3, True
        else:
            p = richardson_order(e1, e2, GRIDS)
            ok = THR["P2_order"][0] <= p <= THR["P2_order"][1]
            qs = q3 + (q3 - q2) * GRIDS[2] ** p / (GRIDS[1] ** p - GRIDS[2] ** p) if ok else q3
        a, f = cmp(qs, q3)
        qstar[w] = qs
        p2[str(w)] = {"order": p, "order_ok": bool(ok), "corr_amp": a, "corr_phi": f, "diff_mid_fine": e2,
                      "q": [[abs(x), float(np.angle(x))] for x in (q1, q2, q3)], "qstar": [abs(qs), float(np.angle(qs))]}
    p2_pass = all(v["order_ok"] and v["corr_amp"] <= THR["P2_amp"] and v["corr_phi"] <= THR["P2_phi"] for v in p2.values())
    print(f"FDTD analysed; RCWA ({time.time() - t0:.0f}s)", flush=True)
    # P5
    t437 = rcwa_ratios(437, 200, WIDTHS)
    t317 = rcwa_ratios(41 if brk else 317, 200, WIDTHS)
    t437g = rcwa_ratios(437, 400, WIDTHS)
    i0 = WIDTHS.index(W_REF)
    p5 = {}
    for i, w in enumerate(WIDTHS):
        if w == W_REF:
            continue
        a1, f1 = cmp(t317[i] / t317[i0], t437[i] / t437[i0])
        a2, f2 = cmp(t437g[i] / t437g[i0], t437[i] / t437[i0])
        p5[str(w)] = {"terms_amp": a1, "terms_phi": f1, "grid_amp": a2, "grid_phi": f2}
    p5_pass = all(v["terms_amp"] <= THR["P5_amp"] and v["terms_phi"] <= THR["P5_phi"] and v["grid_amp"] <= THR["P5_amp"]
                  and v["grid_phi"] <= THR["P5_phi"] for v in p5.values())
    # P3
    tr = rcwa_ratios(437, 400, [w + 0.025 for w in WIDTHS]) if brk else t437g
    p3 = {}
    for i, w in enumerate(WIDTHS):
        if w == W_REF:
            continue
        rr = tr[i] / tr[i0]
        a, f = cmp(rr, qstar[w])
        p3[str(w)] = {"amp_err": a, "phi_err": f, "rcwa": [abs(rr), float(np.angle(rr))], "fdtd": [abs(qstar[w]), float(np.angle(qstar[w]))]}
    p3_pass = all(v["amp_err"] <= THR["P3_amp"] and v["phi_err"] <= THR["P3_phi"] for v in p3.values())
    res.update({"P1": p1, "P1_order": order1, "P1_pass": bool(p1_pass), "P2": p2, "P2_pass": p2_pass, "P5": p5, "P5_pass": p5_pass,
                "P3": p3, "P3_pass": p3_pass, "P3_conclusive": bool(p1_pass and p2_pass and p5_pass), "seconds": time.time() - t0})
    verdicts = {"P1": bool(p1_pass), "P2": p2_pass, "P5": p5_pass, "P3": p3_pass}
    for k, v in verdicts.items():
        extra = f" (order {order1:.2f})" if k == "P1" else ""
        print(f"{k}: {'PASS' if v else 'FAIL'}{extra}  {json.dumps(res[k])}", flush=True)
    ok = (not any(verdicts.values())) if brk else all(verdicts.values())
    print(f"{'break run: every check must fail' if brk else 'normal run: every check must pass'} -> {'OK' if ok else 'NOT OK'}")
    write_json(out / ("metrics_break.json" if brk else "metrics.json"), res)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main("--break" in sys.argv))
