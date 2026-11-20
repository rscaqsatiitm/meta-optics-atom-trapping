"""RCWA against FDTD with converged references: three FDTD grids, both polarisations, RCWA digits.

Usage: python scripts/rcwa_fdtd_convergence.py [--break]

Supersedes the comparison of scripts/rcwa_fdtd_check.py, which had three defects found on
2026-11-20: (i) fdtd.PeriodicBoundary identifies the first and the last cell (E[0] = E[-1],
H[-1] = H[0]), so a grid of Nx cells has a period of Nx - 1 cells; that script used Nx = p / d,
i.e. periods of 0.49 um (10 nm grid) and 0.495 um (5 nm grid) instead of 0.50 um, an error
first order in d; (ii) it centred the ridges on cell centres, so on the 10 nm grid the ridges
had 14 / 24 / 34 cells (0.01 um too narrow; exact on the 5 nm grid); (i) and (ii) plausibly
caused much of its C2 failure; (iii) it covered only the field parallel to the ridges. Here
Nx = p / d + 1 (period exactly p; cells 0 and Nx - 1 both air), the ridge fills whole cells
from a cell boundary, and the grids 25, 12.5 and 6.25 nm put every width, the period and every
layer on whole cells (asserted).

Case as in rcwa_fdtd_check.py: substrate 1.45, ridges 2.0, air, period 0.50 um, height 1.00 um,
1.075 um, normal incidence from the substrate, ridge widths 0.15 / 0.25 / 0.35 um; compared
quantity t(w) / t(0.25 um). TE: E parallel to the ridges (FDTD E_z, fmmax t_y, from C4 of
rcwa_fdtd_check.py); TM: E perpendicular (FDTD E_x, fmmax t_x).

Checks and thresholds (declared before the first run; never changed after a result):
  C1  FDTD uniform layer against Airy, grids 25 / 12.5 / 6.25 nm:  |amp - 1| <= 0.01, |dphi| <= 0.03 rad.
  C2  FDTD sampling, per polarisation and width (Richardson from the three grids, coarse to fine
      q1, q2, q3): observed order p = log2(|q1 - q2| / |q2 - q3|) within [0.7, 3.5], and
      extrapolation correction |q* - q3| <= 0.01 in amplitude and <= 0.02 rad in phase, with
      q* = q3 + (q3 - q2) / (2^p - 1). If |q2 - q3| < 1e-3 the grids are converged and the
      order test is skipped.
  C5  RCWA digits (fmmax pol_fourier): 317 against 437 terms (grid 200) and grid 200 against
      400 (437 terms), every width and both polarisations: |damp| <= 0.002, |dphi| <= 0.005 rad.
  C3  RCWA (437 terms, grid 400) against the extrapolated FDTD q*, both polarisations:
      |damp| <= 0.02, |dphi| <= 0.05 rad. Conclusive only if C1, C2 and C5 pass.
--break: C1 with n = 2.05 in the analytic layer; C2 with the coarsest-grid ratio of another width;
C5 with 41 against 437 terms; C3 with RCWA TM against FDTD TE. Every check must then fail.
Output: results/rcwa_fdtd_convergence/{metrics,metrics_break}.json, FDTD phasors cached in
fdtd_phasors.json (the break run reuses them).

--fine (declared 2026-11-20, after the first run, before this one): the same checks C1, C2, C5
and C3 with the same thresholds, on the three finest grids 12.5 / 6.25 / 3.125 nm (all widths,
the period and the layers on whole cells). Reason: in the first run C1 failed only on the 25 nm
grid (0.0320 rad against 0.03; 0.0079 and 0.0020 rad on the finer grids). That failure stays on
record (metrics.json); this run decides whether the rule is met on a finer grid sequence.
The missing FDTD runs are made four at a time in parallel (the solver uses one core).
Output: metrics_fine.json, metrics_fine_break.json.

--break-c2-mid (added 2026-11-20 after the --fine --break run): the declared C2 break corrupts
only the coarsest grid, which C2 does not read when the two finer grids agree to < 1e-3 (the
converged branch). On the fine grids every case took that branch, so the declared C2 break did
not fire (metrics_fine_break.json: C2 PASS). This control corrupts the middle grid (the ratio of
another width), i.e. the clause C2 actually used, and C2 must fail. Added after the result; the
declared break and its non-failure stay on record. Output: metrics_fine_break_c2_mid.json.
"""

import json
import sys
import time

import numpy as np

from metatrap.runio import REPO_ROOT, write_json

LAM, P, H = 1.075, 0.50, 1.00
N_SUB, N_RIDGE = 1.45, 2.0
WIDTHS = [0.15, 0.25, 0.35]
W_REF = 0.25
GRIDS = [0.025, 0.0125, 0.00625]  # every width, the period and every layer are whole numbers of cells
GRIDS_FINE = [0.0125, 0.00625, 0.003125]
THR = {"C1_amp": 0.01, "C1_phi": 0.03, "C2_order": [0.7, 3.5], "C2_amp": 0.01, "C2_phi": 0.02, "C2_converged": 1e-3,
       "C5_amp": 0.002, "C5_phi": 0.005, "C3_amp": 0.02, "C3_phi": 0.05}
COMP = {"TE": 2, "TM": 0}  # FDTD field component: E_z parallel to the ridges, E_x perpendicular


def airy_t(n1, n2, n3, h, lam):
    t12, t23 = 2 * n1 / (n1 + n2), 2 * n2 / (n2 + n3)
    r12, r23 = (n1 - n2) / (n1 + n2), (n2 - n3) / (n2 + n3)
    d = 2 * np.pi / lam * n2 * h
    return t12 * t23 * np.exp(1j * d) / (1 + r12 * r23 * np.exp(2j * d))


def fdtd_t(width_um, d_um, layer, pol):
    """x-averaged phasor of the transmitted E (component of pol) above the layer, soft plane-wave source."""
    import fdtd

    fdtd.set_backend("numpy")
    d = d_um * 1e-6
    n_per = int(round(LAM / (0.7 * d_um)))
    courant = LAM / (n_per * d_um)
    lengths = {"pml": 0.8, "sub": 1.0, "h": H, "air": 1.0, "src": 0.3, "mon": 0.8, "p": P, "w": width_um}
    cells = {k: int(round(v / d_um)) for k, v in lengths.items()}
    assert all(abs(cells[k] * d_um - v) < 1e-9 for k, v in lengths.items()), "geometry not on the grid"
    npml, nsub, nh, nair, ncell = cells["pml"], cells["sub"], cells["h"], cells["air"], cells["p"]
    nx, ny = ncell + 1, 2 * npml + nsub + nh + nair  # periodic boundary copies cell 0 from cell nx - 1
    eps = np.ones((nx, ny, 1))
    eps[:, : npml + nsub, 0] = N_SUB**2
    y0, y1 = npml + nsub, npml + nsub + nh
    if layer == "slab":
        eps[:, y0:y1, 0] = N_RIDGE**2
    elif layer == "grating":
        # unique cells 1 .. ncell; the ridge fills cells 1 .. w/d exactly (a lateral shift of the
        # grating leaves the zeroth-order coefficient unchanged at normal incidence); cells 0 and
        # nx - 1 (identified by the periodic boundary) stay air.
        eps[1 : 1 + cells["w"], y0:y1, 0] = N_RIDGE**2
    assert np.allclose(eps[0], eps[-1])
    g = fdtd.Grid((nx, ny, 1), grid_spacing=d, permittivity=eps, courant_number=courant)
    g[0, :, :] = fdtd.PeriodicBoundary(name="xper")
    g[:, 0:npml, :] = fdtd.PML(name="pml_lo")
    g[:, -npml:, :] = fdtd.PML(name="pml_hi")
    comp = COMP[pol]
    y_src, y_mon = npml + cells["src"], y1 + cells["mon"]
    omega = 2 * np.pi / n_per
    n_ramp, n_win = 10 * n_per, 20 * n_per
    n_total = n_ramp + 30 * n_per + 2 * n_win
    acc = [np.zeros(nx - 1, complex), np.zeros(nx - 1, complex)]
    for q in range(n_total):
        g.step()
        ramp = np.sin(0.5 * np.pi * min(1.0, q / n_ramp)) ** 2
        g.E[:, y_src, 0, comp] += ramp * np.sin(omega * q)
        k = (q - (n_total - 2 * n_win)) // n_win
        if 0 <= k < 2:
            acc[k] += g.E[1:, y_mon, 0, comp] * np.exp(1j * omega * q)  # unique cells only
    a1, a2 = (2 / n_win) * acc[0].mean(), (2 / n_win) * acc[1].mean()
    return a2, abs(a2 / a1 - 1)


def rcwa_ratios(n_terms, n_grid, period=P, widths=WIDTHS):
    from metatrap.rcwa.solver import solve_eps

    xc = (np.arange(n_grid) + 0.5) / n_grid * period - period / 2
    eps = np.ones((len(widths), n_grid, n_grid))
    for i, w in enumerate(widths):
        eps[i][np.abs(xc) < w / 2 - 1e-9, :] = N_RIDGE**2
    r = solve_eps(eps.astype(complex), period, H, LAM, N_SUB**2, 1.0, n_terms=n_terms, formulation="pol_fourier")
    i0 = widths.index(W_REF)
    return {"TE": {w: r.t_y[i] / r.t_y[i0] for i, w in enumerate(widths)}, "TM": {w: r.t_x[i] / r.t_x[i0] for i, w in enumerate(widths)}}


def cmp(a, b):
    return float(abs(abs(a) - abs(b))), float(abs(np.angle(a / b)))


def _run(job):
    key, w, d, layer, pol = job
    v, steady = fdtd_t(w, d, layer, pol)
    return key, [v.real, v.imag, steady]


def main(brk, fine=False, brk_mid=False):
    global GRIDS
    if fine:
        GRIDS = GRIDS_FINE
    t0 = time.time()
    out = REPO_ROOT / "results" / "rcwa_fdtd_convergence"
    out.mkdir(parents=True, exist_ok=True)
    case = {"lam": LAM, "pitch": P, "height": H, "n_sub": N_SUB, "n_ridge": N_RIDGE, "widths": WIDTHS, "grids": GRIDS, "nx": "p/d+1"}
    base = {k: v for k, v in case.items() if k != "grids"}  # cache entries are per run ("d:name")
    cache_path = out / "fdtd_phasors.json"
    cache = json.load(open(cache_path)) if cache_path.exists() else {}
    old = json.loads(cache["key"]) if "key" in cache else None
    if old is None or {k: v for k, v in old.items() if k != "grids"} != base:
        cache = {}
    cache["key"] = json.dumps(case, sort_keys=True)
    jobs = [(f"{d}:{name}", w, d, layer, pol) for d in GRIDS
            for name, layer, w, pol in [("slab", "slab", 0, "TE"), ("air", "air", 0, "TE")] +
            [(f"{pol}:{w}", "grating", w, pol) for pol in ("TE", "TM") for w in WIDTHS]]
    jobs = [j for j in jobs if j[0] not in cache]
    if jobs:
        from multiprocessing import Pool

        with Pool(4) as pool:
            for key, val in pool.imap_unordered(_run, jobs):
                cache[key] = val
                json.dump(cache, open(cache_path, "w"), indent=1)
                print(f"FDTD {key} done ({time.time() - t0:.0f}s, steady-state change {val[2]:.1e})", flush=True)
    ph = {k: complex(v[0], v[1]) for k, v in cache.items() if k != "key"}
    res = {"thresholds": THR, "case": case, "break": brk,
           "fdtd_max_steady_state_change": max(v[2] for k, v in cache.items() if k != "key")}
    # C1
    ana = airy_t(N_SUB, 2.05 if brk else N_RIDGE, 1.0, H, LAM) / airy_t(N_SUB, 1.0, 1.0, H, LAM)
    c1 = {}
    for d in GRIDS:
        r = ph[f"{d}:slab"] / ph[f"{d}:air"]
        c1[str(d)] = {"amp_err": float(abs(abs(r / ana) - 1)), "phi_err": float(abs(np.angle(r / ana)))}
    c1_pass = all(v["amp_err"] <= THR["C1_amp"] and v["phi_err"] <= THR["C1_phi"] for v in c1.values())
    # C2: Richardson per polarisation and width.
    q = {pol: {d: {w: ph[f"{d}:{pol}:{w}"] / ph[f"{d}:{pol}:{W_REF}"] for w in WIDTHS} for d in GRIDS} for pol in COMP}
    c2, qstar = {}, {pol: {} for pol in COMP}
    for pol in COMP:
        for i, w in enumerate(WIDTHS):
            q1 = q[pol][GRIDS[0]][WIDTHS[(i + 1) % 3] if brk else w]
            q2, q3 = q[pol][GRIDS[1]][WIDTHS[(i + 1) % 3] if brk_mid else w], q[pol][GRIDS[2]][w]
            d1, d2 = abs(q1 - q2), abs(q2 - q3)
            if w == W_REF:
                qstar[pol][w] = 1.0 + 0j
                continue
            if d2 < THR["C2_converged"]:
                p, qs, ok_order = float("nan"), q3, True
            else:
                p = float(np.log2(d1 / d2)) if d1 > 0 else float("nan")
                ok_order = THR["C2_order"][0] <= p <= THR["C2_order"][1]
                qs = q3 + (q3 - q2) / (2**p - 1) if ok_order else q3
            a, f = cmp(qs, q3)
            qstar[pol][w] = qs
            c2[f"{pol}:{w}"] = {"order": p, "order_ok": bool(ok_order), "corr_amp": a, "corr_phi": f,
                                "q1": [abs(q1), float(np.angle(q1))], "q2": [abs(q2), float(np.angle(q2))],
                                "q3": [abs(q3), float(np.angle(q3))], "qstar": [abs(qs), float(np.angle(qs))],
                                "diff_mid_fine_phi": cmp(q2, q3)[1]}
    c2_pass = all(v["order_ok"] and v["corr_amp"] <= THR["C2_amp"] and v["corr_phi"] <= THR["C2_phi"] for v in c2.values())
    print(f"FDTD done; RCWA ({time.time() - t0:.0f}s)", flush=True)
    # C5: RCWA digits.
    r437 = rcwa_ratios(437, 200)
    r317 = rcwa_ratios(41 if brk else 317, 200)
    r437g = rcwa_ratios(437, 400)
    c5 = {}
    for pol in COMP:
        for w in WIDTHS:
            if w == W_REF:
                continue
            a1, f1 = cmp(r317[pol][w], r437[pol][w])
            a2, f2 = cmp(r437g[pol][w], r437[pol][w])
            c5[f"{pol}:{w}"] = {"terms_amp": a1, "terms_phi": f1, "grid_amp": a2, "grid_phi": f2}
    c5_pass = all(v["terms_amp"] <= THR["C5_amp"] and v["terms_phi"] <= THR["C5_phi"] and v["grid_amp"] <= THR["C5_amp"]
                  and v["grid_phi"] <= THR["C5_phi"] for v in c5.values())
    # C3: RCWA (437 terms, grid 400) against FDTD q*.
    c3 = {}
    for pol in COMP:
        for w in WIDTHS:
            if w == W_REF:
                continue
            fd = qstar["TE" if (brk and pol == "TM") else pol][w]
            a, f = cmp(r437g[pol][w], fd)
            c3[f"{pol}:{w}"] = {"amp_err": a, "phi_err": f, "rcwa": [abs(r437g[pol][w]), float(np.angle(r437g[pol][w]))],
                                "fdtd_extrapolated": [abs(fd), float(np.angle(fd))]}
    c3_pass = all(v["amp_err"] <= THR["C3_amp"] and v["phi_err"] <= THR["C3_phi"] for v in c3.values())
    res.update({"C1": c1, "C1_pass": c1_pass, "C2": c2, "C2_pass": c2_pass, "C5": c5, "C5_pass": c5_pass, "C3": c3,
                "C3_pass": c3_pass, "C3_conclusive": bool(c1_pass and c2_pass and c5_pass), "seconds": time.time() - t0})
    verdicts = {"C1": c1_pass, "C2": c2_pass, "C5": c5_pass, "C3": c3_pass}
    for k, v in verdicts.items():
        print(f"{k}: {'PASS' if v else 'FAIL'}  {json.dumps(res[k])}", flush=True)
    ok = (not any(verdicts.values())) if brk else all(verdicts.values())
    print(f"{'break run: every check must fail' if brk else 'normal run: every check must pass'} -> {'OK' if ok else 'NOT OK'}")
    tag = "metrics_fine" if fine else "metrics"
    if brk_mid:  # C2 control on the middle grid only (see the docstring)
        write_json(out / f"{tag}_break_c2_mid.json", {"C2": res["C2"], "C2_pass": res["C2_pass"]})
        print(f"C2 middle-grid control: C2 {'PASS' if c2_pass else 'FAIL'} (must fail)")
        return 0 if not c2_pass else 1
    write_json(out / (f"{tag}_break.json" if brk else f"{tag}.json"), res)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main("--break" in sys.argv, fine="--fine" in sys.argv, brk_mid="--break-c2-mid" in sys.argv))
