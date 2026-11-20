"""Independent route for the RCWA meta-atom response: 2D FDTD of a lamellar grating, plus two
analytic limits.

SUPERSEDED (2026-11-20) by scripts/rcwa_fdtd_convergence.py. The FDTD set-up here has two
geometry defects: fdtd.PeriodicBoundary identifies the first and last cell, so Nx = p / d cells
give a period of p - d (0.49 / 0.495 um on the 10 / 5 nm grids); and ridges centred on cell
centres have 14 / 24 / 34 cells on the 10 nm grid (0.01 um too narrow). Its C2 failure and its
C3 agreement are therefore not evidence either way. Kept, unchanged otherwise, as the record of
results/rcwa_fdtd_check.

Usage: python scripts/rcwa_fdtd_check.py [--break]

Why. fmmax and grcwa are both Fourier modal methods (RCWA) and share the layers where their
most likely errors live (Fourier factorisation and truncation of the permittivity, the same
staircased permittivity grid built by rcwa/cell.py, the same material values); see
docs/rcwa_lineage.md. Their agreement cannot rule out an error in those layers. This script
compares fmmax (the production solver, pol_fourier formulation, as in the libraries) with a
method different in kind: finite-difference time-domain (the third-party `fdtd` package,
pure numpy; no Fourier expansion, no eigenmodes, no S-matrix), on a 1D-periodic lamellar
grating with the electric field parallel to the ridges.

Case (fixed constants, not taken from materials.py): substrate n = 1.45, ridges n = 2.0,
air above and between the ridges, period 0.50 um, height 1.00 um, wavelength 1.075 um,
normal incidence from the substrate, ridge widths 0.15, 0.25 and 0.35 um. The compared
quantity is convention-free: the complex ratio t(w) / t(0.25 um) of zeroth-order transmission
coefficients (amplitude ratio and phase difference).

Checks and thresholds (declared before the first run; never changed after a result):
  C1  FDTD planted truth: uniform layer (n = 2.0, 1.0 um) against the analytic Airy result,
      as the ratio t(slab) / t(air layer):        |amp - 1| <= 0.01,  |dphi| <= 0.03 rad.
  C2  FDTD sampling: grid 10 nm against 5 nm, for the C3 ratios:
                                                  |damp| <= 0.01,     |dphi| <= 0.02 rad.
  C3  RCWA (fmmax, pol_fourier, 197 terms, 2.5 nm grid) against FDTD (5 nm) for the grating
      ratios:                                     |damp| <= 0.02,     |dphi| <= 0.05 rad.
  C4  Polarisation convention of the wrapper: at period 0.05 um (fill 0.5) fmmax must match
      the zeroth-order effective-medium layers, E || ridges with n_TE^2 = f n^2 + (1 - f) and
      E perpendicular with 1 / n_TM^2 = f / n^2 + (1 - f), each within |dphi| <= 0.10 rad
      (analytic Airy coefficient of the effective layer).
  If C1 or C2 fails, C3 is inconclusive (the FDTD route itself is not converged), not a pass.

--break applies one deliberate corruption per check (C1: n = 2.05 in the analytic layer;
C2: the 5 nm ratios compared with the 10 nm ones of a different width; C3: the RCWA
polarisation perpendicular to the ridges; C4: ridge index 2.1 in the effective-medium limit)
and must make every check fail. (A first version swapped the two effective indices for C4;
the break run showed that C4 infers the orientation from the data, so a swap cannot make it
fail, and that the C3 swap was undone by the same inference. Both breaks were replaced.)
Output: results/rcwa_fdtd_check/metrics.json (or metrics_break.json). Exit status 0 only if
every check passes (without --break) or every check fails (with --break).
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
THR = {"C1_amp": 0.01, "C1_phi": 0.03, "C2_amp": 0.01, "C2_phi": 0.02, "C3_amp": 0.02, "C3_phi": 0.05, "C4_phi": 0.10}


def airy_t(n1, n2, n3, h, lam):
    t12, t23 = 2 * n1 / (n1 + n2), 2 * n2 / (n2 + n3)
    r12, r23 = (n1 - n2) / (n1 + n2), (n2 - n3) / (n2 + n3)
    d = 2 * np.pi / lam * n2 * h
    return t12 * t23 * np.exp(1j * d) / (1 + r12 * r23 * np.exp(2j * d))


# ----------------------------------------------------------------------------- FDTD route
def fdtd_t(width_um, d_um, layer="grating", n_layer=N_RIDGE):
    """Complex zeroth-order E_z phasor above the layer for a soft plane-wave source in the substrate.

    layer: 'grating' (ridge of the given width per period), 'slab' (uniform n_layer) or 'air'.
    Ratios of this phasor between two runs with the same d are transmission ratios."""
    import fdtd

    fdtd.set_backend("numpy")
    c = 299792458.0
    d = d_um * 1e-6
    n_per = int(round(LAM / (0.7 * d_um)))  # steps per optical period (integer), Courant = LAM / (n_per d) <= 1 / sqrt(2)
    courant = LAM / (n_per * d_um)
    npml, nsub, nh, nair = int(round(0.8 / d_um)), int(round(1.0 / d_um)), int(round(H / d_um)), int(round(1.0 / d_um))
    nx, ny = int(round(P / d_um)), 2 * npml + nsub + nh + nair
    eps = np.ones((nx, ny, 1))
    eps[:, : npml + nsub, 0] = N_SUB**2
    y0, y1 = npml + nsub, npml + nsub + nh
    if layer == "slab":
        eps[:, y0:y1, 0] = n_layer**2
    elif layer == "grating":
        xc = (np.arange(nx) + 0.5) * d_um - P / 2
        ridge = np.abs(xc) < width_um / 2 - 1e-9
        eps[np.ix_(ridge, np.arange(y0, y1), [0])] = n_layer**2
    g = fdtd.Grid((nx, ny, 1), grid_spacing=d, permittivity=eps, courant_number=courant)
    g[0, :, :] = fdtd.PeriodicBoundary(name="xper")
    g[:, 0:npml, :] = fdtd.PML(name="pml_lo")
    g[:, -npml:, :] = fdtd.PML(name="pml_hi")
    y_src, y_mon = npml + int(round(0.3 / d_um)), y1 + int(round(0.8 / d_um))
    omega = 2 * np.pi / n_per  # per step
    n_ramp, n_win = 10 * n_per, 20 * n_per
    n_total = n_ramp + 30 * n_per + 2 * n_win
    acc = [np.zeros(nx, complex), np.zeros(nx, complex)]
    for q in range(n_total):
        g.step()
        ramp = np.sin(0.5 * np.pi * min(1.0, q / n_ramp)) ** 2
        g.E[:, y_src, 0, 2] += ramp * np.sin(omega * q)  # soft source: transparent to returning waves
        k = (q - (n_total - 2 * n_win)) // n_win
        if 0 <= k < 2:
            acc[k] += g.E[:, y_mon, 0, 2] * np.exp(1j * omega * q)
    a1, a2 = (2 / n_win) * acc[0].mean(), (2 / n_win) * acc[1].mean()
    return a2, abs(a2 / a1 - 1)  # phasor (x-averaged = zeroth order), steady-state change between windows


# ----------------------------------------------------------------------------- RCWA route
def rcwa_t(widths, period=P, n_grid=200, n_terms=197):
    from metatrap.rcwa.solver import solve_eps

    xc = (np.arange(n_grid) + 0.5) / n_grid * period - period / 2
    eps = np.ones((len(widths), n_grid, n_grid))
    for i, w in enumerate(widths):
        eps[i][np.abs(xc) < w / 2 - 1e-9, :] = N_RIDGE**2  # varies along the first grid axis
    r = solve_eps(eps.astype(complex), period, H, LAM, N_SUB**2, 1.0, n_terms=n_terms, formulation="pol_fourier")
    return r.t_x, r.t_y


def ratio_err(a, b):
    """(|amplitude ratio difference|, |phase difference|) between two complex ratios."""
    return float(abs(abs(a) - abs(b))), float(abs(np.angle(a / b)))


def main(brk):
    res = {"thresholds": THR, "case": {"lam": LAM, "pitch": P, "height": H, "n_sub": N_SUB, "n_ridge": N_RIDGE, "widths": WIDTHS},
           "break": brk}
    t0 = time.time()
    # C1: FDTD against the analytic layer.
    # FDTD phasors are cached (the --break run reuses them; the corruptions act on the comparisons).
    out = REPO_ROOT / "results" / "rcwa_fdtd_check"
    out.mkdir(parents=True, exist_ok=True)
    cache_path = out / "fdtd_phasors.json"
    key = json.dumps(res["case"], sort_keys=True)
    cache = json.load(open(cache_path)) if cache_path.exists() else {}
    if cache.get("key") != key:
        cache = {"key": key}
        for d in (0.010, 0.005):
            runs = {"slab": fdtd_t(0, d, "slab"), "air": fdtd_t(0, d, "air")}
            runs.update({str(w): fdtd_t(w, d, "grating") for w in WIDTHS})
            cache[str(d)] = {k: [v[0].real, v[0].imag, v[1]] for k, v in runs.items()}
            print(f"FDTD d = {d} um done ({time.time() - t0:.0f}s)", flush=True)
        json.dump(cache, open(cache_path, "w"), indent=1)
    fd = {}
    for d in (0.010, 0.005):
        ph = {k: complex(v[0], v[1]) for k, v in cache[str(d)].items()}
        fd[d] = {"slab_over_air": ph["slab"] / ph["air"], "steady": max(v[2] for v in cache[str(d)].values()),
                 "ratios": {w: ph[str(w)] / ph[str(W_REF)] for w in WIDTHS}}
    n_an = 2.05 if brk else N_RIDGE
    ana = airy_t(N_SUB, n_an, 1.0, H, LAM) / airy_t(N_SUB, 1.0, 1.0, H, LAM)
    c1 = {}
    for d in fd:
        a, ph = ratio_err(fd[d]["slab_over_air"], ana)
        c1[str(d)] = {"amp_err": float(abs(abs(fd[d]["slab_over_air"] / ana) - 1)), "phi_err": ph}
    c1_pass = all(v["amp_err"] <= THR["C1_amp"] and v["phi_err"] <= THR["C1_phi"] for v in c1.values())
    # C2: FDTD sampling (10 nm against 5 nm).
    c2 = {}
    for w in WIDTHS:
        wc = WIDTHS[(WIDTHS.index(w) + 1) % 3] if brk else w
        a, ph = ratio_err(fd[0.010]["ratios"][w], fd[0.005]["ratios"][wc])
        c2[str(w)] = {"amp_err": a, "phi_err": ph}
    c2_pass = all(v["amp_err"] <= THR["C2_amp"] and v["phi_err"] <= THR["C2_phi"] for v in c2.values())
    # C3: RCWA against FDTD (5 nm). E || ridges: ridges vary along the first grid axis, so the
    # parallel polarisation is the second axis (t_y) if the first axis is x; C4 decides.
    tx, ty = rcwa_t(WIDTHS)
    tx4, ty4 = rcwa_t([0.025], period=0.05)
    n_emt = N_RIDGE + 0.1 if brk else N_RIDGE  # break: wrong ridge index in the analytic limit
    nte = np.sqrt(0.5 * n_emt**2 + 0.5)
    ntm = 1 / np.sqrt(0.5 / n_emt**2 + 0.5)
    a_te, a_tm = airy_t(N_SUB, nte, 1.0, H, LAM), airy_t(N_SUB, ntm, 1.0, H, LAM)
    # C4: which wrapper output is E || ridges.
    # Orientation from the unperturbed effective index (so the C4 break cannot flip it).
    a_te0 = airy_t(N_SUB, np.sqrt(0.5 * N_RIDGE**2 + 0.5), 1.0, H, LAM)
    par_is_y = abs(np.angle(ty4[0] / a_te0)) < abs(np.angle(tx4[0] / a_te0))
    t_par4, t_perp4 = (ty4[0], tx4[0]) if par_is_y else (tx4[0], ty4[0])
    c4 = {"parallel_output": "t_y" if par_is_y else "t_x", "phi_err_TE": float(abs(np.angle(t_par4 / a_te))),
          "phi_err_TM": float(abs(np.angle(t_perp4 / a_tm)))}
    c4_pass = c4["phi_err_TE"] <= THR["C4_phi"] and c4["phi_err_TM"] <= THR["C4_phi"]
    tp = ty if par_is_y else tx
    if brk:
        tp = tx if par_is_y else ty
    i0 = WIDTHS.index(W_REF)
    c3 = {}
    for i, w in enumerate(WIDTHS):
        a, ph = ratio_err(tp[i] / tp[i0], fd[0.005]["ratios"][w])
        c3[str(w)] = {"amp_err": a, "phi_err": ph, "rcwa_ratio": [float(abs(tp[i] / tp[i0])), float(np.angle(tp[i] / tp[i0]))],
                      "fdtd_ratio": [float(abs(fd[0.005]["ratios"][w])), float(np.angle(fd[0.005]["ratios"][w]))]}
    c3_pass = all(v["amp_err"] <= THR["C3_amp"] and v["phi_err"] <= THR["C3_phi"] for v in c3.values())
    res.update({"C1": c1, "C1_pass": c1_pass, "C2": c2, "C2_pass": c2_pass, "C3": c3, "C3_pass": c3_pass,
                "C3_conclusive": bool(c1_pass and c2_pass), "C4": c4, "C4_pass": bool(c4_pass),
                "fdtd_steady_state_change": {str(d): fd[d]["steady"] for d in fd}, "seconds": time.time() - t0})
    verdicts = [c1_pass, c2_pass, c3_pass, bool(c4_pass)]
    ok = (not any(verdicts)) if brk else all(verdicts)
    for name, v in zip(("C1", "C2", "C3", "C4"), verdicts):
        print(f"{name}: {'PASS' if v else 'FAIL'}  {json.dumps(res[name])}")
    print(f"{'break run: every check must fail' if brk else 'normal run: every check must pass'} -> {'OK' if ok else 'NOT OK'}")
    write_json(REPO_ROOT / "results" / "rcwa_fdtd_check" / ("metrics_break.json" if brk else "metrics.json"), res)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main("--break" in sys.argv))
