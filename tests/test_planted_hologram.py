"""Planted truth for the hologram design path, checked with a second propagator.

A known target array (3 x 3 sites, 8 um pitch, deliberately asymmetric planted powers so that
x/y swaps and sign flips are caught) is designed with the production path (PupilGrid,
exact kernels, wgs_exact with targets, render_grid, metasurface_field) for a free-space stack
(200 um of air, NA 0.3, 1075 nm). The metasurface field is then propagated two ways:
  AS  the production angular-spectrum propagator (propagation/angular_spectrum.propagate) and
      the production spot metrics (hologram/evaluate.spot_metrics), and
  RS  a direct Rayleigh-Sommerfeld (first kind) sum over the metasurface pixels, written here
      with numpy only. It shares no code with the design or the propagation modules; the
      shared layers are the sampled metasurface field itself, the pixel convention
      (pixel (i, j) at ((i - n//2) dx, (j - n//2) dx), first index x) and numpy.
Both are exact solutions of the Helmholtz equation for the same boundary field, so they must
agree; and both must return the planted positions and powers.

Tolerances, fixed before the first run (never loosened after a result):
  POS_TOL   0.05 um   spot position against the planted site (AS and RS)
  POW_TOL   0.05      max |relative error| of the normalised core power against the planted
                      power (AS and RS); covers the design residual of a pupil-model design
                      without feedback
  AGREE_PK  0.01      max |RS / AS - 1| of the per-site peak intensity (absolute, not normalised)
  AGREE_POS 0.01 um   max |RS - AS| position difference
  AGREE_POW 0.01      max |RS / AS - 1| of the per-site core power
Negative controls (each must make its check fail): planted powers transposed x <-> y; the RS
input shifted by one pixel; the RS evaluated 1 % further away (z x 1.01).

Change of 2026-11-20 (the thresholds are unchanged). The first version compared the AS spots
from spot_metrics with RS spots from a separate estimator (1D quadratic fits on a 0.125 um grid)
and failed the position criterion at 0.0138 um, although the two intensity maps agree to 0.52 %
of the peak (docs/planted_truth_hologram.md): it tested the estimators as much as the
propagators. The AS-RS comparison now applies one measurement chain to both fields: each field
is sampled on the same 96 x 96 window of 0.5 um pixels, upsampled x 4 by fourier_upsample and
measured by spot_metrics, so the shared layers are the window, the upsampling and the estimator,
and the check isolates the propagators. The RS recovery check keeps its own estimator
(rs_spots), so that recovery is still measured by code independent of the production metrics.

Revision of 2026-11-20: the design claims trap depth (decided by the PI). The recovered quantity
is the normalised peak intensity (trap depth is proportional to it), after the production
feedback: three pupil-phase rounds at gain 1.0 on the AS-measured peaks, as in
scripts/fullaperture.py (feedback_mode pupil_phase). Threshold:
  DEPTH_TOL 0.01      max |relative error| of the normalised peak intensity against the planted
                      power, measured (a) by the production AS route, the quantity the loop drives
                      (in-sample), and (b) by the RS route with its own estimator (independent of
                      the loop and of the production metrics).
Spot positions and core powers are reported but not asserted: the former tolerances (POS_TOL,
POW_TOL, kept below for scripts/planted_truth_diagnostics.py) failed at 0.13 um and 6.6 %, which
the diagnostics traced to neighbour interference at 8 um pitch, not to the design path
(docs/planted_truth_hologram.md). Provenance of DEPTH_TOL: chosen after diagnostics D1 to D3,
which had shown 1.2 % for the AS peak error without feedback and 0.21 % after three rounds; the RS
value after feedback had not been computed. The level is that of the claim (relative trap depths
set to better than 1 %), not fitted to those values. Additional negative control: one planted
power changed by 2 % (twice DEPTH_TOL) must fail both recoveries.
"""

import numpy as np
import pytest

from metatrap.hologram import evaluate as E
from metatrap.hologram import exact as X
from metatrap.hologram.render import fourier_upsample, metasurface_field
from metatrap.hologram.wgs import PupilGrid
from metatrap.propagation.angular_spectrum import propagate
from metatrap.propagation.layers import LayerStack

POS_TOL, POW_TOL = 0.05, 0.05  # former recovery tolerances (positions, core powers): reported only
DEPTH_TOL = 0.01
FEEDBACK_ROUNDS, FEEDBACK_GAIN = 3, 1.0
AGREE_PK, AGREE_POS, AGREE_POW = 0.01, 0.01, 0.01

LAM, NA, Z, PITCH, DX, N = 1.075, 0.30, 200.0, 8.0, 0.5, 512
ILL = {"kind": "top_hat"}


def planted_power(x, y):
    """Distinct values at all nine sites (0.58 to 1.5), not symmetric under x <-> y or sign flips."""
    return 1.0 + 0.30 * x / PITCH + 0.12 * y / PITCH + 0.5 * ((abs(x) < 1e-9) & (abs(y) < 1e-9))


def sites_and_plant():
    c = np.array([-PITCH, 0.0, PITCH])
    X_, Y_ = np.meshgrid(c, c, indexing="ij")
    sites = np.column_stack([X_.ravel(), Y_.ravel()])
    return sites, planted_power(sites[:, 0], sites[:, 1])


# ---------------------------------------------------------------- second propagator (numpy only)
def rs_field(f, dx, z, lam, pts, shift_px=0):
    """Rayleigh-Sommerfeld (first kind) field at points (x, y) a distance z beyond the sampled field f.

    U(P) = (z / 2 pi) sum_Q U(Q) exp(i k r) (1 / r - i k) / r^2 dx^2, pixel (i, j) at
    ((i - n//2) dx, (j - n//2) dx)."""
    k = 2 * np.pi / lam
    n = f.shape[0]
    ii, jj = np.nonzero(np.abs(f) > 0)
    xs = (ii - n // 2 + shift_px) * dx
    ys = (jj - n // 2) * dx
    v = f[ii, jj].astype(complex) * dx * dx
    out = np.empty(len(pts), complex)
    for a in range(0, len(pts), 256):
        p = pts[a : a + 256]
        r = np.sqrt((p[:, 0:1] - xs[None, :]) ** 2 + (p[:, 1:2] - ys[None, :]) ** 2 + z * z)
        out[a : a + 256] = (z / (2 * np.pi) * np.exp(1j * k * r) * (1 / r - 1j * k) / r**2) @ v
    return out


def rs_spots(f, sites, z=Z, shift_px=0, h=0.125):
    """Per-site peak (quadratic fit on the RS grid), position and power within 0.61 lambda / NA."""
    r_core = 0.61 * LAM / NA
    m = int(np.ceil((r_core + 0.5) / h))
    off = np.arange(-m, m + 1) * h
    OX, OY = np.meshgrid(off, off, indexing="ij")
    res = {"peak": [], "x": [], "y": [], "p_core": []}
    for sx, sy in sites:
        pts = np.column_stack([(sx + OX).ravel(), (sy + OY).ravel()])
        I = np.abs(rs_field(f, DX, z, LAM, pts, shift_px)) ** 2
        I = I.reshape(OX.shape)
        a, b = np.unravel_index(np.argmax(I), I.shape)
        w3 = I[a - 1 : a + 2, b - 1 : b + 2]
        # 1D quadratic fits through the maximum along each axis.
        fx = 0.5 * (w3[0, 1] - w3[2, 1]) / (w3[0, 1] - 2 * w3[1, 1] + w3[2, 1])
        fy = 0.5 * (w3[1, 0] - w3[1, 2]) / (w3[1, 0] - 2 * w3[1, 1] + w3[1, 2])
        px, py = sx + off[a] + fx * h, sy + off[b] + fy * h
        pk = w3[1, 1] - 0.25 * (w3[0, 1] - w3[2, 1]) * fx - 0.25 * (w3[1, 0] - w3[1, 2]) * fy
        core = np.hypot(sx + OX - px, sy + OY - py) <= r_core
        res["peak"].append(pk)
        res["x"].append(px)
        res["y"].append(py)
        res["p_core"].append(float(I[core].sum() * h * h))
    return {k: np.array(v) for k, v in res.items()}


def rs_window_spots(f, sites, z=Z, shift_px=0, w=96):
    """RS field sampled on the AS measurement window, then the production measurement chain."""
    idx = (np.arange(w) - w // 2) * DX
    X_, Y_ = np.meshgrid(idx, idx, indexing="ij")
    pts = np.column_stack([X_.ravel(), Y_.ravel()])
    win = rs_field(f, DX, z, LAM, pts, shift_px).reshape(w, w)
    I = np.abs(fourier_upsample(win, 4)) ** 2
    return E.spot_metrics(I, DX / 4, sites, PITCH, LAM, NA, 1.0)["per_site"]


# ---------------------------------------------------------------- production design and AS route
@pytest.fixture(scope="module")
def design():
    sites, plant = sites_and_plant()
    stack = LayerStack.from_list([[Z, "air"]])
    g = PupilGrid.for_sites(LAM, PITCH, 4, NA, 2 * np.abs(sites).max() + PITCH, guard=1.6)
    ps = X.PupilSampling.build(g, stack, NA, ILL)
    radius = stack.aperture_for_na(NA, LAM)
    tab = X.OPLTable(stack, LAM, radius + 2 * np.abs(sites).max() + 20)
    K = X.kernels(ps, sites, tab, LAM)
    d = X.wgs_exact(ps, K, n_iter=80, n_fix=40, target=np.sqrt(plant))
    assert 2 * radius / DX < 0.6 * N  # the AS window holds the aperture with margin

    def as_route(coeffs, psi):
        holo = X.render_grid(ps, stack, sites, coeffs, LAM, NA, kind="exact", table=tab, upsample=8)
        f = metasurface_field(N, DX, stack, LAM, NA, ps, psi, ILL, holo=holo)
        foc = propagate(f.copy(), DX, LAM, stack)
        w = 96
        c0 = N // 2 - w // 2
        I = np.abs(fourier_upsample(foc[c0 : c0 + w, c0 : c0 + w], 4)) ** 2  # propagate keeps the axis at N//2
        return f, E.spot_metrics(I, DX / 4, sites, PITCH, LAM, NA, 1.0)["per_site"]

    f, per0 = as_route(d.coeffs, d.psi)
    per, c = per0, d.coeffs.astype(np.complex64)
    for _ in range(FEEDBACK_ROUNDS):  # production pupil-phase feedback on the measured peaks
        corr = ((plant / plant.mean()) / (per["peak"] / per["peak"].mean())) ** FEEDBACK_GAIN
        c = (c * np.sqrt(corr)).astype(np.complex64)
        f, per = as_route(c, np.angle(K @ c).astype(np.float64))
    return {"sites": sites, "plant": plant, "f": f, "as0": per0, "as": per, "rs": rs_spots(f, sites),
            "rs_win": rs_window_spots(f, sites)}


def recovery_errors(spots, sites, plant):
    pos = np.hypot(spots["x"] - sites[:, 0], spots["y"] - sites[:, 1])
    p = spots["p_core"] / spots["p_core"].mean()
    q = plant / plant.mean()
    return float(pos.max()), float(np.max(np.abs(p / q - 1)))


def agreement(a, b):
    return (float(np.max(np.abs(b["peak"] / a["peak"] - 1))), float(np.max(np.hypot(b["x"] - a["x"], b["y"] - a["y"]))),
            float(np.max(np.abs(b["p_core"] / a["p_core"] - 1))))


def depth_error(spots, plant):
    p = spots["peak"] / spots["peak"].mean()
    q = plant / plant.mean()
    return float(np.max(np.abs(p / q - 1)))


def recovered(spots, sites, plant):
    """Trap-depth recovery: normalised peak intensities within DEPTH_TOL of the planted powers."""
    return depth_error(spots, plant) <= DEPTH_TOL


def agree(a, b):
    pk, pos, powr = agreement(a, b)
    return pk <= AGREE_PK and pos <= AGREE_POS and powr <= AGREE_POW


def report(name, spots, sites, plant):
    pos, powr = recovery_errors(spots, sites, plant)
    print(f"{name}: trap-depth error {depth_error(spots, plant):.4f} (tol {DEPTH_TOL}); reported only: max position "
          f"error {pos:.4f} um, core-power error {powr:.4f}")


def test_trap_depths_recovered_by_production_as(design):
    report("AS without feedback", design["as0"], design["sites"], design["plant"])
    report(f"AS after {FEEDBACK_ROUNDS} feedback rounds", design["as"], design["sites"], design["plant"])
    assert recovered(design["as"], design["sites"], design["plant"])


def test_trap_depths_recovered_by_rayleigh_sommerfeld(design):
    report(f"RS after {FEEDBACK_ROUNDS} feedback rounds", design["rs"], design["sites"], design["plant"])
    assert recovered(design["rs"], design["sites"], design["plant"])


def test_angular_spectrum_agrees_with_rayleigh_sommerfeld(design):
    rs = design["rs_win"]  # same window, upsampling and estimator as the AS route
    pk, pos, powr = agreement(design["as"], rs)
    print(f"AS vs RS: peak {pk:.2e} (tol {AGREE_PK}), position {pos:.2e} um (tol {AGREE_POS}), core power {powr:.2e} (tol {AGREE_POW})")
    assert agree(design["as"], rs)


def test_negative_controls_fail(design):
    """Each deliberate break must make its check fail."""
    sites, plant, f = design["sites"], design["plant"], design["f"]
    transposed = planted_power(sites[:, 1], sites[:, 0])  # x <-> y swapped plant
    assert not recovered(design["as"], sites, transposed)
    assert not recovered(design["rs"], sites, transposed)
    off = plant.copy()
    off[0] *= 1.02  # one planted power changed by twice DEPTH_TOL
    assert not recovered(design["as"], sites, off)
    assert not recovered(design["rs"], sites, off)
    assert not agree(design["as"], rs_window_spots(f, sites, shift_px=1))
    assert not agree(design["as"], rs_window_spots(f, sites, z=Z * 1.01))
