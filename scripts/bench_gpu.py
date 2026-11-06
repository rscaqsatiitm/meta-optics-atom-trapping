"""Benchmark and validate the GPU kernels on the fullaperture_100k geometry.

Usage: python scripts/bench_gpu.py [n_sites_per_side]      (default 316, the 100k design)

Builds the wedge pupil and fundamental sites as scripts/fullaperture.py does, then
(i) compares LazyKGGPU with the CPU LazyKG on 2,000 pupil rows and 500 sites (relative error),
(ii) times one full K @ x and one conj(v) @ K on the GPU and reports kernel evaluations per
second, and (iii) extrapolates the WGS time of the 10^6-site design (configs/fullaperture_1M.yaml:
5.2e12 evaluations per product, 120 iterations of two products).
"""

import sys
import time
import types

import numpy as np

from metatrap.hologram import lazy as L
from metatrap.hologram import symmetric as S
from metatrap.hologram import targets
from metatrap.hologram.exact import OPLTable
from metatrap.hologram.lazy_gpu import LazyKGGPU
from metatrap.hologram.wgs import PupilGrid
from metatrap.propagation.layers import LayerStack

nx = int(sys.argv[1]) if len(sys.argv) > 1 else 316
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
