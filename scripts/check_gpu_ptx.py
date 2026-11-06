"""Compile the GPU kernels of hologram/lazy_gpu.py to PTX and check that the per-pair work is
single precision (consumer GPUs run float64 at 1/32 to 1/64 of the float32 rate).

Usage: python scripts/check_gpu_ptx.py        (needs numba-cuda with the CUDA toolkit; no GPU)

Prints the float64 instructions in each kernel; only the final accumulator update (one
add.f64 per output element and launch) is expected. Exit status 1 if float64 sin/cos or
float64 FMAs appear.
"""

import re
import sys
from collections import Counter

from numba import cuda, types

from metatrap.hologram import lazy_gpu as G

f32, f32a, f32a2, f64a = types.float32, types.float32[:], types.float32[:, :], types.float64[:]
SIG = (f32a, f32a, f32a, f32a2, f32a2, f32a, f32a, f32a, f32, f32, types.int64, types.int64, f64a, f64a)

bad = False
for name in ("_matvec_kernel", "_rmatvec_kernel"):
    ptx, _ = cuda.compile_ptx(getattr(G, name).py_func, SIG, fastmath=True, cc=(8, 6))
    f64 = Counter(re.findall(r"\b([a-z]+(?:\.[a-z0-9]+)*\.f64)\b", ptx))
    trig32 = len(re.findall(r"(?:sin|cos)\.approx[.a-z]*f32", ptx))
    slow = "trig_reduction_slowpathd" in ptx or any(k.startswith("fma") for k in f64)
    print(f"{name}: float64 ops {dict(f64)}; sin/cos.approx.f32: {trig32}; {'FAIL' if slow else 'ok'}")
    bad |= slow
sys.exit(1 if bad else 0)
