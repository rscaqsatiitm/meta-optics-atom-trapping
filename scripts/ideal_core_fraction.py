"""Core (first-dark-ring) power fraction of the aberration-free focus of a metasurface pupil.

Usage: python scripts/ideal_core_fraction.py

The power budget of scripts/fullaperture.py converts the simulated core efficiency into a
utilisation as core_efficiency / core_fraction_ideal, where core_fraction_ideal is the fraction
of the power that an ideal lens with the same pupil puts inside r_core = 0.61 lambda / NA, in the
convention of the full-pixel metric (scalar |u|^2 in the focal plane, total power = incident
power). Here it is computed with the pupil model: the angular spectrum of the ideal focus is
B(s) = a(s) sqrt(cos theta), with a(s) the pupil amplitude of hologram.wgs.pupil_amplitude (which
carries the power-flux factor 1 / sqrt(cos theta)), and the focal field is its FFT. Checked
against the full-pixel bench references (results/fullpixel_bench5mm*: 0.7926 hard top hat, 0.8171
soft edge, NA 0.5). Output: results/ideal_core_fraction.json.
"""

import numpy as np

from metatrap.hologram.wgs import PupilGrid, pupil_amplitude
from metatrap.propagation.layers import LayerStack
from metatrap.runio import REPO_ROOT, write_json


def core_fraction(stack, na, lam, ill, n=4096, m=220):
    ds = na / m
    g = PupilGrid(n=n, ds=ds, lam_um=lam)
    a = pupil_amplitude(g, stack, na, ill)
    s1 = g.s_coords()
    SX, SY = np.meshgrid(s1, s1, indexing="ij")
    B = a * np.sqrt(np.sqrt(np.clip(1 - SX**2 - SY**2, 0, None)))
    I = np.abs(np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(B)))) ** 2
    x = (np.arange(n) - n // 2) * lam / (n * ds)
    R = np.hypot(x[:, None], x[None, :])
    return float(I[R <= 0.61 * lam / na].sum() / I.sum())


def main():
    bench = LayerStack.from_list([[1000, "air"], [1000, "fused_silica"], [2693, "air"]])
    cell = LayerStack.from_list([[2500, "air"], [4000, "fused_silica"], [9000, "air"]])
    res = {
        "bench_na0.5_top_hat": {"pupil_model": core_fraction(bench, 0.5, 1.075, {"kind": "top_hat"}), "full_pixel": 0.7926},
        "bench_na0.5_soft_edge_250um": {"pupil_model": core_fraction(bench, 0.5, 1.075, {"kind": "top_hat", "edge_um": 250.0}), "full_pixel": 0.8171},
        "cell_na0.55_top_hat": {"pupil_model": core_fraction(cell, 0.55, 1.075, {"kind": "top_hat"})},
    }
    print(res)
    write_json(REPO_ROOT / "results" / "ideal_core_fraction.json", res)


if __name__ == "__main__":
    main()
