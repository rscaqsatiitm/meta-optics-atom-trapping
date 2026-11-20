"""Spot metrics on planted spots with known sub-pixel centre and power.

Thresholds, fixed before the first run:
  disc centroid       within 0.05 pixel of the sub-pixel peak (ox along the first axis, oy along
                      the second), for ox = 0.3, oy = -0.2, radius 15 pixels
  Gaussian spot       position within 0.005 um (0.05 pixel), peak and core power within 0.5 %
                      of the analytic values (w = 0.8 um on a 0.1 um grid, core radius 1.2 um)
The disc test is the one the swapped-axis bug in the disc centring must fail (found in the
planted-truth work, docs/planted_truth_hologram.md, finding 3).
"""

import numpy as np

from metatrap.hologram import evaluate as E


def test_disc_centred_on_subpixel_peak():
    rr, ox, oy = 20, 0.3, -0.2
    d = E.disc_distances(rr, ox, oy, 1.0)
    m = d <= 15.0
    i, j = np.nonzero(m)
    ci, cj = i.mean() - rr, j.mean() - rr
    assert abs(ci - ox) <= 0.05 and abs(cj - oy) <= 0.05, (ci, cj)


def test_gaussian_spot_recovered():
    dx, n, w, P = 0.1, 129, 0.8, 2.0
    x0, y0 = 0.03, -0.02  # sub-pixel offsets 0.3 and -0.2 pixel, along the first and second axis
    na = 0.61 / 1.2       # core radius 0.61 lambda / NA = 1.2 um with lambda = 1
    c = n // 2
    ax = (np.arange(n) - c) * dx
    X, Y = np.meshgrid(ax, ax, indexing="ij")
    I = 2 * P / (np.pi * w**2) * np.exp(-2 * ((X - x0) ** 2 + (Y - y0) ** 2) / w**2)
    m = E.spot_metrics(I, dx, np.array([[0.0, 0.0]]), 6.0, 1.0, na, P)["per_site"]
    assert abs(m["x"][0] - x0) <= 0.005 and abs(m["y"][0] - y0) <= 0.005, (m["x"][0], m["y"][0])
    assert abs(m["peak"][0] / (2 * P / (np.pi * w**2)) - 1) <= 0.005, m["peak"][0]
    assert abs(m["p_core"][0] / (P * (1 - np.exp(-2 * 1.2**2 / w**2))) - 1) <= 0.005, m["p_core"][0]
