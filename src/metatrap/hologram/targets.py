"""Target site patterns for the static tweezer array."""

from __future__ import annotations

import numpy as np


def square_grid(nx: int, ny: int, pitch_um: float, center=(0.0, 0.0)) -> np.ndarray:
    """(N, 2) site coordinates in um of an nx x ny square grid centred on ``center``."""
    ix = (np.arange(nx) - (nx - 1) / 2) * pitch_um + center[0]
    iy = (np.arange(ny) - (ny - 1) / 2) * pitch_um + center[1]
    X, Y = np.meshgrid(ix, iy, indexing="ij")
    return np.column_stack([X.ravel(), Y.ravel()])


def triangular_grid(nx: int, ny: int, pitch_um: float) -> np.ndarray:
    rows = []
    for j in range(ny):
        y = (j - (ny - 1) / 2) * pitch_um * np.sqrt(3) / 2
        off = 0.5 * pitch_um * (j % 2)
        x = (np.arange(nx) - (nx - 1) / 2) * pitch_um + off
        rows.append(np.column_stack([x, np.full(nx, y)]))
    return np.vstack(rows)


def make(kind: str, **kw) -> np.ndarray:
    if kind == "square":
        return square_grid(kw["nx"], kw["ny"], kw["pitch_um"], tuple(kw.get("center", (0.0, 0.0))))
    if kind == "triangular":
        return triangular_grid(kw["nx"], kw["ny"], kw["pitch_um"])
    raise ValueError(kind)
