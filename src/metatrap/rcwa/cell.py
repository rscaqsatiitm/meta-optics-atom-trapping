"""Unit-cell geometry for polarisation-independent pillar meta-atoms."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class PillarCell:
    """Square lattice of pillars of height h on a substrate.

    Here, pitch and height are in um; ``shape`` is 'square' or 'circle'; for a square
    pillar ``width`` is the side, for a circle it is the diameter.
    """

    pitch: float
    height: float
    shape: str = "square"

    def mask(self, widths, n_grid: int, subpixel: int = 4) -> np.ndarray:
        """Fill fraction of each grid pixel, shape (len(widths), n_grid, n_grid).

        Sub-pixel averaging (``subpixel`` x ``subpixel`` samples per pixel) gives a
        smooth dependence of the permittivity on width, which removes staircase
        artefacts in the library.
        """
        widths = np.atleast_1d(np.asarray(widths, dtype=float))
        ns = n_grid * subpixel
        x = (np.arange(ns) + 0.5) / ns * self.pitch - self.pitch / 2
        X, Y = np.meshgrid(x, x, indexing="ij")
        out = np.empty((widths.size, n_grid, n_grid))
        for i, w in enumerate(widths):
            if self.shape == "square":
                m = (np.abs(X) <= w / 2) & (np.abs(Y) <= w / 2)
            elif self.shape == "circle":
                m = X**2 + Y**2 <= (w / 2) ** 2
            else:
                raise ValueError(f"unknown shape {self.shape}")
            out[i] = m.reshape(n_grid, subpixel, n_grid, subpixel).mean(axis=(1, 3))
        return out

    def permittivity(self, widths, eps_pillar: complex, eps_clad: complex, n_grid: int) -> np.ndarray:
        """Linear permittivity mixing on the sub-pixel fill fraction."""
        f = self.mask(widths, n_grid)
        return eps_clad + (eps_pillar - eps_clad) * f

    def aspect_ratio(self, widths) -> np.ndarray:
        return self.height / np.asarray(widths, dtype=float)
