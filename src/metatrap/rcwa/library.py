"""Meta-atom library: complex transmission of square pillars over (material, pitch, height, width, wavelength).

Storage: one HDF5 file per run with datasets
    t[material, pitch, height, wavelength, width]  complex128  zeroth-order co-pol transmission
    T, R, T0                                      float64     total T, total R, zeroth-order T
    t_cross                                       float64     |cross-pol| (symmetry check)
and coordinate datasets 'materials', 'pitch_um', 'height_um', 'wavelength_um',
'width_frac' (width / pitch). The sweep is resumable: a boolean 'done' array marks
finished (material, pitch, height, wavelength) batches.
"""

from __future__ import annotations

import time
from dataclasses import dataclass
from pathlib import Path

import h5py
import numpy as np

from .. import materials
from .cell import PillarCell
from .solver import solve_pillars


@dataclass
class LibrarySpec:
    materials: list[str]
    substrate: str
    pitch_um: list[float]
    height_um: list[float]
    wavelength_um: list[float]
    width_frac: list[float]
    shape: str = "square"
    n_terms: int = 197
    n_grid: int = 256
    formulation: str = "jones_direct_fourier"

    @classmethod
    def from_config(cls, cfg: dict) -> "LibrarySpec":
        wf = cfg["width_frac"]
        width_frac = list(np.linspace(wf["min"], wf["max"], wf["n"])) if isinstance(wf, dict) else list(wf)
        h = cfg["height_um"]
        heights = list(np.round(np.arange(h["min"], h["max"] + 1e-9, h["step"]), 6)) if isinstance(h, dict) else list(h)
        return cls(
            materials=list(cfg["materials"]),
            substrate=cfg["substrate"],
            pitch_um=list(cfg["pitch_um"]),
            height_um=heights,
            wavelength_um=list(cfg["wavelength_um"]),
            width_frac=width_frac,
            shape=cfg.get("shape", "square"),
            n_terms=cfg.get("n_terms", 197),
            n_grid=cfg.get("n_grid", 256),
            formulation=cfg.get("formulation", "jones_direct_fourier"),
        )

    @property
    def shape_tuple(self):
        return (len(self.materials), len(self.pitch_um), len(self.height_um), len(self.wavelength_um), len(self.width_frac))


def _create(path: Path, spec: LibrarySpec):
    with h5py.File(path, "w") as f:
        f.create_dataset("materials", data=np.array(spec.materials, dtype="S"))
        for k in ("pitch_um", "height_um", "wavelength_um", "width_frac"):
            f.create_dataset(k, data=np.asarray(getattr(spec, k), float))
        shp = spec.shape_tuple
        f.create_dataset("t", shape=shp, dtype=np.complex128, fillvalue=complex(np.nan, np.nan))
        for k in ("T", "R", "T0", "t_cross"):
            f.create_dataset(k, shape=shp, dtype=np.float64, fillvalue=np.nan)
        f.create_dataset("done", shape=shp[:4], dtype=bool, fillvalue=False)
        f.attrs.update(
            substrate=spec.substrate, shape=spec.shape, n_terms=spec.n_terms, n_grid=spec.n_grid,
            formulation=spec.formulation, incidence="from substrate, exit into air, normal incidence",
        )


def run_sweep(spec: LibrarySpec, path: str | Path, log=print, max_batches: int | None = None) -> Path:
    """Run (or resume) the sweep. ``max_batches`` stops after that many new batches.

    JAX's CPU compilation cache grows over a long sweep until LLVM fails to allocate
    memory (observed after a few hundred batches on 15 GB), so long sweeps should run
    as a loop of fresh processes with a batch limit (scripts/run_library.py --max-batches).
    """
    path = Path(path)
    n_new = 0
    if not path.exists():
        _create(path, spec)
    # Wavelength is the outermost loop so that the first wavelength listed (the
    # design wavelength) is complete for every geometry before the others start.
    for il, lam in enumerate(spec.wavelength_um):
        for im, mat in enumerate(spec.materials):
            for ip, p in enumerate(spec.pitch_um):
                widths = np.asarray(spec.width_frac) * p
                for ih, h in enumerate(spec.height_um):
                    cell = PillarCell(p, h, spec.shape)
                    with h5py.File(path, "r") as f:
                        if f["done"][im, ip, ih, il]:
                            continue
                    t0 = time.time()
                    n_pil = materials.get(mat).nk(lam)
                    n_sub = materials.get(spec.substrate).nk(lam)
                    r = solve_pillars(cell, widths, lam, n_pil**2, n_sub**2, n_terms=spec.n_terms, n_grid=spec.n_grid, formulation=spec.formulation)
                    with h5py.File(path, "a") as f:
                        f["t"][im, ip, ih, il, :] = r.t_x
                        f["T"][im, ip, ih, il, :] = r.T
                        f["R"][im, ip, ih, il, :] = r.R
                        f["T0"][im, ip, ih, il, :] = r.T0
                        f["t_cross"][im, ip, ih, il, :] = r.t_cross
                        f["done"][im, ip, ih, il] = True
                    log(f"{mat} p={p} h={h} lam={lam}: {time.time()-t0:.1f}s, T0 range {r.T0.min():.3f}-{r.T0.max():.3f}")
                    n_new += 1
                    if max_batches is not None and n_new >= max_batches:
                        return path
    return path


@dataclass
class Library:
    materials: list[str]
    pitch_um: np.ndarray
    height_um: np.ndarray
    wavelength_um: np.ndarray
    width_frac: np.ndarray
    t: np.ndarray
    T0: np.ndarray
    T: np.ndarray
    R: np.ndarray
    attrs: dict

    @classmethod
    def load(cls, path: str | Path) -> "Library":
        with h5py.File(path, "r") as f:
            return cls(
                materials=[m.decode() for m in f["materials"][:]],
                pitch_um=f["pitch_um"][:], height_um=f["height_um"][:],
                wavelength_um=f["wavelength_um"][:], width_frac=f["width_frac"][:],
                t=f["t"][:], T0=f["T0"][:], T=f["T"][:], R=f["R"][:], attrs=dict(f.attrs),
            )

    def index(self, material: str, pitch: float, height: float, lam: float):
        im = self.materials.index(material)
        ip = int(np.argmin(np.abs(self.pitch_um - pitch)))
        ih = int(np.argmin(np.abs(self.height_um - height)))
        il = int(np.argmin(np.abs(self.wavelength_um - lam)))
        return im, ip, ih, il

    def curve(self, material: str, pitch: float, height: float, lam: float):
        """(widths_um, complex t normalised to unit max |t| reference-free, T0) for one geometry."""
        im, ip, ih, il = self.index(material, pitch, height, lam)
        return self.width_frac * self.pitch_um[ip], self.t[im, ip, ih, il], self.T0[im, ip, ih, il]


def is_complete(path: str | Path) -> bool:
    with h5py.File(path, "r") as f:
        return bool(np.all(f["done"][:]))
