"""Refractive indices of the films and substrates used in the meta-optics.

Sources (retrieved and read 2026-10-31 from the refractiveindex.info database,
https://github.com/polyanskiy/refractiveindex.info-database, CC0):

* ``sin_lpcvd``: K. Luke et al., Opt. Lett. 40, 4823 (2015), DOI 10.1364/OL.40.004823.
  Stoichiometric LPCVD Si3N4, formula 1, valid 0.310-5.504 um.
* ``fused_silica``: I. H. Malitson, J. Opt. Soc. Am. 55, 1205 (1965),
  DOI 10.1364/JOSA.55.001205. Valid 0.21-3.71 um, 20 C.
* ``sapphire_o``: I. H. Malitson and M. J. Dodge, J. Opt. Soc. Am. 62, 1405 (1972)
  (abstract); M. J. Dodge, CRC Handbook of Laser Science and Technology IV (1986).
  Ordinary ray, 20 C, valid 0.20-5.0 um.
* ``sin_si_rich``: constant n = 2.3, the PECVD silicon-rich nitride of Wang et al.,
  Chin. Phys. Lett. 43, 010606 (2026) [proposal ref. 14], as quoted in the proposal
  (Section 4, B1). Not retrieved this session: [VALUE?: n = 2.3 at 852 nm, verify
  against ref. 14 and measure dispersion at 1075 nm].
* ``sin_si_rich_est``: literature-based estimate for a Si-rich PECVD film with n = 2.30
  at 852 nm (the ref. 14 value quoted in the proposal). Dispersion: single-oscillator
  (Wemple-DiDomenico) model n^2 - 1 = E0 Ed / (E0^2 - E^2), with E0 extrapolated
  linearly in n(852 nm) from fits to five measured SiN films (Luke et al. 2015; Vogt,
  PhD thesis 2015, films #1-#3; Beliaev et al., Thin Solid Films 763, 139568 (2022),
  DOI 10.1016/j.tsf.2022.139568; all via the refractiveindex.info database), giving
  E0 = 5.834 eV, Ed = 23.469 eV, n(1075 nm) = 2.277. Absorption: k = 6.7e-5, from the
  photothermal absorption coefficient 7.94 +/- 0.50 cm^-1 at 1064 nm of Si-rich
  low-stress SiNx membranes (Kumar et al., ACS Photonics 12, 6381 (2025),
  DOI 10.1021/acsphotonics.5c01886), consistent with 6.9 +/- 0.7 cm^-1 of Steinlechner et
  al., Phys. Rev. D 96, 022007 (2017). Those films are LPCVD; PECVD films contain more
  hydrogen, so k for our films is uncertain by a factor of a few (Provisional; band
  3e-5 to 3e-4 used in the sensitivity runs). k is taken as wavelength independent near
  1.07 um.
* LPCVD stoichiometric SiN absorption: k = 1.8e-7 at 1064 nm (Kumar et al. 2025, 1.5-2.9e-2
  cm^-1; Land et al. value quoted there).
* ``pyrex``: constant n = 1.47 [VALUE?: approx. 1.47 for borosilicate 3.3 near
  850-1075 nm, verify]. Used only to reproduce the glass cell of proposal ref. 4.

Measured films are loaded from CSV files (columns: wavelength_um, n, k) with
:func:`from_csv`. Extinction coefficient k is zero for all literature models; film
absorption at 1075 nm is a free parameter of the thermal model (hypothesis H2a).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

import numpy as np


def _sellmeier_formula1(coeffs: tuple[float, ...]) -> Callable[[np.ndarray], np.ndarray]:
    """refractiveindex.info 'formula 1': n^2 - 1 = C0 + sum B_i l^2 / (l^2 - C_i^2)."""
    c0 = coeffs[0]
    pairs = list(zip(coeffs[1::2], coeffs[2::2]))

    def n_of(lam_um: np.ndarray) -> np.ndarray:
        l2 = np.asarray(lam_um, dtype=float) ** 2
        n2 = 1.0 + c0
        for b, c in pairs:
            n2 = n2 + b * l2 / (l2 - c**2)
        return np.sqrt(n2)

    return n_of


@dataclass(frozen=True)
class Material:
    """A dispersive, isotropic material."""

    name: str
    n_func: Callable[[np.ndarray], np.ndarray]
    k_func: Callable[[np.ndarray], np.ndarray] = field(default=lambda lam: np.zeros_like(np.asarray(lam, float)))
    valid_um: tuple[float, float] = (0.2, 5.0)
    source: str = ""

    def n(self, lam_um) -> np.ndarray:
        lam = np.asarray(lam_um, dtype=float)
        if np.any(lam < self.valid_um[0]) or np.any(lam > self.valid_um[1]):
            raise ValueError(f"{self.name}: wavelength outside validity range {self.valid_um} um")
        return self.n_func(lam)

    def k(self, lam_um) -> np.ndarray:
        return self.k_func(np.asarray(lam_um, dtype=float))

    def nk(self, lam_um) -> np.ndarray:
        """Complex index n + i k."""
        return self.n(lam_um) + 1j * self.k(lam_um)

    def eps(self, lam_um) -> np.ndarray:
        """Complex relative permittivity (n + i k)^2."""
        return self.nk(lam_um) ** 2


def _const(value: float) -> Callable[[np.ndarray], np.ndarray]:
    return lambda lam: np.full_like(np.asarray(lam, dtype=float), value)


SIN_LPCVD = Material(
    name="sin_lpcvd",
    n_func=_sellmeier_formula1((0.0, 3.0249, 0.1353406, 40314.0, 1239.842)),
    valid_um=(0.310, 5.504),
    source="Luke et al., Opt. Lett. 40, 4823 (2015), DOI 10.1364/OL.40.004823",
)

# Literature absorption estimates (see module docstring). Kept separate from the
# lossless entries so that earlier results remain reproducible.
K_SIN_LPCVD_1064 = 1.8e-7
K_SIN_SI_RICH_1064 = 6.7e-5


def _wdd(e0: float, ed: float) -> Callable[[np.ndarray], np.ndarray]:
    """Single-oscillator (Wemple-DiDomenico) dispersion; energies in eV, lambda in um."""

    def n_of(lam_um):
        e = 1.23984 / np.asarray(lam_um, dtype=float)
        return np.sqrt(1.0 + e0 * ed / (e0**2 - e**2))

    return n_of


SIN_LPCVD_K = Material(
    name="sin_lpcvd_k",
    n_func=SIN_LPCVD.n_func,
    k_func=_const(K_SIN_LPCVD_1064),
    valid_um=(0.310, 5.504),
    source=SIN_LPCVD.source + "; k = 1.8e-7 (Kumar et al., ACS Photonics 2025)",
)


def sin_si_rich_estimate(n852: float = 2.30, k: float = K_SIN_SI_RICH_1064) -> Material:
    """Si-rich SiN with n(852 nm) = n852; E0 from the linear E0(n852) trend of five films."""
    e0 = -8.38462007 * n852 + 25.1182883
    e852 = 1.23984 / 0.852
    ed = (n852**2 - 1.0) * (e0**2 - e852**2) / e0
    return Material(
        name=f"sin_si_rich_est_n{n852:.2f}_k{k:.1e}",
        n_func=_wdd(e0, ed),
        k_func=_const(k),
        valid_um=(0.5, 1.7),
        source=f"literature-based estimate (E0 = {e0:.3f} eV, Ed = {ed:.3f} eV, k = {k:.2e}); see module docstring",
    )


SIN_SI_RICH_EST = sin_si_rich_estimate()
SIN_SI_RICH_EST = Material(name="sin_si_rich_est", n_func=SIN_SI_RICH_EST.n_func, k_func=SIN_SI_RICH_EST.k_func,
                           valid_um=SIN_SI_RICH_EST.valid_um, source=SIN_SI_RICH_EST.source)

FUSED_SILICA = Material(
    name="fused_silica",
    n_func=_sellmeier_formula1((0.0, 0.6961663, 0.0684043, 0.4079426, 0.1162414, 0.8974794, 9.896161)),
    valid_um=(0.21, 6.7),
    source="Malitson, J. Opt. Soc. Am. 55, 1205 (1965), DOI 10.1364/JOSA.55.001205",
)

SAPPHIRE_O = Material(
    name="sapphire_o",
    n_func=_sellmeier_formula1((0.0, 1.4313493, 0.0726631, 0.65054713, 0.1193242, 5.3414021, 18.028251)),
    valid_um=(0.20, 5.0),
    source="Malitson and Dodge, JOSA 62, 1405 (1972); Dodge, CRC Handbook IV (1986); ordinary ray",
)

SIN_SI_RICH = Material(
    name="sin_si_rich",
    n_func=_const(2.3),
    valid_um=(0.4, 2.0),
    source="[VALUE?: n = 2.3, PECVD Si-rich SiN as quoted from proposal ref. 14; verify, dispersion neglected]",
)

PYREX = Material(
    name="pyrex",
    n_func=_const(1.47),
    valid_um=(0.4, 2.0),
    source="[VALUE?: approx. 1.47 for borosilicate 3.3, verify]",
)

GAN_O = Material(
    name="gan_o",
    n_func=_sellmeier_formula1((2.60, 1.75, 0.256, 4.1, 17.86)),
    valid_um=(0.35, 10.0),
    source="Barker and Ilegems, Phys. Rev. B 7, 743 (1973), DOI 10.1103/PhysRevB.7.743; wurtzite GaN, ordinary ray "
    "(n_e lower by ~0.02 to 0.04 at 0.85 to 1.06 um; birefringence neglected); via refractiveindex.info",
)

AIR = Material(name="air", n_func=_const(1.0), valid_um=(0.1, 100.0), source="n = 1 (vacuum/air, 3e-4 neglected)")

REGISTRY: dict[str, Material] = {
    m.name: m for m in (SIN_LPCVD, SIN_LPCVD_K, SIN_SI_RICH_EST, FUSED_SILICA, SAPPHIRE_O, SIN_SI_RICH, PYREX, GAN_O, AIR)
}

_DATA = Path(__file__).resolve().parents[2] / "data" / "materials"
# Tabulated literature films (CSV with provenance in the header line).
TABULATED = {
    "tio2_ald": _DATA / "tio2_ald_siefke2016.csv",  # Siefke et al., Adv. Opt. Mater. 4, 1780 (2016)
    "si_crystalline": _DATA / "si_crystalline_schinke2015.csv",  # Schinke et al., AIP Adv. 5, 067168 (2015)
}


def get(name: str) -> Material:
    """Look up a material by name, or load ``csv:<path>`` as a measured film."""
    if name.startswith("csv:"):
        return from_csv(name[4:])
    if name.startswith("sin_si_rich_est:"):
        # e.g. "sin_si_rich_est:2.25:3e-4" -> n(852) = 2.25, k = 3e-4
        parts = name.split(":")
        return sin_si_rich_estimate(float(parts[1]), float(parts[2]) if len(parts) > 2 else K_SIN_SI_RICH_1064)
    if name in TABULATED:
        return from_csv(TABULATED[name], name)
    try:
        return REGISTRY[name]
    except KeyError as exc:
        raise KeyError(f"unknown material '{name}', known: {sorted(REGISTRY)}") from exc


def from_csv(path: str | Path, name: str | None = None) -> Material:
    """Load measured (wavelength_um, n, k) data, linearly interpolated.

    The file needs a header line; columns are wavelength in um, n, k.
    """
    path = Path(path)
    data = np.loadtxt(path, delimiter=",", skiprows=1)
    lam, n, k = data[:, 0], data[:, 1], data[:, 2]
    order = np.argsort(lam)
    lam, n, k = lam[order], n[order], k[order]
    return Material(
        name=name or path.stem,
        n_func=lambda x: np.interp(x, lam, n),
        k_func=lambda x: np.interp(x, lam, k),
        valid_um=(float(lam[0]), float(lam[-1])),
        source=f"measured, {path}",
    )
