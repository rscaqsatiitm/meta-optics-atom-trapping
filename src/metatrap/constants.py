"""Physical constants and atomic line data.

Every value carries its provenance. Atomic data are from D. A. Steck, "Cesium D
Line Data" and "Rubidium 87 D Line Data", revision 2.3.4 (8 August 2025),
http://steck.us/alkalidata, retrieved and read on 2026-10-31.

Steck's reduced dipole matrix elements use the convention
    Gamma = omega^3 / (3 pi eps0 hbar c^3) * (2J+1)/(2J'+1) * |<J||er||J'>|^2,
so the oscillator-sum polarisability of a J = 1/2 ground state reads
    alpha(omega) = (2/3) * sum_k |<J||er||J'_k>|^2 * omega_k / (omega_k^2 - omega^2)
in atomic units (see metatrap.atoms.alkali).
"""

from __future__ import annotations

from dataclasses import dataclass

# CODATA 2018 exact or recommended values (SI).
C = 299_792_458.0  # m/s, exact
H = 6.626_070_15e-34  # J s, exact
HBAR = H / (2.0 * 3.141_592_653_589_793)
KB = 1.380_649e-23  # J/K, exact
EPS0 = 8.854_187_8128e-12  # F/m
A0 = 0.529_177_210_903e-10  # m, Bohr radius (also Steck Table 1)
HARTREE_CM = 219_474.631_363_20  # cm^-1, Hartree energy in wavenumbers
AMU = 1.660_539_066_60e-27  # kg

# Conversion: polarisability in atomic units -> SI (C m^2 / V).
ALPHA_AU_TO_SI = 4.0 * 3.141_592_653_589_793 * EPS0 * A0**3
# Conversion: Steck's alpha/h in Hz/(V/cm)^2 -> atomic units.
# alpha[cm^3] = 5.955_214_861_985_943e-22 * (alpha/h)[Hz/(V/cm)^2]  (Steck sec. 3.2)
_HZ_PER_VCM2_TO_CM3 = 5.955_214_861_985_943e-22
_A0_CM3 = (A0 * 100.0) ** 3
HZ_PER_VCM2_TO_AU = _HZ_PER_VCM2_TO_CM3 / _A0_CM3


@dataclass(frozen=True)
class AlkaliLines:
    """D-line data of an alkali atom (Steck convention for matrix elements)."""

    name: str
    mass_amu: float
    lambda_d1_nm: float  # vacuum wavelength
    lambda_d2_nm: float
    d1_au: float  # <J=1/2||er||J'=1/2>, units of e a0, Steck convention
    d2_au: float  # <J=1/2||er||J'=3/2>
    alpha0_hz_per_vcm2: float  # static ground-state polarisability alpha/h
    source: str


CESIUM = AlkaliLines(
    name="Cs-133",
    mass_amu=132.905_451_931,  # Steck Cs Table 1
    lambda_d1_nm=894.592_959_86,  # Steck Cs Table 4
    lambda_d2_nm=852.347_275_82,  # Steck Cs Table 3
    d1_au=3.1869,  # Steck Cs Table 7
    d2_au=4.4837,  # Steck Cs Table 7
    alpha0_hz_per_vcm2=0.099_77,  # Steck Cs Table 6
    source="Steck, Cesium D Line Data, rev. 2.3.4 (2025), Tables 3, 4, 6, 7",
)

RUBIDIUM87 = AlkaliLines(
    name="Rb-87",
    mass_amu=86.909_180_520,  # Steck Rb87 Table 1
    lambda_d1_nm=794.978_851_156,  # Steck Rb87 D1 table
    lambda_d2_nm=780.241_209_686,  # Steck Rb87 D2 table
    d1_au=2.9931,  # Steck Rb87 Table (dipole matrix elements)
    d2_au=4.227_52,  # Steck Rb87 Table (dipole matrix elements)
    alpha0_hz_per_vcm2=0.0794,  # Steck Rb87 Table (polarizabilities)
    source="Steck, Rubidium 87 D Line Data, rev. 2.3.4 (2025)",
)
