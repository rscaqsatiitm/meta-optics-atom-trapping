"""Array-level power budget: number of sites a laser supports at a given depth.

Model of proposal Section 4 (ii): a fraction ``utilisation`` of the laser power
reaches the tweezers, and each tweezer is an Airy-like focus of effective NA
``na_eff`` (I0 = pi NA_eff^2 P_site / lambda^2). Default parameters are those the
proposal takes from ref. 4 (Wang et al., arXiv:2606.02715): 12.2 W of effective
tweezer power from 33 W incident (utilisation 0.37), and a measured first dark
radius of 1.11 um at 850 nm, i.e. NA_eff = 0.61 x 0.85 / 1.11 = 0.467 (DERIVED).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..constants import CESIUM, AlkaliLines
from .alkali import airy_peak_intensity, depth_K, intensity_for_depth

# Ref. 4 numbers (SOURCE: arXiv:2606.02715 HTML, retrieved 2026-10-31).
REF4_P_INCIDENT_W = 33.0
REF4_P_EFFECTIVE_W = 12.2
REF4_N_TWEEZERS = 18_225
REF4_DARK_RADIUS_UM = 1.11
REF4_LAMBDA_NM = 850.0
REF4_DEPTH_MEASURED_MK = 0.305


@dataclass(frozen=True)
class PowerBudget:
    lam_nm: float = 1075.0
    utilisation: float = REF4_P_EFFECTIVE_W / REF4_P_INCIDENT_W
    na_eff: float = 0.47
    atom: AlkaliLines = CESIUM

    def power_per_site_W(self, depth_K_: float) -> float:
        i0 = intensity_for_depth(depth_K_, self.lam_nm, self.atom)
        lam = self.lam_nm * 1e-9
        return float(i0 * lam**2 / (np.pi * self.na_eff**2))

    def sites(self, laser_power_W, depth_K_: float):
        return self.utilisation * np.asarray(laser_power_W) / self.power_per_site_W(depth_K_)

    def depth_K(self, laser_power_W, n_sites):
        p_site = self.utilisation * np.asarray(laser_power_W) / np.asarray(n_sites)
        i0 = airy_peak_intensity(p_site, self.lam_nm * 1e-9, self.na_eff)
        return depth_K(i0, self.lam_nm, self.atom)
