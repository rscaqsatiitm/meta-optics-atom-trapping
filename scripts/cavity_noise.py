"""Noise requirements of the cavity lattice: parametric and position-noise heating.
Usage:
    python scripts/cavity_noise.py configs/cavity_noise.yaml

Heating model of Savard, O'Hara and Thomas, Phys. Rev. A 56, R1095 (1997) (text read in session):
  intensity noise (their Eq. 12):  1 / T_I = pi^2 nu^2 S_eps(2 nu)
  position noise  (their Eq. 16):  1 / T_x = pi^2 nu^2 S_x(nu) / <x^2>
with S one-sided power spectra per Hz, nu the trap frequency in Hz, T_I the energy e-folding time and
T_x the energy-doubling time; <x^2> = k_B T / (m (2 pi nu)^2) for atoms at temperature T.
Couplings in the self-intersecting ring (DERIVED, docs/cavity_concept.md Section 6c):
  site shift = a dphi / (2 pi), dphi = k L1 + ..., L1 the path between the two crossings (X -> X),
  - laser frequency change dnu (field following the laser within the cavity linewidth):
        shift = a L1 dnu / c
  - normal displacement d of one mirror at incidence AOI, laser locked to the round trip L_rt
    (the lock re-tunes k):  shift = a (2 d cos AOI / lambda) (1 - L1 / L_rt)  [= a d cos AOI / lambda for L1 = L_rt / 2]
  - crossing-angle change dtheta: pitch change cot(theta) dtheta, edge-site shift (S/2) cot(theta) dtheta
  - residual laser-cavity detuning dnu_r (Gaussian, flat spectrum over bandwidth B, linewidth kappa FWHM):
        eps = -(2 dnu_r / kappa)^2, S_eps ~ 32 S_dnu^2 B / kappa^4 (zero offset);
        with a static offset d0: S_eps = 64 d0^2 S_dnu / kappa^4.
Writes noise.json (requirements, as one-sided amplitude spectral densities at the relevant frequency).
"""

import sys

import numpy as np

from metatrap.cavity import budget as B
from metatrap.constants import AMU, CESIUM, KB
from metatrap.runio import load_config, run_dir, write_json, write_params

C = 299792458.0
H = 6.62607015e-34


def heating_intensity_limit(nu, T_s):
    """Largest S_eps(2 nu) (1/Hz) for an energy e-folding time T_s."""
    return 1.0 / (T_s * np.pi**2 * nu**2)


def heating_position_limit(nu, T_s, x2):
    """Largest S_x(nu) (m^2/Hz) for an energy-doubling time T_s."""
    return x2 / (T_s * np.pi**2 * nu**2)


def main(cfg_path):
    cfg = load_config(cfg_path)
    out = run_dir(cfg["run_id"])
    write_params(out, cfg)
    m = CESIUM.mass_amu * AMU
    lam = cfg["wavelength_nm"] * 1e-9
    # Check against the numerical examples of Savard et al. (100 s, 10 kHz).
    check = {"intensity_sqrtS_for_100s_10kHz": float(np.sqrt(heating_intensity_limit(10e3, 100.0))),
             "position_sqrtS_um_for_100s_10kHz_1um": float(np.sqrt(heating_position_limit(10e3, 100.0, 1e-12)) * 1e6)}
    rows = []
    for ax in cfg["axes"]:
        a = ax["period_um"] * 1e-6
        theta = B.crossing_half_angle(ax["period_um"], cfg["wavelength_nm"] * 1e-3)
        for depth_mK in cfg["depths_mK"]:
            U = depth_mK * 1e-3 * ax.get("depth_scale", 1.0)
            nu = float(B.lattice_omega(U, a) / (2 * np.pi))
            x2 = KB * cfg["atom_temperature_K"] / (m * (2 * np.pi * nu) ** 2)
            for T_s in cfg["heating_times_s"]:
                s_eps = heating_intensity_limit(nu, T_s)
                s_x = heating_position_limit(nu, T_s, x2)
                lever_mirror = a * 2 * np.cos(np.radians(ax["typical_aoi_deg"])) / lam * (1 - ax["L1_m"] / ax["Lrt_m"])
                lever_freq = a * ax["L1_m"] / C
                lever_angle = cfg["array_side_um"] * 1e-6 / 2 / np.tan(theta)
                row = {"axis": ax["name"], "depth_mK": depth_mK, "trap_freq_kHz": nu / 1e3, "x_rms_um": float(np.sqrt(x2) * 1e6),
                       "heating_time_s": T_s,
                       "heating_rate_uK_per_s": float(cfg["atom_temperature_K"] / T_s * 1e6),
                       "heating_rate_quanta_per_s": float(KB * cfg["atom_temperature_K"] / T_s / (H * nu)),
                       "RIN_at_2nu_dBc_per_Hz": float(10 * np.log10(s_eps)),
                       "site_position_pm_per_rtHz": float(np.sqrt(s_x) * 1e12),
                       "mirror_displacement_pm_per_rtHz": float(np.sqrt(s_x) / lever_mirror * 1e12),
                       "laser_frequency_Hz_per_rtHz": float(np.sqrt(s_x) / lever_freq),
                       "crossing_angle_nrad_per_rtHz_edge_sites": float(np.sqrt(s_x) / lever_angle * 1e9),
                       "levers": {"site_shift_per_mirror_displacement": lever_mirror, "site_shift_m_per_Hz": lever_freq,
                                  "edge_shift_m_per_rad": lever_angle}}
                for g in cfg["enhancements"]:
                    loss = 0.75 / g  # over-coupled design T = 3L gives G = 0.75 / L (Section 3)
                    kappa = ax["fsr_Hz"] / float(B.finesse(3 * loss, loss))  # FWHM linewidth
                    bw = cfg["lock_noise_bandwidth_Hz"]
                    s_dnu_zero = np.sqrt(s_eps * kappa**4 / (32 * bw))
                    s_dnu_off = s_eps * kappa**4 / (64 * (cfg["lock_offset_fraction"] * kappa) ** 2)
                    row[f"G{g:.0f}"] = {"linewidth_kHz": kappa / 1e3,
                                         "lock_residual_Hz_per_rtHz_zero_offset": float(np.sqrt(s_dnu_zero)),
                                         "lock_residual_rms_kHz_zero_offset": float(np.sqrt(s_dnu_zero * bw) / 1e3),
                                         "lock_residual_Hz_per_rtHz_with_offset": float(np.sqrt(s_dnu_off))}
                rows.append(row)
    write_json(out / "noise.json", {"savard_check": check, "rows": rows})
    print(check)
    for r in rows:
        print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in r.items() if k != "levers"})


if __name__ == "__main__":
    main(sys.argv[1])
