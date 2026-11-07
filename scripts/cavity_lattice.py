"""S2: site shape, trap frequencies and registration sensitivities of the large-period lattice.
Usage:
    python scripts/cavity_lattice.py configs/cavity_lattice.yaml

Plane-wave model of the three mutually incoherent lattice axes at the array centre (each axis
two beams of equal intensity, barrier U per axis), compared with the meta-hologram tweezer of
the same barrier depth. Numerical Hessian at the site gives the trap frequencies. Also: the
saddle (barrier) between neighbours, the potential maximum, and the site displacement at the
array edge per unit error of crossing angle and wavelength. Writes lattice.json.
"""

import sys

import numpy as np

from metatrap.cavity import budget as B
from metatrap.cavity import lattice as LA
from metatrap.constants import AMU, CESIUM, KB
from metatrap.runio import load_config, run_dir, write_json, write_params


def main(cfg_path):
    cfg = load_config(cfg_path)
    out = run_dir(cfg["run_id"])
    write_params(out, cfg)
    lam = cfg["wavelength_nm"] * 1e-3  # um
    k = 2 * np.pi / lam
    a, az = cfg["period_um"], cfg["vertical_period_um"]
    th, thz = B.crossing_half_angle(a, lam), B.crossing_half_angle(az, lam)
    kappa = B.depth_per_intensity_K(cfg["wavelength_nm"])
    m = CESIUM.mass_amu * AMU
    rep = {"theta_deg": np.degrees(th), "theta_z_deg": np.degrees(thz), "cases": []}
    c, s, cz, sz = np.cos(th), np.sin(th), np.cos(thz), np.sin(thz)
    d = np.array([1, 1, 0]) / np.sqrt(2)
    dp = np.array([-1, 1, 0]) / np.sqrt(2)
    for U_mK in cfg["depths_mK"]:
        U = U_mK * 1e-3
        i1 = U / (4 * kappa)  # each beam's intensity for barrier U (equal beams, contrast 1)
        amp = np.sqrt(i1)
        fam_y = [(amp, k * np.array([c, s, 0])), (amp, k * np.array([c, -s, 0]))]  # fringes along y
        fam_x = [(amp, k * np.array([s, c, 0])), (amp, k * np.array([-s, c, 0]))]  # fringes along x
        az_amp = np.sqrt(cfg["vertical_ratio"]) * amp
        fam_z = [(az_amp, k * (cz * d + sz * np.array([0, 0, 1]))), (az_amp, k * (cz * d - sz * np.array([0, 0, 1])))]
        fams = [fam_x, fam_y, fam_z]

        def V(p):
            return -kappa * KB * LA.plane_wave_intensity(np.asarray(p, float), fams)  # J

        h = 1e-3  # um
        hess = np.zeros((3, 3))
        e = np.eye(3) * h
        for i in range(3):
            for j in range(3):
                hess[i, j] = (V(e[i] + e[j]) - V(e[i] - e[j]) - V(-e[i] + e[j]) + V(-e[i] - e[j])) / (4 * h * h)
        w2 = np.linalg.eigvalsh(hess) / 1e-12 / m  # per um^2 -> per m^2
        f = np.sqrt(np.abs(w2)) / (2 * np.pi)
        v0 = V([0, 0, 0])
        saddle_x = V([a / 2, 0, 0]) - v0
        saddle_z = V([0, 0, az / 2]) - v0
        corner = V([a / 2, a / 2, 0]) - v0
        w0 = cfg["tweezer_waist_um"] * 1e-6
        f_tw_r = B.tweezer_omega_radial(U, w0) / (2 * np.pi)
        zr = np.pi * w0**2 / (lam * 1e-6)
        f_tw_z = np.sqrt(2 * U * KB / (m * zr**2)) / (2 * np.pi)
        rep["cases"].append({
            "barrier_mK": U_mK, "beam_intensity_W_cm2": i1 * 1e-4,
            "trap_freq_kHz_sorted": (f / 1e3).tolist(),
            "analytic_inplane_kHz": float(B.lattice_omega(U, a * 1e-6) / (2 * np.pi) / 1e3),
            "analytic_vertical_kHz": float(B.lattice_omega(cfg["vertical_ratio"] * U, az * 1e-6) / (2 * np.pi) / 1e3),
            "barrier_x_mK": saddle_x / KB * 1e3, "barrier_z_mK": saddle_z / KB * 1e3, "cell_corner_mK": corner / KB * 1e3,
            "site_light_shift_mK": -v0 / KB * 1e3,
            "tweezer_radial_kHz": f_tw_r / 1e3, "tweezer_axial_kHz": f_tw_z / 1e3,
        })
    side = (np.sqrt(cfg["n_sites"]) - 1) * a
    edge = side / 2
    rep["registration"] = {
        "array_side_um": side,
        "site_shift_at_edge_per_urad_um": edge / np.tan(th) * 1e-6,
        "site_shift_at_edge_per_ppm_wavelength_um": edge * 1e-6,
        "site_shift_per_nm_differential_wavefront_um": 1e-3 / (2 * s),
        "laser_detuning_100MHz_ppm": 100e6 / (299792458.0 / (lam * 1e-6)) * 1e6,
    }
    write_json(out / "lattice.json", rep)
    for r in rep["cases"]:
        print(r)
    print(rep["registration"])


if __name__ == "__main__":
    main(sys.argv[1])
