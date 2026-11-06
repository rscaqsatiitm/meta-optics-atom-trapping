"""S1: power budget of the cavity-enhanced large-period lattice against the tweezer routes.
Usage:
    python scripts/cavity_budget.py configs/cavity_budget.yaml

Lattice: two in-plane axes, each one self-intersecting ring whose two crossing passes form a
4 um lattice with flat-top x Gaussian sheets; one vertical axis, a ring whose two passes cross
in a vertical plane (propagation along the array diagonal) with beams flat-topped in both
directions. Each axis provides the same barrier depth (vertical_ratio scales the vertical).
All equations in src/metatrap/cavity/budget.py and lattice.py. Writes budget.json.
"""

import sys

import numpy as np

from metatrap.cavity import budget as B
from metatrap.cavity import lattice as LA
from metatrap.runio import load_config, run_dir, write_json, write_params


def lattice_geometry(cfg, n_sites):
    lam_um = cfg["wavelength_nm"] * 1e-3
    a = cfg["period_um"]
    th = B.crossing_half_angle(a, lam_um)
    side = (np.sqrt(n_sites) - 1) * a
    n, tol = cfg["flat_order"], cfg["flat_tol"]
    h = LA.array_halfrange(side, th)
    w = LA.supergauss_w_for_flat(h, n, tol)
    weq = LA.supergauss_equivalent_width(w, n)
    # Vertical axis: beams propagate along the diagonal; crossing half-angle theta_z in a vertical plane.
    thz = B.crossing_half_angle(cfg["vertical_period_um"], lam_um)
    hv = side / np.sqrt(2)
    wv = LA.supergauss_w_for_flat(hv, n, tol)
    weq_v = LA.supergauss_equivalent_width(wv, n)
    hz = hv * np.tan(thz)  # vertical walk-off of each beam's centre at the array's diagonal ends
    wz_flat = LA.supergauss_w_for_flat(max(hz, 5.0), n, tol)
    heq_v = LA.supergauss_equivalent_width(wz_flat, n)
    return {"side_um": side, "theta_rad": th, "theta_z_rad": thz, "flat_halfrange_um": h, "w_um": w, "W_eq_um": weq,
            "vertical": {"flat_halfrange_um": hv, "W_eq_um": weq_v, "vertical_flat_half_um": hz, "H_eq_um": heq_v}}


def circulating_powers(cfg, geo, depth_K, contrast_inplane=1.0, vertical_ratio=None):
    lam = cfg["wavelength_nm"]
    vr = cfg["vertical_ratio"] if vertical_ratio is None else vertical_ratio
    i1 = B.beam_intensity_for_barrier(depth_K, lam, contrast_inplane)
    p_axis = float(i1 * B.sheet_peak_intensity(1.0, geo["W_eq_um"] * 1e-6, cfg["w_z_um"] * 1e-6) ** -1)
    i1v = B.beam_intensity_for_barrier(vr * depth_K, lam)
    v = geo["vertical"]
    p_vert = float(i1v * v["W_eq_um"] * 1e-6 * v["H_eq_um"] * 1e-6)
    return {"inplane_per_axis_W": p_axis, "vertical_W": p_vert, "total_circulating_W": 2 * p_axis + p_vert,
            "beam_intensity_W_cm2": float(i1) * 1e-4}


def loss_budgets(cfg):
    out = {}
    nref, ngr = cfg["reflections_per_roundtrip"], cfg["graded_per_roundtrip"]
    tc, bw = cfg["through_cell"], cfg["brewster"]
    bulk = tc["bulk_absorption_per_cm"] * tc["wall_cm"]
    for lm in cfg["mirror_loss"]:
        for lg in cfg["graded_extra_loss"]:
            base = {"mirrors": (nref, lm), "graded_extra": (ngr, lg)}
            out[f"vacuum_m{lm:g}_g{lg:g}"] = {"placement": "in_vacuum", "mirror_loss": lm, "graded_extra": lg, "L": B.roundtrip_loss(base)}
            for la in tc["ar_loss"]:
                d = dict(base, ar=(tc["ar_surfaces_per_roundtrip"], la), bulk=(tc["bulk_passes_per_roundtrip"], bulk))
                out[f"cell_m{lm:g}_g{lg:g}_ar{la:g}"] = {"placement": "through_cell_AR", "mirror_loss": lm, "graded_extra": lg, "ar_loss": la, "L": B.roundtrip_loss(d)}
            d = dict(base, brewster=(bw["surfaces_per_roundtrip"], bw["loss"]), bulk=(tc["bulk_passes_per_roundtrip"], bulk))
            out[f"brewster_m{lm:g}_g{lg:g}"] = {"placement": "through_cell_brewster", "mirror_loss": lm, "graded_extra": lg, "L": B.roundtrip_loss(d)}
    for k, v in out.items():
        L = v["L"]
        v["G_max"] = float(B.g_max(L))
        v["G_design_T3L"] = float(B.enhancement(3 * L, L))  # over-coupled x3: robust choice
        v["finesse_T3L"] = float(B.finesse(3 * L, L))
    return out


def main(cfg_path):
    cfg = load_config(cfg_path)
    out = run_dir(cfg["run_id"])
    write_params(out, cfg)
    lam = cfg["wavelength_nm"]
    w0 = cfg["tweezer_waist_um"] * 1e-6
    a = cfg["period_um"] * 1e-6
    th = B.crossing_half_angle(cfg["period_um"], lam * 1e-3)
    losses = loss_budgets(cfg)
    rep = {"constants": {"depth_per_intensity_K_per_W_m2": B.depth_per_intensity_K(lam), "theta_deg": np.degrees(th),
                          "equal_omega_factor_inplane": float(B.depth_for_equal_omega(1.0, w0, a))},
           "losses": losses, "cases": []}
    zr = np.pi * w0**2 / (lam * 1e-9)
    az = cfg["vertical_period_um"] * 1e-6
    rep["constants"]["equal_omega_factor_vertical"] = float((az / (np.pi * zr)) ** 2)
    sel = {"vacuum": "vacuum_m5e-05_g0", "vacuum_meta": "vacuum_m5e-05_g0.01", "cell_AR": "cell_m5e-05_g0_ar0.001",
           "cell_AR_good": "cell_m5e-05_g0_ar0.00025", "brewster": "brewster_m5e-05_g0"}
    for n in cfg["n_sites"] + [x for x in cfg["scaling_n_sites"] if x not in cfg["n_sites"]]:
        geo = lattice_geometry(cfg, n)
        for d in cfg["depths_mK"]:
            for match in ("equal_depth", "equal_omega"):
                if match == "equal_depth":
                    u_in, vr = d * 1e-3, cfg["vertical_ratio"]
                else:
                    u_in = d * 1e-3 * rep["constants"]["equal_omega_factor_inplane"]
                    vr = rep["constants"]["equal_omega_factor_vertical"] / rep["constants"]["equal_omega_factor_inplane"]
                pc = circulating_powers(cfg, geo, u_in, vertical_ratio=vr)
                pcb = circulating_powers(cfg, geo, u_in, contrast_inplane=np.cos(2 * th), vertical_ratio=vr)
                row = {"n_sites": n, "depth_mK": d, "match": match, "geometry": geo, "circulating": pc,
                       "tweezer_hologram_W": float(B.tweezer_power_W(n, d)),
                       "tweezer_lenslet_W": float(B.tweezer_power_W(n, d, B.LENSLET_SQUARE_MK_SITE_PER_W)),
                       "free_space_lattice_W": pc["total_circulating_W"] / cfg["mode_matching"],
                       "cavity_input_W": {}}
                for lab, key in sel.items():
                    pt = pcb["total_circulating_W"] if lab == "brewster" else pc["total_circulating_W"]
                    G = losses[key]["G_design_T3L"]
                    row["cavity_input_W"][lab] = {"G": G, "P_in_W": pt / (G * cfg["mode_matching"])}
                rep["cases"].append(row)
    # Bounds for tweezers through a cavity or with retro-reflection.
    eh, eo, er = cfg["hologram_efficiency"], cfg["objective_efficiency"], cfg["refractive_objective_efficiency"]
    rep["tweezer_bounds"] = {
        "cavity_G_max_hologram_2pass": float(B.intracavity_element_bound(eh, 2)),
        "cavity_G_max_hologram_plus_meta_objective_4pass": float(1 / (1 - eh**2 * eo**2)),
        "cavity_G_max_meta_objective_pair_4pass": float(B.intracavity_element_bound(eo, 4)),
        "cavity_G_max_refractive_pair_4pass": float(B.intracavity_element_bound(er, 4)),
        "retro_gain_meta_objective": float(B.retro_tweezer_gain(eo**2 * 0.99)),
        "retro_gain_refractive_objective": float(B.retro_tweezer_gain(er**2 * 0.99)),
        "retro_gain_limit": 4.0,
    }
    write_json(out / "budget.json", rep)
    for r in rep["cases"]:
        if r["n_sites"] in cfg["n_sites"]:
            c = r["cavity_input_W"]
            print(f"N={r['n_sites']:.0e} U={r['depth_mK']} mK {r['match']:11s} holo {r['tweezer_hologram_W']:7.0f} W | "
                  f"free lattice {r['free_space_lattice_W']:6.0f} W (axis {r['circulating']['inplane_per_axis_W']:.0f}, vert {r['circulating']['vertical_W']:.0f}) | "
                  + " ".join(f"{k} G{v['G']:.0f}: {v['P_in_W']:.2f} W" for k, v in c.items()), flush=True)
    print(rep["tweezer_bounds"])
    for k in sel.values():
        print(k, losses[k])


if __name__ == "__main__":
    main(sys.argv[1])
