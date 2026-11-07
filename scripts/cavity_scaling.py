"""How large a step is 10^6 sites? The cavity-lattice analysis at 10^4, 10^5 and 10^6 sites.
Usage:
    python scripts/cavity_scaling.py configs/cavity_scaling.yaml

Collects, for each array size, the results of the same models run at that size (Fox-Li designs and
correction, vertical ring, vertical-ring mirror intensities; configs/*_n1e4.yaml, *_n1e5.yaml, and
the 10^6 runs of the main study) and adds DERIVED quantities:
  - power: tweezers (hologram), lattice in free space, cavity input at G = 300 and 1500
    (P_circ / (G x mode matching)), per ring; and an in-plane variant whose sheet thickness w_z
    scales as sqrt(side) (constant 1 % roll-off along the beam, Rayleigh range ~ side);
  - in-plane mirror intensities: I1 x (horizontal peak ratio from Fox-Li) x (w_z at the atoms /
    w_z on the mirror), the vertical design of the 70 mm layout (f_z = 40 mm) for the telescope
    ring and free Gaussian propagation over the half trip for the plain ring; absorbing-defect
    limit 50 K kappa / I (criterion of scripts/cavity_thermal.py);
  - noise: only the crossing-angle lever (S/2) cot(theta) depends on the array size;
  - loading: atoms needed for a density of two atoms per site volume over the array (Gaussian
    cloud, peak 3x to 5x that density, matched at the array corners) and the number delivered
    into the in-plane sheets (2 w_z thick).
Writes scaling.json and scaling.png.
"""

import json
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from metatrap.cavity import budget as B  # noqa: E402
from metatrap.runio import REPO_ROOT, load_config, run_dir, write_json, write_params  # noqa: E402
from metatrap.thermal.axisym_fd import SUBSTRATES  # noqa: E402

SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]  # validated categorical slots 1-4 (light surface)
INK, INK2 = "#0b0b0b", "#52514e"


def jload(rel):
    p = REPO_ROOT / rel
    return json.loads(p.read_text()) if p.exists() else None


def final_depth(e):
    if e is None:
        return None
    if "final" in e:
        return {"depth_rms": e["final"]["depth_rms"], "g_ratio": e["final"].get("g_ratio")}
    return {"depth_rms": e["depth_rms"], "g_ratio": e.get("g_ratio")}


def correction_summary(corr, variant, cases, seeds=(1, 2)):
    out = {}
    if corr is None:
        return out
    for case in cases:
        for label, ckey in (("uncorrected", "uncorrected"), ("perfect_metrology", "openloop_noise0.0nm"), ("atom_loop", "closedloop_noise0.01")):
            vals = [final_depth(corr.get(f"{variant}_{case}_seed{s}_{ckey}")) for s in seeds]
            vals = [v for v in vals if v]
            if vals:
                out[f"{case}_{label}"] = vals
    return out


def main(cfg_path):
    cfg = load_config(cfg_path)
    out = run_dir(cfg["run_id"])
    write_params(out, cfg)
    bud = jload(cfg["budget_json"])
    noise = jload(cfg["noise_json"])
    lay = jload(cfg["layout_json"])
    vd = [r for r in lay["vertical_designs"] if r.get("f_z_m") == cfg["vertical_f_z_m"]][0]
    lam_m = cfg["wavelength_nm"] * 1e-9
    theta = B.crossing_half_angle(cfg["period_um"], cfg["wavelength_nm"] * 1e-3)
    wz0 = cfg["w_z_um"] * 1e-6
    zr0 = np.pi * wz0**2 / lam_m
    noise_row = [r for r in noise["rows"] if r["axis"] == "in-plane" and r["depth_mK"] == 0.2 and r["heating_time_s"] == 10.0][0]
    side0 = cfg["reference_side_um"]
    rep = {"sizes": {}}
    for sz in cfg["sizes"]:
        n = sz["n_sites"]
        row = {"n_sites": n}
        case = [c for c in bud["cases"] if c["n_sites"] == n and c["depth_mK"] == cfg["depth_mK"] and c["match"] == "equal_depth"][0]
        geo, pc = case["geometry"], dict(case["circulating"])
        # A vertical profile other than the budget's flat-top (e.g. a Gaussian at 10^4): same peak intensity
        # at the atoms, circulating power scaled by the ratio of equivalent heights.
        vfac = 1.0
        if sz.get("scale_vertical_power_to_profile"):
            vfac = jload(sz["vertical_json"])["z_design"]["H_eq_um"] / geo["vertical"]["H_eq_um"]
            pc["vertical_W"] *= vfac
            pc["total_circulating_W"] = 2 * pc["inplane_per_axis_W"] + pc["vertical_W"]
        row["vertical_power_factor"] = vfac
        side = geo["side_um"]
        mm = cfg["mode_matching"]
        wz_scaled = wz0 * np.sqrt(side / side0)
        p_inplane_scaled = pc["inplane_per_axis_W"] * wz_scaled / wz0
        row["geometry"] = {"side_mm": side * 1e-3, "W_eq_mm": geo["W_eq_um"] * 1e-3, "vertical_W_eq_mm": geo["vertical"]["W_eq_um"] * 1e-3,
                           "vertical_H_eq_um": geo["vertical"]["H_eq_um"], "vertical_flat_half_um": geo["vertical"]["vertical_flat_half_um"],
                           "w_z_scaled_um": wz_scaled * 1e6}
        tot = pc["total_circulating_W"]
        row["power"] = {
            "tweezer_hologram_W": case["tweezer_hologram_W"], "tweezer_lenslet_W": case["tweezer_lenslet_W"],
            "free_space_lattice_W": tot / mm,  # = budget value when the vertical profile is the flat-top
            "circulating_inplane_per_axis_W": pc["inplane_per_axis_W"], "circulating_vertical_W": pc["vertical_W"], "circulating_total_W": tot,
            "vertical_share": pc["vertical_W"] / tot,
            "input_total_G300_W": tot / (300 * mm), "input_total_G1500_W": tot / (1500 * mm),
            "input_inplane_ring_G300_W": pc["inplane_per_axis_W"] / (300 * mm), "input_vertical_ring_G300_W": pc["vertical_W"] / (300 * mm),
            "inplane_per_axis_scaled_wz_W": p_inplane_scaled,
            "input_total_G300_scaled_wz_W": (2 * p_inplane_scaled + pc["vertical_W"]) / (300 * mm),
            "tweezer_over_cavity_G300": case["tweezer_hologram_W"] / (tot / (300 * mm)),
        }
        # In-plane rings: Fox-Li designs and correction.
        fox = jload(sz["foxli_json"])
        corr_main = jload(sz["correction_json"])
        corr_plain = jload(sz.get("correction_plain_json", sz["correction_json"]))
        i1 = pc["beam_intensity_W_cm2"] * 1e4
        rings = {}
        for vname, label in (("telK_M0.2_d70", "telescope"), ("plane150", "plain")):
            d = fox.get(f"design_{vname}") if fox else None
            if d is None:
                continue
            spec = d["spectrum"]
            comp = spec["lowest_loss_modes"][0] if spec["lowest_loss_modes"] else None
            r_g = d["graded_mirror_peak_over_flat_intensity"]
            els = d["elements"]
            if label == "telescope":
                r_p2 = [v for k, v in els.items() if k.startswith("element1")][0]["peak_over_flat_intensity"]
                inten = {"P1": i1 * wz0 / vd["w_z_at_primary_m"], "P2": i1 * r_p2 * wz0 / (0.5 * (vd["w_z_at_primary_m"] + vd["w_z_at_graded_m"])),
                         "G": i1 * r_g * wz0 / vd["w_z_at_graded_m"]}
            else:
                p_m = d["p_um"] * 1e-6
                wz_g = wz0 * np.sqrt(1 + (p_m / zr0) ** 2)
                wz_f = wz0 * np.sqrt(1 + (0.5 * p_m / zr0) ** 2)
                inten = {"fold": i1 * wz0 / wz_f, "G": i1 * r_g * wz0 / wz_g}
            imax = max(inten.values())
            kfs, ksa = SUBSTRATES["fused_silica"].kappa, SUBSTRATES["sapphire"].kappa
            cv = sz["variant_keys"][label]
            rings[label] = {
                "fresnel_number_plain": d["fresnel_number_plain"], "ideal_depth_rms": d["array"]["depth_rms"],
                "fundamental_loss": spec["fundamental_loss"], "nearest_mode": comp,
                "graded_mirror_peak_over_flat": r_g, "graded_mirror_pv_rad": d["graded_mirror"].get("pv_rad"),
                "peak_intensity_kW_cm2": {k: v * 1e-7 for k, v in inten.items()},
                "max_defect_um_fused_silica": cfg["defect_dT_K"] * kfs / imax * 1e6,
                "max_defect_um_sapphire": cfg["defect_dT_K"] * ksa / imax * 1e6,
                "errors": correction_summary(corr_main if label == "telescope" else corr_plain, cv, cfg["cases"]),
            }
        row["inplane_rings"] = rings
        # Vertical ring.
        row["adopted_inplane"] = sz["adopted_inplane"]
        row["adopted_vertical_variant"] = sz["adopted_vertical_variant"]
        vs = jload(sz["vertical_scaled_design_json"]) if sz.get("vertical_scaled_design_json") else None
        if vs:
            c = vs["cases"]
            row["vertical_ring_scaled_10to6_design"] = {
                "ideal": c["ideal"], "spectrum_x": vs.get("spectrum_x"),
                "x_uncorrected": [c[k]["depth_rms"] for k in c if "_x_" in k and k.endswith("uncorrected")],
                "x_perfect": [c[k]["depth_rms"] for k in c if "_x_" in k and k.endswith("openloop_perfect")],
                "closed_loop": [{"start": c[k]["history"][0]["depth_rms"], "final": c[k]["history"][-1]["depth_rms"]} for k in c if k.endswith("closedloop")]}
        vert = jload(sz["vertical_json"])
        vm = jload(sz["vertical_mirrors_json"])
        if vert:
            cs = vert["cases"]
            row["vertical_ring"] = {
                "ideal": cs["ideal"], "spectrum_x": vert.get("spectrum_x"), "spectrum_z": vert.get("spectrum_z"),
                "x_uncorrected": [cs[k]["depth_rms"] for k in cs if "_x_" in k and k.endswith("uncorrected")],
                "x_perfect": [cs[k]["depth_rms"] for k in cs if "_x_" in k and k.endswith("openloop_perfect")],
                "z_uncorrected": [cs[k]["depth_rms"] for k in cs if "_z_" in k and k.endswith("uncorrected")],
                "z_perfect": [cs[k]["depth_rms"] for k in cs if "_z_" in k and k.endswith("openloop_perfect")],
                "closed_loop": [{"start": cs[k]["history"][0]["depth_rms"], "final": cs[k]["history"][-1]["depth_rms"],
                                 "g_x": cs[k]["history"][-1]["g_ratio_x"], "g_z": cs[k]["history"][-1]["g_ratio_z"]} for k in cs if k.endswith("closedloop")],
            }
        if vm:
            dk = str(cfg["depth_mK"])
            mir = lambda d: {m: x for m, x in d.items() if not m.startswith("atoms")}  # mirrors only, not the crossing
            row["vertical_mirrors"] = {name: {"peak_kW_cm2_max": vfac * max(x[dk] for x in mir(v["peak_intensity_kW_cm2"]).values()),
                                              "peak_kW_cm2": {m: vfac * x[str(cfg["depth_mK"])] for m, x in v["peak_intensity_kW_cm2"].items()},
                                              "min_defect_um_fused_silica": min(x[dk] for x in mir(v["max_defect_radius_um"]["fused_silica"]).values()) / vfac,
                                              "min_defect_um_sapphire": min(x[dk] for x in mir(v["max_defect_radius_um"]["sapphire"]).values()) / vfac,
                                              "ideal_depth_rms": v["ideal"]["depth_rms"], "g_ratio_z": v["ideal"]["g_ratio_z"],
                                              "z_uncorrected": [v[k]["depth_rms"] for k in v if k.endswith("uncorrected")],
                                              "z_perfect": [v[k]["depth_rms"] for k in v if k.endswith("openloop_perfect")],
                                              "z_nearest_mode": v["spectrum_z"]["lowest_loss_modes"][0] if v["spectrum_z"]["lowest_loss_modes"] else None}
                                       for name, v in vm["variants"].items()}
        # Noise: crossing-angle requirement scales as 1 / side; everything else is size-independent.
        row["noise_10s_0p2mK"] = {"RIN_dBc_per_Hz": noise_row["RIN_at_2nu_dBc_per_Hz"], "mirror_pm_per_rtHz": noise_row["mirror_displacement_pm_per_rtHz"],
                                  "crossing_angle_nrad_per_rtHz": noise_row["crossing_angle_nrad_per_rtHz_edge_sites"] * side0 / side,
                                  "edge_lever_mm_per_rad": side * 1e-3 / 2 / np.tan(theta)}
        # Loading.
        site_cm3 = (cfg["period_um"] * 1e-4) ** 2 * cfg["vertical_period_um"] * 1e-4
        n0 = cfg["atoms_per_site_needed"] / site_cm3
        rc = side * 1e-4 / np.sqrt(2)
        cloud = {}
        for k in cfg["cloud_peak_over_needed"]:
            sig = rc / np.sqrt(2 * np.log(k))
            cloud[f"peak_{k:g}x"] = float(k * n0 * (2 * np.pi) ** 1.5 * sig**3)
        sheet = n0 * (side * 1e-4) ** 2 * 2 * wz0 * 1e2
        row["loading"] = {"density_needed_cm3": n0, "mot_atoms_cloud": cloud, "atoms_into_sheets": sheet,
                          "planes_in_sheets": 2 * wz0 / (cfg["vertical_period_um"] * 1e-6)}
        rep["sizes"][f"{n:.0e}"] = row
    write_json(out / "scaling.json", rep)
    plot(rep, out / "scaling.png")
    for k, r in rep["sizes"].items():
        p = r["power"]
        print(k, {"tweezer": round(p["tweezer_hologram_W"], 1), "free": round(p["free_space_lattice_W"], 1), "G300": round(p["input_total_G300_W"], 3),
                  "vert_share": round(p["vertical_share"], 2)})
        for lab, ring in r["inplane_rings"].items():
            print("  ", lab, {"F": round(ring["fresnel_number_plain"], 2), "ideal": round(ring["ideal_depth_rms"] * 100, 3), "loss": ring["fundamental_loss"],
                              "mode": ring["nearest_mode"], "I": {a: round(b, 1) for a, b in ring["peak_intensity_kW_cm2"].items()},
                              "err": {a: [round(x["depth_rms"] * 100, 2) for x in b] for a, b in ring["errors"].items()}})


def plot(rep, path):
    keys = list(rep["sizes"])
    n = np.array([rep["sizes"][k]["n_sites"] for k in keys])
    fig, ax = plt.subplots(1, 3, figsize=(13, 4.0))
    for a in ax:
        a.set_xscale("log")
        a.set_xticks(n)
        a.set_xticklabels([f"$10^{{{int(np.log10(x))}}}$" for x in n])
        a.set_xlabel("Number of sites", color=INK2)
        a.grid(True, which="major", color="#e6e5e1", lw=0.8)
        a.tick_params(colors=INK2)
        for s in ("top", "right"):
            a.spines[s].set_visible(False)
    p = [rep["sizes"][k]["power"] for k in keys]
    series = [("Hologram tweezers", [x["tweezer_hologram_W"] for x in p], "o"), ("Lattice, free space", [x["free_space_lattice_W"] for x in p], "s"),
              ("Cavity lattice, G = 300", [x["input_total_G300_W"] for x in p], "^"), ("Cavity lattice, G = 1500", [x["input_total_G1500_W"] for x in p], "D")]
    for (lab, y, mk), c in zip(series, SERIES):
        ax[0].plot(n, y, color=c, lw=2, marker=mk, ms=8, label=lab)
    ax[0].set_yscale("log")
    ax[0].set_ylabel("Laser power at the optics input (W)", color=INK2)
    ax[0].set_title("(a) Power, 0.2 mK barrier", color=INK, loc="left")
    ax[0].legend(frameon=False, fontsize=8)
    def err(k, ring, kind):
        e = rep["sizes"][k]["inplane_rings"].get(ring, {}).get("errors", {}).get(kind)
        return 100 * float(np.median([v["depth_rms"] for v in e])) if e else None

    lines_b = [("Adopted design, uncorrected", lambda k: err(k, rep["sizes"][k]["adopted_inplane"], "fig1000_uncorrected"), "^"),
               ("Adopted design, atom loop", lambda k: err(k, rep["sizes"][k]["adopted_inplane"], "fig1000_atom_loop"), "v"),
               ("0.2x telescope ring at every size, uncorrected", lambda k: err(k, "telescope", "fig1000_uncorrected"), "s"),
               ("Plain ring at every size, uncorrected", lambda k: err(k, "plain", "fig1000_uncorrected"), "D")]
    order = [2, 3, 0, 1]  # comparisons (dashed) underneath, adopted design on top
    for i in order:
        (lab, f, mk), c = lines_b[i], SERIES[i]
        pts = [(x, f(k)) for x, k in zip(n, keys) if f(k) is not None]
        if pts:
            ax[1].plot(*zip(*pts), color=c, lw=2 if i < 2 else 1.5, ls="-" if i < 2 else "--", marker=mk, ms=8 if i < 2 else 10,
                       mfc=c if i < 2 else "none", mew=2, label=lab, zorder=3 if i < 2 else 2)
    h, l = ax[1].get_legend_handles_labels()
    ax[1].legend([h[l.index(x[0])] for x in lines_b if x[0] in l], [x[0] for x in lines_b if x[0] in l], frameon=False, fontsize=8)
    ax[1].set_yscale("log")
    ax[1].set_ylabel("Depth non-uniformity, rms (%)", color=INK2)
    ax[1].set_title("(b) $\\lambda$/1000 figure error on every mirror", color=INK, loc="left")
    for (lab, getter, mk), c in zip([("In-plane, 0.2x telescope ring (the 10^6 design)", lambda r: max(r["inplane_rings"]["telescope"]["peak_intensity_kW_cm2"].values()), "^"),
                                     ("In-plane, plain ring", lambda r: max(r["inplane_rings"]["plain"]["peak_intensity_kW_cm2"].values()), "s"),
                                     ("Vertical ring, adopted design", lambda r: r["vertical_mirrors"][r["adopted_vertical_variant"]]["peak_kW_cm2_max"], "o"),
                                     ("Vertical ring without vertical telescope", lambda r: r["vertical_mirrors"]["no_z_telescope"]["peak_kW_cm2_max"], "D")], SERIES):
        xs, ys = [], []
        for x, k in zip(n, keys):
            try:
                ys.append(getter(rep["sizes"][k]))
                xs.append(x)
            except (KeyError, TypeError):
                pass
        ax[2].plot(xs, ys, color=c, lw=2, marker=mk, ms=8, label=lab, zorder=4 if "adopted" in lab else 2)
    ax[2].set_yscale("log")
    ax[2].set_ylabel("Peak intensity on the hottest mirror (kW/cm$^2$)", color=INK2)
    ax[2].set_title("(c) Mirror intensity, 0.2 mK", color=INK, loc="left")
    ax[2].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=200, facecolor="#fcfcfb")


if __name__ == "__main__":
    main(sys.argv[1])
