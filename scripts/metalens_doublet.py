"""Air-spaced aplanatic metasurface doublet (side-B objective): design, field, chromatic and
tolerance study in ray optics, plus an axisymmetric wave-optics check.

Usage: python scripts/metalens_doublet.py configs/metalens_doublet.yaml

Results (results/<run_id>/doublet.json, checkpointed per section): design data, spacing scan,
field headroom, chromatic behaviour of one phase map and the phase pairs a multi-wavelength
library must supply, mechanical tolerances, and the wave-optics focal check.
"""

import json
import sys
import time

import numpy as np

from metatrap.metalens import doublet as D
from metatrap.runio import load_config, run_dir, write_json, write_params


def spec_for(cfg, spacing_um=None):
    sp = [list(l) for l in cfg["spacer"]]
    if spacing_um is not None:
        sp[0][0] = spacing_um - sum(float(l[0]) for l in sp[1:])
    return D.DoubletSpec(cfg["na"], cfg["cell"], sp)


def summary_fields(rows):
    return [{k: (round(v, 5) if isinstance(v, float) else v) for k, v in r.items()} for r in rows]


def main(cfg_path):
    cfg = load_config(cfg_path)
    out = run_dir(cfg["run_id"])
    write_params(out, cfg)
    path = out / "doublet.json"
    res = json.load(open(path)) if path.exists() else {}

    def save():
        write_json(path, res)

    lams = cfg["wavelengths_um"]
    lam0 = lams[0]
    spec = spec_for(cfg)
    d0 = D.make(spec, lam0)
    R1 = d0.info["R1_um"]
    f = d0.info["f_um"]
    designs = {lam: (d0 if lam == lam0 else D.make(spec, lam, R1_um=R1)) for lam in lams}
    heights = cfg["field_heights_um"]

    # 1. Design data.
    if "design" not in res:
        res["design"] = {str(l): designs[l].info for l in lams}
        save()
        print("design", res["design"], flush=True)

    # 2. Spacing scan: common-plane Strehl over the array field, per wavelength (own design).
    if "spacing_scan" not in res:
        rows = []
        for sp in cfg["spacing_scan_um"]:
            s2 = spec_for(cfg, sp)
            for lam in lams:
                dd = D.make(s2, lam, R1_um=R1)
                sysm = D.system(dd, lam)
                dz, q = D.common_plane(sysm, R1, lam, heights, f)
                rows.append({"spacing_um": sp, "lambda_um": lam, "dz_common_um": dz, "strehl_common": [r["strehl"] for r in q],
                             "min_strehl": min(r["strehl"] for r in q), "max_deflection_corrector": dd.info["max_deflection_corrector"]})
                print("spacing", rows[-1], flush=True)
        res["spacing_scan"] = rows
        save()

    # 3. Field headroom (free focus per field) and the common plane over the array field.
    if "headroom" not in res:
        hr = {}
        for lam in lams:
            sysm = D.system(designs[lam], lam)
            free = D.field_scan(sysm, R1, lam, cfg["headroom_heights_um"], f)
            dz, q = D.common_plane(sysm, R1, lam, heights, f)
            hr[str(lam)] = {"free_focus": summary_fields(free), "common_plane_dz_um": dz, "common_plane": summary_fields(q),
                            "marginal_na": D.marginal_na(sysm, R1, lam)}
            print("headroom", lam, [(r["height_um"], round(r["strehl"], 4), round(r["dz_um"], 2)) for r in free], "common", dz, [round(r["strehl"], 4) for r in q], flush=True)
        res["headroom"] = hr
        save()

    # 4. Chromatic behaviour of one phase map (designed at lam0), and multi-wavelength phase pairs.
    if "chromatic" not in res:
        ch = {"detuning": []}
        for dl_nm in cfg["chromatic"]["detunings_nm"] + [round((l - lam0) * 1e3, 2) for l in lams[1:]]:
            lam = lam0 + dl_nm * 1e-3
            sysm = D.system(d0, lam)
            fixed = D.quality(sysm, R1, lam, 0.0, None, 0.0)
            free = D.quality(sysm, R1, lam, 0.0)
            row = {"detuning_nm": dl_nm, "lambda_um": lam, "strehl_fixed_plane": fixed["strehl"], "free_focus_dz_um": free["dz_um"],
                   "strehl_free_focus": free["strehl"]}
            try:
                v = D.vergence_for_focus(sysm, R1, lam, 0.0, v_max=cfg["chromatic"]["v_max_per_um"])
                L = None if v == 0 else 1.0 / v
                fr = D.field_scan(sysm, R1, lam, [0.0, 200.0, 400.0], f, L, dz_common=0.0)
                row |= {"source_distance_um": (1.0 / v if v else None), "strehl_vergence_corrected": [r["strehl"] for r in fr],
                        "rms_waves_vergence_corrected": fr[0]["rms_waves"], "image_x_at_400_um": fr[2]["x_um"], "marginal_na": D.marginal_na(sysm, R1, lam, L)}
            except RuntimeError as e:
                row["vergence_error"] = str(e)
            ch["detuning"].append(row)
            print("chromatic", {k: (round(v, 5) if isinstance(v, float) else v) for k, v in row.items()}, flush=True)
        pairs = {}
        for i, la in enumerate(lams):
            for lb in lams[i + 1 :]:
                da, db = designs[la], designs[lb]
                mn, mx, H = D.phase_pair_uniformity(da.m2, db.m2, la, lb, da.info["R2_um"])
                mn1, mx1, _ = D.phase_pair_uniformity(da.m1, db.m1, la, lb, R1)
                pairs[f"{la}/{lb}"] = {"focusing_min_over_uniform": mn, "focusing_max_over_uniform": mx, "corrector_min_over_uniform": mn1,
                                       "corrector_max_over_uniform": mx1}
                np.save(out / f"phase_pairs_{la}_{lb}.npy", H)
                # Phase difference between the two maps, in waves at the edge (how often it wraps).
                r = np.array([da.info["R2_um"]])
                pairs[f"{la}/{lb}"]["focusing_phase_difference_waves_at_edge"] = float((da.m2.phase(r, la) - db.m2.phase(r, lb))[0] / (2 * np.pi))
        ch["phase_pairs"] = pairs
        print("phase pairs", pairs, flush=True)
        res["chromatic"] = ch
        save()

    # 5. Mechanical tolerances (at lam0): refocus allowed (common plane re-optimised).
    if "tolerances" not in res:
        tol = []
        hx = cfg["tolerance_heights_um"]
        sysn = D.system(d0, lam0)
        nom = {h: D.quality(sysn, R1, lam0, D.theta_for_height(sysn, R1, lam0, h, f), None, None)["x_um"] for h in hx}
        for case in cfg["tolerances"]:
            pert = D.Perturbation(**{k: v for k, v in case.items() if k not in ("name", "compensate")})
            if case.get("compensate"):
                # Alignment by coma nulling: the focusing-plate decentre that maximises the on-axis
                # Strehl ratio (a tilt eps acts like a decentre of ~ spacing x eps).
                from scipy.optimize import minimize_scalar

                c0 = float(d0.spec.spacer[0][0]) * pert.m2_tilt_rad

                def neg(dec):
                    p2 = D.Perturbation(**(vars(pert) | {"m2_decentre_um": dec}))
                    return -D.best_plane(D.system(d0, lam0, p2), R1, lam0, [0.0])[1][0]["strehl"]

                rr = minimize_scalar(neg, bounds=(c0 - 20, c0 + 20), method="bounded", options={"xatol": 0.01})
                pert.m2_decentre_um = float(rr.x)
                case = case | {"compensating_decentre_um": float(rr.x)}
            sysm = D.system(d0, lam0, pert)
            # Same input angles as the nominal system: displacement and magnification are measured.
            # (Tilts rotate the downstream system about the plate vertex, so the absolute image
            # shift of a tilt case mostly reflects that pivot choice; the Strehl ratios do not.)
            thetas = [D.theta_for_height(sysn, R1, lam0, h, f) for h in hx]
            dz, q = D.best_plane(sysm, R1, lam0, thetas)
            r = type("R", (), {"x": dz})
            row = {"case": case, "dz_common_um": float(r.x), "strehl": [x["strehl"] for x in q], "min_strehl": min(x["strehl"] for x in q),
                   "image_shift_um": [x["x_um"] - nom[h] for x, h in zip(q, hx)], "heights_um": hx}
            tol.append(row)
            print("tol", case["name"], round(row["dz_common_um"], 3), np.round(row["strehl"], 4), np.round(row["image_shift_um"], 3), flush=True)
        res["tolerances"] = tol
        save()

    # 7. Baselines: hyperbolic singlet (stop at the lens) and the aplanat on a 2 mm silica spacer.
    if cfg.get("baselines") and "baselines" not in res:
        from metatrap.metalens.aplanat import SplineMeta
        from metatrap.metalens.raytrace import Medium, System
        from metatrap.propagation.layers import LayerStack

        bl = {}
        hb = cfg["baselines"]["heights_um"]
        st = LayerStack.from_list(cfg["cell"])
        R = st.aperture_for_na(cfg["na"], lam0)
        r = np.linspace(0, 1.3 * R, 6001)
        hyp = SplineMeta(r, st.lens_phase(r, lam0), lam0, "hyperbolic")
        sysh = System([hyp] + D._media(cfg["cell"], lam0))
        fh = R / cfg["na"]
        free = D.field_scan(sysh, R, lam0, hb, fh)
        bl["hyperbolic_singlet"] = {"free_focus": summary_fields(free)}
        sub = D.DoubletSpec(cfg["na"], cfg["cell"], [[cfg["baselines"]["silica_spacer_um"], "fused_silica"]])
        ds = D.make(sub, lam0)
        syss = D.system(ds, lam0)
        free2 = D.field_scan(syss, ds.info["R1_um"], lam0, hb, ds.info["f_um"])
        bl["silica_spaced_aplanat"] = {"free_focus": summary_fields(free2), "info": ds.info}
        for k, v in bl.items():
            print("baseline", k, [(x["height_um"], round(x["strehl"], 4)) for x in v["free_focus"]], flush=True)
        res["baselines"] = bl
        save()

    # 6. Axisymmetric wave-optics check of the on-axis focus.
    if cfg.get("wave_axial") and "wave_axial" not in res:
        from metatrap.metalens import wave_axial as W

        wa = {}
        for lam in lams:
            t0 = time.time()
            r = W.on_axis_focus(designs[lam], lam, **cfg["wave_axial"])
            wa[str(lam)] = r | {"seconds": time.time() - t0}
            print("wave", lam, {k: v for k, v in r.items() if not isinstance(v, list)}, flush=True)
        res["wave_axial"] = wa
        save()

    # 8. Off-axis wave check (superposition from the focusing plate), optionally apodised by the
    #    meta-atom efficiency versus deflection (results/<run_id>/efficiency.json).
    wo = cfg.get("wave_offaxis")
    if wo:
        from metatrap.metalens import wave_offaxis as WO

        res.setdefault("wave_offaxis", {})
        eff_path = out / "efficiency.json"
        eff = json.load(open(eff_path)) if eff_path.exists() else {}
        for lam in wo["wavelengths_um"]:
            dd = designs[lam]
            sysm = D.system(dd, lam)
            for h in wo["heights_um"]:
                th = D.theta_for_height(sysm, R1, lam, h, f)
                q = D.quality(sysm, R1, lam, th)
                cases = [("ideal", None)]
                if h in wo.get("apodised_heights_um", []):
                    for mat, e in eff.get(str(lam), {}).items():
                        cases.append((mat, (np.array(e["sin_theta"]), np.array(e["eta_xy_mean"]))))
                for tag, eta in cases:
                    key = f"{lam}:{h}:{tag}"
                    if key in res["wave_offaxis"]:
                        continue
                    t0 = time.time()
                    r = WO.strehl_offaxis(dd, lam, th, q["x_um"], q["dz_um"], eta=eta, n_pupil=wo["n_pupil"])
                    res["wave_offaxis"][key] = r | {"ray_strehl": q["strehl"], "ray_x_um": q["x_um"], "ray_dz_um": q["dz_um"], "seconds": time.time() - t0}
                    print("wave_offaxis", key, {k: (round(v, 4) if isinstance(v, float) else v) for k, v in res["wave_offaxis"][key].items()}, flush=True)
                    save()


if __name__ == "__main__":
    main(sys.argv[1])
