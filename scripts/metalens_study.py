"""Ray-optics design study of wide-field metalens architectures. Usage:
    python scripts/metalens_study.py configs/metalens_study.yaml
For each wavelength, field target and architecture: optimise the radial phase polynomial(s)
for the worst-field Strehl ratio with a common image plane at the target NA; write
metalens_study.json and Strehl-versus-field plots.
"""

import sys
import time

import numpy as np

from metatrap.metalens import design as D
from metatrap.runio import load_config, run_dir, write_json, write_params


def main(cfg_path):
    cfg = load_config(cfg_path)
    out = run_dir(cfg["run_id"])
    write_params(out, cfg)
    na, nt = cfg["na"], cfg["n_terms"]
    res = []
    for lam in cfg["wavelengths_um"]:
        geom = D.Geometry(lam_design_um=lam, cell=cfg["cell"], substrate_um=cfg["substrate_um"])
        st_probe = D.hyperbolic_start(geom, lam, 1000.0, nt)[1]
        R0 = st_probe.aperture_for_na(na, lam)
        L_eff = sum(t / st_probe.indices(lam)[i] for i, (t, m) in enumerate(cfg["cell"]))
        for h in cfg["field_heights_um"]:
            th_max = np.arctan(h / L_eff)
            thetas = th_max * np.array(cfg["field_fractions"])
            for a in cfg["architectures"]:
                arch, extra = a["kind"], dict(a.get("extra", {}))
                t0 = time.time()
                if arch == "singlet":
                    R = R0
                    p0, _ = D.hyperbolic_start(geom, lam, R, nt)
                elif arch == "stop_shifted":
                    R = R0
                    extra["lens_R"] = R + extra["stop_distance_um"] * np.tan(th_max) * 1.2
                    p0, _ = D.hyperbolic_start(geom, lam, extra["lens_R"], nt)
                else:
                    R = R0
                    extra["lens_R"] = R * 1.02 + extra["doublet_substrate_um"] * np.tan(th_max)
                    pf, _ = D.hyperbolic_start(geom, lam, extra["lens_R"], nt)
                    p0 = np.concatenate([np.zeros(nt), pf])
                q0 = D.evaluate(arch, p0, geom, lam, R, nt, extra, thetas, n_rays=36)
                popt, dz, opt = D.optimise(arch, p0, geom, lam, R, nt, extra, thetas, na, n_rays=cfg["n_rays_opt"], maxiter=cfg["maxiter"])
                fine = th_max * np.linspace(0, 1.25, 11)
                q = D.evaluate(arch, popt, geom, lam, R, nt, extra, fine, n_rays=48, dz_common=dz)
                row = {"lambda_um": lam, "field_height_um": h, "arch": arch, "extra": extra, "stop_radius_um": R,
                       "theta_max_deg": float(np.degrees(th_max)),
                       "start_strehl": [qq["strehl"] for qq in q0],
                       "field_fraction": list(np.linspace(0, 1.25, 11)), "strehl": [qq["strehl"] for qq in q],
                       "image_x_um": [qq["x_um"] for qq in q], "dz_common_um": float(dz),
                       "na_achieved": D.marginal_na(arch, popt, geom, lam, R, nt, extra),
                       "min_strehl_to_target": float(min(qq["strehl"] for qq in q[:9])), "params": popt, "seconds": time.time() - t0}
                res.append(row)
                print(f"lam {lam} h {h} {arch} {a.get('extra', {})}: start {np.round(row['start_strehl'], 3)} -> min Strehl to target {row['min_strehl_to_target']:.3f}, NA {row['na_achieved']:.3f}, {row['seconds']:.0f}s", flush=True)
                write_json(out / "metalens_study.json", res)
    plot(res, out)


def plot(res, out):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    keys = sorted({(r["lambda_um"], r["field_height_um"]) for r in res})
    fig, axes = plt.subplots(1, len(keys), figsize=(5 * len(keys), 4), squeeze=False)
    for ax, (lam, h) in zip(axes[0], keys):
        for r in res:
            if (r["lambda_um"], r["field_height_um"]) != (lam, h):
                continue
            lbl = r["arch"] + "".join(f" {k}={v/1000:.0f}mm" for k, v in r["extra"].items() if k in ("stop_distance_um", "doublet_substrate_um"))
            ax.plot(np.array(r["field_fraction"]) * h, r["strehl"], "o-", ms=3, label=lbl)
        ax.axhline(0.8, color="k", lw=0.5)
        ax.set_xlabel("field height at the atoms (um)")
        ax.set_ylabel("Strehl ratio (common focal plane)")
        ax.set_title(f"{lam*1e3:.0f} nm, NA 0.5, target +/-{h:.0f} um")
        ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(out / "strehl_vs_field.png", dpi=120)
    plt.close(fig)


if __name__ == "__main__":
    main(sys.argv[1])
