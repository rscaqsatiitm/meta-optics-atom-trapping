"""Two-wavelength efficiency of the 852 / 1064 nm meta-atom screen.

Usage: python scripts/analyse_dual_library.py configs/dual_library.yaml

For each material, pitch and height: the ring family, the two-pillar family and their union
(same lattice and height, so they can share one metasurface), with the aspect-ratio limit 10
and 15. The figure of merit is molecules.dual_efficiency: the focusing efficiency at each
wavelength when every phase pair on the torus must be realised (the aplanatic doublet needs a
nearly uniform torus, docs/metalens_report.md Section 5). The single-wavelength optimum of
the same atoms is given for reference. Output: results/<run_id>/dual_efficiency.json and a
scatter plot of the achievable phase pairs for the best case.
"""

import json
import sys

import numpy as np

from metatrap.rcwa.molecules import dual_efficiency
from metatrap.runio import REPO_ROOT, load_config, run_dir, write_json


def min_feature(fam, params):
    p = np.asarray(params, float)
    if fam == "ring":
        return np.where(p[:, 1] > 0, (p[:, 0] - p[:, 1]) / 2, p[:, 0])
    w = np.where(p > 0, p, np.inf)
    return w.min(axis=1)


def load(od, fam, mat, p, h, lams):
    out = []
    for lam in lams:
        f = od / f"{fam}_{mat}_p{p:.2f}_h{h:.1f}_l{lam:.4f}.npz"
        if not f.exists():
            return None
        out.append(np.load(f))
    return out


def plates_eval(cfg, od, out, la, lb):
    """Efficiency of the best libraries for the actual phase-pair distributions of the two
    doublet plates (configs/metalens_doublet.yaml designs; 852.3 / 1064.4 nm maps evaluated
    against the 852.0 / 1064.4 nm library)."""
    from metatrap.metalens import doublet as D
    from metatrap.rcwa.molecules import plate_efficiency

    mc = load_config(REPO_ROOT / "configs/metalens_doublet.yaml")
    spec = D.DoubletSpec(mc["na"], mc["cell"], [list(l) for l in mc["spacer"]])
    d1 = D.make(spec, 0.8523)
    d2 = D.make(spec, 1.0644, R1_um=d1.info["R1_um"])
    R1 = d1.info["R1_um"]
    r = np.linspace(0, R1, 20001)
    from metatrap.propagation.layers import LayerStack

    h2 = LayerStack.from_list(mc["cell"]).ray_height(r / d1.info["f_um"], 0.8523)  # zone followed by the ray of h1
    plates = {"corrector": (d1.m1.phase(r, 0.8523), d2.m1.phase(r, 1.0644)),
              "focusing": (d1.m2.phase(h2, 0.8523), d2.m2.phase(h2, 1.0644))}
    res = []
    for mat in cfg["materials"]:
        for h in cfg["height_um"]:
            for ar_max in (10.0, 15.0):
                ta, tb, nra, nrb = [], [], None, None
                for fam in cfg["families"]:
                    z = load(od, fam, mat, cfg["pitch_um"][-1], h, (la, lb))
                    if z is None:
                        continue
                    ok = h / min_feature(fam, z[0]["params"]) <= ar_max + 1e-9
                    ta.append(z[0]["t"][ok])
                    tb.append(z[1]["t"][ok])
                    nra, nrb = 1 / float(z[0]["n_in"]), 1 / float(z[1]["n_in"])
                if not ta:
                    continue
                ta, tb = np.concatenate(ta), np.concatenate(tb)
                row = {"material": mat, "height_um": h, "max_aspect": ar_max, "n_atoms": int(len(ta))}
                for name, (pa, pb) in plates.items():
                    ea, eb, oa, ob = plate_efficiency(ta, tb, pa, pb, r, nra, nrb)
                    row[f"{name}_eta_852"], row[f"{name}_eta_1064"] = ea, eb
                row["doublet_lpa_852"] = row["corrector_eta_852"] * row["focusing_eta_852"]
                row["doublet_lpa_1064"] = row["corrector_eta_1064"] * row["focusing_eta_1064"]
                res.append(row)
                print("plates", {k: (round(v, 3) if isinstance(v, float) else v) for k, v in row.items()}, flush=True)
    write_json(out / "plate_efficiency.json", res)


def main(cfg_path):
    cfg = load_config(cfg_path)
    out = run_dir(cfg["run_id"])
    od = REPO_ROOT / cfg["output_dir"]
    la, lb = cfg["wavelengths_um"]
    rows = []
    best = None
    for mat in cfg["materials"]:
        for p in cfg["pitch_um"]:
            for h in cfg["height_um"]:
                fams = {}
                for fam in cfg["families"]:
                    z = load(od, fam, mat, p, h, (la, lb))
                    if z is not None:
                        fams[fam] = (z[0]["t"], z[1]["t"], h / min_feature(fam, z[0]["params"]), 1 / float(z[0]["n_in"]), 1 / float(z[1]["n_in"]),
                                     np.minimum(z[0]["T"], 1.0), z[0]["params"])
                if not fams:
                    continue
                sets = {k: [k] for k in fams}
                if len(fams) > 1:
                    sets["union"] = list(fams)
                for name, members in sets.items():
                    ta = np.concatenate([fams[m][0] for m in members])
                    tb = np.concatenate([fams[m][1] for m in members])
                    ar = np.concatenate([fams[m][2] for m in members])
                    nra, nrb = fams[members[0]][3], fams[members[0]][4]
                    for ar_max in (10.0, 15.0):
                        ok = ar <= ar_max + 1e-9
                        if ok.sum() < 3:
                            continue
                        ea, eb, sa, sb, _ = dual_efficiency(ta[ok], tb[ok], nra, nrb)
                        sa1 = dual_efficiency(ta[ok], tb[ok], nra, nrb, wa=1.0, wb=1e-6)[0]
                        sb1 = dual_efficiency(ta[ok], tb[ok], nra, nrb, wa=1e-6, wb=1.0)[1]
                        row = {"material": mat, "pitch_um": p, "height_um": h, "set": name, "max_aspect": ar_max, "n_atoms": int(ok.sum()),
                               "eta_852": ea, "eta_1064": eb, "rms_phase_852": sa, "rms_phase_1064": sb,
                               "single_eta_852": sa1, "single_eta_1064": sb1}
                        rows.append(row)
                        print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in row.items()}, flush=True)
                        score = min(ea, eb)
                        if best is None or score > best[0]:
                            best = (score, row, ta[ok], tb[ok], nra)
    write_json(out / "dual_efficiency.json", rows)
    plates_eval(cfg, od, out, la, lb)
    if best is not None:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        _, row, ta, tb, nra = best
        fig, ax = plt.subplots(figsize=(4.6, 4.2))
        sc = ax.scatter(np.mod(np.angle(ta), 2 * np.pi), np.mod(np.angle(tb), 2 * np.pi), c=np.minimum(np.abs(ta) ** 2 * nra, 1),
                        s=10, cmap="Blues", vmin=0, vmax=1, edgecolors="none")
        ax.set_xlabel("transmission phase at 852 nm (rad, mod 2 pi)")
        ax.set_ylabel("transmission phase at 1064 nm (rad, mod 2 pi)")
        ax.set_xlim(0, 2 * np.pi)
        ax.set_ylim(0, 2 * np.pi)
        ax.set_title(f"{row['material']}, p {row['pitch_um']} um, h {row['height_um']} um, {row['set']}, AR <= {row['max_aspect']:g}", fontsize=8)
        fig.colorbar(sc, label="zeroth-order transmittance at 852 nm")
        fig.tight_layout()
        fig.savefig(out / "phase_pairs_best.png", dpi=130)


if __name__ == "__main__":
    main(sys.argv[1])
