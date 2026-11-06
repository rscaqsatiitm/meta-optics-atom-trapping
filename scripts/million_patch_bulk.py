"""Bulk interference value for a large lenslet-array trap patch (shift-invariant model).

Usage: python scripts/million_patch_bulk.py configs/million.yaml

scripts/million_patch.py propagates every trap of a 7 x 7 patch exactly. That patch is too small
to separate the boundary from the bulk. Near the axis the 25 mm meta-objective is shift-invariant:
each trap field is a translated copy of the on-axis amplitude PSF U0(r), which is rotationally
symmetric (axisymmetric doublet, circular stop, uniform fill). U0(r) is computed once on a radial
line with the same superposition method (FocalSum, exact stack paths) and the patch field
E(x) = sum_j exp(i psi_j) U0(|x - x_j|) is formed by FFT convolution on a 0.1 um grid.

Steps (results/million/patch_bulk.json):
(i) validation against the exact 7 x 7 results of patch.json (same pistons and peak search);
(ii) N x N patches, N up to 51: depth non-uniformity of the central 7 x 7, depth versus distance
     from the patch edge, bulk depth relative to a single trap, brightest interstitial point in
     the bulk; random pistons over several seeds;
(iii) shift invariance: exact trap fields at 0, 50 and 100 um from the axis compared with U0
     (peak position, complex overlap, and the residual piston that would change the relative
     phase of neighbours).
"""

import sys

import numpy as np
import scipy.fft as sfft

from metatrap.metalens import doublet as D
from metatrap.metalens import wave_offaxis as WO
from metatrap.propagation.layers import LayerStack
from metatrap.runio import REPO_ROOT, load_config, write_json


def build(pc):
    lam, na = pc["wavelength_um"], pc["na"]
    cell = [[pc["air_gap_um"], "air"], [4000, "fused_silica"], [9000, "air"]]
    f0 = LayerStack.from_list(cell).aperture_for_na(na, lam) / na
    d = D.make(D.DoubletSpec(na, cell, [[f0, "air"]], n_pts=2001, ext=min(1.6, 0.95 / na)), lam)
    s = D.system(d, lam)
    dz, _ = D.common_plane(s, d.info["R1_um"], lam, list(np.linspace(0, 2828, 6)), d.info["f_um"], n_rays=28)
    fs = WO.FocalSum(d.spec.cell, lam, d.info["R2_um"] * 1.6 + 3000)
    return d, s, dz, fs


def radial_psf(d, fs, lam, dz, n_pupil, r_max, dr=0.01):
    pf = WO.plate_field(d, lam, 0.0, n_pupil)
    r = np.arange(0, r_max + dr, dr)
    u = fs.field(pf, pf["amp"], np.column_stack([r, np.zeros_like(r)]), dz)
    return r, u / u[0]


def lattice_intensity(r, u, N, a, pistons, h, r_t):
    """|sum_j psi_j U0(|x - x_j|)|^2 on a grid of spacing h; sites at multiples of a (a/h integer)."""
    step = int(round(a / h))
    half = int(np.ceil(r_t / h))
    m = (N - 1) * step + 2 * half + 1
    L = sfft.next_fast_len(m + 2 * half)
    lat = np.zeros((L, L), complex)
    idx = half + np.arange(N) * step
    lat[np.ix_(idx, idx)] = np.exp(1j * pistons.reshape(N, N))
    k = np.arange(-half, half + 1) * h
    R = np.hypot(k[:, None], k[None, :])
    ker = np.where(R <= r_t, np.interp(R, r, u.real) + 1j * np.interp(R, r, u.imag), 0)
    # Kernel centred at index 0 (circular convolution; L leaves room for the full support).
    K = np.zeros((L, L), complex)
    K[: half + 1, : half + 1] = ker[half:, half:]
    K[: half + 1, -half:] = ker[half:, :half]
    K[-half:, : half + 1] = ker[:half, half:]
    K[-half:, -half:] = ker[:half, :half]
    E = sfft.ifft2(sfft.fft2(lat) * sfft.fft2(K))
    return np.abs(E) ** 2, idx, step


def site_peaks(I, idx, h, search_um=0.6):
    w = int(round(search_um / h))
    pk = np.empty((len(idx), len(idx)))
    for a_, i in enumerate(idx):
        for b_, j in enumerate(idx):
            pk[a_, b_] = I[i - w : i + w + 1, j - w : j + w + 1].max()
    return pk


def interstitial(I, idx, step, w):
    c = idx[:-1] + step // 2
    return np.array([[I[i - w : i + w + 1, j - w : j + w + 1].max() for j in c] for i in c])


def pistons_for(kind, N, rng=None):
    if kind == "equal":
        return np.zeros(N * N)
    if kind == "alternating":
        i, j = np.meshgrid(np.arange(N), np.arange(N), indexing="ij")
        return (np.pi * ((i + j) % 2)).ravel()
    return rng.uniform(0, 2 * np.pi, N * N)


def main(cfg_path):
    cfg = load_config(cfg_path)
    pc = cfg["patch"]
    lam, a = pc["wavelength_um"], pc["pitch_um"]
    d, s, dz, fs = build(pc)
    r_t = 60.0
    r, u = radial_psf(d, fs, lam, dz, pc["n_pupil"], r_t + 5)
    res = {"method": "shift-invariant convolution of the exact on-axis PSF", "kernel_radius_um": r_t, "pitch_um": a,
           "psf_first_zero_um": float(r[np.argmax(np.abs(u) < 0.02 * 1)])}

    # (i) validation against the exact 7 x 7 (0.2 um grid, grid-maximum peaks as in million_patch.py)
    rng = np.random.default_rng(5)
    rnd = [rng.uniform(0, 2 * np.pi, 49), rng.uniform(0, 2 * np.pi, 49)]
    val = {}
    for name, ps in (("equal", pistons_for("equal", 7)), ("alternating", pistons_for("alternating", 7)), ("random_draw1", rnd[0]),
                     ("random_draw2", rnd[1])):
        I, idx, step = lattice_intensity(r, u, 7, a, ps, 0.2, r_t)
        pk = site_peaks(I, idx, 0.2)[1:-1, 1:-1]
        val[name] = float(pk.std() / pk.mean())
    res["validation_7x7"] = val
    print("validation 7x7 (exact: 0.0519, 0.0400, 0.1106, 0.1369):", val, flush=True)

    # Single-trap peak on the 0.1 um grid (normalisation of the bulk depth).
    I1, idx1, _ = lattice_intensity(r, u, 1, a, np.zeros(1), 0.1, r_t)
    p1 = site_peaks(I1, idx1, 0.1)[0, 0]

    # (ii) large patches
    h = 0.1
    big = {}
    for N in (7, 15, 31, 51):
        for kind in ("equal", "alternating"):
            I, idx, step = lattice_intensity(r, u, N, a, pistons_for(kind, N), h, r_t)
            pk = site_peaks(I, idx, h)
            c = N // 2
            core = pk[c - 3 : c + 4, c - 3 : c + 4]
            ring = [float(np.mean(np.concatenate([pk[e, e:N - e], pk[N - 1 - e, e:N - e], pk[e:N - e, e], pk[e:N - e, N - 1 - e]])) / core.mean())
                    for e in range(min(8, c + 1))]
            it = interstitial(I, idx, step, int(round(0.8 / h)))
            ic = it[c - 3 : c + 3, c - 3 : c + 3] if N >= 7 else it
            big[f"{kind}:{N}"] = {"core7_nonuniformity_rms": float(core.std() / core.mean()), "core_depth_over_single": float(core.mean() / p1),
                                  "inner_nonuniformity_rms": float(pk[1:-1, 1:-1].std() / pk[1:-1, 1:-1].mean()),
                                  "ring_depth_over_core": ring, "core_interstitial_max_over_peak": float(ic.max() / core.mean())}
            print(f"{kind}:{N}", {k: (np.round(v, 4).tolist() if isinstance(v, list) else round(v, 4)) for k, v in big[f'{kind}:{N}'].items()}, flush=True)
        # random pistons: bulk statistics over seeds
        vals = []
        for seed in range(4 if N <= 31 else 2):
            I, idx, step = lattice_intensity(r, u, N, a, pistons_for("random", N, np.random.default_rng(100 + seed)), h, r_t)
            pk = site_peaks(I, idx, h)
            c = N // 2
            core = pk[max(0, c - 7) : c + 8, max(0, c - 7) : c + 8]
            vals.append((core.std() / core.mean(), core.mean() / p1, core.min() / core.mean()))
        v = np.array(vals)
        big[f"random:{N}"] = {"core_nonuniformity_rms": float(v[:, 0].mean()), "core_depth_over_single": float(v[:, 1].mean()),
                              "core_min_over_mean": float(v[:, 2].mean()), "n_seeds": len(vals)}
        print(f"random:{N}", {k: round(v_, 4) for k, v_ in big[f'random:{N}'].items()}, flush=True)
    # Kernel truncation check (bulk of the 51 x 51, equal pistons).
    trunc = {}
    for rt in (30.0, 45.0, 60.0):
        I, idx, step = lattice_intensity(r, u, 51, a, np.zeros(51 * 51), h, rt)
        pk = site_peaks(I, idx, h)
        trunc[str(rt)] = float(pk[22:29, 22:29].mean() / p1)
    big["truncation_core_depth_over_single_equal_51"] = trunc
    print("truncation", trunc, flush=True)
    res["patches"] = big

    # (iii) shift invariance: exact fields of traps at 0, 50, 100 um against U0
    f, R1 = d.info["f_um"], d.info["R1_um"]
    g = np.arange(-4, 4.0001, 0.1)
    GX, GY = np.meshgrid(g, g, indexing="ij")
    ref = np.interp(np.hypot(GX, GY), r, u.real) + 1j * np.interp(np.hypot(GX, GY), r, u.imag)
    inv = {}
    for rho in (0.0, 50.0, 100.0):
        th = D.theta_for_height(s, R1, lam, rho, f) if rho > 0 else 0.0
        x_img = D.quality(s, R1, lam, th, None, dz)["x_um"] if rho > 0 else 0.0
        pf = WO.plate_field(d, lam, th, pc["n_pupil"])
        U = fs.field(pf, pf["amp"], np.column_stack([(x_img + GX).ravel(), GY.ravel()]), dz).reshape(GX.shape)
        ov = np.vdot(ref, U)
        inv[str(rho)] = {"x_img_um": float(x_img), "overlap_abs": float(abs(ov) / np.sqrt(np.vdot(ref, ref).real * np.vdot(U, U).real)),
                         "piston_rad": float(np.angle(ov))}
        print("shift invariance", rho, inv[str(rho)], flush=True)
    p = [inv[k]["piston_rad"] for k in ("0.0", "50.0", "100.0")]
    p = np.unwrap(np.array(p) - p[0])
    # Quadratic piston phi = alpha rho^2 -> neighbour phase difference 2 alpha rho a at rho.
    alpha = p[2] / 100.0**2
    inv["quadratic_piston_rad_per_um2"] = float(alpha)
    inv["neighbour_phase_step_at_100um_rad"] = float(2 * alpha * 100 * a)
    res["shift_invariance"] = inv
    write_json(REPO_ROOT / "results" / "million" / "patch_bulk.json", res)


if __name__ == "__main__":
    main(sys.argv[1])
