"""Lenslet-array relay at Fresnel number ~1: objective-pupil field of a square lenslet.

Usage: python scripts/million_lenslet_fresnel.py

The earlier pupil-fill model (scripts/million_pupil_fill.py) took the objective-pupil amplitude
to be the uniformly lit lenslet aperture (geometric mapping). That needs a Fresnel number
N = b^2 / (lambda f_L) >> 1 (b: lenslet half-width), whereas the relay has
N_S = a NA / (2 lambda) = 0.93 at a = 4 um, NA 0.5, 1075 nm. Here the pupil field is the angular
spectrum of the lenslet field, propagated to the lenslet focal plane (which the relay images onto
the atoms). Per axis, with the aperture coordinate u in units of b and the pupil coordinate sigma
in units of the stop radius (sigma = 1 is the objective NA):

    F(sigma) = exp(-i pi N sigma'^2) int_{-1}^{1} exp(-i pi (N u^2 + 2 N_S sigma u)) du,

where N = N_S (theta_L / theta_S) is the lenslet Fresnel number (theta_L = b / f_L, the
geometric edge angle; theta_S the stop radius in lenslet space; N = N_S for a lenslet
circumscribing the stop) and the leading quadratic phase (sigma' = N_S sigma / N) is the
propagation from the lenslet plane to its focal plane. The square-lenslet pupil field is
F(sx) F(sy), normalised to unit power; the circular stop is |sigma| <= 1. In atom space the
pupil is refocused by exp(i k z sqrt(1 - NA^2 sigma^2)), z optimised (the relay can be focused
on the best plane).

Single trap: peak per watt relative to a uniformly filled disc, |int_disc L e^{i phi}|^2 / pi.
Infinite lattice of pitch a: the pupil field is L sampled on the reciprocal lattice
sigma = (m + q) Delta (Delta = lambda / (a NA) = 1 / (2 N_S); q = 0 for equal pistons,
(1/2, 1/2) for alternating 0 / pi (checkerboard), (1/2, 0) for alternating rows), and the bulk
site field per unit lenslet power is Delta^2 sum_G L(G) e^{i phi(G)} over G in the stop. The
uniform-disc limit of this formula is (n_G / n_P)^2 (docs/million_sites.md, infinite-lattice
check). Finite 51 x 51 patches (FFT convolution of the 2D single-lenslet PSF, central 15 x 15)
check the lattice values and give the random-piston spread. Output:
results/million/lenslet_fresnel.json.
"""

import numpy as np

from metatrap.runio import REPO_ROOT, write_json

LAM, NA = 1.075, 0.5


def F1(sig, N, NS, nu=4001):
    """Per-axis pupil field (unnormalised) at the lenslet focal plane, sigma in stop units."""
    u = np.linspace(-1, 1, nu)
    du = u[1] - u[0]
    w = np.full(nu, du)
    w[0] = w[-1] = du / 2
    sig = np.atleast_1d(sig)
    out = np.empty(sig.shape, complex)
    for k, s in enumerate(sig.ravel()):
        out.ravel()[k] = np.sum(w * np.exp(-1j * np.pi * (N * u**2 + 2 * NS * s * u)))
    # Propagation to the lenslet focus: angular-spectrum defocus by f_L (paraxial in lenslet space).
    sp = NS * sig / N if N > 0 else 0 * sig
    return out * np.exp(-1j * np.pi * N * sp**2)


def norm_1d(N, NS):
    """Integral of |F1|^2 over all sigma (Parseval: (1 / (NS)) x int |g|^2 du = 2 / NS)."""
    return 2.0 / NS


def refocus(sig2, z):
    return np.exp(1j * 2 * np.pi / LAM * z * (np.sqrt(np.clip(1 - NA**2 * sig2, 0, None)) - 1))


def single(N, NS, z_grid, n=801):
    s = np.linspace(-1, 1, n)
    ds = s[1] - s[0]
    f = F1(s, N, NS) / np.sqrt(norm_1d(N, NS))
    SX, SY = np.meshgrid(s, s, indexing="ij")
    disc = SX**2 + SY**2 <= 1
    L = np.outer(f, f)
    trans = float(np.sum(np.abs(L[disc]) ** 2) * ds**2)
    pk = [abs(np.sum((L * refocus(SX**2 + SY**2, z))[disc]) * ds**2) ** 2 / np.pi for z in z_grid]
    i = int(np.argmax(pk))
    return {"stop_transmission": trans, "peak_per_W_focal_plane": float(pk[int(np.argmin(np.abs(z_grid)))]),
            "peak_per_W_best": float(pk[i]), "best_z_um": float(z_grid[i])}


def lattice(N, NS, q, z_grid):
    D = 1.0 / (2 * NS)
    m = np.arange(-int(2 / D) - 2, int(2 / D) + 3)
    gx, gy = (m + q[0]) * D, (m + q[1]) * D
    fx = F1(gx, N, NS) / np.sqrt(norm_1d(N, NS))
    fy = F1(gy, N, NS) / np.sqrt(norm_1d(N, NS))
    GX, GY = np.meshgrid(gx, gy, indexing="ij")
    L = np.outer(fx, fy)
    inside = GX**2 + GY**2 <= 1
    trans = float(np.sum(np.abs(L[inside]) ** 2) * D**2)
    pk = [abs(D**2 * np.sum((L * refocus(GX**2 + GY**2, z))[inside])) ** 2 / np.pi for z in z_grid]
    i = int(np.argmax(pk))
    # Intensity over one unit cell at the best plane: is the site the maximum, and how bright is
    # the brightest point away from it (|x| or |y| >= a / 4)?
    a = LAM / (NA * D)
    t = np.linspace(-0.5, 0.5, 81) * a
    TX, TY = np.meshgrid(t, t, indexing="ij")
    amp = (L * refocus(GX**2 + GY**2, z_grid[i]))[inside]
    kx, ky = 2 * np.pi * NA / LAM * GX[inside], 2 * np.pi * NA / LAM * GY[inside]
    E = D**2 * np.einsum("g,gxy->xy", amp, np.exp(1j * (kx[:, None, None] * TX + ky[:, None, None] * TY)))
    I = np.abs(E) ** 2 / np.pi
    away = np.maximum(np.abs(TX), np.abs(TY)) >= a / 4
    return {"n_orders_in_stop": int(inside.sum()), "stop_transmission": trans,
            "bulk_peak_per_W": float(pk[i]), "best_z_um": float(z_grid[i]), "bulk_peak_per_W_z0": float(pk[int(np.argmin(np.abs(z_grid)))]),
            "cell_max_over_site": float(I.max() / pk[i]), "brightest_away_from_site_over_site": float(I[away].max() / pk[i]),
            "cell_centre_over_site": float(I[0, 0] / pk[i]), "edge_midpoint_over_site": float(max(I[0, 40], I[40, 0]) / pk[i])}


def uniform_lattice(NS, q):
    D = 1.0 / (2 * NS)
    m = np.arange(-10, 11)
    GX, GY = np.meshgrid((m + q[0]) * D, (m + q[1]) * D, indexing="ij")
    nG = int(np.sum(GX**2 + GY**2 <= 1))
    return nG, (nG * D**2 / np.pi) ** 2


def psf_2d(N, NS, z, dx=0.1, n_sig=256, r_max=40.0):
    """Amplitude PSF of one lenslet (stop-limited Fresnel pupil, refocus z) on a grid of spacing dx."""
    from scipy import fft as sfft

    ds = 2.0 / n_sig
    m = n_sig // 2 + 1
    sig = (np.arange(2 * m + 1) - m) * ds  # sigma = 0 on a sample (no spurious tilt of the PSF)
    f = F1(sig, N, NS) / np.sqrt(norm_1d(N, NS))
    SX, SY = np.meshgrid(sig, sig, indexing="ij")
    P = np.where(SX**2 + SY**2 <= 1, np.outer(f, f) * refocus(SX**2 + SY**2, z), 0)
    nfft = int(round(LAM / (NA * dx) / ds))
    buf = np.zeros((nfft, nfft), complex)
    buf[: 2 * m + 1, : 2 * m + 1] = P
    buf = np.roll(buf, (-m, -m), axis=(0, 1))  # sigma = 0 at index 0
    U = sfft.fftshift(sfft.ifft2(buf)) * nfft**2 * ds**2  # sum over the pupil x dsigma^2
    c = nfft // 2
    k = int(r_max / dx)
    U = U[c - k : c + k + 1, c - k : c + k + 1]
    i0 = np.unravel_index(np.argmax(np.abs(U)), U.shape)
    return U, i0


def patch(U, i0, Nside, a, pist, dx=0.1, core=15, w_um=0.6):
    """Peaks of |sum_j psi_j U(x - x_j)|^2 at the sites of an Nside^2 lattice (FFT convolution)."""
    from scipy.signal import fftconvolve

    step = int(round(a / dx))
    lat = np.zeros(((Nside - 1) * step + 1,) * 2, complex)
    lat[::step, ::step] = np.exp(1j * pist.reshape(Nside, Nside))
    E = fftconvolve(lat, U, mode="full")
    off = np.array(i0)  # site j lands at its lattice index + i0 (the PSF maximum)
    w = int(round(w_um / dx))
    pk = np.empty((Nside, Nside))
    for p_ in range(Nside):
        for q in range(Nside):
            i, j = p_ * step + off[0], q * step + off[1]
            pk[p_, q] = np.max(np.abs(E[i - w : i + w + 1, j - w : j + w + 1]) ** 2)
    c = Nside // 2
    h = core // 2
    cpk = pk[c - h : c + h + 1, c - h : c + h + 1]
    single = np.abs(U[i0]) ** 2
    return {"core_nonuniformity_rms": float(cpk.std() / cpk.mean()), "core_min_over_mean": float(cpk.min() / cpk.mean()),
            "core_peak_per_W": float(cpk.mean() / np.pi), "core_over_single": float(cpk.mean() / single)}


def main():
    z_grid = np.linspace(-3, 3, 121)
    res = {"wavelength_um": LAM, "na": NA, "model": "square lenslet, Fourier relay, circular stop; scalar, aberration-free relay"}
    pist = {"equal": (0.0, 0.0), "alternating_checkerboard": (0.5, 0.5), "alternating_rows": (0.5, 0.0)}
    for a in (4.0, 3.0):
        NS = a * NA / (2 * LAM)
        r = {"N_S": NS, "reciprocal_spacing_stop_units": 1 / (2 * NS)}
        r["circumscribed"] = {"single": single(NS, NS, z_grid)}
        for k, q in pist.items():
            lat = lattice(NS, NS, q, z_grid)
            nG, unif = uniform_lattice(NS, q)
            lat["uniform_pupil_n_G"], lat["uniform_pupil_bulk_over_single"] = nG, unif
            lat["bulk_over_single_best"] = lat["bulk_peak_per_W"] / r["circumscribed"]["single"]["peak_per_W_best"]
            r["circumscribed"][k] = lat
        # Lenslet focal length scan: theta_L / theta_S from 0 (no lens phase) to 2.
        scan = []
        for ratio in np.round(np.arange(0.0, 2.01, 0.1), 2):
            N = NS * ratio
            row = {"thetaL_over_thetaS": float(ratio), "N": float(N), "single_peak_per_W_best": single(N, NS, z_grid, n=401)["peak_per_W_best"]}
            for k, q in pist.items():
                lat = lattice(N, NS, q, z_grid)
                row[f"{k}_bulk_peak_per_W"], row[f"{k}_stop_transmission"] = lat["bulk_peak_per_W"], lat["stop_transmission"]
            scan.append(row)
        r["focal_length_scan"] = scan
        for k in pist:
            best = max(scan, key=lambda x: x[f"{k}_bulk_peak_per_W"])
            r[f"best_{k}"] = {"thetaL_over_thetaS": best["thetaL_over_thetaS"], "bulk_peak_per_W": best[f"{k}_bulk_peak_per_W"]}
        # Finite patches with the Fresnel pupil (51 x 51, central 15 x 15), at the single-trap
        # best plane: periodic pistons against the lattice values, and random pistons.
        z1 = r["circumscribed"]["single"]["best_z_um"]
        U, i0 = psf_2d(NS, NS, z1)
        r["single_peak_per_W_psf_check"] = float(np.abs(U[i0]) ** 2 / np.pi)
        Ns, ii, jj = 51, *np.meshgrid(np.arange(51), np.arange(51), indexing="ij")
        pp = {"equal": np.zeros(Ns * Ns), "alternating_checkerboard": (np.pi * ((ii + jj) % 2)).ravel(),
              "alternating_rows": (np.pi * (ii % 2)).ravel()}
        r["patch51"] = {k: patch(U, i0, Ns, a, v) for k, v in pp.items()}
        rnd = [patch(U, i0, Ns, a, np.random.default_rng(200 + sd).uniform(0, 2 * np.pi, Ns * Ns)) for sd in range(3)]
        r["patch51"]["random_3_seeds"] = {k: float(np.mean([x[k] for x in rnd])) for k in rnd[0]}
        r["patch51"]["plane_z_um"] = z1
        print(f"  patch51 at z = {z1}: {r['patch51']}")
        res[f"a_{a:g}um"] = r
        c = r["circumscribed"]
        print(f"a = {a} um, N_S = {NS:.3f}: single {c['single']}")
        for k in pist:
            print(f"  {k}: {c[k]}; best over f_L: {r['best_' + k]}")
    write_json(REPO_ROOT / "results" / "million" / "lenslet_fresnel.json", res)


if __name__ == "__main__":
    main()
