"""2D (horizontal x vertical) Fox-Li check of the in-plane telescope ring: separability of the
design, effect of non-separable figure errors, and what a correction in x alone can recover.
Usage:
    python scripts/cavity_2d_check.py configs/cavity_2d_check.yaml

The 1D strip model (cavity_foxli.py) assumes a separable mode. Here the full transverse field
E(z, x) is propagated by the 2D angular spectrum. Mirrors: cylindrical telescope in x, vertical
focal length f_z on the primary (S4 design), graded mirrors designed in 2D by conjugation through
the actual half trips. Figure errors: isotropic smooth 2D maps on every mirror (correlation length
referred to the crossing width and scaled by the horizontal beam size on each mirror, the same
length vertically), so on the small mirrors they vary across the vertical beam (non-separable).
Corrections: (i) 2D phase on G from perfect metrology; (ii) the same projected onto x only
(intensity-weighted average over z), i.e. a cylindrical corrector; (iii) closed loop from the
atom-plane depth map with x-only target updates, nominal model. Writes check2d.json.
"""

import sys
import time

import numpy as np

from metatrap.cavity import budget as B
from metatrap.cavity import foxli as FL
from metatrap.cavity import lattice as LA
from metatrap.runio import REPO_ROOT, load_config, run_dir, write_json, write_params

sys.path.insert(0, str(REPO_ROOT / "scripts"))
import cavity_correction as C  # noqa: E402


class Grid2D:
    """Transverse grid (z rows, x columns); .x is the horizontal coordinate array used for apertures."""

    def __init__(self, nx, width_x, nz, width_z, lam):
        self.nx, self.nz, self.lam_um = nx, nz, lam
        self.dx, self.dz = width_x / nx, width_z / nz
        self.x1 = (np.arange(nx) - nx // 2) * self.dx
        self.z1 = (np.arange(nz) - nz // 2) * self.dz
        self.Z, self.X = np.meshgrid(self.z1, self.x1, indexing="ij")
        self.x = self.X
        self.n = nx * nz
        self._h = {}

    def transfer(self, z):
        key = round(float(z), 6)
        if key not in self._h:
            k = 2 * np.pi / self.lam_um
            kx = 2 * np.pi * np.fft.fftfreq(self.nx, d=self.dx)
            kz = 2 * np.pi * np.fft.fftfreq(self.nz, d=self.dz)
            k2 = k**2 - kz[:, None] ** 2 - kx[None, :] ** 2
            self._h[key] = np.where(k2 > 0, np.exp(1j * (np.sqrt(np.clip(k2, 0, None)) - k) * z), 0.0)
        return self._h[key]

    def prop(self, u, z):
        if z == 0:
            return u
        return np.fft.ifft2(np.fft.fft2(u) * self.transfer(z))


class Flat:
    """Adapter so that the 1D lock/solve routines see a flat vector space."""

    def __init__(self, ring):
        self.r, self.shape = ring, ring.grid.X.shape
        self.grid = type("G", (), {"n": ring.grid.n})()

    def roundtrip(self, v):
        return self.r.roundtrip(np.asarray(v).reshape(self.shape)).ravel()


def design_2d(ring, t_in, t_out, which, floor=1e-4):
    a = ring.to_mirror(t_in, which)
    b = ring.from_mirror_adjoint(t_out, which)
    ph = np.angle(b * np.conj(a))
    return np.where(np.abs(a) >= floor * np.abs(a).max(), ph, 0.0)


def screen_2d(g, rms, corr, rng, hw_x, hw_z):
    white = rng.standard_normal((g.nz, g.nx))
    kx = 2 * np.pi * np.fft.fftfreq(g.nx, d=g.dx)
    kz = 2 * np.pi * np.fft.fftfreq(g.nz, d=g.dz)
    s = np.real(np.fft.ifft2(np.fft.fft2(white) * np.exp(-((kz[:, None] * corr) ** 2 + (kx[None, :] * corr) ** 2) / 8.0)))
    sel = (np.abs(g.X) <= hw_x) & (np.abs(g.Z) <= hw_z)
    s = s - s[sel].mean()
    return s * rms / s[sel].std()


def beam_sizes(ring, F):
    """(2 x rms half-width in x, in z) of the unperturbed beam on each stop of half trip 1 and on G."""
    g = ring.grid

    def size(u):
        i = np.abs(u) ** 2
        return 2 * np.sqrt(np.sum(g.X**2 * i) / i.sum()), 2 * np.sqrt(np.sum(g.Z**2 * i) / i.sum())

    eG = ring.to_mirror(F, 1)
    after = eG * ring._ap_g * np.exp(1j * ring.phase_g1)
    out = [size(ring._segment(F, 0.0, z - 1e-6, 1) if z < ring.p_um else ring._segment(after, ring.p_um, z - 1e-6, 1)) for z in ring.stop_positions()]
    return out, size(eG)


def lock_recentre(fr, e_in, t_in, loss, psi, iters=6, x0=None):
    """Lock to the round-trip phase of the dominant circulating mode by re-centring only (each step a
    warm-started driven solve); cheaper than the golden-section refinement used in 1D, and equal to
    it within a small fraction of a linewidth when one mode dominates."""
    r = np.sqrt((1 - t_in) * (1 - loss))
    lw = 2 * (1 - r) / np.sqrt(r)
    e = x0
    for _ in range(iters):
        e, _ = FL.driven(fr, e_in, t_in, loss, psi, x0=e, tol=1e-7)
        lam_e = np.vdot(e, fr.roundtrip(e)) / np.vdot(e, e)
        new = psi + np.angle(np.exp(1j * (-np.angle(lam_e) - psi)))
        done = abs(new - psi) < 0.02 * lw
        psi = new
        if done:
            break
    e, _ = FL.driven(fr, e_in, t_in, loss, psi, x0=e, tol=1e-9)
    return {"field": e, "enhancement": np.vdot(e, e).real / np.vdot(e_in, e_in).real}


def evaluate(cfg, ctx, ring):
    t_in = B.coupler_for_enhancement(cfg["g_target"], cfg["loss_other"])
    F = ctx["F"]
    fr = Flat(ring)
    lam0 = np.vdot(F.ravel(), fr.roundtrip(F.ravel())) / np.vdot(F.ravel(), F.ravel())
    res = lock_recentre(fr, F.ravel(), t_in, cfg["loss_other"], -np.angle(lam0))
    e = res["field"].reshape(F.shape)
    eb = ring.half(e, 1)
    g = ring.grid
    planes = LA.map_planes(cfg["array_side_um"], ctx["theta"])
    iz = g.nz // 2
    st_a = np.array([g.prop(e, s)[iz] for s in planes])
    st_b = np.array([g.prop(eb, s)[iz] for s in planes])
    m = LA.site_map_from_stacks(st_a, st_b, g.x1, planes, ctx["lam"], ctx["theta"], cfg["array_side_um"], cfg["map_sites_per_side"], flip_b=True)
    met = LA.map_metrics(m)
    return {"depth_rms": met["depth_rms"], "min_over_mean": met["depth_min_over_mean"], "dy_affine_res_nm": met["dy_affine_residual_rms_um"] * 1e3,
            "g_ratio": res["enhancement"] / float(B.enhancement(t_in, cfg["loss_other"]))}, m


def main(cfg_path):
    cfg = load_config(cfg_path)
    out = run_dir(cfg["run_id"])
    write_params(out, cfg)
    lam = cfg["wavelength_um"]
    k = 2 * np.pi / lam
    theta = B.crossing_half_angle(cfg["period_um"], lam)
    h = LA.array_halfrange(cfg["array_side_um"], theta)
    w = LA.supergauss_w_for_flat(h, cfg["order"], cfg["flat_tol"])
    gd = cfg["grid"]
    g = Grid2D(gd["nx"], gd["width_x_um"], gd["nz"], gd["width_z_um"], lam)
    t = cfg["telescope"]
    f2 = t["mag"] * t["f1_um"]
    px1 = -k * g.X**2 / (2 * t["f1_um"])
    px2 = -k * g.X**2 / (2 * f2)
    pz1 = -k * g.Z**2 / (2 * cfg["f_z_um"])
    els = [(t["d0_um"], t["ap1_um"], px1 + pz1), (t["d0_um"] + t["f1_um"] + f2, t["ap2_um"], px2)]
    p = t["d0_um"] + t["f1_um"] + f2 + t["L_um"]
    zero = np.zeros(g.X.shape)
    ring0 = FL.RingStrip(g, p, 0.0, zero, zero, ap_g_um=t["ap_g_um"], elements=els)
    F = (LA.supergauss(g.X, w, cfg["order"]) * np.exp(-(g.Z / cfg["w_z_um"]) ** 2)).astype(complex)
    ph0 = design_2d(ring0, F, F, 1)
    ring0.phase_g1 = ring0.phase_g2 = ph0
    ctx = {"F": F, "lam": lam, "theta": theta}
    rep = {"grid": gd, "cases": {}}

    def run(name, ring):
        t0 = time.time()
        ev, m = evaluate(cfg, ctx, ring)
        ev["seconds"] = time.time() - t0
        rep["cases"][name] = ev
        write_json(out / "check2d.json", rep)
        print(name, ev, flush=True)
        return ev, m

    run("ideal_2d", ring0)
    sizes, (gx, gz) = beam_sizes(ring0, F)
    rep["beam_sizes_um"] = {"stops": sizes, "G": [gx, gz]}
    W_x = 2 * np.sqrt(np.sum(g.X**2 * np.abs(F) ** 2) / np.sum(np.abs(F) ** 2))
    for seed in cfg["seeds"]:
        rng = np.random.default_rng(seed)
        rms = 4 * np.pi * cfg["fig_waves"]
        changes = [(wh, i, screen_2d(g, rms, cfg["corr_um"] * sx / W_x, rng, sx, sz)) for wh in (1, 2) for i, (sx, sz) in enumerate(sizes)]
        g_err = {wh: screen_2d(g, rms, cfg["corr_um"] * gx / W_x, rng, gx, gz) for wh in (1, 2)}
        true = ring0.perturbed(changes)

        def with_g(p1, p2):
            r = true.perturbed()
            r.phase_g1, r.phase_g2 = p1 + g_err[1], p2 + g_err[2]
            return r

        run(f"seed{seed}_uncorrected", with_g(ph0, ph0))
        # (i) perfect metrology, 2D phase on G (fabricated surface absorbs the known G error).
        q1, q2 = design_2d(true, F, F, 1), design_2d(true, F, F, 2)
        run(f"seed{seed}_openloop_2d", with_g(q1 - g_err[1], q2 - g_err[2]))
        # (ii) the same correction projected onto x only (cylindrical corrector).
        iG = np.abs(ring0.to_mirror(F, 1)) ** 2
        proj = []
        for q, e in ((q1, g_err[1]), (q2, g_err[2])):
            dlt = np.angle(np.exp(1j * (q - e - ph0)))
            dx_only = np.sum(dlt * iG, axis=0) / np.maximum(np.sum(iG, axis=0), 1e-30)
            proj.append(ph0 + dx_only[None, :])
        run(f"seed{seed}_openloop_x_only", with_g(*proj))
        # (iii) closed loop from the atom-plane depth map, x-only target updates, nominal model.
        TA, TB = F.copy(), F.copy()
        p1, p2 = ph0, ph0
        gauss = np.exp(-(g.Z / cfg["w_z_um"]) ** 2)
        for rnd in range(cfg["closed_loop_rounds"] + 1):
            ev, m = run(f"seed{seed}_closedloop_round{rnd}", with_g(p1, p2))
            if rnd == cfg["closed_loop_rounds"]:
                break
            knots, al, be = C.fit_profiles(m, theta, h * 1.02, 48, 1.0)
            a = C.profile_on_grid(g.x1, knots, al, h * 1.02)
            b = C.profile_on_grid(g.x1, knots, be, h * 1.02)
            TA = TA * (a[None, :] ** -0.5)
            TB = TB * (b[None, :] ** -0.5)
            TA = np.abs(TA[g.nz // 2])[None, :] * gauss
            TB = np.abs(TB[g.nz // 2])[None, :] * gauss
            p1, p2 = design_2d(ring0, TA.astype(complex), TB.astype(complex), 1), design_2d(ring0, TB.astype(complex), TA.astype(complex), 2)
    write_json(out / "check2d.json", rep)


if __name__ == "__main__":
    main(sys.argv[1])
