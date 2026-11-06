"""Memory-bounded render grid and tiled focal metrics against the untiled versions."""

import importlib.util
import types
from pathlib import Path

import numpy as np

from metatrap.hologram import symmetric as S
from metatrap.hologram import targets
from metatrap.hologram.exact import OPLTable
from metatrap.hologram.lazy import LazyKG, render_grid_symmetric_lazy
from metatrap.hologram.tiled import focal_metrics_tiled, fourier_interp_2d, tiled_render_grid
from metatrap.hologram.wgs import PupilGrid
from metatrap.propagation.layers import LayerStack

_spec = importlib.util.spec_from_file_location("fa", Path(__file__).parents[1] / "scripts" / "fullaperture.py")
FA = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(FA)


def _setup(nx=8):
    stack = LayerStack.from_list([[2500, "air"], [4000, "fused_silica"], [9000, "air"]])
    sst = LayerStack.from_list([[2500, "air"], [4000, "fused_silica"], [9050, "air"]])
    R = stack.aperture_for_na(0.55, 1.075)
    sites = targets.make(kind="square", nx=nx, ny=nx, pitch_um=4.0)
    fs, fstab, _ = S.wedge_sites(sites)
    rmax = np.abs(sites).max() * np.sqrt(2)
    g = PupilGrid.for_sites(1.075, 4.0, 4, 0.55, 2 * rmax, guard=1.25)
    tab = OPLTable(sst, 1.075, R + 2 * rmax + 20)
    rng = np.random.default_rng(3)
    coeffs = np.exp(2j * np.pi * rng.random(len(fs))).astype(np.complex64)
    return stack, g, fs, fstab, tab, coeffs


def test_fourier_interp_alignment():
    x = np.exp(2j * np.pi * (3 * np.arange(32)[:, None] + 5 * np.arange(32)[None, :]) / 32).astype(np.complex64)
    y = fourier_interp_2d(x, 4)
    assert np.allclose(y[::4, ::4], x, atol=1e-4)
    k = np.arange(128) / 4
    ref = np.exp(2j * np.pi * (3 * k[:, None] + 5 * k[None, :]) / 32)
    assert np.abs(y - ref).max() < 1e-3


def _exact_unit_field(stack, fs, fstab, tab, coeffs, sx, sy, lam=1.075):
    """Exact pupil field (lazy kernel sum) at arbitrary direction cosines, unit modulus."""
    s = np.hypot(sx, sy)
    rho = stack.ray_height(s, lam)
    ux = np.where(s > 0, sx / np.maximum(s, 1e-30), 0.0)
    uy = np.where(s > 0, sy / np.maximum(s, 1e-30), 0.0)
    wp = types.SimpleNamespace(r=np.column_stack([-ux * rho, -uy * rho]), opl0=stack.opl(s, lam))
    q = LazyKG(wp, fs, tab, lam, float(rho.max()) + 1.5 * float(np.abs(fs).max()) + 50) @ (coeffs / fstab)
    return q / np.abs(q), np.abs(q) ** 2


def test_tiled_render_against_exact_field():
    # 16 x 16 sites in the focal plane (no array defocus), so that the coarse pupil grid
    # samples the field as in the production designs (array frequencies up to 1/guard of Nyquist).
    stack, g, fs, fstab, _, coeffs = _setup(16)
    tab = OPLTable(stack, 1.075, stack.aperture_for_na(0.55, 1.075) + 200)
    ill = {"kind": "top_hat"}
    U = 5
    ref = render_grid_symmetric_lazy(g, stack, 0.55, ill, fs, fstab, coeffs, tab, 1.075, upsample=U)
    til = tiled_render_grid(g, stack, 0.55, ill, fs, fstab, coeffs, tab, 1.075, upsample=U)
    assert til.c0 == ref.c0 and np.isclose(til.ds, ref.ds)
    rng = np.random.default_rng(0)
    m = ref.field.shape[0]
    fi = rng.integers(0, m, 20000)
    fj = rng.integers(0, m, 20000)
    sx, sy = (fi - ref.c0) * ref.ds, (fj - ref.c0) * ref.ds
    ok = np.hypot(sx, sy) <= 0.55
    truth, w = _exact_unit_field(stack, fs, fstab, tab, coeffs, sx[ok], sy[ok])
    a = ref.field[fi[ok], fj[ok]]  # reference: FFT-upsampled grid, exact at its nodes
    b = til.sample(fi[ok].astype(float), fj[ok].astype(float))
    ea = np.abs(np.angle(a * np.conj(truth)))
    eb = np.abs(np.angle(b * np.conj(truth)))
    rms_a = np.sqrt(np.sum(w * ea**2) / w.sum())
    rms_b = np.sqrt(np.sum(w * eb**2) / w.sum())
    assert rms_b < 0.02 and rms_b <= rms_a  # 0.013 and 0.032 rad when written
    # Off the fine nodes (where the reference would interpolate bilinearly) as well.
    fi2, fj2 = fi[ok] + rng.random(ok.sum()), fj[ok] + rng.random(ok.sum())
    truth2, w2 = _exact_unit_field(stack, fs, fstab, tab, coeffs, (fi2 - ref.c0) * ref.ds, (fj2 - ref.c0) * ref.ds)
    e2 = np.abs(np.angle(til.sample(fi2, fj2) * np.conj(truth2)))
    assert np.sqrt(np.sum(w2 * e2**2) / w2.sum()) < 0.02


def test_tiled_focal_metrics_match():
    nf, dx, pitch, lam, na = 256, 0.5, 4.0, 1.075, 0.55
    sites = targets.make(kind="square", nx=10, ny=10, pitch_um=pitch)
    rng = np.random.default_rng(1)
    c = (np.arange(nf) - nf // 2) * dx
    X, Y = np.meshgrid(c, c, indexing="ij")
    F = np.zeros((nf, nf), np.complex128)
    for sx, sy in sites:
        sx2, sy2 = sx + 0.1 * rng.normal(), sy + 0.1 * rng.normal()
        F += (1 + 0.1 * rng.normal()) * np.exp(2j * np.pi * rng.random()) * np.exp(-((X - sx2) ** 2 + (Y - sy2) ** 2) / (2 * 0.6**2))
    F += 0.05 * (rng.normal(size=F.shape) + 1j * rng.normal(size=F.shape))
    # Band-limit to the pupil (|k| <= NA / lambda), as for a focal field.
    k = np.fft.fftfreq(nf, dx)
    F = np.fft.ifft2(np.fft.fft2(F) * (np.hypot(k[:, None], k[None, :]) <= na / lam))
    F = np.fft.ifftshift(F)  # focal_metrics applies fftshift after the (identity) transfer function
    cfg = {"site_pitch_um": pitch, "upsample": 4}
    H = np.ones((nf, nf))
    m_ref, per_ref, _ = FA.focal_metrics(cfg, F, dx, lam, na, sites, 1.0, H=H)
    m_til, per_til, _ = focal_metrics_tiled(cfg, F, dx, lam, na, sites, 1.0, H=H, tile_px=16)
    assert np.allclose(per_til["peak"], per_ref["peak"], rtol=2e-3)
    assert np.allclose(per_til["x"], per_ref["x"], atol=5e-3)  # um
    assert np.allclose(per_til["p_cell"], per_ref["p_cell"], rtol=2e-3)
    for k in ("depth_nonuniformity_rms", "array_efficiency", "core_efficiency", "waist_mean_um", "axis_intensity_over_mean_site_peak"):
        assert np.isclose(m_til[k], m_ref[k], rtol=5e-3, atol=1e-4), k
    assert np.isclose(m_til["interstitial_max"], m_ref["interstitial_max"], rtol=2e-2, atol=1e-3)
