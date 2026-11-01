"""Anchors for the air-spaced aplanatic doublet, the frame changes and the QDHT."""

import numpy as np

from metatrap.metalens import doublet as D
from metatrap.metalens.raytrace import Frame, Medium, System, point_source_rays
from metatrap.metalens.wave_axial import QDHT

SPEC = D.DoubletSpec(0.5, [[2500, "air"], [4000, "fused_silica"], [9000, "air"]], [[16000, "air"]], n_pts=2001)


def test_frame_identity_and_roundtrip():
    rng = np.random.default_rng(1)
    x, y = rng.uniform(-100, 100, (2, 50))
    px, py = rng.uniform(-0.3, 0.3, (2, 50))
    base = System([Medium(1.0, 1000.0)]).trace(x, y, px, py, 0.85)
    for fr in ((0.01, 0.0), (0.0, 3.0)):  # pure tilt and pure decentre are exact round trips
        tilt = System([Frame(fr[0], fr[1], 0.0), Frame(-fr[0], -fr[1], 0.0), Medium(1.0, 1000.0)]).trace(x, y, px, py, 0.85)
        for a, b in zip(base[:5], tilt[:5]):
            assert np.allclose(a, b, atol=1e-9)


def test_point_source_collimated_limit():
    x = np.linspace(-1000, 1000, 11)
    y = np.zeros_like(x)
    px, py, _ = point_source_rays(x, y, 0.02, 1e12)
    assert np.allclose(px, np.sin(0.02), atol=1e-8) and np.allclose(py, 0)


def test_doublet_on_axis_and_field():
    d = D.make(SPEC, 0.8523)
    s = D.system(d, 0.8523)
    q0 = D.quality(s, d.info["R1_um"], 0.8523, 0.0)
    assert q0["strehl"] > 0.9999 and abs(q0["dz_um"]) < 0.01
    dz, q = D.common_plane(s, d.info["R1_um"], 0.8523, [0, 400], d.info["f_um"])
    assert min(r["strehl"] for r in q) > 0.995


def test_wall_shift_follows_focus():
    d = D.make(SPEC, 0.8523)
    s = D.system(d, 0.8523, D.Perturbation(wall_delta_um=50))
    q = D.quality(s, d.info["R1_um"], 0.8523, 0.0)
    # Relative to the extended evaluation plane: -dt cos t / sqrt(n^2 - sin^2 t) at the rms ray.
    assert -36 < q["dz_um"] < -30


def test_qdht_gaussian():
    q = QDHT(50.0, 2.0)
    w = 5.0
    F = q.forward(np.exp(-q.r**2 / w**2))
    assert np.allclose(F, np.pi * w**2 * np.exp(-(np.pi * w * q.nu) ** 2), atol=1e-12)
    assert abs(np.sum(F / (np.pi * q.R**2 * q.J1**2)) - 1.0) < 1e-10
