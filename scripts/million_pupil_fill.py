"""Pupil fill of a lenslet aperture imaged onto the circular objective pupil (Fourier mapping).

Usage: python scripts/million_pupil_fill.py

In the lenslet relay each lenslet aperture (uniformly lit) is imaged into the objective pupil,
so its shape is the pupil-plane amplitude: a phase cannot turn a square into a disc. For a
uniform aperture S scaled to the stop radius 1, the transmitted fraction is |S n D| / |S| and the
focal peak per watt relative to a uniformly filled disc is |S n D|^2 / (|S| pi) (D: unit disc).
Squares (half-width h) and regular hexagons (inradius r) are evaluated, circumscribed and
optimised over size. Output: results/million/pupil_fill.json.
"""

import numpy as np

from metatrap.runio import REPO_ROOT, write_json

n = 4001
c = (np.arange(n) + 0.5) / n * 4 - 2
X, Y = np.meshgrid(c, c, indexing="ij")
dA = (c[1] - c[0]) ** 2
disc = X**2 + Y**2 <= 1


def metrics(mask):
    s, inter = mask.sum() * dA, (mask & disc).sum() * dA
    return inter**2 / (s * np.pi), inter / s


def square(h):
    return (np.abs(X) <= h) & (np.abs(Y) <= h)


def hexagon(r):  # flat sides at distance r, three orientations of the side normals
    m = np.ones_like(X, bool)
    for a in (0, np.pi / 3, 2 * np.pi / 3):
        m &= np.abs(X * np.cos(a) + Y * np.sin(a)) <= r
    return m


res = {}
pk, tr = metrics(square(1.0))
res["square_circumscribed"] = {"peak_per_W": pk, "transmitted": tr}
pk, tr = metrics(square(1 / np.sqrt(2)))
res["square_inscribed"] = {"peak_per_W": pk, "transmitted": tr}
hs = np.round(np.arange(0.70, 1.001, 0.01), 2)
v = [metrics(square(h)) for h in hs]
i = int(np.argmax([a for a, _ in v]))
res["square_optimum"] = {"peak_per_W": v[i][0], "half_width": float(hs[i]), "transmitted": v[i][1]}
pk, tr = metrics(hexagon(1.0))
res["hex_circumscribed"] = {"peak_per_W": pk, "transmitted": tr}
rs = np.round(np.arange(0.80, 1.001, 0.01), 2)
v = [metrics(hexagon(r)) for r in rs]
i = int(np.argmax([a for a, _ in v]))
res["hex_optimum"] = {"peak_per_W": v[i][0], "inradius": float(rs[i]), "transmitted": v[i][1]}
res = {k: {a: float(b) for a, b in d.items()} for k, d in res.items()}
print(res)
write_json(REPO_ROOT / "results" / "million" / "pupil_fill.json", res)
