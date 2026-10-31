import numpy as np
import pytest

from metatrap import materials as m


def test_fused_silica_1064_matches_corning_datasheet():
    # Corning 7980 product sheet: n = 1.449633 at 1064.00 nm (retrieved 2026-10-31).
    assert m.FUSED_SILICA.n(1.064) == pytest.approx(1.449633, abs=2e-5)


def test_sin_lpcvd_reasonable():
    n = m.SIN_LPCVD.n(np.array([0.852, 1.075]))
    assert np.all((n > 1.99) & (n < 2.03))
    assert n[0] > n[1]  # normal dispersion


def test_validity_range_enforced():
    with pytest.raises(ValueError):
        m.SIN_LPCVD.n(0.2)


def test_csv_roundtrip(tmp_path):
    f = tmp_path / "film.csv"
    f.write_text("wavelength_um,n,k\n0.8,2.05,1e-5\n1.1,2.00,2e-5\n")
    mat = m.get(f"csv:{f}")
    assert mat.n(0.95) == pytest.approx(2.025)
    assert mat.k(0.95) == pytest.approx(1.5e-5)


def test_sapphire_ordinary_1064():
    # Vendor sheet quoted in search (n_o = 1.75449 at 1.06 um); Dodge formula agrees.
    assert m.SAPPHIRE_O.n(1.06) == pytest.approx(1.75449, abs=3e-4)
