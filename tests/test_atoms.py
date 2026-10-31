"""Anchors for the polarisability and power-budget model (proposal Sections 3, 4)."""

import numpy as np
import pytest

from metatrap.atoms.alkali import alpha_au, alpha_static_measured_au
from metatrap.atoms.power_budget import (
    REF4_DARK_RADIUS_UM,
    REF4_DEPTH_MEASURED_MK,
    REF4_LAMBDA_NM,
    REF4_N_TWEEZERS,
    REF4_P_INCIDENT_W,
    PowerBudget,
)
from metatrap.constants import RUBIDIUM87


def test_static_limit_recovers_measured_alpha0():
    assert alpha_au(1.0e9) == pytest.approx(alpha_static_measured_au(), rel=1e-6)
    assert alpha_static_measured_au() == pytest.approx(401.0, abs=0.5)  # Steck Table 6


def test_cs_alpha_1064_matches_proposal_value():
    # Proposal Section 3: ~1,160 a.u. at 1064 nm (ref. 36); model within 3 %.
    assert alpha_au(1064.0) == pytest.approx(1160.0, rel=0.03)
    # Proposal Section 4 (ii): ~1,120 a.u. at 1075 nm.
    assert alpha_au(1075.0) == pytest.approx(1120.0, rel=0.03)


def test_ref4_depth_reproduced():
    na_eff = 0.61 * REF4_LAMBDA_NM * 1e-3 / REF4_DARK_RADIUS_UM
    pb = PowerBudget(lam_nm=REF4_LAMBDA_NM, na_eff=na_eff, atom=RUBIDIUM87)
    depth_mk = pb.depth_K(REF4_P_INCIDENT_W, REF4_N_TWEEZERS) * 1e3
    assert depth_mk == pytest.approx(REF4_DEPTH_MEASURED_MK, rel=0.02)


def test_site_counts_section4():
    pb = PowerBudget()  # Cs, 1075 nm, utilisation 0.37, NA_eff 0.47
    assert pb.sites(50.0, 1.0e-3) == pytest.approx(2800, rel=0.02)
    assert pb.sites(50.0, 0.305e-3) == pytest.approx(9200, rel=0.02)


def test_transport_tweezer_power():
    p = PowerBudget(lam_nm=1064.435, na_eff=0.5).power_per_site_W(1.2e-3)
    assert p * 1e3 == pytest.approx(6.6, rel=0.02)
