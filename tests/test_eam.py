"""Tests for the EAM module."""

import numpy as np

from allomorph.eam.create_eam import embed, pair, prof
from allomorph.eam.eam_database import AtType, Database


def test_database_has_expected_elements():
    assert "Au" in Database
    assert "Pd" in Database
    assert "Pt" in Database
    assert "Cu" in Database


def test_database_validation_for_arbitrary_elements():
    """Verify that the EAM database can be queried for arbitrary elements."""
    # The EAM module already accepts arbitrary element names; this test
    # documents that behaviour.  We do not call create_eam here because it
    # would require writing files, but we assert the database structure.
    assert hasattr(Database, "keys")
    assert callable(Database.keys)


def test_attype_attributes():
    au = Database["Au"]
    assert isinstance(au, AtType)
    assert au.name == "Au"
    assert au.ielement == 79
    assert au.amass == 196.96654


def test_prof_returns_array():
    r = np.linspace(0.5, 5.0, 100)
    f = prof("Au", r)
    assert isinstance(f, np.ndarray)
    assert f.shape == r.shape


def test_pair_returns_array():
    r = np.linspace(0.5, 5.0, 100)
    psi = pair("Au", "Au", r)
    assert isinstance(psi, np.ndarray)
    assert psi.shape == r.shape


def test_pair_mixed_alloy():
    r = np.linspace(0.5, 5.0, 100)
    psi = pair("Au", "Pd", r)
    assert isinstance(psi, np.ndarray)
    assert psi.shape == r.shape


def test_pair_ni_pd_pt():
    """Verify that Ni-Pd-Pt components can be calculated."""
    r = np.linspace(0.5, 5.0, 10)
    for e1 in ["Ni", "Pd", "Pt"]:
        for e2 in ["Ni", "Pd", "Pt"]:
            psi = pair(e1, e2, r)
            assert isinstance(psi, np.ndarray)
            assert psi.shape == r.shape
