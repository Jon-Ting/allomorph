"""Tests for the EAM module."""

import numpy as np

from tnp_gen.eam.eam_database import AtType, Database
from tnp_gen.eam.create_eam import prof, pair, embed


def test_database_has_expected_elements():
    assert "Au" in Database
    assert "Pd" in Database
    assert "Pt" in Database
    assert "Cu" in Database


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


def test_embed_returns_array():
    rho = np.linspace(0.0, 50.0, 100)
    emb = embed("Au", rho)
    assert isinstance(emb, np.ndarray)
    assert emb.shape == rho.shape
