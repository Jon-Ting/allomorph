"""Tests for the init_struct module."""

import pytest

from tnp_gen.constants import ELE_DICT
from tnp_gen.init_struct.gen_bnp_al import gen_bnp
from tnp_gen.init_struct.gen_mnp import gen_mnp
from tnp_gen.init_struct.gen_tnp_al import gen_tnp


def test_gen_mnp_creates_atoms():
    """Test that gen_mnp returns an ASE Atoms object."""
    element = "Au"
    diameter = 30
    shape = "OT"
    lat_const = ELE_DICT[element]["lc"]["FCC"]
    atoms = gen_mnp(shape, diameter, element, lat_const)
    assert atoms is not None
    assert len(atoms) > 0
    assert all(s == element for s in atoms.symbols)


def test_gen_mnp_cube():
    element = "Pd"
    diameter = 30
    shape = "CU"
    lat_const = ELE_DICT[element]["lc"]["FCC"]
    atoms = gen_mnp(shape, diameter, element, lat_const)
    assert len(atoms) > 0


def test_gen_mnp_icosahedron():
    element = "Pt"
    diameter = 30
    shape = "IC"
    lat_const = ELE_DICT[element]["lc"]["FCC"]
    atoms = gen_mnp(shape, diameter, element, lat_const)
    assert len(atoms) > 0


@pytest.mark.parametrize("shape", ["TH", "RD", "TO", "CO"])
def test_gen_mnp_new_shapes(shape):
    """Test that newly added shapes generate valid nanoparticles."""
    element = "Au"
    diameter = 20
    lat_const = ELE_DICT[element]["lc"]["FCC"]
    atoms = gen_mnp(shape, diameter, element, lat_const)
    assert len(atoms) > 0
    assert all(s == element for s in atoms.symbols)


def test_gen_bnp_ral():
    """Test that gen_bnp converts atoms to an alloy."""
    element = "Au"
    diameter = 30
    shape = "OT"
    lat_const = ELE_DICT[element]["lc"]["FCC"]
    mnp = gen_mnp(shape, diameter, element, lat_const)
    bnp = gen_bnp(mnp.copy(), "Pd", shape, 40, "RAL", 0)
    symbols = list(bnp.symbols)
    assert "Au" in symbols
    assert "Pd" in symbols


def test_gen_tnp_ral():
    """Test that gen_tnp converts atoms to a trimetallic alloy."""
    element = "Au"
    diameter = 30
    shape = "OT"
    lat_const = ELE_DICT[element]["lc"]["FCC"]
    mnp = gen_mnp(shape, diameter, element, lat_const)
    tnp = gen_tnp(mnp, "Au", "Pd", "Pt", 40, 30, 30, "RAL", "RAL", 0)
    symbols = list(tnp.symbols)
    assert "Au" in symbols
    assert "Pd" in symbols
    assert "Pt" in symbols


def test_gen_tnp_l10_l10():
    """Test the fully ordered LL10 distribution."""
    from tnp_gen.constants import VACUUM_THICKNESS
    element = "Au"
    diameter = 20
    shape = "OT"
    lat_const = ELE_DICT[element]["lc"]["FCC"]
    mnp = gen_mnp(shape, diameter, element, lat_const)
    # The L10 logic assumes atoms have been translated by the vacuum offset
    # (as write_mnp does before saving to file).
    mnp.translate([VACUUM_THICKNESS / 2] * 3)
    tnp = gen_tnp(mnp, "Au", "Pd", "Pt", 33, 33, 34, "L10", "L10", 0)
    symbols = list(tnp.symbols)
    assert "Au" in symbols
    assert "Pd" in symbols
    assert "Pt" in symbols


def test_tnp_distribution_mapping():
    """Test that all BNP-pair combinations map to a valid TNP distribution."""
    from tnp_gen.constants import TNP_DISTRIB_LIST
    from tnp_gen.init_struct.gen_tnp_al import _map_tnp_distribution

    # Every valid (distrib1, distrib2) pair should map to a TNP name
    valid_pairs = [
        ("L10", "RAL"), ("L10", "L10"), ("L10", "RCS"),
        ("RAL", "RAL"), ("RAL", "RCS"),
        ("RCS", "RAL"), ("RCS", "L10"), ("RCS", "RCS"),
    ]
    for d1, d2 in valid_pairs:
        name = _map_tnp_distribution(d1, d2)
        assert name in TNP_DISTRIB_LIST


def test_tnp_distribution_invalid_pair():
    """Test that invalid distribution pairs raise ValueError."""
    from tnp_gen.init_struct.gen_tnp_al import _map_tnp_distribution
    with pytest.raises(ValueError):
        _map_tnp_distribution("RAL", "L10")
