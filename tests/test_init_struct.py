"""Tests for the init_struct module."""

import pytest

from tnp_gen.init_struct.gen_mnp import gen_mnp, write_mnp
from tnp_gen.init_struct.gen_bnp_al import gen_bnp, write_bnp
from tnp_gen.init_struct.gen_tnp_al import gen_tnp, write_tnp
from tnp_gen.constants import ELE_DICT, DIAMETER_LIST, SHAPE_LIST


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
