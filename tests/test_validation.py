"""Rigorous structure validation tests for generated nanoparticles."""

import numpy as np
import pytest

from tnp_gen.constants import ELE_DICT
from tnp_gen.init_struct.gen_bnp_cs import gen_hard_core_shell
from tnp_gen.init_struct.gen_mnp import gen_mnp


def get_min_dist(atoms):
    """Calculate minimum distance between any two atoms."""
    dists = atoms.get_all_distances()
    # Mask the diagonal (distance to self)
    masked_dists = np.where(dists > 0, dists, np.inf)
    return np.min(masked_dists)


@pytest.mark.parametrize("element", ["Au", "Pd", "Pt"])
@pytest.mark.parametrize("shape", ["OT", "SP", "IC", "CU", "TH", "RD", "TO", "CO"])
def test_mnp_physical_validity(element, shape):
    """Verify that generated MNPs have reasonable interatomic distances."""
    diameter = 20
    lat_const = ELE_DICT[element]["lc"]["FCC"]
    atoms = gen_mnp(shape, diameter, element, lat_const)

    # Check minimum distance between atoms
    expected_nn = lat_const / np.sqrt(2)
    min_dist = get_min_dist(atoms)

    assert min_dist > expected_nn * 0.8
    assert min_dist < expected_nn * 1.2


def test_hard_core_shell_validity():
    """Verify that hard core-shell BNPs have no overlapping atoms."""
    diameter_shell = 25
    diameter_core = 15
    shape = "OT"

    shell_atoms = gen_mnp(shape, diameter_shell, "Au", ELE_DICT["Au"]["lc"]["FCC"])
    core_atoms = gen_mnp(shape, diameter_core, "Pd", ELE_DICT["Pd"]["lc"]["FCC"])

    del_cutoff = (ELE_DICT["Au"]["radius"] + ELE_DICT["Pd"]["radius"]) / 2

    bnp = gen_hard_core_shell(shell_atoms, core_atoms, del_cutoff)

    symbols = bnp.get_chemical_symbols()
    assert "Au" in symbols
    assert "Pd" in symbols

    min_dist = get_min_dist(bnp)
    assert min_dist > del_cutoff * 0.7


def test_atom_counts_consistency():
    """Verify that atom counts are reasonable for the given diameter."""
    element = "Au"
    lat_const = ELE_DICT[element]["lc"]["FCC"]

    # Volume of a sphere: 4/3 * pi * r^3
    # Volume per atom in FCC: lat_const^3 / 4
    # Expected N = (4/3 * pi * r^3) / (lat_const^3 / 4)

    diameter = 30
    radius = diameter / 2
    vol_sphere = (4 / 3) * np.pi * (radius ** 3)
    vol_per_atom = (lat_const ** 3) / 4
    expected_n = vol_sphere / vol_per_atom

    atoms = gen_mnp("SP", diameter, element, lat_const)
    actual_n = len(atoms)

    # Spherical approximation should be within 20% for these sizes
    assert actual_n > expected_n * 0.8
    assert actual_n < expected_n * 1.2


def test_new_shapes_generate_atoms():
    """Verify that newly added shapes (TH, RD, TO, CO) generate valid atoms."""
    element = "Pd"
    diameter = 20
    lat_const = ELE_DICT[element]["lc"]["FCC"]
    for shape in ["TH", "RD", "TO", "CO"]:
        atoms = gen_mnp(shape, diameter, element, lat_const)
        assert len(atoms) > 0
        assert all(s == element for s in atoms.symbols)
