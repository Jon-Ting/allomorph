# Purpose: Store variables for generation of NPs using ASE
# Author: Jonathan Yik Chang Ting
# Date: 19/20/2020

from math import sqrt

import numpy as np

LMP_DATA_DIR = "."
MNP_DIR = "MNP"
BNP_DIR = "BNP"
TNP_DIR = "TNP"
CS_DIR = "CS"
GOLDEN_RATIO = (1 + sqrt(5)) / 2
VACUUM_THICKNESS = 40.0
RANDOM_DISTRIB_NO = 3

# Elements of interest, their lattice parameters, metallic radii, and physical properties.
# - Pd, Pt, Au values were obtained from N. W. Ashcroft and N. D. Mermin,
#   Solid State Physics (Holt, Rinehart, and Winston, New York, 1976.
# - The lattice constants are 3.859, 3.912, and 4.065 Angstroms, for the
#   respective FCC metals at 300 K according to W. P. Davey,
#   "Precision Measurements of the Lattice Constants of Twelve Common Metals,"
#   Physical Review, vol. 25, (6), pp. 753-761, 1925.
# - Metallic radii taken from Greenwood, Norman N.; Earnshaw, Alan (1997).
#   Chemistry of the Elements (2nd ed.).
# - Density (rho) in kg/m^3, molar mass (m) in kg/mol, and bulk cohesive
#   energy (bulkE) in eV/atom are used by the feature-extraction pipeline.
ELE_DICT = {
    "Pd": {"lc": {"FCC": 3.89}, "radius": 1.37, "mass": 106.42,
           "rho": 12020, "m": 0.10642, "bulkE": 3.89},
    "Pt": {"lc": {"FCC": 3.92}, "radius": 1.39, "mass": 195.08,
           "rho": 21450, "m": 0.195084, "bulkE": 5.84},
    "Au": {"lc": {"FCC": 4.09}, "radius": 1.44, "mass": 196.97,
           "rho": 19320, "m": 0.196967, "bulkE": 3.81},
}

DIAMETER_LIST = [10, 15, 20, 25, 30]  # NP diameters of interest (Angstrom)
SHAPE_LIST = ["OT", "SP", "IC", "CU", "DH", "TH", "RD", "TO", "CO"]  # Shapes of interest
BNP_DISTRIB_LIST = ["L10", "RAL", "RCS"]  # BNP distributions of interest
TNP_DISTRIB_LIST = [
    "L10R", "CS", "CL10S", "CRALS", "RRAL", "CSRAL", "CSL10", "CRSR", "LL10"
]  # TNP distributions of interest
RATIO_LIST = [20, 40, 60, 80]  # Ratios of interest


def validate_ele_dict(ele_dict):
    """Validate a user-supplied element dictionary.

    Args:
        ele_dict: Dictionary mapping element symbols to property dicts.

    Returns:
        The validated dictionary.

    Raises:
        ValueError: If required keys are missing or invalid.
    """
    required_keys = {"lc", "radius", "mass"}
    for element, props in ele_dict.items():
        missing = required_keys - set(props.keys())
        if missing:
            raise ValueError(
                f"Element '{element}' is missing required keys: {missing}"
            )
        if "FCC" not in props.get("lc", {}):
            raise ValueError(
                f"Element '{element}' must have lattice constant 'lc' with an 'FCC' entry"
            )
    return ele_dict


def load_ele_dict_from_file(path):
    """Load an element dictionary from a JSON file.

    Args:
        path: Path to a JSON file containing the element dictionary.

    Returns:
        Validated element dictionary.
    """
    import json

    with open(path, "r") as f:
        ele_dict = json.load(f)
    return validate_ele_dict(ele_dict)


def dist_1d(coord1, coord2, dim):
    """Compute distance between 2 points in one of their real space coordinates."""
    return round(np.sqrt(np.sum((coord2[dim] - coord1[dim]) ** 2)), 3)


def dist_3d(coord1, coord2):
    """Compute real space distance between 2 points."""
    return round(np.sqrt(np.sum((coord2 - coord1) ** 2)), 3)


def calc_rcs_prob(obj, shape):
    """Compute the core-shell probability for each atom in *obj*.

    Returns a list where higher values correspond to surface-like positions.

    Args:
        obj: ASE Atoms object.
        shape: Nanoparticle shape string (e.g. 'IC', 'DH', 'OT').

    Returns:
        List of probabilities, one per atom.
    """
    prob_list = []
    if shape == 'IC':
        mass_center = obj.get_center_of_mass()
        radius = (obj.cell[0][0] - VACUUM_THICKNESS) / 2
        for atom in obj:
            prob_list.append(dist_3d(mass_center, atom.position) / radius)
    else:
        x_slices = {round(atom.position[0], 3) for atom in obj}
        y_slices = {round(atom.position[1], 3) for atom in obj}
        z_slices = {round(atom.position[2], 3) for atom in obj}

        z_thread = {(x, y): {'max': 0, 'min': 0, 'mid': []} for x in x_slices for y in y_slices}
        x_thread = {(y, z): {'max': 0, 'min': 0, 'mid': []} for y in y_slices for z in z_slices}
        y_thread = {(z, x): {'max': 0, 'min': 0, 'mid': []} for z in z_slices for x in x_slices}
        for atom in obj:
            x, y, z = round(atom.position[0], 3), round(atom.position[1], 3), round(atom.position[2], 3)
            z_thread[(x, y)]['mid'].append(z)
            x_thread[(y, z)]['mid'].append(x)
            y_thread[(z, x)]['mid'].append(y)
        for d in (z_thread, x_thread, y_thread):
            empty = [k for k, v in d.items() if len(v['mid']) == 0]
            for k in empty:
                del d[k]
            for k, v in d.items():
                v['max'] = max(v['mid'])
                v['min'] = min(v['mid'])
                v['mid'] = (v['max'] + v['min']) / 2

        for atom in obj:
            x, y, z = round(atom.position[0], 3), round(atom.position[1], 3), round(atom.position[2], 3)
            z_half = (z_thread[(x, y)]['max'] - z_thread[(x, y)]['min']) / 2
            x_half = (x_thread[(y, z)]['max'] - x_thread[(y, z)]['min']) / 2
            y_half = (y_thread[(z, x)]['max'] - y_thread[(z, x)]['min']) / 2
            z_rel = abs(z - z_thread[(x, y)]['mid']) / z_half if round(z_half, 3) != 0.0 else abs(z - z_thread[(x, y)]['mid'])
            x_rel = abs(x - x_thread[(y, z)]['mid']) / x_half if round(x_half, 3) != 0.0 else abs(x - x_thread[(y, z)]['mid'])
            y_rel = abs(y - y_thread[(z, x)]['mid']) / y_half if round(y_half, 3) != 0.0 else abs(y - y_thread[(z, x)]['mid'])
            if shape == 'DH':
                prob = 1.0 if (round(z_rel, 3) == 1.0) or (round(z_half, 3) == 0.0) else z_rel
            else:
                if (round(z_rel, 3) == 1.0) or (round(x_rel, 3) == 1.0) or (round(y_rel, 3) == 1.0) or (round(z_half, 3) == 0.0) or (round(x_half, 3) == 0.0) or (round(y_half, 3) == 0.0):
                    prob = 1.0
                else:
                    prob = (z_rel + x_rel + y_rel) / 3
            prob_list.append(prob)
    return prob_list


def parse_ele_comb(ele_comb, ele_dict=None):
    """Parse a concatenated element combination string into individual symbols.

    Args:
        ele_comb: String like 'AuPtPd' or 'CuAgAu'.
        ele_dict: Optional element dictionary to use for matching. If None,
            the built-in ELE_DICT is used.

    Returns:
        List of element symbols, e.g. ['Au', 'Pt', 'Pd'].

    Raises:
        ValueError: If the string cannot be parsed.
    """
    if ele_dict is None:
        ele_dict = ELE_DICT
    symbols = sorted(ele_dict.keys(), key=len, reverse=True)
    elements = []
    i = 0
    while i < len(ele_comb):
        matched = False
        for sym in symbols:
            if ele_comb[i:].startswith(sym):
                elements.append(sym)
                i += len(sym)
                matched = True
                break
        if not matched:
            raise ValueError(
                f"Could not parse element combination '{ele_comb}' at position {i}. "
                f"Known symbols: {list(ele_dict.keys())}"
            )
    return elements
