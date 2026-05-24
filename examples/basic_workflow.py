"""
Basic workflow example for allomorph.

This script demonstrates how to use the allomorph library to generate
monometallic, bimetallic, and trimetallic nanoparticles with arbitrary elements.
"""

from allomorph.constants import ELE_DICT, update_constants, validate_ele_dict
from allomorph.init_struct.gen_bnp_al import gen_bnp
from allomorph.init_struct.gen_mnp import gen_mnp
from allomorph.init_struct.gen_tnp_al import gen_tnp


def main():
    """Main function to demonstrate the workflow."""
    print("=== Custom Configuration ===")
    # You can define custom elements and update the global configuration
    custom_ele_dict = {
        "Cu": {"lc": {"FCC": 3.61}, "radius": 1.28, "mass": 63.55,
               "rho": 8960, "m": 0.063546, "bulkE": 3.49},
        "Ag": {"lc": {"FCC": 4.09}, "radius": 1.44, "mass": 107.87,
               "rho": 10490, "m": 0.107868, "bulkE": 2.95},
        "Au": {"lc": {"FCC": 4.09}, "radius": 1.44, "mass": 196.97,
               "rho": 19320, "m": 0.196967, "bulkE": 3.81},
    }
    # Validate before updating
    validate_ele_dict(custom_ele_dict)

    # Update global constants (this affects all modules using ELE_DICT)
    update_constants({"ELE_DICT": custom_ele_dict, "RANDOM_DISTRIB_NO": 2})

    print("=== Monometallic nanoparticle ===")
    # Generate an octahedral Cu nanoparticle with diameter 30 Angstrom
    cu_np = gen_mnp(
        shape='OT',
        diameter=30,
        element='Cu',
        latConst=ELE_DICT['Cu']['lc']['FCC']
    )
    print(f"Cu MNP: {len(cu_np)} atoms")

    print("\n=== Bimetallic alloy (RAL) ===")
    # Convert 30% of atoms to Ag randomly
    bnp = gen_bnp(
        obj=cu_np.copy(),
        element2='Ag',
        shape='OT',
        ratio2=30,
        distrib='RAL',
        rseed=42
    )
    print(f"Cu-Ag BNP formula: {bnp.get_chemical_formula()}")

    print("\n=== Trimetallic alloy (RAL + RAL) ===")
    # Convert some atoms to Au to make it trimetallic
    tnp = gen_tnp(
        obj=bnp.copy(),
        element1='Cu',
        element2='Ag',
        element3='Au',
        ele1Ratio=40,
        ele2Ratio=30,
        ele3Ratio=30,
        distrib1='RAL',
        distrib2='RAL',
        rseed=0,
    )
    print(f"Cu-Ag-Au TNP formula: {tnp.get_chemical_formula()}")

    # In a real scenario you would use write_mnp / write_bnp / write_tnp
    # to save .lmp files.  Uncomment below to visualise with ASE GUI:
    # from ase.visualize import view
    # view(tnp)


if __name__ == "__main__":
    main()
