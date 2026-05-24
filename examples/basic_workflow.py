"""
Basic workflow example for tnp-gen.

This script demonstrates how to use the tnp-gen library to generate
monometallic and bimetallic nanoparticles with arbitrary elements.
"""

from tnp_gen.init_struct.gen_mnp import gen_mnp
from tnp_gen.init_struct.gen_bnp_al import gen_bnp
from tnp_gen.constants import ELE_DICT
from ase.visualize import view


def main():
    print("=== Monometallic nanoparticle ===")
    # Generate an octahedral Pd nanoparticle with diameter 30 Angstrom
    pd_np = gen_mnp(
        shape='OT',
        diameter=30,
        element='Pd',
        latConst=ELE_DICT['Pd']['lc']['FCC']
    )
    print(f"Pd MNP: {len(pd_np)} atoms")

    print("\n=== Bimetallic alloy (RAL) ===")
    # Convert 30% of atoms to Pt randomly
    bnp = gen_bnp(
        obj=pd_np.copy(),
        element2='Pt',
        shape='OT',
        ratio2=30,
        distrib='RAL',
        rseed=42
    )
    print(f"Pd-Pt BNP formula: {bnp.get_chemical_formula()}")

    print("\n=== Trimetallic alloy (L10 + RAL) ===")
    from tnp_gen.init_struct.gen_tnp_al import gen_tnp
    tnp = gen_tnp(
        obj=bnp.copy(),
        element1='Pd',
        element2='Pt',
        element3='Au',
        ele1Ratio=40,
        ele2Ratio=30,
        ele3Ratio=30,
        distrib1='L10',
        distrib2='RAL',
        rseed=0,
    )
    print(f"Pd-Pt-Au TNP formula: {tnp.get_chemical_formula()}")

    # In a real scenario you would use write_mnp / write_bnp / write_tnp
    # to save .lmp files.  Uncomment below to visualise with ASE GUI:
    # view(tnp)


if __name__ == "__main__":
    main()
