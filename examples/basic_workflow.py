"""
Basic workflow example for tnp-gen.
This script demonstrates how to use the tnp-gen library to generate nanoparticles.
"""

from tnp_gen.init_struct.gen_mnp import gen_mnp
from tnp_gen.init_struct.gen_bnp_al import gen_bnp
from ase.visualize import view

def main():
    print("Generating a monometallic Palladium nanoparticle...")
    # Generate an octahedral Pd nanoparticle with diameter 30 Angstrom
    # Lattice constant for Pd is ~3.89
    pd_np = gen_mnp(shape='OT', diameter=30, element='Pd', latConst=3.89)
    print(f"Number of atoms: {len(pd_np)}")
    
    print("Converting it to a Randomly distributed Alloy (RAL) with Platinum...")
    # Convert 30% of atoms to Pt randomly
    bnp = gen_bnp(obj=pd_np.copy(), element2='Pt', shape='OT', ratio2=30, distrib='RAL', rseed=42)
    print(f"BNP Formula: {bnp.get_chemical_formula()}")
    
    # In a real scenario, you would use write_mnp and write_bnp to save files.
    # view(bnp) # Uncomment to visualize if ASE GUI is available

if __name__ == "__main__":
    main()
