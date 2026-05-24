"""Generate hard core-shell BNP structures."""

import numpy as np
from pathlib import Path
from ase.io.lammpsdata import read_lammps_data, write_lammps_data
from ase.visualize import view
import logging

from tnp_gen.constants import (
    LMP_DATA_DIR, MNP_DIR, BNP_DIR, CS_DIR,
    ELE_DICT, DIAMETER_LIST, SHAPE_LIST
)

logger = logging.getLogger(__name__)

def gen_hard_core_shell(shell_atoms, core_atoms, del_cutoff):
    """
    Generate a hard core-shell structure by combining shell and core atoms.
    Shell atoms overlapping with core atoms within del_cutoff are deleted.
    """
    # Center core relative to shell
    shell_center = shell_atoms.get_center_of_mass()
    core_center = core_atoms.get_center_of_mass()
    core_atoms.translate(shell_center - core_center)
    
    # Identify shell atoms to delete
    shell_pos = shell_atoms.get_positions()
    core_pos = core_atoms.get_positions()
    
    to_delete = []
    # Using broadcasting for efficiency if atoms count is not too large
    # For very large NPs, a NeighborList or KDTree would be better.
    # Given the max diameter in DIAMETER_LIST is 30A, atom count is ~thousands.
    for i, s_p in enumerate(shell_pos):
        dists = np.linalg.norm(core_pos - s_p, axis=1)
        if np.any(dists < del_cutoff):
            to_delete.append(i)
            
    del shell_atoms[to_delete]
    
    # Combine
    # Note: we need to update atom types correctly
    # Shell atoms are type 1, core atoms are type 2 (as in LAMMPS script)
    combined = shell_atoms + core_atoms
    return combined

def write_hard_core_shell(replace=False, vis=False):
    """
    Orchestrates the generation of hard core-shell BNPs.
    Iterates through elements, sizes, and shapes.
    """
    output_base_dir = Path(LMP_DATA_DIR) / BNP_DIR / CS_DIR
    mnp_dir = Path(LMP_DATA_DIR) / MNP_DIR
    output_base_dir.mkdir(parents=True, exist_ok=True)

    logger.info("Generating hard core-shell bimetallic nanoparticles:")
    
    for element1, prop1 in ELE_DICT.items():
        radius1 = prop1['radius']
        for element2, prop2 in ELE_DICT.items():
            if element1 == element2:
                continue
            
            radius2 = prop2['radius']
            del_cutoff = (radius1 + radius2) / 2
            
            logger.info(f"  {element2}@{element1} (core@shell)")
            
            # Sizes: k is shell, l is core (shell must be larger)
            for k in range(1, len(DIAMETER_LIST)):
                size1 = DIAMETER_LIST[k]
                for l in range(k):
                    size2 = DIAMETER_LIST[l]
                    
                    for shape1 in SHAPE_LIST:
                        file_name1 = f"{element1}{size1}{shape1}.lmp"
                        path1 = mnp_dir / file_name1
                        if not path1.exists():
                            continue
                            
                        shell_atoms = read_lammps_data(str(path1), style='atomic', units='metal')
                        shell_atoms.set_chemical_symbols([element1] * len(shell_atoms))
                        
                        for shape2 in SHAPE_LIST:
                            file_name2 = f"{element2}{size2}{shape2}.lmp"
                            path2 = mnp_dir / file_name2
                            if not path2.exists():
                                continue
                                
                            core_atoms = read_lammps_data(str(path2), style='atomic', units='metal')
                            core_atoms.set_chemical_symbols([element2] * len(core_atoms))
                            
                            file_out = f"{element1}{size1}{shape1}{element2}{size2}{shape2}CS.lmp"
                            output_path = output_base_dir / file_out
                            
                            if not replace and output_path.exists():
                                logger.info(f"    {file_out} already exists, skipping...")
                                continue
                                
                            logger.info(f"    Generating {file_out}...")
                            
                            bnp = gen_hard_core_shell(shell_atoms.copy(), core_atoms.copy(), del_cutoff)
                            write_lammps_data(str(output_path), atoms=bnp, units='metal', atom_style='atomic')
                            
                            if vis:
                                view(bnp)

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    write_hard_core_shell()
