# Goal: Generate initial NP structures for MD simulations
# Authors: Jonathan Yik Chang Ting, Kaihan Lu, Haotai Peng
# Date: 19/10/2020
"""
Note:
- PBC was not used during structure generation
- Abbreviations:
    - CU = cube
    - TH = tetrahedron
    - RD = rhombic dodecahedron
    - OT = octahedron
    - TO = truncated octahedron (regular)
    - CO = cuboctahedron
    - DH = decahedron (pentagonal bipyramid)
    - IC = icosahedron (regular convex)
    - SP = sphere
- To do:
    - Warn if box size is abnormal
    - Implement shapes for hcp
    - Could make use of Octahedron(alloy=True) to generate L1_2 alloys
"""

from pathlib import Path
from math import sqrt, radians, cos, sin
from os.path import isfile, isdir
from os import mkdir
import numpy as np
from ase.build import add_vacuum
from ase.cluster.cubic import FaceCenteredCubic
from ase.cluster import Octahedron, Decahedron, Icosahedron
from ase.io.lammpsdata import write_lammps_data
from ase.io.xyz import write_xyz
from ase.visualize import view

from tnp_gen.constants import (
    LMP_DATA_DIR,
    MNP_DIR,
    GOLDEN_RATIO,
    VACUUM_THICKNESS,
    ELE_DICT,
    DIAMETER_LIST,
    SHAPE_LIST,
)


def gen_mnp(shape, diameter, element, latConst):
    """Generate a monometallic nanoparticle using ASE."""
    if shape in ['CU', 'RD']:
        estLatNo = diameter / latConst
        latNoCU = round(estLatNo)
        if shape == 'CU':
            return FaceCenteredCubic(
                symbols=element,
                surfaces=[(1, 0, 0)],
                layers=[latNoCU],
                latticeconstant=latConst,
            )
        elif shape == 'RD':
            if latNoCU % 2 == 0:
                diagLayerNoRD = latNoCU
            else:
                diagLayerNoRD = latNoCU - 1 if (round(estLatNo + 0.5) == latNoCU) else latNoCU + 1
            return FaceCenteredCubic(
                symbols=element,
                surfaces=[(1, 0, 0), (1, 1, 0)],
                layers=[diagLayerNoRD, diagLayerNoRD],
                latticeconstant=latConst,
            )
    elif shape == 'TH':
        edgeAtomNoTH = round(1 + 2 * diameter / latConst)
        diagLayerNoTH = edgeAtomNoTH + 2
        return FaceCenteredCubic(
            symbols=element,
            surfaces=[
                (1, 0, 0),
                (1, 1, 1),
                (1, 1, 1),
                (1, 1, -1),
                (1, -1, 1),
                (-1, 1, 1),
            ],
            layers=[diagLayerNoTH - 1, diagLayerNoTH, -1, -1, -1, -1],
            latticeconstant=latConst,
        )
    elif shape in ['OT', 'TO', 'CO']:
        edgeLengthOT = diameter / sqrt(2)
        edgeAtomNoOT = 1 + round(edgeLengthOT / (latConst / sqrt(2)))
        if shape == 'OT':
            return Octahedron(
                symbol=element,
                length=edgeAtomNoOT,
                cutoff=0,
                latticeconstant=latConst,
                alloy=False,
            )
        elif shape == 'TO':
            compLayerNo = diameter // latConst
            if compLayerNo % 2 == 0:
                cutLayerNoTO = compLayerNo / 2
            else:
                cutLayerNoTO = (compLayerNo + 1) / 2 if (int(diameter / latConst) == compLayerNo) else (compLayerNo - 1) / 2
            edgeAtomNoTO = 3 * cutLayerNoTO + 1
            return Octahedron(
                symbol=element,
                length=edgeAtomNoTO,
                cutoff=cutLayerNoTO,
                latticeconstant=latConst,
                alloy=False,
            )
        elif shape == 'CO':
            cutLayerNoCO = round(diameter / latConst)
            edgeAtomNoCO = 2 * cutLayerNoCO + 1
            return Octahedron(
                symbol=element,
                length=edgeAtomNoCO,
                cutoff=cutLayerNoCO,
                latticeconstant=latConst,
                alloy=False,
            )
    elif shape == 'DH':
        edgeLengthDH = 2 * diameter * cos(radians(72))
        edgeAtomDistDH = (
            latConst / sqrt(2) * cos(radians(30))
            * (1 + sin(radians(54)))
            / (GOLDEN_RATIO * sin(radians(72)))
        )
        edgeAtomNoDH = 1 + round(edgeLengthDH / edgeAtomDistDH)
        return Decahedron(
            symbol=element,
            p=edgeAtomNoDH,
            q=1,
            r=0,
            latticeconstant=latConst,
        )
    elif shape == 'IC':
        circRadIC = sqrt(diameter ** 2 * (1 + GOLDEN_RATIO ** 2) / (4 * GOLDEN_RATIO ** 2))
        shellNoIC = 1 + round(circRadIC / (latConst / sqrt(2)))
        return Icosahedron(
            symbol=element,
            noshells=shellNoIC,
            latticeconstant=latConst,
        )
    elif shape == 'SP':
        # Generate a large cube and cut it into a sphere
        mnp = gen_mnp('CU', diameter, element, latConst)
        center = mnp.get_center_of_mass()
        to_delete = [atom.index for atom in mnp if np.linalg.norm(center - atom.position) > diameter / 2]
        del mnp[to_delete]
        return mnp


def write_mnp(element, diameter, lat_const, shape, replace=False, vis=False):
    """Write a monometallic nanoparticle to LAMMPS data and XYZ files."""
    file_name_lmp = f"{element}{diameter}{shape}.lmp"
    output_dir = Path(LMP_DATA_DIR) / MNP_DIR
    output_path_lmp = output_dir / file_name_lmp
    
    if not replace and output_path_lmp.exists():
        print(f"      {file_name_lmp} already exists, skipping...")
        return

    # Generate and add vacuum
    mnp = gen_mnp(shape, diameter, element, lat_const)
    box_size = [dim[i] + VACUUM_THICKNESS for (i, dim) in enumerate(mnp.cell)]
    mnp.set_cell(box_size)
    mnp.translate([VACUUM_THICKNESS / 2] * 3)

    output_dir.mkdir(parents=True, exist_ok=True)
    
    write_lammps_data(
        str(output_path_lmp),
        atoms=mnp,
        units='metal',
        atom_style='atomic',
    )
    
    file_name_xyz = f"{element}{diameter}{shape}.xyz"
    output_path_xyz = output_dir / file_name_xyz
    with open(output_path_xyz, 'w') as f:
        write_xyz(f, [mnp])

    print(f"      Generated {file_name_lmp}, diameter: {box_size[0]:.1f} A, size: {mnp.get_global_number_of_atoms()} atoms")
    if vis:
        view(mnp)


def main(replace=False, vis=False):
    """Generate all monometallic nanoparticles defined in constants."""
    print(f"Generating NPs with {VACUUM_THICKNESS} Angstrom of vacuum on each dimension:")
    for diameter in DIAMETER_LIST:
        print(f"\n  Size {diameter} Angstrom for:")
        for element in ELE_DICT:
            print(f"    Element {element}:")
            lat_const = ELE_DICT[element]['lc']['FCC']
            for shape in SHAPE_LIST:
                write_mnp(element, diameter, lat_const, shape, replace=replace, vis=vis)


if __name__ == '__main__':
    main(replace=False, vis=False)
    print('ALL DONE!')
