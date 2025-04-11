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
from constants import LMP_DATA_DIR, MNP_DIR, GOLDEN_RATIO, VACUUM_THICKNESS, eleDict, diameterList, shapeList


def genMNP(shape, diameter, element, latConst):
    if shape in ['CU', 'RD']:
        estLatNo = diameter / latConst
        latNoCU = round(estLatNo)
        if shape == 'CU':
            return FaceCenteredCubic(symbols=element, 
                                     surfaces=[(1,0,0)], 
                                     layers=[latNoCU], 
                                     latticeconstant=latConst)
        elif shape == 'RD':
            if (latNoCU % 2 == 0): diagLayerNoRD = latNoCU
            else:
                diagLayerNoRD = latNoCU-1 if (round(estLatNo+0.5) == latNoCU) else latNoCU+1
            return FaceCenteredCubic(symbols=element, 
                                     surfaces=[(1,0,0), (1,1,0)], 
                                     layers=[diagLayerNoRD, diagLayerNoRD], 
                                     latticeconstant=latConst)
    elif shape == 'TH':
        edgeAtomNoTH = round(1 + 2*diameter/latConst)
        diagLayerNoTH = edgeAtomNoTH + 2
        return FaceCenteredCubic(symbols=element, 
                                 surfaces=[(1,0,0), (1,1,1), (1,1,1), (1,1,-1), (1,-1,1), (-1,1,1)], 
                                 layers=[diagLayerNoTH - 1, diagLayerNoTH, -1, -1, -1, -1], 
                                 latticeconstant=latConst)
    elif shape in ['OT', 'TO', 'CO']:
        edgeLengthOT = diameter / sqrt(2)
        edgeAtomNoOT = 1 + round(edgeLengthOT / (latConst/sqrt(2)))
        if shape == 'OT':
            return Octahedron(symbol=element, 
                              length=edgeAtomNoOT, 
                              cutoff=0, 
                              latticeconstant=latConst, 
                              alloy=False)
        elif shape == 'TO':
            compLayerNo = diameter // latConst
            if (compLayerNo % 2 == 0): cutLayerNoTO = compLayerNo / 2
            else:
                cutLayerNoTO = (compLayerNo+1)/2 if (int(diameter/latConst) == compLayerNo) else (compLayerNo-1)/2
            edgeAtomNoTO = 3*cutLayerNoTO + 1
            return Octahedron(symbol=element, 
                              length=edgeAtomNoTO, 
                              cutoff=cutLayerNoTO, 
                              latticeconstant=latConst, 
                              alloy=False)
        elif shape == 'CO':
            cutLayerNoCO = round(diameter / latConst)
            edgeAtomNoCO = 2*cutLayerNoCO + 1
            return Octahedron(symbol=element, 
                              length=edgeAtomNoCO, 
                              cutoff=cutLayerNoCO, 
                              latticeconstant=latConst, 
                              alloy=False)
    elif shape == 'DH':
        edgeLengthDH = 2 * diameter * cos(radians(72))
        edgeAtomDistDH = latConst/sqrt(2)*cos(radians(30)) * (1+sin(radians(54))) / (GOLDEN_RATIO*sin(radians(72)))
        edgeAtomNoDH = 1 + round(edgeLengthDH / edgeAtomDistDH)
        return Decahedron(symbol=element, 
                          p=edgeAtomNoDH, 
                          q=1, 
                          r=0, 
                          latticeconstant=latConst)
    elif shape == 'IC':
        circRadIC = sqrt(diameter**2*(1+GOLDEN_RATIO**2) / (4*GOLDEN_RATIO**2))
        shellNoIC = 1 + round(circRadIC / (latConst/sqrt(2)))
        return Icosahedron(symbol=element, 
                           noshells=shellNoIC, 
                           latticeconstant=latConst)


def writeMNP(element, diameter, latConst, shape, replace=False, vis=False):
    fileName = f"{element}{diameter}{shape}.lmp"
    if not replace:
        if isfile(f"{LMP_DATA_DIR}/{MNP_DIR}/{fileName}"): 
            print(f"      {fileName} already exist, skipping...")
            return
    
    # Add vacuum
    mnp = genMNP(shape, diameter, element, latConst)
    boxSize = [dim[i] + VACUUM_THICKNESS for (i, dim) in enumerate(mnp.cell)]
    mnp.set_cell(boxSize)
    mnp.translate([VACUUM_THICKNESS / 2] * 3)

    if not isdir(f"{LMP_DATA_DIR}/{MNP_DIR}"): mkdir(f"{LMP_DATA_DIR}/{MNP_DIR}")
    write_lammps_data(f"{LMP_DATA_DIR}/{MNP_DIR}/{fileName}", atoms=mnp, units='metal', atom_style='atomic')
    fileName = f"{element}{diameter}{shape}.xyz"
    xyzFile = open(f"{LMP_DATA_DIR}/{MNP_DIR}/{fileName}", 'w')
    write_xyz(xyzFile, [mnp])
    xyzFile.close()

    print(f"      Generated {fileName}, diameter: {boxSize[0]:.1f} A^3, size: {mnp.get_global_number_of_atoms()} atoms")
    if vis: view(mnp)


def writeMNP_sphere(element, diameter, latConst, replace=False, vis=False):
    fileName = f"{element}{diameter}SP.lmp"
    if not replace:
        if isfile(f"{LMP_DATA_DIR}/{MNP_DIR}/{fileName}"):
            print(f"      {fileName} already exist, skipping...")
            return

    # Add vacuum
    mnp = genMNP('CU', diameter, element, latConst)

    # Delete atoms that are out of range
    center = mnp.get_center_of_mass()
    del mnp[[atom.index for atom in mnp if np.linalg.norm(center - atom.position) > diameter / 2]]

    boxSize = [dim[i] + VACUUM_THICKNESS for (i, dim) in enumerate(mnp.cell)]  # * np.array([diameter / 10, diameter / 10, diameter / 10])
    mnp.set_cell(boxSize)
    mnp.translate([VACUUM_THICKNESS / 2] * 3)

    # Write atoms structure to a lmp file and create file if not exists
    if not isdir(f"{LMP_DATA_DIR}/{MNP_DIR}"): mkdir(f"{LMP_DATA_DIR}/{MNP_DIR}")
    write_lammps_data(f"{LMP_DATA_DIR}/{MNP_DIR}/{fileName}", atoms=mnp, units='metal', atom_style='atomic')
    fileName = f"{element}{diameter}SP.xyz"
    xyzFile = open(f"{LMP_DATA_DIR}/{MNP_DIR}/{fileName}", 'w')
    write_xyz(xyzFile, [mnp])
    xyzFile.close()

    print(f"      Generated {fileName}, diameter: {boxSize[0]:.1f} A^3, size: {mnp.get_global_number_of_atoms()} atoms")
    if vis: view(mnp)


def main(replace=False, vis=False):
    if not isdir(LMP_DATA_DIR): mkdir(LMP_DATA_DIR)
    print(f"Generating NPs with {VACUUM_THICKNESS} Angstrom of vacuum on each dimension:")
    for diameter in diameterList:
        print(f"\n  Size {diameter} Angstrom for:")
        for element in eleDict:
            if element != 'Pd': 
                continue
            print(f"    Element {element}:")
            latConst = eleDict[element]['lc']['FCC']
            for shape in shapeList:
                writeMNP(element, diameter, latConst, shape, replace=replace, vis=vis)
                # writeMNP_sphere(element, diameter, latConst, replace=replace, vis=vis)
            

if __name__ == '__main__':
    main(replace=False, vis=False)
    print('ALL DONE!')


