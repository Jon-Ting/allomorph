# Goal: Generate initial alloy TNP structures for MD simulations
# Author: Jonathan Yik Chang Ting
# Date: 22/10/2020
'''
Note:
- Abbreviations:
    - RAL = randomly distributed alloy
    - RCS = randomly distributed core-shell-like alloy
    - (R)L10 = L1_0 intermetallic alloy (with/without random component)
    - (R)L12 = L1_2 intermetallic alloy (with/without random component)
- To do:
    - FCC is currently being hard-coded for lattice constant retrieval, might need to be flexible
    - Perhaps could add a parameter to control core thickness
'''

import numpy as np
from os.path import isfile, isdir
from os import mkdir
from numpy.random import seed, rand, RandomState
from ase.io.lammpsdata import read_lammps_data, write_lammps_data
from ase.visualize import view
from constants import LMP_DATA_DIR, MNP_DIR, BNP_DIR, TNP_DIR, RANDOM_DISTRIB_NO, VACUUM_THICKNESS, ELE_DICT, DIAMETER_LIST, SHAPE_LIST, BNP_DISTRIB_LIST, TNP_DISTRIB_LIST, RATIO_LIST


def dist1D(coord1, coord2, dim):
    """Compute distance between 2 points in one of their real space coordinates"""
    return round(np.sqrt(np.sum((coord2[dim]-coord1[dim]) ** 2)), 3)


def dist3D(coord1, coord2):
    """Compute real space distance between 2 points"""
    return round(np.sqrt(np.sum((coord2 - coord1) ** 2)), 3)


def randConv(obj, element1, element2, element3, ele1Ratio, ele2Ratio, ele3Ratio, rseed, prob):
    """Randomly convert elements of atoms until specified ratio is reached"""
    ele1Arr, ele2Arr, ele3Arr = obj.symbols.search(element1), obj.symbols.search(element2), obj.symbols.search(element3)
    ele1IdealNum, ele2IdealNum, ele3IdealNum = round(ele1Ratio / 100 * len(obj)), round(ele2Ratio / 100 * len(obj)), round(ele3Ratio / 100 * len(obj))
    diff1, diff2, diff3 = len(ele1Arr) - ele1IdealNum, len(ele2Arr) - ele2IdealNum, len(ele3Arr) - ele3IdealNum

    randGen = RandomState(rseed)
    # Put excessive element into temporary array
    idxArr = np.array([], dtype=int)
    if diff1 > 0: idxArr = np.concatenate((idxArr, randGen.choice(a=ele1Arr, size=abs(diff1), replace=False, p=None)))
    if diff2 > 0: idxArr = np.concatenate((idxArr, randGen.choice(a=ele2Arr, size=abs(diff2), replace=False, p=None)))
    if diff3 > 0: idxArr = np.concatenate((idxArr, randGen.choice(a=ele3Arr, size=abs(diff3), replace=False, p=None)))
    # Assign atoms from this array to element with insufficient amount
    if diff1 < 0:
        if abs(diff1) > len(idxArr): diff1 = len(idxArr)
        idxArr1 = randGen.choice(a=idxArr, size=abs(diff1), replace=False, p=None)
        for idx in idxArr1: obj[idx].symbol = element1
        idxArr = np.array([idx for idx in idxArr if idx not in idxArr1])
    if diff2 < 0:
        if abs(diff2) > len(idxArr): diff2 = len(idxArr)
        idxArr2 = randGen.choice(a=idxArr, size=abs(diff2), replace=False, p=None)
        for idx in idxArr2: obj[idx].symbol = element2
        idxArr = np.array([idx for idx in idxArr if idx not in idxArr2])
    if diff3 < 0:
        if abs(diff3) > len(idxArr): diff3 = len(idxArr)
        idxArr3 = randGen.choice(a=idxArr, size=abs(diff3), replace=False, p=None)
        for idx in idxArr3: obj[idx].symbol = element3
    return obj


def genTNP(obj, element1, element2, element3, ele1Ratio, ele2Ratio, ele3Ratio, distrib1, distrib2, rseed=0):
    probList = []
    if distrib2 == 'RAL': pass
    elif distrib1 == 'L10' and distrib2 == 'L10':
        lc = ELE_DICT[element1]['lc']['FCC']
        vacOffset = VACUUM_THICKNESS / 2
        for (i, atom) in enumerate(obj):
            yModulo = round((round(obj.positions[i][1], 3) - vacOffset) % (lc + lc / 2), 3)
            if (yModulo == lc): atom.symbol = element2
            if (yModulo == 0.0) | (yModulo == lc + lc / 2): atom.symbol = element3
    elif distrib2 == 'L10':
        lc = ELE_DICT[obj[0].symbol]['lc']['FCC']
        vacOffset = VACUUM_THICKNESS / 2
        for (i, atom) in enumerate(obj):
            yModulo = round((round(obj.positions[i][1], 3) - vacOffset) % lc, 3)
            if (yModulo == 0.0) | (yModulo == lc): atom.symbol = element3
    else: raise Exception('Specified distribution type unrecognised!')

    if 'R' in distrib2: obj = randConv(obj, element1, element2, element3, ele1Ratio, ele2Ratio, ele3Ratio, rseed, probList)
    return obj


def writeTNP(element1, element2, element3, diameter, shape, ele1Ratio, ele2Ratio, ele3Ratio, distrib1, distrib2, rep1=0, replace=False, vis=False):
    if not isdir(LMP_DATA_DIR): mkdir(LMP_DATA_DIR)
    if not isdir(f"{LMP_DATA_DIR}/{TNP_DIR}"): mkdir(f"{LMP_DATA_DIR}/{TNP_DIR}")

    # Get input file name and read data from previous generated file
    if distrib1 == 'L10': bnpRatio1, bnpRatio2, rep1, dirName = 50, 50, '', BNP_DISTRIB_LIST[0]
    else: bnpRatio1, bnpRatio2, dirName = 100 - ele2Ratio, ele2Ratio, BNP_DISTRIB_LIST[1]
    if distrib1 == 'L10' and distrib2 == 'L10':
        fileNameMNP = f"{element1}{diameter}{shape}.lmp"
        bnp = read_lammps_data(f"{LMP_DATA_DIR}/{MNP_DIR}/{fileNameMNP}", atom_style='atomic', units='metal')
        bnp.set_chemical_symbols(symbols=[element1] * len(bnp))
    else:
        fileNameBNP = f"{element1}{element2}{diameter}{shape}{bnpRatio1}{bnpRatio2}{distrib1}{rep1}.lmp"
        bnp = read_lammps_data(f"{LMP_DATA_DIR}/{BNP_DIR}/{dirName}/{fileNameBNP}", atom_style='atomic', units='metal')
        bnp.set_chemical_symbols(symbols=[element1 if bnp.arrays['type'][i] == 1 else element2 for i in range(bnp.arrays['type'].size)])

    # Generate the new file name
    if distrib1 == 'L10' and distrib2 == 'RAL': dirName = TNP_DISTRIB_LIST[0]
    elif distrib1 == 'RAL' and distrib2 == 'RAL': dirName = TNP_DISTRIB_LIST[1]
    for rep2 in range(RANDOM_DISTRIB_NO):
        if distrib1 == 'L10' and distrib2 == 'L10':
            ele1Ratio, ele2Ratio, ele3Ratio, rep2, dirName = '', '', '', '', TNP_DISTRIB_LIST[2]
        fileNameTNP = f"{element1}{element2}{element3}{diameter}{shape}{ele1Ratio}{ele2Ratio}{ele3Ratio}{distrib1}{rep1}{distrib2}{rep2}.lmp"
        if not replace:
            if isfile(f"{LMP_DATA_DIR}/{TNP_DIR}/{dirName}/{fileNameTNP}"):
                print(f"          {fileNameTNP} already exist, skipping...")
                return

        tnp = genTNP(bnp, element1, element2, element3, ele1Ratio, ele2Ratio, ele3Ratio, distrib1, distrib2, rep2)
        if not isdir(f"{LMP_DATA_DIR}/{TNP_DIR}/{dirName}"): mkdir(f"{LMP_DATA_DIR}/{TNP_DIR}/{dirName}")
        write_lammps_data(f"{LMP_DATA_DIR}/{TNP_DIR}/{dirName}/{fileNameTNP}", atoms=tnp, units='metal', atom_style='atomic')
        print(f"          Generated {fileNameTNP}, formula: {tnp.get_chemical_formula()}")
        if vis: view(tnp)


def main(replace=False, vis=False):
    print('Generating TNP alloys of:')
    for diameter in DIAMETER_LIST:
        if diameter != 30: continue
        print(f"\n  Size {diameter} Angstrom for:")
        for element1 in ELE_DICT:
            for element2 in ELE_DICT:
                if element1 is element2: continue
                for element3 in ELE_DICT:
                    if element3 is element1 or element3 is element2: continue
                    print(f"    Element 1: {element1}, Element 2: {element2}, Element 3: {element3}")
                    for shape in SHAPE_LIST:
                        for ratio1 in RATIO_LIST:
                            for ratio2 in RATIO_LIST:
                                ratio3 = 100 - ratio1 - ratio2
                                if ratio3 <= 0 or sum((ratio1, ratio2, ratio3)) != 100: continue
                                print(f"      Ratio1: {ratio1}, Ratio 2: {ratio2}, Ratio 3: {ratio3}")
                                for distrib1 in BNP_DISTRIB_LIST:
                                    for distrib2 in BNP_DISTRIB_LIST:
                                        for rep1 in range(RANDOM_DISTRIB_NO):
                                            # if distrib1 == 'L10' and distrib2 == 'L10': continue
                                            if distrib1 == 'RAL' and distrib2 == 'L10': continue
                                            print(f"        Distrib 1: {distrib1}, Distrib 2: {distrib2}")
                                            writeTNP(
                                                element1=element1,
                                                element2=element2,
                                                element3=element3,
                                                diameter=diameter,
                                                shape=shape,
                                                ele1Ratio=ratio1,
                                                ele2Ratio=ratio2,
                                                ele3Ratio=ratio3,
                                                distrib1=distrib1,
                                                distrib2=distrib2,
                                                rep1=rep1,
                                                replace=replace,
                                                vis=vis
                                            )


if __name__ == '__main__':
    main(replace=False, vis=False)
    print('ALL DONE!')
