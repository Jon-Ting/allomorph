# Goal: Generate initial alloy BNP structures for MD simulations
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
from constants import LMP_DATA_DIR, MNP_DIR, BNP_DIR, RANDOM_DISTRIB_NO, VACUUM_THICKNESS, ELE_DICT, DIAMETER_LIST, SHAPE_LIST, BNP_DISTRIB_LIST, RATIO_LIST


def dist1D(coord1, coord2, dim):
    """Compute distance between 2 points in one of their real space coordinates"""
    return round(np.sqrt(np.sum((coord2[dim]-coord1[dim]) ** 2)), 3)


def dist3D(coord1, coord2):
    """Compute real space distance between 2 points"""
    return round(np.sqrt(np.sum((coord2 - coord1) ** 2)), 3)


def randConv(obj, element2, ele2Ratio, rseed, prob):
    """Randomly convert elements of atoms until specified ratio is reached"""
    elements = np.array(list(obj.symbols.species()))
    element1 = elements[elements != element2][0]
    ele1Arr, ele2Arr = obj.symbols.search(element1), obj.symbols.search(element2)
    ele2IdealNum = round(ele2Ratio / 100 * len(obj))
    diff = len(ele2Arr) - ele2IdealNum
    (convEleArr, targetEle) = (ele1Arr, element2) if diff < 0 else (ele2Arr, element1)
    randGen = RandomState(rseed)
    probArr = np.array(prob)[convEleArr] / np.array(prob)[convEleArr].sum() if len(prob) > 0 else None
    idxArr = randGen.choice(a=convEleArr, size=abs(diff), replace=False, p=probArr)
    for idx in idxArr: obj[idx].symbol = targetEle
    return obj


def genBNP(obj, element2, shape, ratio2, distrib, rseed):
    probList = []
    if distrib == 'RAL':
        seed(rseed)
        randList = rand(len(obj))  # Uniform distribution
        for (i, atom) in enumerate(obj):
            if randList[i] > (100 - ratio2) / 100: atom.symbol = element2

    elif distrib == 'RCS':
        if shape == 'IC':
            massCenter = obj.get_center_of_mass()
            radius = (obj.cell[0][0]-VACUUM_THICKNESS)/2
        else:
            xSlices = set([round(atom[0], 3) for atom in obj.positions])
            ySlices = set([round(atom[1], 3) for atom in obj.positions])
            zSlices = set([round(atom[2], 3) for atom in obj.positions])
            
            # Find the center of each line in the box
            zThreadDict = {(x, y): {'max': 0, 'min': 0, 'mid': []} for x in xSlices for y in ySlices}
            xThreadDict = {(y, z): {'max': 0, 'min': 0, 'mid': []} for y in ySlices for z in zSlices}
            yThreadDict = {(z, x): {'max': 0, 'min': 0, 'mid': []} for z in zSlices for x in xSlices}
            for (atomIdx, atom) in enumerate(obj):
                xCoord, yCoord, zCoord = round(atom.position[0], 3), round(atom.position[1], 3), round(atom.position[2], 3)
                zThreadDict[(xCoord, yCoord)]['mid'].append(zCoord)
                xThreadDict[(yCoord, zCoord)]['mid'].append(xCoord)
                yThreadDict[(zCoord, xCoord)]['mid'].append(yCoord)
            threadDicts = [zThreadDict, xThreadDict, yThreadDict]
            for threadDict in threadDicts:
                emptyDictKeys = []
                for (dim1, dim2) in threadDict.keys():
                    if len(threadDict[(dim1, dim2)]['mid']) == 0: emptyDictKeys.append((dim1, dim2)); continue
                    threadDict[(dim1, dim2)]['max'] = max(threadDict[(dim1, dim2)]['mid'])
                    threadDict[(dim1, dim2)]['min'] = min(threadDict[(dim1, dim2)]['mid'])
                    threadDict[(dim1, dim2)]['mid'] = (threadDict[(dim1, dim2)]['max']+threadDict[(dim1, dim2)]['min']) / 2
                for key in emptyDictKeys: del threadDict[key]
     
        # Alchemical change
        seed(rseed)
        randList = rand(len(obj))
        for (atomIdx, atom) in enumerate(obj):
            if shape == 'IC':
                prob = dist3D(massCenter, atom.position) / radius
            else:
                xCoord, yCoord, zCoord = round(atom.position[0], 3), round(atom.position[1], 3), round(atom.position[2], 3)
                zHalfLen = (zThreadDict[(xCoord, yCoord)]['max']-zThreadDict[(xCoord, yCoord)]['min']) / 2
                xHalfLen = (xThreadDict[(yCoord, zCoord)]['max']-xThreadDict[(yCoord, zCoord)]['min']) / 2
                yHalfLen = (yThreadDict[(zCoord, xCoord)]['max']-yThreadDict[(zCoord, xCoord)]['min']) / 2
                zRelPos = abs(zCoord-zThreadDict[(xCoord, yCoord)]['mid']) / zHalfLen if round(zHalfLen, 3) != 0.0 else abs(zCoord-zThreadDict[(xCoord, yCoord)]['mid'])
                xRelPos = abs(xCoord-xThreadDict[(yCoord, zCoord)]['mid']) / xHalfLen if round(xHalfLen, 3) != 0.0 else abs(xCoord-xThreadDict[(yCoord, zCoord)]['mid'])
                yRelPos = abs(yCoord-yThreadDict[(zCoord, xCoord)]['mid']) / yHalfLen if round(yHalfLen, 3) != 0.0 else abs(yCoord-yThreadDict[(zCoord, xCoord)]['mid'])
                if shape == 'DH':
                    prob = 1 if (round(zRelPos, 3) == 1.0) | (round(zHalfLen, 3) == 0.0) else zRelPos
                else:
                    if (round(zRelPos, 3) == 1.0) | (round(xRelPos, 3) == 1.0) | (round(yRelPos, 3) == 1.0) | (round(zHalfLen, 3) == 0.0) | (round(xHalfLen, 3) == 0.0) | (round(yHalfLen, 3) == 0.0): prob = 1
                    else: prob = (zRelPos+xRelPos+yRelPos) / 3
            if randList[atomIdx] < prob: obj[atomIdx].symbol = element2
            probList.append(prob)

    elif distrib in ['L10', 'L12', 'RL10', 'RL12']:
        lc = ELE_DICT[obj[0].symbol]['lc']['FCC']
        vacOffset = VACUUM_THICKNESS / 2
        for (i, atom) in enumerate(obj):
            yModulo = round((round(obj.positions[i][1], 3) - vacOffset) % lc, 3)
            if (yModulo == 0.0) | (yModulo == lc): atom.symbol = element2
            if distrib == 'L12':
                xModulo = round((round(obj.positions[i][0], 3) - vacOffset) % lc, 3)
                if (xModulo == 0.0) | (xModulo == lc): atom.symbol = element2

    else:
        raise Exception('Specified distribution type unrecognised!')

    if 'R' in distrib: obj = randConv(obj=obj, element2=element2, ele2Ratio=ratio2, rseed=rseed, prob=probList)
    return obj


def writeBNP(element1, diameter, shape, ratio2, distrib, replace=False, vis=False):
    if not isdir(LMP_DATA_DIR): mkdir(LMP_DATA_DIR)
    if not isdir(f"{LMP_DATA_DIR}/{BNP_DIR}"): mkdir(f"{LMP_DATA_DIR}/{BNP_DIR}")

    for element2 in ELE_DICT:
        if element2 is element1: continue
        fileNameMNP = f"{element1}{diameter}{shape}.lmp"
        mnp = read_lammps_data(f"{LMP_DATA_DIR}/{MNP_DIR}/{fileNameMNP}", style='atomic', units='metal')
        mnp.set_chemical_symbols(symbols=[element1] * len(mnp))
        ratio1 = 100 - ratio2
        for rep in range(RANDOM_DISTRIB_NO):
            if 'R' not in distrib: ratio1, ratio2, rep = 50, 50, ''
            fileNameBNP = f"{element1}{element2}{diameter}{shape}{ratio1}{ratio2}{distrib}{rep}.lmp"
            if not isdir(f"{LMP_DATA_DIR}/{BNP_DIR}/{distrib}"): mkdir(f"{LMP_DATA_DIR}{BNP_DIR}/{distrib}")
            if not replace:
                if isfile(f"{LMP_DATA_DIR}/{BNP_DIR}/{distrib}/{fileNameBNP}"):
                    print(f"      {fileNameBNP} already exist, skipping...")
                    continue

            bnp = genBNP(obj=mnp.copy(), element2=element2, shape=shape, ratio2=ratio2, distrib=distrib, rseed=rep)
            write_lammps_data(f"{LMP_DATA_DIR}/{BNP_DIR}/{distrib}/{fileNameBNP}", atoms=bnp, units='metal', atom_style='atomic')
            print(f"      Generated {fileNameBNP}, formula: {bnp.get_chemical_formula()}")
            if vis: view(bnp)
            if 'R' not in distrib: break


def main(replace=False, vis=False):
    print('Generating BNP alloys of:')
    for diameter in DIAMETER_LIST:
        # if diameter != 30: continue  # DEBUG
        print(f"\n  Size {diameter} Angstrom for:")
        for element in ELE_DICT:
            # if element != 'Pt':  continue  # DEBUG
            print(f"    Element: {element}")
            for shape in SHAPE_LIST:
                # if shape not in ['DH']:  continue  # DEBUG
                print(f"    Shape: {shape}")
                for ratio2 in RATIO_LIST:
                    # if ratio2 != 10:  continue  # DEBUG
                    print(f"    Element 2 Ratio: {ratio2}")
                    for distrib in BNP_DISTRIB_LIST:
                        # if (shape in ['DH', 'IC']) & ('L1' in distrib): continue
                        # if distrib != 'RCS':  continue  # DEBUG
                        writeBNP(element1=element, diameter=diameter, shape=shape, ratio2=ratio2, distrib=distrib, replace=replace, vis=vis)
    

if __name__ == '__main__':
    main(replace=False, vis=False)
    print('ALL DONE!')


