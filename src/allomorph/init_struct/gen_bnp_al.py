# Goal: Generate initial alloy BNP structures for MD simulations
# Author: Jonathan Yik Chang Ting
# Date: 22/10/2020
"""
Note:
- Abbreviations:
    - RAL = randomly distributed alloy
    - RCS = randomly distributed core-shell-like alloy
    - (R)L10 = L1_0 intermetallic alloy (with/without random component)
    - (R)L12 = L1_2 intermetallic alloy (with/without random component)
- To do:
    - FCC is currently being hard-coded for lattice constant retrieval, might need to be flexible
    - Perhaps could add a parameter to control core thickness
"""

from multiprocessing import Pool
from pathlib import Path

import numpy as np
from ase.io.lammpsdata import read_lammps_data, write_lammps_data
from ase.visualize import view
from numpy.random import RandomState, rand, seed

from allomorph.constants import (
    BNP_DIR,
    BNP_DISTRIB_LIST,
    DIAMETER_LIST,
    ELE_DICT,
    LMP_DATA_DIR,
    MNP_DIR,
    RANDOM_DISTRIB_NO,
    RATIO_LIST,
    SHAPE_LIST,
    VACUUM_THICKNESS,
    calc_rcs_prob,
)


def rand_conv(obj, element1, element2, ele2Ratio, rseed, prob):
    """Randomly convert elements of atoms until specified ratio is reached"""
    ele1Arr, ele2Arr = obj.symbols.search(element1), obj.symbols.search(element2)
    ele2IdealNum = round(ele2Ratio / 100 * len(obj))
    diff = len(ele2Arr) - ele2IdealNum
    (convEleArr, targetEle) = (ele1Arr, element2) if diff < 0 else (ele2Arr, element1)
    randGen = RandomState(rseed)
    probArr = None
    if len(prob) > 0:
        weights = np.array(prob)[convEleArr]
        total_weight = weights.sum()
        if total_weight > 0:
            probArr = weights / total_weight
    idxArr = randGen.choice(a=convEleArr, size=abs(diff), replace=False, p=probArr)
    for idx in idxArr:
        obj[idx].symbol = targetEle
    return obj


def gen_bnp(obj, element1, element2, shape, ratio2, distrib, rseed, ele_dict=None):
    if ele_dict is None:
        ele_dict = ELE_DICT
    probList = []
    if distrib == 'RAL':
        # Handled by rand_conv if 'R' in distrib
        pass

    elif distrib == 'RCS':
        probList = calc_rcs_prob(obj, shape)
        # Handled by rand_conv if 'R' in distrib

    elif distrib in ['L10', 'L12', 'RL10', 'RL12']:
        lc = ele_dict[element1]['lc']['FCC']
        vacOffset = VACUUM_THICKNESS / 2
        for (i, atom) in enumerate(obj):
            yModulo = round((round(obj.positions[i][1], 3) - vacOffset) % lc, 3)
            if (yModulo == 0.0) | (yModulo == lc):
                atom.symbol = element2
            if distrib == 'L12':
                xModulo = round((round(obj.positions[i][0], 3) - vacOffset) % lc, 3)
                if (xModulo == 0.0) | (xModulo == lc):
                    atom.symbol = element2

    else:
        raise Exception('Specified distribution type unrecognised!')

    if 'R' in distrib:
        obj = rand_conv(obj=obj, element1=element1, element2=element2, ele2Ratio=ratio2, rseed=rseed, prob=probList)
    return obj


def write_bnp(element1, element2, diameter, shape, ratio2, distrib, replace=False, vis=False, ele_dict=None):
    """Generates and writes a specific BNP alloy structure."""
    if ele_dict is None:
        ele_dict = ELE_DICT
    output_base_dir = Path(LMP_DATA_DIR) / BNP_DIR
    mnp_dir = Path(LMP_DATA_DIR) / MNP_DIR

    file_name_mnp = f"{element1}{diameter}{shape}.lmp"
    mnp_path = mnp_dir / file_name_mnp

    if not mnp_path.exists():
        # print(f"      MNP file {mnp_path} not found, skipping...")
        return

    mnp = read_lammps_data(str(mnp_path), atom_style='atomic', units='metal')
    mnp.set_chemical_symbols(symbols=[element1] * len(mnp))
    ratio1 = 100 - ratio2

    distrib_dir = output_base_dir / distrib
    distrib_dir.mkdir(parents=True, exist_ok=True)

    for rep in range(RANDOM_DISTRIB_NO):
        if 'R' not in distrib:
            r1, r2, rep_suffix = 50, 50, ''
        else:
            r1, r2, rep_suffix = ratio1, ratio2, str(rep)

        file_name_bnp = f"{element1}{element2}{diameter}{shape}{r1}{r2}{distrib}{rep_suffix}.lmp"
        output_path = distrib_dir / file_name_bnp

        if not replace and output_path.exists():
            continue

        bnp = gen_bnp(obj=mnp.copy(), element1=element1, element2=element2, shape=shape, ratio2=ratio2, distrib=distrib, rseed=rep if 'R' in distrib else 0, ele_dict=ele_dict)
        write_lammps_data(str(output_path), atoms=bnp, units='metal', atom_style='atomic')
        print(f"      Generated {file_name_bnp}, formula: {bnp.get_chemical_formula()}")

        if vis:
            view(bnp)
        if 'R' not in distrib:
            break


def main(replace=False, vis=False, ele_dict=None):
    """Main entry point for BNP generation."""
    if ele_dict is None:
        ele_dict = ELE_DICT
    print('Generating BNP alloys:')

    work_items = []
    for diameter in DIAMETER_LIST:
        for element in ele_dict:
            for element2 in ele_dict:
                if element == element2:
                    continue
                for shape in SHAPE_LIST:
                    for ratio2 in RATIO_LIST:
                        for distrib in BNP_DISTRIB_LIST:
                            work_items.append((element, element2, diameter, shape, ratio2, distrib, replace, vis, ele_dict))

    # Run in parallel unless visualization is requested
    if vis:
        print("  Visualization enabled: running serially...")
        for item in work_items:
            write_bnp(*item)
    elif len(work_items) > 1:
        with Pool() as p:
            p.starmap(write_bnp, work_items)
    elif work_items:
        write_bnp(*work_items[0])


if __name__ == '__main__':
    main(replace=False, vis=False)
    print('ALL DONE!')
