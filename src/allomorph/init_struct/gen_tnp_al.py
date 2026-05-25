# Goal: Generate initial alloy TNP structures for MD simulations
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
    TNP_DIR,
    VACUUM_THICKNESS,
    calc_rcs_prob,
)


def rand_conv(obj, element1, element2, element3, ele1Ratio, ele2Ratio, ele3Ratio, rseed, prob):
    """Randomly convert elements of atoms until specified ratio is reached."""
    ele1Arr, ele2Arr, ele3Arr = obj.symbols.search(element1), obj.symbols.search(element2), obj.symbols.search(element3)
    ele1IdealNum, ele2IdealNum, ele3IdealNum = round(ele1Ratio / 100 * len(obj)), round(ele2Ratio / 100 * len(obj)), round(ele3Ratio / 100 * len(obj))
    diff1, diff2, diff3 = len(ele1Arr) - ele1IdealNum, len(ele2Arr) - ele2IdealNum, len(ele3Arr) - ele3IdealNum

    randGen = RandomState(rseed)
    # Put excessive element into temporary array
    idxArr = np.array([], dtype=int)
    if diff1 > 0:
        p = _prob_array(prob, ele1Arr)
        idxArr = np.concatenate((idxArr, randGen.choice(a=ele1Arr, size=abs(diff1), replace=False, p=p)))
    if diff2 > 0:
        p = _prob_array(prob, ele2Arr)
        idxArr = np.concatenate((idxArr, randGen.choice(a=ele2Arr, size=abs(diff2), replace=False, p=p)))
    if diff3 > 0:
        p = _prob_array(prob, ele3Arr)
        idxArr = np.concatenate((idxArr, randGen.choice(a=ele3Arr, size=abs(diff3), replace=False, p=p)))
    # Assign atoms from this array to element with insufficient amount
    if diff1 < 0:
        if abs(diff1) > len(idxArr):
            diff1 = len(idxArr)
        p = _prob_array(prob, idxArr)
        idxArr1 = randGen.choice(a=idxArr, size=abs(diff1), replace=False, p=p)
        for idx in idxArr1:
            obj[idx].symbol = element1
        idxArr = np.array([idx for idx in idxArr if idx not in idxArr1])
    if diff2 < 0:
        if abs(diff2) > len(idxArr):
            diff2 = len(idxArr)
        p = _prob_array(prob, idxArr)
        idxArr2 = randGen.choice(a=idxArr, size=abs(diff2), replace=False, p=p)
        for idx in idxArr2:
            obj[idx].symbol = element2
        idxArr = np.array([idx for idx in idxArr if idx not in idxArr2])
    if diff3 < 0:
        if abs(diff3) > len(idxArr):
            diff3 = len(idxArr)
        p = _prob_array(prob, idxArr)
        idxArr3 = randGen.choice(a=idxArr, size=abs(diff3), replace=False, p=p)
        for idx in idxArr3:
            obj[idx].symbol = element3
    return obj


def _prob_array(prob, indices):
    """Return a normalised probability array for *indices* or None."""
    if prob is None or len(prob) == 0:
        return None
    arr = np.array(prob)[indices]
    s = arr.sum()
    if s > 0 and np.isfinite(s):
        p = arr / s
        if np.all(np.isfinite(p)):
            return p
    return None


def gen_tnp(obj, element1, element2, element3, ele1Ratio, ele2Ratio, ele3Ratio, distrib1, distrib2, rseed=0, shape=None, ele_dict=None):
    if ele_dict is None:
        ele_dict = ELE_DICT
    prob_list = []
    if distrib2 == 'RAL':
        pass
    elif distrib1 == 'L10' and distrib2 == 'L10':
        lc = ele_dict[element1]['lc']['FCC']
        vac_offset = VACUUM_THICKNESS / 2
        for i, atom in enumerate(obj):
            y_modulo = round((round(obj.positions[i][1], 3) - vac_offset) % (lc + lc / 2), 3)
            if y_modulo == lc:
                atom.symbol = element2
            if (y_modulo == 0.0) or (y_modulo == lc + lc / 2):
                atom.symbol = element3
    elif distrib2 == 'L10':
        lc = ele_dict.get(obj[0].symbol, ele_dict[element1])['lc']['FCC']
        vac_offset = VACUUM_THICKNESS / 2
        for i, atom in enumerate(obj):
            y_modulo = round((round(obj.positions[i][1], 3) - vac_offset) % lc, 3)
            if (y_modulo == 0.0) or (y_modulo == lc):
                atom.symbol = element3
    elif distrib2 == 'RCS':
        if shape is None:
            raise ValueError("shape must be provided when distrib2 is 'RCS'")
        prob_list = calc_rcs_prob(obj, shape)
        seed(rseed)
        rand_list = rand(len(obj))
        for atom_idx in range(len(obj)):
            if rand_list[atom_idx] < prob_list[atom_idx]:
                obj[atom_idx].symbol = element3
    else:
        raise ValueError(f"Specified distribution type '{distrib2}' unrecognised!")

    if 'R' in distrib2:
        obj = rand_conv(obj, element1, element2, element3, ele1Ratio, ele2Ratio, ele3Ratio, rseed, prob_list)
    return obj


def _map_tnp_distribution(distrib1, distrib2):
    """Map a (distrib1, distrib2) pair to the canonical TNP distribution name."""
    mapping = {
        ('L10', 'RAL'): 'L10R',
        ('L10', 'L10'): 'LL10',
        ('L10', 'RCS'): 'CL10S',
        ('RAL', 'RAL'): 'RRAL',
        ('RAL', 'RCS'): 'CRALS',
        ('RAL', 'L10'): 'CRSR',
        ('RCS', 'RAL'): 'CSRAL',
        ('RCS', 'L10'): 'CSL10',
        ('RCS', 'RCS'): 'CS',
    }
    key = (distrib1, distrib2)
    if key in mapping:
        return mapping[key]
    raise ValueError(
        f"No TNP distribution mapping for ({distrib1}, {distrib2})."
    )


def write_tnp(element1, element2, element3, diameter, shape, ele1Ratio, ele2Ratio, ele3Ratio, distrib1, distrib2, rep1=0, replace=False, vis=False, ele_dict=None):
    """Generates and writes TNP alloy structures."""
    if ele_dict is None:
        ele_dict = ELE_DICT
    output_base_dir = Path(LMP_DATA_DIR) / TNP_DIR
    mnp_dir = Path(LMP_DATA_DIR) / MNP_DIR
    bnp_dir = Path(LMP_DATA_DIR) / BNP_DIR

    # Get input file name and read data from previous generated file
    if distrib1 == 'L10':
        bnpRatio1, bnpRatio2, rep1_val, bnp_dist_name = 50, 50, '', 'L10'
    else:
        bnpRatio1, bnpRatio2, rep1_val, bnp_dist_name = 100 - ele2Ratio, ele2Ratio, rep1, distrib1

    if distrib1 == 'L10' and distrib2 == 'L10':
        file_name_mnp = f"{element1}{diameter}{shape}.lmp"
        mnp_path = mnp_dir / file_name_mnp
        if not mnp_path.exists():
            return
        bnp = read_lammps_data(str(mnp_path), atom_style='atomic', units='metal')
        bnp.set_chemical_symbols(symbols=[element1] * len(bnp))
    else:
        file_name_bnp = f"{element1}{element2}{diameter}{shape}{bnpRatio1}{bnpRatio2}{distrib1}{rep1_val}.lmp"
        bnp_path = bnp_dir / bnp_dist_name / file_name_bnp
        if not bnp_path.exists():
            return
        bnp = read_lammps_data(str(bnp_path), atom_style='atomic', units='metal')
        # Map atom types back to symbols for processing
        symbols = []
        for i in range(len(bnp)):
            atype = bnp.arrays['type'][i]
            symbols.append(element1 if atype == 1 else element2)
        bnp.set_chemical_symbols(symbols)

    tnp_dist_name = _map_tnp_distribution(distrib1, distrib2)
    distrib_dir = output_base_dir / tnp_dist_name
    distrib_dir.mkdir(parents=True, exist_ok=True)

    for rep2 in range(RANDOM_DISTRIB_NO):
        if tnp_dist_name == 'LL10':
            ele1_r, ele2_r, ele3_r, rep1_s, rep2_s = '', '', '', '', ''
        else:
            ele1_r, ele2_r, ele3_r, rep1_s, rep2_s = ele1Ratio, ele2Ratio, ele3Ratio, rep1_val, rep2

        file_name_tnp = f"{element1}{element2}{element3}{diameter}{shape}{ele1_r}{ele2_r}{ele3_r}{distrib1}{rep1_s}{distrib2}{rep2_s}.lmp"
        output_path = distrib_dir / file_name_tnp

        if not replace and output_path.exists():
            continue

        tnp = gen_tnp(bnp.copy(), element1, element2, element3, ele1Ratio, ele2Ratio, ele3Ratio, distrib1, distrib2, rep2, shape=shape, ele_dict=ele_dict)
        write_lammps_data(str(output_path), atoms=tnp, units='metal', atom_style='atomic')
        print(f"          Generated {file_name_tnp}, formula: {tnp.get_chemical_formula()}")

        if vis:
            view(tnp, block=True)
        if distrib1 == 'L10' and distrib2 == 'L10':
            break


def main(replace=False, vis=False, ele_dict=None):
    """Main entry point for TNP generation."""
    if ele_dict is None:
        ele_dict = ELE_DICT
    print('Generating TNP alloys:')
    
    work_items = []
    for diameter in DIAMETER_LIST:
        for element1 in ele_dict:
            for element2 in ele_dict:
                if element1 == element2:
                    continue
                for element3 in ele_dict:
                    if element3 == element1 or element3 == element2:
                        continue
                    for shape in SHAPE_LIST:
                        for ratio1 in RATIO_LIST:
                            for ratio2 in RATIO_LIST:
                                ratio3 = 100 - ratio1 - ratio2
                                if ratio3 <= 0 or sum((ratio1, ratio2, ratio3)) != 100:
                                    continue
                                for distrib1 in BNP_DISTRIB_LIST:
                                    for distrib2 in BNP_DISTRIB_LIST:
                                        for rep1 in range(RANDOM_DISTRIB_NO):
                                            work_items.append((element1, element2, element3, diameter, shape, ratio1, ratio2, ratio3, distrib1, distrib2, rep1, replace, vis, ele_dict))

    # Run in parallel unless visualization is requested
    if vis:
        print("  Visualization enabled: running serially...")
        for item in work_items:
            write_tnp(*item)
    elif len(work_items) > 1:
        with Pool() as p:
            p.starmap(write_tnp, work_items)
    elif work_items:
        write_tnp(*work_items[0])


if __name__ == '__main__':
    main(replace=False, vis=False)
    print('ALL DONE!')
