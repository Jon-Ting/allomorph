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

from math import cos, radians, sin, sqrt
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from ase.cluster import Decahedron, Icosahedron, Octahedron
from ase.cluster.cubic import FaceCenteredCubic
from ase.io.lammpsdata import write_lammps_data
from ase.io.xyz import write_xyz
from ase.visualize import view

from np_gen.constants import (
    DIAMETER_LIST,
    ELE_DICT,
    GOLDEN_RATIO,
    LMP_DATA_DIR,
    MNP_DIR,
    SHAPE_LIST,
    VACUUM_THICKNESS,
)


def _build_cu(diameter, element, latConst):
    estLatNo = diameter / latConst
    latNoCU = round(estLatNo)
    return FaceCenteredCubic(
        symbols=element,
        surfaces=[(1, 0, 0)],
        layers=[latNoCU],
        latticeconstant=latConst,
    )


def _build_rd(diameter, element, latConst):
    estLatNo = diameter / latConst
    latNoCU = round(estLatNo)
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


def _build_th(diameter, element, latConst):
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


def _build_ot(diameter, element, latConst):
    edgeLengthOT = diameter / sqrt(2)
    edgeAtomNoOT = 1 + round(edgeLengthOT / (latConst / sqrt(2)))
    return Octahedron(
        symbol=element,
        length=edgeAtomNoOT,
        cutoff=0,
        latticeconstant=latConst,
        alloy=False,
    )


def _build_to(diameter, element, latConst):
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


def _build_co(diameter, element, latConst):
    cutLayerNoCO = round(diameter / latConst)
    edgeAtomNoCO = 2 * cutLayerNoCO + 1
    return Octahedron(
        symbol=element,
        length=edgeAtomNoCO,
        cutoff=cutLayerNoCO,
        latticeconstant=latConst,
        alloy=False,
    )


def _build_dh(diameter, element, latConst):
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


def _build_ic(diameter, element, latConst):
    circRadIC = sqrt(diameter ** 2 * (1 + GOLDEN_RATIO ** 2) / (4 * GOLDEN_RATIO ** 2))
    shellNoIC = 1 + round(circRadIC / (latConst / sqrt(2)))
    return Icosahedron(
        symbol=element,
        noshells=shellNoIC,
        latticeconstant=latConst,
    )


def _build_sp(diameter, element, latConst):
    # Generate a large cube and cut it into a sphere
    mnp = _build_cu(diameter, element, latConst)
    center = mnp.get_center_of_mass()
    to_delete = [atom.index for atom in mnp if np.linalg.norm(center - atom.position) > diameter / 2]
    del mnp[to_delete]
    return mnp


SHAPE_BUILDERS = {
    'CU': _build_cu,
    'RD': _build_rd,
    'TH': _build_th,
    'OT': _build_ot,
    'TO': _build_to,
    'CO': _build_co,
    'DH': _build_dh,
    'IC': _build_ic,
    'SP': _build_sp,
}


def gen_mnp(shape, diameter, element, latConst, custom_shapes=None):
    """Generate a monometallic nanoparticle using ASE."""
    if custom_shapes and shape in custom_shapes:
        config = custom_shapes[shape]
        builder_name = config.get("builder", "FaceCenteredCubic")
        params = {k: v for k, v in config.items() if k != "builder"}
        params["symbols"] = element
        params["latticeconstant"] = latConst
        if builder_name == "FaceCenteredCubic":
            return FaceCenteredCubic(**params)
        # Add more builders if needed
        raise ValueError(f"Unknown builder '{builder_name}' for custom shape '{shape}'")

    if shape in SHAPE_BUILDERS:
        return SHAPE_BUILDERS[shape](diameter, element, latConst)

    raise ValueError(f"Unknown shape '{shape}'")


def write_mnp(element, diameter, lat_const, shape, replace=False, vis=False, custom_shapes=None):
    """Write a monometallic nanoparticle to LAMMPS data and XYZ files."""
    file_name_lmp = f"{element}{diameter}{shape}.lmp"
    output_dir = Path(LMP_DATA_DIR) / MNP_DIR
    output_path_lmp = output_dir / file_name_lmp

    if not replace and output_path_lmp.exists():
        # Using print here for consistency with original script, though logging might be better
        # print(f"      {file_name_lmp} already exists, skipping...")
        return

    # Generate and add vacuum
    try:
        mnp = gen_mnp(shape, diameter, element, lat_const, custom_shapes=custom_shapes)
    except Exception as e:
        print(f"      Error generating {file_name_lmp}: {e}")
        return

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


def main(replace=False, vis=False, ele_dict=None):
    """Generate all monometallic nanoparticles defined in constants."""
    if ele_dict is None:
        ele_dict = ELE_DICT

    # Load custom shapes from constants if present
    from np_gen.constants import CUSTOM_SHAPES
    custom_shapes = CUSTOM_SHAPES

    print(f"Generating NPs with {VACUUM_THICKNESS} Angstrom of vacuum on each dimension:")

    work_items = []
    for diameter in DIAMETER_LIST:
        for element in ele_dict:
            lat_const = ele_dict[element]['lc']['FCC']
            for shape in SHAPE_LIST:
                work_items.append((element, diameter, lat_const, shape, replace, vis, custom_shapes))

    # Run in parallel
    if len(work_items) > 1:
        with Pool() as p:
            p.starmap(write_mnp, work_items)
    elif work_items:
        write_mnp(*work_items[0])


if __name__ == '__main__':
    main(replace=False, vis=False)
    print('ALL DONE!')
