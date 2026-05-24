"""Nanoparticle generation toolkit.

A package for generating monometallic to trimetallic nanoparticle structural
datasets for machine learning applications, including:

- Initial structure generation (monometallic, bimetallic, trimetallic)
- EAM potential file creation
- Molecular dynamics simulation pipeline
- Feature extraction and dataset merging
"""

__version__ = "0.1.0"

from np_gen.constants import (
    BNP_DIR,
    BNP_DISTRIB_LIST,
    DIAMETER_LIST,
    ELE_DICT,
    GOLDEN_RATIO,
    LMP_DATA_DIR,
    MNP_DIR,
    RANDOM_DISTRIB_NO,
    RATIO_LIST,
    SHAPE_LIST,
    TNP_DIR,
    TNP_DISTRIB_LIST,
    VACUUM_THICKNESS,
    load_ele_dict_from_file,
    parse_ele_comb,
    validate_ele_dict,
)
from np_gen.eam.create_eam import create_eam
from np_gen.eam.eam_database import AtType, Database
from np_gen.init_struct.gen_bnp_al import gen_bnp, write_bnp
from np_gen.init_struct.gen_bnp_cs import gen_hard_core_shell, write_hard_core_shell
from np_gen.init_struct.gen_mnp import gen_mnp, write_mnp
from np_gen.init_struct.gen_tnp_al import gen_tnp, write_tnp

__all__ = [
    "__version__",
    "AtType",
    "Database",
    "create_eam",
    "LMP_DATA_DIR",
    "MNP_DIR",
    "BNP_DIR",
    "TNP_DIR",
    "GOLDEN_RATIO",
    "VACUUM_THICKNESS",
    "RANDOM_DISTRIB_NO",
    "ELE_DICT",
    "DIAMETER_LIST",
    "SHAPE_LIST",
    "BNP_DISTRIB_LIST",
    "TNP_DISTRIB_LIST",
    "RATIO_LIST",
    "validate_ele_dict",
    "load_ele_dict_from_file",
    "parse_ele_comb",
    "gen_mnp",
    "write_mnp",
    "gen_bnp",
    "write_bnp",
    "gen_tnp",
    "write_tnp",
    "gen_hard_core_shell",
    "write_hard_core_shell",
]
