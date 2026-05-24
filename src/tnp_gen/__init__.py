"""TNP generation toolkit.

A package for generating trimetallic nanoparticle (TNP) structural datasets
for machine learning applications, including:

- Initial structure generation (monometallic, bimetallic, trimetallic)
- EAM potential file creation
- Feature extraction and dataset merging
"""

__version__ = "0.1.0"

from tnp_gen.eam.eam_database import AtType, Database
from tnp_gen.eam.create_eam import create_eam
from tnp_gen.constants import (
    LMP_DATA_DIR,
    MNP_DIR,
    BNP_DIR,
    TNP_DIR,
    GOLDEN_RATIO,
    VACUUM_THICKNESS,
    RANDOM_DISTRIB_NO,
    ELE_DICT,
    DIAMETER_LIST,
    SHAPE_LIST,
    BNP_DISTRIB_LIST,
    TNP_DISTRIB_LIST,
    RATIO_LIST,
)
from tnp_gen.init_struct.gen_mnp import gen_mnp, write_mnp
from tnp_gen.init_struct.gen_bnp_al import gen_bnp, write_bnp
from tnp_gen.init_struct.gen_tnp_al import gen_tnp, write_tnp
from tnp_gen.init_struct.gen_bnp_cs import gen_hard_core_shell, write_hard_core_shell

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
    "gen_mnp",
    "write_mnp",
    "gen_bnp",
    "write_bnp",
    "gen_tnp",
    "write_tnp",
    "gen_hard_core_shell",
    "write_hard_core_shell",
]
