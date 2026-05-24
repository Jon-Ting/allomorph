# Purpose: Store variables for generation of NPs using ASE
# Author: Jonathan Yik Chang Ting
# Date: 19/20/2020

from math import sqrt

LMP_DATA_DIR = "."
MNP_DIR = "MNP"
BNP_DIR = "BNP"
TNP_DIR = "TNP"
CS_DIR = "CS"
GOLDEN_RATIO = (1 + sqrt(5)) / 2
VACUUM_THICKNESS = 40.0
RANDOM_DISTRIB_NO = 3

# Elements of interest, their lattice parameters and metallic radii
# - Pd, Pt, Au values were obtained from N. W. Ashcroft and N. D. Mermin,
#   Solid State Physics (Holt, Rinehart, and Winston, New York, 1976.
# - The lattice constants are 3.859, 3.912, and 4.065 Angstroms, for the
#   respective FCC metals at 300 K according to W. P. Davey,
#   "Precision Measurements of the Lattice Constants of Twelve Common Metals,"
#   Physical Review, vol. 25, (6), pp. 753-761, 1925.
# - Metallic radii taken from Greenwood, Norman N.; Earnshaw, Alan (1997). 
#   Chemistry of the Elements (2nd ed.).
ELE_DICT = {
    "Pd": {"lc": {"FCC": 3.89}, "radius": 1.37, "mass": 106.42},
    "Pt": {"lc": {"FCC": 3.92}, "radius": 1.39, "mass": 195.08},
    "Au": {"lc": {"FCC": 4.09}, "radius": 1.44, "mass": 196.97},
}

DIAMETER_LIST = [10, 15, 20, 25, 30]  # NP diameters of interest (Angstrom)
SHAPE_LIST = ["OT", "SP", "IC", "CU", "DH"]  # Shapes of interest
BNP_DISTRIB_LIST = ["L10", "RAL", "RCS"]  # BNP distributions of interest
TNP_DISTRIB_LIST = ["L10R", "RRAL", "LL10"]  # TNP distributions of interest
RATIO_LIST = [20, 40, 60, 80]  # Ratios of interest
