# AlloMorph

[![ci-cd](https://github.com/Jon-Ting/allomorph/actions/workflows/ci-cd.yml/badge.svg)](https://github.com/Jon-Ting/allomorph/actions/workflows/ci-cd.yml)

---

## Description

AlloMorph is a Python toolkit for generating initial structural models of monometallic, bimetallic, and trimetallic nanoparticles. It is designed to create large-scale datasets for atomistic simulations and machine learning applications.

![](https://github.com/Jon-Ting/allomorph/blob/main/docs/figs/TNPs_generation_flowchart.png)

---

## Installation

This package uses [uv](https://docs.astral.sh/uv/) for dependency management and packaging. For detailed information on external tool requirements (LAMMPS, NCPac, etc.), see [DEPENDENCIES.md](DEPENDENCIES.md).

```bash
# Clone the repository
git clone https://github.com/jonathan-ting/allomorph.git
cd allomorph

# Create a virtual environment and install dependencies
uv venv
uv pip install -e ".[dev]"
```

Alternatively, you can install directly with pip:

```bash
pip install allomorph
```

---

## Usage

The core `allomorph` command provides a streamlined interface for structure generation:

```bash
# Generate a full suite of nanoparticles (MNP, BNP, TNP)
allomorph init-struct --stage all

# Generate only monometallic nanoparticles
allomorph init-struct --stage mnp

# Generate bimetallic nanoparticles with visual check (ASE GUI)
allomorph init-struct --stage bnp --vis
```

Check out the [basic demonstration](https://github.com/Jon-Ting/allomorph/blob/main/docs/demo.ipynb) notebook for further explanations and demonstrations!

---

## Features

### Core Package Structure

- `allomorph.init_struct` — Generation logic for monometallic, bimetallic, and trimetallic structures.
- `allomorph.constants` — Physical constants, element properties, and global configuration.

### Project Extras & Research Tools

While the core `allomorph` package focuses on structure generation, this repository includes additional tools and legacy resources used in the original research pipeline:

#### 1. Research Extras (`extras/`)
These are standalone Python modules for managing the full simulation lifecycle:
- `extras/eam/` — Creation of EAM alloy potential files.
- `extras/md_sim/` — Management of LAMMPS simulations and HPC job submission.
- `extras/feat_ext_eng/` — Integration with NCPac for structural feature extraction.

#### 2. Physical Resources
- **InitStruct**: LAMMPS and Bash scripts for alternative structure generation workflows.
- **EAM**: Reference databases and tools for interatomic potentials.
- **MDsim**: LAMMPS templates and scripts for simulation orchestration.
- **FeatExtEng**: Source code and configuration for NCPac.

### Degrees of Freedom

AlloMorph supports a wide range of configurable parameters:

- **Elemental composition:** Supports 1–3 arbitrary metallic elements (e.g., Au, Pt, Pd, Cu, Ni, Ag).
- **Size:** Configurable diameter range (default: 10–30 Å).
- **Shape:** CU (Cube), TH (Tetrahedron), RD (Rhombic Dodecahedron), OT (Octahedron), TO (Truncated Octahedron), CO (Cuboctahedron), DH (Decahedron), IC (Icosahedron), SP (Sphere).
- **Ratio:** Arbitrary stoichiometric ratios (e.g., 20:40:40).
- **Atomic ordering:**
    - **BNP**: L10 (Ordered), RAL (Random), RCS (Random Core-Shell).
    - **TNP**: L10R, CS, CL10S, CRALS, RRAL, CSRAL, CSL10, CRSR, LL10.

### Future Enhancements

- **Support >3 metals:** Extend logic to quadrimetallic and beyond.
- **Non-FCC lattice support:** Support for BCC, HCP, and other lattices.
- **Modern Formats:** Support for HDF5 and ASE trajectory outputs.

## Documentation

Detailed [documentations](https://allomorph.readthedocs.io/en/latest/) are hosted by `Read the Docs`.

## Contributing

`Sphractal` appreciates your enthusiasm and welcomes your expertise! 

Please check out the [contributing guidelines](https://github.com/Jon-Ting/allomorph/blob/main/CONTRIBUTING.md) and 
[code of conduct](https://github.com/Jon-Ting/allomorph/blob/main/CONDUCT.md). 
By contributing to this project, you agree to abide by its terms.

## License

The project is distributed under an [MIT License](https://github.com/Jon-Ting/allomorph/blob/main/LICENSE).
