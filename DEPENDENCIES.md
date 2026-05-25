# Dependencies and External Tools

This document outlines the requirements for various parts of the AlloMorph project.

## 1. Core Structure Generation
The main `allomorph` package requires only a standard Python environment.
- **Python:** >= 3.9
- **Libraries:** `ase`, `numpy`, `scipy`, `pandas`
- **Installation:** `pip install allomorph` or `uv pip install -e .`

*Note: LAMMPS is NOT required to generate initial structures.*

## 2. Molecular Dynamics Simulations (`extras/md_sim` & `MDsim/`)
To run the generated models through MD simulations, you must have:
- **LAMMPS:** Recommended version `29Sep2021` or later.
- **MPI:** For parallel execution (`mpirun`).
- **Potentials:** Embedded Atom Method (EAM) files for relevant alloy systems.

## 3. High-Performance Computing (`extras/md_sim/submission.py`)
The automated job submission logic is designed for HPC environments using:
- **PBS (Portable Batch System):** e.g., Gadi (NCI) or similar clusters.
- **Commands:** `qsub`, `qstat`, `qselect`.

## 4. Feature Extraction (`extras/feat_ext_eng` & `FeatExtEng/`)
The feature extraction pipeline relies on:
- **NCPac:** A Fortran-based tool for structural analysis. Requires a Fortran compiler (e.g., `gfortran`) to build from source in `FeatExtEng/NCPac/`.
- **Pandas:** For merging and managing CSV datasets.
