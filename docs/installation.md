# Installation

To install `tnp-gen`, follow these steps:

## Prerequisites

- Python 3.9 or higher
- [ASE (Atomic Simulation Environment)](https://wiki.fysik.dtu.dk/ase/)
- [LAMMPS](https://www.lammps.org/) (for MD simulations)
- [NCPac](https://github.com/ncpac/ncpac) (for feature extraction)

## Local Installation

Clone the repository and install it in editable mode:

```bash
git clone https://github.com/jonathanyct/tnp-gen.git
cd tnp-gen
pip install -e .
```

To install development dependencies (for testing and documentation):

```bash
pip install -e ".[dev]"
```

## HPC (NCI Gadi) Setup

If you are using NCI Gadi, you can load the necessary modules and set up a virtual environment:

```bash
module load python3/3.9.2
module load lammps
# Set up your venv as usual
```
