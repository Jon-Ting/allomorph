# Welcome to tnp-gen

`tnp-gen` is a toolkit for generating trimetallic nanoparticle (TNP) structural datasets for machine learning.

## Key Features

- **Initial Structure Generation**: Generate monometallic, bimetallic, and trimetallic nanoparticles with various shapes and distributions.
- **MD Simulation Orchestration**: Automatically set up and manage LAMMPS simulations for annealing nanoparticles.
- **Feature Extraction**: Extract structural features from nanoparticles using NCPac and MD outputs.
- **Modular Design**: Clean Python API and CLI for easy integration into workflows.

## Quick Start

```bash
# Generate structures
tnp-gen init-struct --stage all

# Setup MD simulations
tnp-gen md-sim setup --init-dir InitStruct/TNP --target-dir MDsim_runs
```

For more details, see the [Installation](installation.md) and [Usage](usage.md) guides.
