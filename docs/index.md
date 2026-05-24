# Welcome to tnp-gen

`tnp-gen` is a toolkit for generating monometallic to trimetallic nanoparticle structural datasets for machine learning.

## Key Features

- **Initial Structure Generation**: Generate monometallic (MNP), bimetallic (BNP), and trimetallic (TNP) nanoparticles with configurable shapes and atomic distributions.
- **Configurable Elements**: Support arbitrary metallic elements via an overridable element-property dictionary.
- **MD Simulation Orchestration**: Automatically set up and manage LAMMPS simulations for annealing nanoparticles.
- **Feature Extraction**: Extract structural features from nanoparticles using NCPac and MD outputs.
- **Modular Design**: Clean Python API and CLI for easy integration into workflows.

## Quick Start

```bash
# Generate structures (all stages: MNP, BNP, CS, TNP)
tnp-gen init-struct --stage all

# Or generate only monometallic nanoparticles
tnp-gen init-struct --stage mnp

# Setup MD simulations
tnp-gen md-sim setup --init-dir InitStruct/TNP --target-dir MDsim_runs
```

For more details, see the [Installation](installation.md) and [Usage](usage.md) guides.
