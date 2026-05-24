# Usage

`tnp-gen` provides both a command-line interface (CLI) and a Python API.

## Command-Line Interface

The CLI is the easiest way to run the full workflow.

### Initial Structure Generation

Generate monometallic, bimetallic (alloy and core-shell), and trimetallic nanoparticles:

```bash
tnp-gen init-struct --stage all
```

You can also specify a single stage:
- `mnp`: Monometallic
- `bnp`: Bimetallic alloy
- `cs`: Bimetallic core-shell
- `tnp`: Trimetallic alloy

### EAM Potential Generation

Create alloy potential files for LAMMPS:

```bash
tnp-gen eam --names Au Pt Pd
```

### MD Simulation Management

Set up and manage LAMMPS simulations:

```bash
# Setup directories
tnp-gen md-sim setup --init-dir InitStruct/TNP --target-dir MDsim_runs

# Generate job list
tnp-gen md-sim gen-joblist --stage 1 --sim-dir MDsim_runs
```

## Python API

You can also use the library directly in your Python scripts.

```python
from tnp_gen.init_struct.gen_mnp import gen_mnp
from tnp_gen.constants import ELE_DICT

# Generate an Au octahedron
atoms = gen_mnp(shape='OT', diameter=30, element='Au', latConst=ELE_DICT['Au']['lc']['FCC'])
print(f"Number of atoms: {len(atoms)}")
```

For more detailed API information, see the [API Reference](api.md).
