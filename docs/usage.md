# Usage

`tnp-gen` provides both a command-line interface (CLI) and a Python API.

## Command-Line Interface

The CLI is the easiest way to run the full workflow.

### Initial Structure Generation

Generate monometallic, bimetallic (alloy and core-shell), and trimetallic nanoparticles:

```bash
# Generate all stages (MNP → BNP → CS → TNP)
tnp-gen init-struct --stage all

# Or run individual stages
tnp-gen init-struct --stage mnp   # monometallic
tnp-gen init-struct --stage bnp   # bimetallic alloy
tnp-gen init-struct --stage cs    # bimetallic core-shell
tnp-gen init-struct --stage tnp   # trimetallic alloy
```

The default element set is **Au, Pt, Pd**. You can override this via the Python API (see below).

### EAM Potential Generation

Create alloy potential files for LAMMPS. The EAM module already accepts arbitrary elements:

```bash
# Trimetallic
tnp-gen eam --names Au Pt Pd

# Bimetallic
tnp-gen eam --names Cu Ag
```

### MD Simulation Management

Set up and manage LAMMPS simulations:

```bash
# Setup directories
tnp-gen md-sim setup --init-dir InitStruct/TNP --target-dir MDsim_runs

# Generate LAMMPS input files
tnp-gen md-sim gen-input --stage 0 --sim-dir MDsim --init-dir InitStruct --eam-dir EAM --template-dir templates

# Generate job list
tnp-gen md-sim gen-joblist --stage 1 --sim-dir MDsim_runs
```

### Feature Extraction

Set up and run NCPac feature extraction:

```bash
# Setup NCPac directories
tnp-gen feat-ext setup --sim-dir MDsim --target-dir NCPac --exe path/to/NCPac.exe --inp path/to/NCPac.inp

# Run NCPac in parallel
tnp-gen feat-ext run --target-dir NCPac --final-dir Features

# Merge features with MD outputs
tnp-gen feat-ext merge --md-out MDout.csv --feat-dir Features --output-dir Merged --ele-comb AuPtPd
```

## Python API

You can also use the library directly in your Python scripts.

### Monometallic nanoparticle

```python
from tnp_gen.init_struct.gen_mnp import gen_mnp
from tnp_gen.constants import ELE_DICT

# Generate an Au octahedron
atoms = gen_mnp(shape='OT', diameter=30, element='Au', latConst=ELE_DICT['Au']['lc']['FCC'])
print(f"Number of atoms: {len(atoms)}")
```

### Using a custom element dictionary

```python
from tnp_gen.init_struct.gen_mnp import main as gen_mnp_main
from tnp_gen.constants import validate_ele_dict

# Define your own elements
custom_ele_dict = {
    "Cu": {"lc": {"FCC": 3.61}, "radius": 1.28, "mass": 63.55,
           "rho": 8960, "m": 0.06355, "bulkE": 3.49},
    "Ag": {"lc": {"FCC": 4.09}, "radius": 1.44, "mass": 107.87,
           "rho": 10490, "m": 0.10787, "bulkE": 2.95},
    "Au": {"lc": {"FCC": 4.09}, "radius": 1.44, "mass": 196.97,
           "rho": 19320, "m": 0.196967, "bulkE": 3.81},
}
validate_ele_dict(custom_ele_dict)
gen_mnp_main(ele_dict=custom_ele_dict)
```

### Parsing element-combination strings

```python
from tnp_gen.constants import parse_ele_comb

elements = parse_ele_comb("AuPtPd")
# → ['Au', 'Pt', 'Pd']

elements = parse_ele_comb("CuAgAu")
# → ['Cu', 'Ag', 'Au']
```

For more detailed API information, see the [API Reference](api.md).
