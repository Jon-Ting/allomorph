# Configuration

`allomorph` allows you to override default constants such as element properties, nanoparticle diameters, shapes, and atomic distributions.

## Configuration File

You can provide a configuration file in **JSON**, **YAML**, or **TOML** format to the CLI using the `--config` argument.

### YAML Example (`config.yaml`)

```yaml
# Global constants
VACUUM_THICKNESS: 50.0
RANDOM_DISTRIB_NO: 5

# Nanoparticle parameters
DIAMETER_LIST: [20, 30, 40]
SHAPE_LIST: ["OT", "TO", "SP", "MY_CUSTOM_CUBE"]
RATIO_LIST: [25, 50, 75]

# Custom shape definitions
CUSTOM_SHAPES:
  MY_CUSTOM_CUBE:
    builder: "FaceCenteredCubic"
    surfaces: [[1, 0, 0]]
    layers: [5]

# Element property dictionary
ELE_DICT:
  Cu:
    lc: {FCC: 3.61}
    radius: 1.28
    mass: 63.55
    rho: 8960
    m: 0.063546
    bulkE: 3.49
  Ag:
    lc: {FCC: 4.09}
    radius: 1.44
    mass: 107.87
    rho: 10490
    m: 0.107868
    bulkE: 2.95
```

### JSON Example (`config.json`)

```json
{
  "VACUUM_THICKNESS": 50.0,
  "RANDOM_DISTRIB_NO": 5,
  "DIAMETER_LIST": [20, 30, 40],
  "ELE_DICT": {
    "Cu": {
      "lc": {"FCC": 3.61},
      "radius": 1.28,
      "mass": 63.55,
      "rho": 8960,
      "m": 0.063546,
      "bulkE": 3.49
    }
  }
}
```

## CLI Usage

Pass the configuration file to any `allomorph` command:

```bash
allomorph --config config.yaml init-struct --stage mnp
```

## Python API

In Python, you can use `update_constants` or pass a custom `ele_dict` to the generation functions.

```python
from allomorph.constants import update_constants, load_config

# Load and update globally
config = load_config("config.yaml")
update_constants(config)

# Or pass to specific functions
from allomorph.init_struct.gen_mnp import main as gen_mnp_main
gen_mnp_main(ele_dict=config["ELE_DICT"])
```

## Custom Shapes

You can define custom nanoparticle shapes by adding them to the `CUSTOM_SHAPES` dictionary in your configuration and including their name in `SHAPE_LIST`.

Currently, custom shapes are built using the `ase.cluster.cubic.FaceCenteredCubic` builder. You can specify `surfaces` and `layers`.

```yaml
CUSTOM_SHAPES:
  LARGE_CUBE:
    builder: "FaceCenteredCubic"
    surfaces: [[1, 0, 0]]
    layers: [10]
```

## Parallel Generation

Structure generation (MNP, BNP, TNP) now automatically utilizes multiple CPU cores via the `multiprocessing` module to speed up the creation of large datasets. No extra configuration is needed.

## Required Element Properties

When defining a custom element in `ELE_DICT`, the following keys are required:

| Key | Description | Unit |
|-----|-------------|------|
| `lc` | Lattice constants (must contain `FCC` entry) | Å |
| `radius` | Metallic radius | Å |
| `mass` | Atomic mass | g/mol |
| `rho` | Density (required for feature extraction) | kg/m³ |
| `m` | Molar mass (required for feature extraction) | kg/mol |
| `bulkE` | Bulk cohesive energy (required for feature extraction) | eV/atom |
