# delphes-to-hdf5

## Installation

### 0. Install micromamba

see https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html

### 1. Install dependencies with micromamba

Create and activate the conda environment:
```fish
micromamba create -y -f environment.yaml
```

### 2. Setup the environment

- bash: `source setup.sh`
- fish: `source setup.fish`

### 3. Build the project

The project uses Meson build system with a `just` command runner for convenience:

```fish
# Configure the build
just setup

# Build and install
just build
```

To reconfigure the build (e.g., after changing build options):

```fish
just resetup
```

### Available commands

- `just setup` - Configure the build system
- `just resetup` - Reconfigure the build
- `just build` - Compile and install
- `just format` - Format source code with clang-format
- `just help` - Show all available commands

## Usage

After building, the `delphes-to-hdf5` executable will be installed in the
`install/bin/` directory, which is registered in PATH by setup scripts.

```bash
delphes-to-hdf5 -i <input.root> -o <output.h5>
```

## Output Format
see https://zenodo.org/records/8113516

## License

See the project repository for license information.
