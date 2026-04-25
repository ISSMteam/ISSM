# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What is ISSM

ISSM (Ice-sheet and Sea-level System Model) is a large-scale thermo-mechanical 2D/3D parallelized multi-purpose finite-element software for ice sheet and sea-level modeling. It is written in C++ (computational core) with MATLAB, Python, and JavaScript interfaces.

## Environment Setup

Before building or running anything, set `ISSM_DIR` in your shell profile (`.bashrc`/`.zshrc`) pointing to the repository root, then source the environment script:

```sh
export ISSM_DIR=/path/to/ISSM
source $ISSM_DIR/etc/environment.sh
```

## Build

External packages must be installed before configuring ISSM. The minimum required set for a basic build:

```sh
cd $ISSM_DIR/externalpackages/triangle && ./install-linux.sh    # or install-mac.sh
cd $ISSM_DIR/externalpackages/m1qn3    && ./install-linux.sh
cd $ISSM_DIR/externalpackages/petsc    && ./install-3.22-linux.sh
```

Then configure and build:

```sh
source $ISSM_DIR/etc/environment.sh
autoreconf -ivf
./configure.sh        # local config script (adjust paths as needed), or run ./configure directly
make -j$(nproc)
make install
```

Key `configure` flags: `--prefix=$ISSM_DIR`, `--with-matlab-dir`, `--with-python`, `--with-petsc-dir`, `--with-triangle-dir`, `--with-m1qn3-dir`, `--enable-debugging`.

## Running Tests

Tests live in `test/NightlyRun/`. Each test is a numbered script (`test101.m` / `test101.py`).

**Python** (from `test/NightlyRun/`):
```sh
export PYTHONPATH="$ISSM_DIR/src/m/dev:$PYTHONPATH"
export PYTHONSTARTUP="$ISSM_DIR/src/m/dev/devpath.py"
cd $ISSM_DIR/test/NightlyRun
./runme.py                        # run all nightly tests
./runme.py -i 101 102             # run specific tests by ID
./runme.py -i SquareShelf         # run by (partial) name
./runme.py --benchmark nightly    # benchmark filter
```

**MATLAB** (from within MATLAB):
```matlab
cd([getenv('ISSM_DIR') '/test/NightlyRun'])
addpath([getenv('ISSM_DIR') '/src/m/dev']); devpath;
runme                          % run all
runme('id', [101 102])         % run specific tests
runme('id', 102, 'procedure', 'update')  % update reference archive (developers only)
```

To update a test's reference archive (after an intentional result change), use `procedure='update'` (MATLAB) or `--procedure update` (Python).

## Code Architecture

### Dual-layer design

ISSM has two layers that work together:

1. **High-level interfaces** (`src/m/`) — MATLAB (`.m`), Python (`.py`), and JavaScript (`.js`) code for building and parameterizing models, running simulations, and post-processing results. The key object is `model` (defined in `src/m/classes/model.m` / `model.py`), which holds all simulation fields as properties (mesh, geometry, materials, boundary conditions, solver settings, results, etc.).

2. **C++ computational core** (`src/c/`) — compiled finite-element parallel engine that does the actual solving.

The high-level interface generates input files (`.bin`, `.queue`, and `.toolkits`) that are read by the computational core. In turn, the results from the simulation are saved in an `.outbin` file that is read by the High-level interface and added to the model (saved in `md.results`)

### Typical model workflow

```
triangle/mesh → setmask → parameterize → setflowequation → solve → results
```

Each step corresponds to functions in `src/m/parameterization/` and `src/m/solve/`. `parameterize()` runs a user-supplied `.par` file that fills the `model` object fields. `solve()` marshals model data to binary, calls the C++ executable (`bin/issm.exe`), and loads results back into `md.results`.

### C++ core layout (`src/c/`)

- **`analyses/`** — One `*Analysis` class per physics type (e.g., `StressbalanceAnalysis`, `ThermalAnalysis`). Each implements the abstract `Analysis` interface: element matrix/vector assembly (`CreateKMatrix`, `CreatePVector`), node/constraint creation, solution update.
- **`cores/`** — Top-level solution entry points (e.g., `stressbalance_core.cpp`, `transient_core.cpp`). These orchestrate which analyses to run and in what order.
- **`solutionsequences/`** — Linear/nonlinear/adjoint solvers that call PETSc (via `toolkits/`).
- **`classes/`** — C++ representations of FEM objects: `Elements/` (Tria, Penta, etc.), `Nodes/`, `Constraints/`, `Loads/`, `Inputs/`, `Params/`, `Materials/`.
- **`modules/`** — Compiled callable modules exposed as MEX/Python wrappers (e.g., mesh generation, interpolation, partitioning).
- **`toolkits/`** — Abstraction layer over PETSc, MPI, MUMPS, METIS for linear algebra and distributed computing.
- **`datastructures/`** — `DataSet` container and `Object` base class used throughout the core.

### Wrappers (`src/wrappers/`)

Glue code that compiles C++ modules as shared libraries loadable from MATLAB (`*_matlab.la`) and Python (`*_python.la`). The `io/` subdirectory handles binary serialization of the `model` object (marshalling) for communication between the interface and the executable.

### External packages (`externalpackages/`)

Each subdirectory has its own `install-linux.sh` / `install-mac.sh` / etc. ISSM only needs a handful of external packages installed depending on the desired configuration. The key dependencies are: 
- **PETSc** (includes MPI/MPICH, BLAS/LAPACK, MUMPS, METIS/ParMETIS, ScaLAPACK)
- **Triangle** (mesh generation),
- **m1qn3** (L-BFGS optimizer for inversions).

Some optional that can be useful depending on the application:
- **Dakota** (UQ/sampling),
- **CoDiPack** (automatic differentiation)

### Path setup for interfaces

- **Python**: `src/m/dev/devpath.py` walks `src/m/` and adds all directories containing `.py` files to `sys.path`, plus `$ISSM_DIR/lib` and `$ISSM_DIR/src/wrappers/python/.libs`.
- **MATLAB**: `src/m/dev/devpath.m` does the equivalent using `addpath` recursively.
