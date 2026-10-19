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

Key `configure` flags: `--prefix=$ISSM_DIR`, `--with-matlab-dir`, `--with-python`, `--with-petsc-dir`, `--with-triangle-dir`, `--enable-debugging`.

### Fast single-file syntax check

A full `make` is slow. To verify one file still compiles after an edit, use the flags from the
generated `src/c/Makefile` directly:

```sh
cd $ISSM_DIR/src/c
mpicxx -std=c++11 -D_DO_NOT_LOAD_GLOBALS_ -Wno-deprecated-register -Wno-return-type \
       -I. -I.. -I../.. -DHAVE_CONFIG_H -fsyntax-only classes/FemModel.cpp
```

Add `-Wunused-variable -Wunused-but-set-variable` to hunt dead locals. Two caveats:

- `classes/Params/EmulatorParam.cpp` needs pybind11 include paths and always fails this way — not a real error.
- This does **not** check code inside optional-dependency guards (`_HAVE_AD_`, `_HAVE_ADOLC_`,
  `_HAVE_CODIPACK_`, `_HAVE_DAKOTA_`, `PETSC_VERSION_*`), since they compile out by default.

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
./runme.py -b nightly             # benchmark filter (-b/--benchmark: all/nightly/ismip/eismint/thermal/mesh/slc/qmu/...)
./runme.py -e Dakota              # exclude tests by ID or name (-e/--exclude)
./runme.py -p update              # update reference archive (-p/--procedure: check/update)
```

**MATLAB** (from within MATLAB):
```matlab
cd([getenv('ISSM_DIR') '/test/NightlyRun'])
addpath([getenv('ISSM_DIR') '/src/m/dev']); devpath;
runme                          % run all
runme('id', [101 102])         % run specific tests
runme('id', 102, 'procedure', 'update')  % update reference archive (developers only)
```

To update a test's reference archive (after an intentional result change), use `procedure='update'` (MATLAB) or `-p update` (Python).

### Inspecting a test archive

Reference values live in `test/Archives/ArchiveNNN.arch`. Field *k* is named `ArchiveNNN_fieldK`,
in the same order as `field_names` in `testNNN.m`. To read one without MATLAB:

```sh
cd $ISSM_DIR
python3 -c "
import sys; sys.path.insert(0,'src/m/archive')
from arch import archdisp, archread
archdisp('test/Archives/Archive517.arch')                        # list fields + sizes
print(archread('test/Archives/Archive517.arch','Archive517_field2'))
"
```

Data comes back flat, row-major — reshape using the size `archdisp` reports.

Check the archive's own vintage (`git log -- test/Archives/ArchiveNNN.arch`) before assuming a
mismatch is a regression: many archives predate large refactors. Archives are also
platform-sensitive, and tolerances at ~1e-11 and below can legitimately fail across arm64 macOS
vs x86 Linux for identical code.

## Code Architecture

### Dual-layer design

ISSM has two layers that work together:

1. **High-level interfaces** (`src/m/`) — MATLAB (`.m`), Python (`.py`), and JavaScript (`.js`) code for building and parameterizing models, running simulations, and post-processing results. The key object is `model` (defined in `src/m/classes/model.m` / `model.py` / `model.js`), which holds all simulation fields as properties (mesh, geometry, materials, boundary conditions, solver settings, results, etc.).

2. **C++ computational core** (`src/c/`) — compiled finite-element parallel engine that does the actual solving.

3. **JavaScript interface** (`src/m/js/`, `src/m/classes/model.js`) — browser/WebAssembly interface built via Emscripten, with wrappers in `src/wrappers/javascript/`.

The high-level interface generates input files (`.bin`, `.queue`, and `.toolkits`) that are read by the computational core. In turn, the results from the simulation are saved in an `.outbin` file that is read by the High-level interface and added to the model (saved in `md.results`)

### Typical model workflow

```
triangle/mesh → setmask → parameterize → setflowequation → solve → results
```

Each step corresponds to functions in `src/m/parameterization/` and `src/m/solve/`. `parameterize()` runs a user-supplied `.par` file that fills the `model` object fields. `solve()` marshals model data to binary, calls the C++ executable (`bin/issm.exe` for ice dynamics, `bin/issm_slc.exe` for sea-level change), and loads results back into `md.results`.

## Gotchas

**IDE/clangd diagnostics are unreliable here.** With no `compile_commands.json` in the tree, clangd
reports spurious errors — most commonly `"Cannot compile with HAVE_CONFIG_H symbol!"` (it does not
pass `-DHAVE_CONFIG_H`), plus bogus type errors such as `DataSet` vs `Nodes`/`Elements`/`Constraints`.
Trust the `mpicxx ... -fsyntax-only` command above, not the editor squiggles.

**Blanket aggregator headers are the idiom.** Most `.cpp` files include `classes.h`, `shared.h`,
`toolkits.h`, `modules.h`, each pulling in a large tree. clangd's "included header ... is not used
directly" hints on these are not actionable — removing them breaks the build. They are also the main
reason compiles are slow; fixing that properly means replacing them with precise per-file includes,
which is a large refactor, not a line deletion.

**Variables that look unused are often `#ifdef`-gated.** Before deleting an "unused" local, check
whether it is consumed inside `#ifdef _HAVE_AD_` / `_HAVE_ADOLC_` / `_HAVE_CODIPACK_` /
`MELTPERTURBATION` / `PETSC_VERSION_LT(...)`. The default build compiles those out, so the compiler
flags the variable even though AD- or PETSc-specific builds need it. Move the declaration inside the
guard instead of removing it (examples: `classes/IoModel.cpp`, `cores/ad_core.cpp`,
`cores/controlvalidation_core.cpp`, `cores/controltao_core.cpp`).

**A dead-looking variable can mark a dropped computation.** Several "set but not used" locals are
values a sibling function does use in a real formula (or that a commented-out line once consumed).
Compare against the sibling before deleting — e.g. `Tria::SealevelchangeGeometrySubElementKernel`
vs `SealevelchangeGeometryInitial`, or `Penpair::PenaltyCreateKMatrixStressbalanceFS` (missing the
`_assert_(numdof==numdof2)` its SSA/HO counterpart has).
