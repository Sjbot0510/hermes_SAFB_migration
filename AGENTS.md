# AGENTS.md

## Role

You are a Hermes coding agent assigned to gradually rewrite the verified Python project
`/sandbox/Sg_init` (Symmetry-Adapted Fourier Basis — SAFB) into a clean C language project
at `/sandbox/hermes_SAFB_migration`.

## Source of Truth

- **Python project**: `/sandbox/Sg_init` — this is the behavioral reference
- **Target C project**: `/sandbox/hermes_SAFB_migration` — incremental C rewrite
- **GitHub**: `https://github.com/Sjbot0510/hermes_SAFB_migration`

## Build Environment

Before ANY compilation or testing, ALWAYS run:

```bash
source /sandbox/setup_build.sh
```

This activates the conda toolchain (`$CC` = `x86_64-conda-linux-gnu-gcc`,
`$CXX` = `x86_64-conda-linux-gnu-g++`).

## Core Rules

1. **Never rewrite the whole project at once** — work in small, reviewable, validated steps.
2. **Preserve Python logic exactly** — algorithms, data structures, numerical behavior.
3. **Use standard C** (C11 preferred) — keep code simple, explicit, readable.
4. **Validate C output against Python** — write C tests that compare numerically.
5. **Do not delete or overwrite the Python project.**
6. **Do not rely only on memory** — always update repository documentation before ending session.
7. **Do not leave the repository in an unexplained state.**
8. **Use `$CC` and `$CXX`** for compilers; never call `gcc`/`g++` directly.
9. **Conda env name is `build`** (not `build-env`). Use `conda install -c conda-forge`.
10. **Conda safety**: `CONDA_NO_PLUGINS=true`, `CONDA_OVERRIDE_CUDA=""`, `solver classic`.

## Project Overview: SAFB (Symmetry-Adapted Fourier Basis)

The Python project generates scalar fields (chemical potentials, density distributions)
that strictly adhere to crystallographic space group symmetry constraints (3D or 2D).
Used for initializing Self-Consistent Field Theory (SCFT) simulations.

### Key Concepts

- **Space Group**: A group of symmetry operations (rotations + translations) describing crystal symmetry.
- **Star**: A set of symmetry-equivalent reciprocal lattice vectors (h,k,l).
- **SAFBBasis**: Collection of valid Fourier modes (stars) for a given space group.
- **Initialization**: Assign amplitudes to stars → solve phase constraints → generate real-space field via iFFT.
- **VTK Export**: Write field to `.vts` format for ParaView visualization.

### Python Modules and C Equivalents

| Python File | Purpose | C Target |
|---|---|---|
| `domain.py` | Data structures (LatticeInfo, Star, SAFBBasis, InitializationResult) | `include/domain.h`, `src/domain.c` |
| `symmetry_ops.py` | Symmetry utilities: ops parsing, star generation, coefficient solver, metric math | `src/symmetry_ops.c`, `include/symmetry_ops.h` |
| `space_group_plane_family.py` | Build SAFB basis from ops and lattice | `src/basis.c`, `include/basis.h` |
| `engine.py` | High-level API class (SpaceGroupInitializationEngine) | `src/engine.c`, `include/engine.h` |
| `initializers.py` | Amplitude assignment (Random, Manual, File-based) | `src/initializers.c`, `include/initializers.h` |
| `field.py` | Real-space field generation (iFFT) + VTK export | `src/field.c`, `include/field.h` |
| `Analytic.py` | Analytical calculations, square norm, symbolic expressions | `src/analytic.c`, `include/analytic.h` |

### Dependencies (Python → C)

| Python | C Equivalent |
|---|---|
| NumPy (arrays, linalg, fft, random) | Manual arrays + BLAS/LAPACK (or custom impl) + FFTW + custom RNG |
| SciPy (ifftn) | FFTW (`fftw3.h`) |
| SymPy (symbolic math, nsimplify, integrate) | Not needed in C — use numerical approximations |
| VTK (`.vts` export) | Custom VTK XML writer |
| fractions.Fraction | `typedef struct { int num, den; } Fraction;` |

## Maintenance Cycle

Every session:
1. Read all migration docs
2. Inspect git status
3. Pick ONE small task
4. Implement, validate, commit
5. Update all docs

## Documentation Files

- `AGENTS.md` — this file
- `MIGRATION_PLAN.md` — phased plan with tasks and dependencies
- `MIGRATION_STATUS.md` — current progress, blockers, next task
- `PYTHON_REFERENCE_MAP.md` — Python→C mapping for every function
- `VALIDATION_PLAN.md` — test strategy, numerical tolerances, run commands
- `CHANGELOG_AGENT.md` — session-by-session change log
- `MIGRATION_REPORTS.txt` — 4-hour status reports

## Git Commit Format

- `migration: translate <module/function>`
- `validation: add test for <feature>`
- `docs: update migration status`
- `fix: correct <specific issue>`
