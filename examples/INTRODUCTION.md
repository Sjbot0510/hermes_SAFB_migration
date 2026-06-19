# SAFB C Library — User Guide

## What Is This?

This library generates **real-space scalar fields** that strictly obey **crystallographic
space group symmetry** constraints. It is a complete C port of a verified Python project
(`Sg_init`), translated function-by-function for bit-exact numerical agreement.

### Typical Use Case

You are running a **Self-Consistent Field Theory (SCFT)** simulation of a polymer or
block-copolymer material. You need a symmetric initial guess for the order parameter
(chemical potential, density distribution) — one that already respects the target crystal
symmetry. This library does exactly that.

### Example Applications

| Space Group | Crystal Structure | Example Use Case |
|---|---|---|
| **Ia-3d** | Gyroid (bicontinuous cubic) | Gyroid polymer phases, zeolites |
| **Pm-3m** | Simple cubic | BCC-like morphologies |
| **P42/mmc** | Tetragonal | Tetragonal block-copolymer phases |
| **p4** | 2D square | Quasi-2D confined polymer films |
| **Fm-3m** | Face-centered cubic | FCC-like morphologies |

---

## Getting Started

### 1. Build Environment

The project uses a conda-managed toolchain. **Always** initialize it before building:

```bash
source /sandbox/setup_build.sh   # sets $CC, $CXX, loads conda env
```

> **Never** call `gcc` or `g++` directly. Always use `$CC` and `$CXX` so the conda
> compiler toolchain is used (which bundles FFTW and the correct standard library).

### 2. Build

```bash
make                            # builds all test executables
source /sandbox/setup_build.sh && make test   # build + run all tests
```

Tests compile one executable per module (e.g., `./engine`, `./e2e`). The `make test`
target runs all of them and reports passing/failing counts.

### 3. Dependencies

| Dependency | Why | Install |
|---|---|---|
| **FFTW3** | 3D inverse FFT for field generation | `conda install -c conda-forge fftw` |
| **libm** | Math library (`sqrt`, `cos`, etc.) | Built into glibc |
| **C11 compiler** | Language standard | `conda install -c conda-forge gcc` |

---

## Architecture Overview

```
┌──────────────────────────────────────────────────────────┐
│                     Engine API (engine.h)                 │
│  ─ One-stop interface for the full workflow              │
│                                                          │
│  engine_create()              Load space group + lattice  │
│  engine_build_basis()          Generate Fourier modes     │
│  engine_random_init()          Random amplitude init      │
│  engine_manual_init()          Specify amplitudes by hand │
│  engine_output_field()         iFFT → VTK file           │
│  engine_full_pipeline()        All of the above in one call│
└──────────────────┬───────────────────────────────────────┘
                   │
  ┌────────────────┼────────────────┐
  ▼                ▼                ▼
┌──────┐      ┌──────────┐    ┌────────┐
│basis │      │initializers│   │ field │
│.c .h │      │.c .h     │    │.c .h  │
│      │      │          │    │       │
│Build │      │Solve     │    │iFFT   │
│SAFB  │      │phase     │    │+VTK   │
│basis │      │constraints│   │write  │
└──────┘      └──────────┘    └────────┘
       ▲                ▲              ▲
       │                │              │
  ┌──────────┐    ┌──────────┐   ┌─────────┐
  │symmetry_ │    │ domain   │   │analytic │
  │ops.c .h  │    │ .c .h    │   │.c .h    │
  │          │    │          │   │         │
  │Parse ops │    │Data      │   │Square   │
  │Star gen  │    │structures│   │norm     │
  │Phase solver│  │          │   │Point eval│
  └──────────┘    └──────────┘   └─────────┘
```

### Module Responsibilities

| Module | Header / Source | What It Does |
|---|---|---|
| **domain** | `include/domain.h` / `src/domain.c` | Data structures: `LatticeInfo`, `Star`, `SAFBBasis`, `FullInitializationResult` |
| **symmetry_ops** | `include/symmetry_ops.h` / `src/symmetry_ops.c` | Parse space group ops files, generate stars, solve phase constraints, family planes |
| **basis** | `include/basis.h` / `src/basis.c` | Build the SAFB basis from ops + lattice + mode count |
| **initializers** | `include/initializers.h` / `src/initializers.c` | Assign amplitudes (random, manual, or from scattering data), solve phase constraints |
| **field** | `include/field.h` / `src/field.c` | 3D inverse FFT via FFTW, VTK XML writer, coefficient placement on grid |
| **engine** | `include/engine.h` / `src/engine.c` | High-level API — orchestrates all lower-level modules |
| **analytic** | `include/analytic.h` / `src/analytic.c` | Analytical field math: square norm integration, point evaluation, star function derivation |

---

## Core Workflow (Step by Step)

Here is the typical end-to-end workflow, shown with the modular API:

```c
#include "engine.h"

/* Step 1: Define your crystal lattice */
LatticeInfo lattice = lattice_info_new(
    4.0, 4.0, 4.0,        /* a, b, c (Å) */
    90.0, 90.0, 90.0,     /* alpha, beta, gamma (degrees) */
    3                      /* dimension */
);

/* Step 2: Create the engine — load space group ops file */
EngineContext ctx;
engine_create(&ctx, "Ia-3d",
    "examples/space_groups_3d/Cubic/I_a_-3_d.txt",
    &lattice,
    5        /* N: generate top 5 modes by magnitude */
);

/* Step 3: Choose initialization method */

/* Option A: Random (deterministic via seed) */
FullInitializationResult result;
engine_random_init(&ctx, &result);

/* Option B: Manual — specify amplitudes for each mode */
const char *keys[] = {"{222}", "{400}", "{420}", "{440}", "{531}"};
double amps[]    = {1.0, 1.0, 1.0, 1.0, 1.0};
engine_manual_init(&ctx, (const char *const *)keys, amps, 5, &result);

/* Option C: From experimental scattering data */
// engine_file_init(&ctx, "examples/scattering_data/scattering_C14.txt", &lattice, NULL, &result);

/* Step 4: Generate real-space field and write VTK */
engine_output_field(
    "output.vts",        /* output file */
    "psi",               /* field name */
    0,                   /* no tiling */
    1, 1, 1,
    &result,             /* initialization result */
    1.0,                 /* resolution (lattice_dim / resol = grid points) */
    0                    /* no lattice transform */
);

/* Step 5: Clean up */
engine_free(&ctx);
```

### One-Liner: `engine_full_pipeline()`

For the common case (random init → field output), there is a single function:

```c
engine_full_pipeline(
    "Ia-3d",                           /* space group symbol */
    "examples/space_groups_3d/Cubic/I_a_-3_d.txt",  /* ops file */
    &(LatticeInfo){4.0, 4.0, 4.0, 90.0, 90.0, 90.0, 3},
    5,                                 /* N: number of modes */
    DIST_UNIFORM,                      /* random distribution type */
    (DistParams){0.0, 1.0},            /* distribution params */
    42,                                /* RNG seed */
    "gyroid_field.vts",                /* output file */
    "psi",                             /* field name */
    1.0,                               /* resolution */
    0, 1, 1, 1                         /* no tiling */
);
```

---

## Data Types Quick Reference

### `LatticeInfo`

```c
typedef struct {
    double a, b, c;       /* unit cell edge lengths */
    double alpha, beta, gamma;  /* angles in degrees */
    int    dim;            /* 2 or 3 */
} LatticeInfo;
```

### `SAFBBasis`

```c
typedef struct {
    char space_group[32];
    int  centrosymmetric_group;   /* is -I in point group? */
    int  has_inversion_at_origin; /* inversion through (0,0,0)? */
    LatticeInfo lattice;
    Star modes[MAX_MODES];         /* valid Fourier modes */
    int modes_count;
} SAFBBasis;
```

Each `Star` has:
- `family_key`: canonical label like `"{222}"`
- `q2`: squared magnitude of reciprocal vectors
- `star_vectors[N][3]`: all symmetry-equivalent (h,k,l) vectors
- `multiplicity`: number of vectors in the star
- `star_close`: whether the star contains its own inverse (-G)

### `FullInitializationResult`

Holds everything after amplitude assignment:
- Space group name and lattice
- Coefficients for each family key (hkl → complex value)
- Amplitudes and keys
- Discretization grid dimensions (`Na`, `Nb`, `Nc`)
- Optional lattice transform coefficients (for cross-lattice work)

---

## Common Pitfalls

### 1. Forgetting to `source setup_build.sh`

If `$CC` is not set, `make` will fail with "gcc: command not found". Always run:

```bash
source /sandbox/setup_build.sh
```

### 2. Wrong Resolution

The `resol` parameter controls grid density: `N = round(lattice_dim / resol)`.
- `resol = 1.0` → ~1 grid point per lattice unit (coarse, fine enough for most visualizations)
- `resol = 0.5` → ~2 points per unit cell (smoother)
- `resol = 2.0` → ~0.5 points per unit cell (too coarse, may alias)

If you get a grid size of 0 or 1, pick a smaller `resol`.

### 3. Family Key Format

When using `engine_manual_init()`, keys must use **curly braces**: `"{222}"`, not `"222"`.
The C code produces keys in `{hkl}` format via `format_family_key()`.

### 4. Space Group File Paths

The examples directory mirrors the Python project structure:
```
examples/space_groups_3d/{Crystal System}/{Space Group Name}.txt
examples/space_groups_2d/{Lattice Type}/{Space Group Name}.txt
```

If you have a custom ops file, place it anywhere and pass the full path.

### 5. Tiling

Tiling copies the field in each dimension. Useful for simulating larger supercells:
```c
engine_output_field("tiled.vts", "psi", 1, 2, 2, 2, &result, 1.0, 0);
// Produces a 2x2x2 tiling of the base field
```

---

## Available Space Group Files

### 3D Cubic

| File | Space Group | Ops |
|---|---|---|
| `Cubic/I_a_-3_d.txt` | Ia-3d (gyroid) | 96 |
| `Cubic/F_m_-3_m.txt` | Fm-3m | 48 |
| `Cubic/P_m_-3_m.txt` | Pm-3m | 48 |
| `Cubic/I_m_-3_m.txt` | Im-3m | 96 |
| `Cubic/I_41_3_2.txt` | I4₁32 | 96 |

### 3D Tetragonal

| File | Space Group | Ops |
|---|---|---|
| `Tetragonal/P_42_m_m_c.txt` | P42/mmc | 16 |
| `Tetragonal/P_42_m_n_m.txt` | P42/mnm | 16 |

### 2D Square

| File | Space Group | Ops |
|---|---|---|
| `square/p_4.txt` | p4 | 4 |
| `square/p_4_m_m.txt` | p4mm | 8 |
| `square/p_4_g_m.txt` | p4gm | 8 |

### 2D Hexagonal

| File | Space Group | Ops |
|---|---|---|
| `Hexagonal/p_6_m_m.txt` | p6mm | 12 |
| `Hexagonal/p_6.txt` | p6 | 6 |

### 2D Other

| File | Space Group | Ops |
|---|---|---|
| `Oblique/p_2.txt` | p2 | 2 |
| `Rectangular_oc/c_m.txt` | cm | 2 |

---

## Verifying Your Results

The project includes a comprehensive test suite:

```bash
source /sandbox/setup_build.sh && make test
```

- **9 test suites**, **119 tests**, all passing
- `e2e` — end-to-end validation against Python reference (Ia-3d, Pm-3m, 2D p4)
- `python_compare` — numerical comparison with the Python implementation
- `analytic` — square norm verification on isolated Fourier terms
- `engine` — all engine API entry points

Each test suite exercises specific functionality. Run an individual suite:
```bash
./engine          # engine layer tests
./e2e             # end-to-end tests
./analytic        # analytical field tests
```

---

## Analytical Utilities

Beyond field generation, the library provides analytical tools in `analytic.h`:

```c
#include "analytic.h"

/* Compute <|f|²> = mean squared value of a Fourier series */
double norm = calculate_square_norm(&field, 64);  /* 64³ quadrature grid */

/* Evaluate f(X,Y,Z) at a single point */
double val;
evaluate_analytic_real(&field, 0.5, 0.5, 0.5, &val);

/* Derive analytical fields from initialization results */
AnalyticBasis ab;
derive_analytical_star_function(&result, 64, &ab);
/* ab.stars[i].sq_norm gives the squared norm of each star */

/* Numerical coefficient cleanup */
AnalyticField cleaned;
extract_basis(&noisy_field, &cleaned);
```

These are useful for:
- Computing field energetics (square norm is proportional to free energy)
- Evaluating order parameters at specific points
- Validating that your initialization produces the expected field energy

---

## VTK Output

All field output is in **VTK XML Structured Grid (.vts)** format, directly readable by:

- **ParaView** — import → apply → adjust color range → visualize gyroid/Double-Diamond/etc.
- **VisIt** — same, drag and drop
- **Custom parsers** — the output is plain-text XML, easy to read programmatically

The field values are **normalized to [0, 1]** via a tanh transform. The field name
defaults to `"psi"` but can be customized.

---

## Performance Notes

Benchmarks from the validation suite (256³ grid, 3 runs, average):

| Space Group | Modes | Total Time | RSS |
|---|---|---|---|
| P4 (2D) | 5 | ~2 ms | 10 MB |
| Pm-3m | 3 | ~1.5 ms | 10 MB |
| Ia-3d | 5 | ~186 ms | 10 MB |
| Fm-3m | 5 | ~98 ms | 10 MB |
| Ia-3d | 15 | ~894 ms | 10 MB |

The dominant cost is **basis construction** for high-symmetry groups (star generation
and phase constraint solving). FFTW field generation is consistently <0.4ms regardless
of space group. Memory usage stays flat at ~10MB for all tested configurations.

---

## Further Reading

- `MIGRATION_PLAN.md` — phased implementation plan
- `MIGRATION_STATUS.md` — detailed progress tracker
- `PYTHON_REFERENCE_MAP.md` — Python-to-C function mapping
- `VALIDATION_PLAN.md` — test strategy and numerical tolerances
- `CHANGELOG_AGENT.md` — session-by-session change history

## Examples

See the `examples/` directory for compilable demos:

| Demo | What It Shows |
|---|---|
| `demo_manual_init.c` | Manual amplitudes for Ia-3d gyroid structure |
| `demo_random_init.c` | One-liner random init → VTK output |
| `demo_2d.c` | 2D square lattice (p4) — quasi-2D polymers |
| `demo_analytic.c` | Analytical field math: square norm, point evaluation |
