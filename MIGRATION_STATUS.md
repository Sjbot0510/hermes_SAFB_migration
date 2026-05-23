# MIGRATION_STATUS.md

## Overview

Migration of `/sandbox/Sg_init` (Python SAFB) to `/sandbox/hermes_SAFB_migration` (C).

## Progress

### Phase 1: Foundation — In Progress

| Task | Status | Notes |
|---|---|---|
| 1.1 Project skeleton | ✅ COMPLETE | `src/`, `include/`, `tests/`, `examples/`, `Makefile`, `.gitignore` created |
| 1.2 `domain.h` / `domain.c` | ✅ COMPLETE | Data structures: LatticeInfo, Star, SAFBBasis, InitializationResult, ScatteringProfile |
| 1.3 `symmetry_ops.c` `read_spacegroup_ops()` | ✅ COMPLETE | Parse PSCF-style text files; Fraction math; lattice basis |
| 1.4 Fraction math utilities | ✅ COMPLETE | Fraction struct, frac_mod1, is_zero_mod1_vec, frac_to_float, equal_int_mat |
| 1.5 Lattice basis functions | ✅ COMPLETE | direct_basis_from_lattice_info, lattice_params_from_basis |
| 2.1 `unique_rotations()`, `star_from_hkl()`, `get_family_key_lexicographical()` | ✅ COMPLETE | Star generation, canonical key selection |
| 2.2 `star_is_closed()`, `point_group_has_neg_identity()`, `find_inversion_ops()` | ✅ COMPLETE | Star closure and inversion checks |
| 2.3 `relationships_in_star()` | ✅ COMPLETE | Star phase constraint graph |
| 2.4 `solve_star_coeffs()` | ✅ COMPLETE | BFS coefficient solver on star graph |
| 2.5 `family_planes_info()` | ✅ COMPLETE | Generate valid Miller index modes (heap-allocated) |
| 3.1 `basis_build()` | ✅ COMPLETE | Assembles SAFBBasis from ops, lattice, mode count |
| 4.1 `build_result_from_coeffs()` | ✅ COMPLETE | Solve phase constraints → coefficients |
| 4.2 `random_initialization()` | ✅ COMPLETE | Random amplitude assignment with RNG |
| 4.3 `manual_initialization()` | ✅ COMPLETE | Manual amplitude assignment |
| 4.4 `file_initialization()` + `read_scattering_data()` | ✅ COMPLETE | File-based amplitude from scattering data |
| 5.1 Reciprocal grid construction (`freqf_int`) | ✅ COMPLETE | fftfreq equivalent — fixed threshold bug |
| 5.2 Coefficient placement on grid | ✅ COMPLETE | Map FamilyCoeffs onto 3D complex grid |
| 5.3 iFFT (Cooley-Tukey + DFT fallback) | ✅ COMPLETE | 3D inverse FFT, tanh normalization |
| 5.4 VTK `.vts` XML writer | ✅ COMPLETE | ASCII VTK StructuredGrid export |

### Phase 6-7: Not Started

All remaining phases listed in MIGRATION_PLAN.md are pending.

## Current Blockers

None.

## Next Recommended Task

**Task 6.1**: Translate `engine.c` — `SpaceGroupInitializationEngine` as procedural functions.

This is the high-level API class that ties all modules together (basis building, initialization, field generation, VTK output).

## Current Session Context

- Session: Initial setup (Session 1)
- Python project fully inspected
- All 6 core modules analyzed
- Documentation created: AGENTS.md, MIGRATION_PLAN.md, PYTHON_REFERENCE_MAP.md, VALIDATION_PLAN.md
- C project skeleton created with Makefile build system
