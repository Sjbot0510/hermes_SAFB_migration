# MIGRATION_PLAN.md

## Phase 1: Foundation
- [x] Task 1.1: Create project skeleton (src/, include/, tests/, Makefile)
- [ ] Task 1.2: Translate `domain.h` / `domain.c` — all data structures
- [ ] Task 1.3: Translate `symmetry_ops.h` — `read_spacegroup_ops_txt()`
- [ ] Task 1.4: Translate `symmetry_ops.h` — fraction math utilities (Fraction struct, frac_mod1, is_zero_mod1_vec, frac_to_float)
- [ ] Task 1.5: Translate `symmetry_ops.h` — `direct_basis_from_lattice_info()` and `lattice_params_from_basis()`

## Phase 2: Symmetry Core
- [ ] Task 2.1: Translate `symmetry_ops.h` — `unique_rotations()`, `star_from_hkl()`, `get_family_key_lexicographical()`
- [ ] Task 2.2: Translate `symmetry_ops.h` — `star_is_closed()`, `point_group_has_neg_identity()`, `find_inversion_ops()`
- [ ] Task 2.3: Translate `symmetry_ops.h` — `relationships_in_star()` — star phase constraint graph
- [ ] Task 2.4: Translate `symmetry_ops.h` — `solve_star_coeffs()` — BFS coefficient solver on star graph
- [ ] Task 2.5: Translate `symmetry_ops.h` — `family_planes_info()` — generate valid Miller index modes

## Phase 3: Basis Construction
- [ ] Task 3.1: Translate `space_group_plane_family.h` — `build_basis()` into `basis.c`
- [ ] Task 3.2: Test basis construction against Python output for simple space group (e.g., P1, P2_1)

## Phase 4: Initialization
- [ ] Task 4.1: Translate `initializers.h` — `BaseInitializer._build_result()` → `build_initialization_result()`
- [ ] Task 4.2: Translate `initializers.h` — `RandomInitializer` → random amplitude assignment
- [ ] Task 4.3: Translate `initializers.h` — `ManualInitializer` → manual amplitude assignment
- [ ] Task 4.4: Translate `initializers.h` — `FileInitializer` + `read_scattering_data()` → file-based amplitude

## Phase 5: Field Generation
- [ ] Task 5.1: Translate `field.h` — reciprocal grid construction (fftfreq, meshgrid equivalents)
- [ ] Task 5.2: Translate `field.h` — coefficient placement on reciprocal grid
- [ ] Task 5.3: Translate `field.h` — iFFT via FFTW, field normalization (tanh normalization)
- [ ] Task 5.4: Translate `field.h` — VTK `.vts` XML file writer

## Phase 6: Engine & Analysis
- [x] Task 6.1: Translate `engine.c` — `SpaceGroupInitializationEngine` as procedural functions
- [x] Task 6.2: Translate `Analytic.py` — `calculate_square_norm()` numerical version
- [x] Task 6.3: Translate `Analytic.py` — `derive_analytical_star_function()` numerical version
- [x] Task 6.4: Translate `Analytic.py` — `extract_basis()` helper + `analyze_star()` helper functions

## Phase 7: Polish
- [x] Task 7.1: Full end-to-end test: generate field for Ia-3d (gyroid) and compare with Python VTK output
- [x] Task 7.2: Benchmark comparison (runtime, memory)
- [x] Task 7.3: Documentation, examples, Makefile refinement

## Dependencies

```
Phase 1 (foundation)
  └─> Phase 2 (symmetry core)
        └─> Phase 3 (basis construction)
              └─> Phase 4 (initialization)
                    └─> Phase 5 (field generation)
                          └─> Phase 6 (engine + analysis)
                                └─> Phase 7 (polish)
```

Each phase must pass validation against Python before proceeding.
