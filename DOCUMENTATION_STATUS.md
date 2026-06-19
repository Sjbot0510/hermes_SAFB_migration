# SAFB C Library — Documentation Status

## Papers Read

- [ ] **No formal papers yet** — The SAFB project is derived from the PSCF (Self-Consistent Field Theory) framework for block copolymer simulations. Key references include:
  - Glotzer, S. C., & Solomon, M. J. (2007). "Anisotropy of building blocks and their assembly into complex structures." *Chemical Reviews*, 107(8), 2891-2950.
  - Matsen, M. W. (2002). "The standard thermodynamic model of BCP self-assembly." *Journal of Physics: Condensed Matter*, 14(2), R21.
  - Wang, Z. G., & Edwards, S. F. (1993). "Phase transitions in block copolymer melts with cubic symmetry." *Macromolecules*, 26(22), 5993-6000. (Ia-3d gyroid paper)
  - International Tables for Crystallography (volume A: Space-Group Symmetry) — reference for space group operations and symmetry operations parsing.

## Code Modules Inspected

- [x] `domain.py` — Data structures (LatticeInfo, Star, SAFBBasis, InitializationResult, ScatteringProfile)
- [x] `symmetry_ops.py` — Space group operations parsing, star generation, coefficient solver, metric math
- [x] `space_group_plane_family.py` — Basis construction from symmetry operations
- [x] `initializers.py` — Amplitude assignment (Random, Manual, File-based)
- [x] `field.py` — Real-space field generation (iFFT via FFTW) + VTK export
- [x] `engine.py` — High-level API (SpaceGroupInitializationEngine)
- [x] `Analytic.py` — Analytical calculations (square norm, star function derivation, basis extraction)

## Theory-to-Code Mapping

### Symmetry-Adapted Fourier Basis

The SAFB project generates scalar fields that obey crystallographic space group symmetry. The core mathematical concept is:

```
ψ(r) = Σ_n A_n · exp(i · q_n · r)
```

where q_n are reciprocal lattice vectors that form "stars" (sets of symmetry-equivalent vectors), and A_n are complex coefficients constrained by symmetry phase relationships.

| Concept | Theory | Python Function | C Function |
|---------|--------|-----------------|------------|
| Space group ops | International Tables for Crystallography | `read_spacegroup_ops_txt()` | `read_spacegroup_ops()` |
| Star generation | Symmetry-equivalent q-vectors | `star_from_hkl()` | `generate_star()` |
| Star closure | Star contains all symmetry equivalents | `star_is_closed()` | `star_is_closed()` |
| Phase constraints | c_q = exp(i·t·q) · c_{Rq} | `relationships_in_star()` | `compute_star_relationships()` |
| Coefficient solver | BFS on star relationship graph | `solve_star_coeffs()` | `solve_star_coeffs()` |
| Basis assembly | Collect all valid stars | `build_basis()` | `build_basis()` |
| Random init | Random amplitudes + phase constraints | `random_initialization()` | `random_initialization()` |
| Field generation | iFFT of symmetry-adapted coefficients | `build_field()` | `generate_field()` |
| Square norm | <|ψ|²> = (1/V) ∫ |ψ|² dr | `calculate_square_norm()` | `calculate_square_norm()` |

### Star Generation Algorithm

```
1. Start from a Miller index (h, k, l)
2. Apply all space group rotation matrices R to get {R·(h,k,l)}
3. Include negatives: {-(h,k,l)} and {-R·(h,k,l)}
4. Sort by canonical key (lexicographic)
5. Verify closure: all symmetry-equivalent vectors present
6. Compute phase relationships from translation parts of ops
```

Python: `star_from_hkl(hkl, rotations, Ginv)` → returns Star with vectors, q², multiplicity, relationships
C: `generate_star(hkl, Rs, Ginv, dim)` → returns Star struct

### Coefficient Phase Constraints

For a star with vectors {q₁, q₂, ...}, symmetry requires:
```
c_{Rq} = exp(2πi · t · q) · c_q
```
where R is a rotation and t is the fractional translation from the space group operation.

This creates a system of linear equations solved via BFS on the star's relationship graph.

Python: `solve_star_coeffs(star, rels, ref_real=1.0)` → Dict[Tuple, complex]
C: `solve_star_coeffs(star, rels, ref_real, &ncoeffs)` → complex*

### Field Generation Pipeline

```
1. Build basis (stars + coefficients) from space group + lattice
2. Place coefficients on reciprocal grid (freqf_int)
3. Apply 3D inverse FFT (FFTW)
4. Normalize with tanh to get real-space field in [0, 1]
5. Export to VTK .vts for visualization
```

Python: `engine.init_random()` → `engine.build_field()` → `engine.Output_field()`
C: `engine_init_random()` → `engine_build_field()` → `engine_output_field()`

### Square Norm Calculation

```
<|ψ|²> = (1 / (2π)³) · ∫₀²π ∫₀²π ∫₀²π |ψ(X,Y,Z)|² dX dY dZ
```

Python: SymPy symbolic integration (`sp.integrate`)
C: Numerical quadrature on Nx×Ny×Nz grid (default 64³)

### Analytical Star Function

For closed stars (star contains all negatives):
```
f(X,Y,Z) = Σ_j c_j · cos(h_j·X + k_j·Y + l_j·Z)
```

For open stars:
```
f_even(X,Y,Z) = Σ_j |c_j| · cos(φ_j + h_j·X + k_j·Y + l_j·Z)
f_odd(X,Y,Z)  = Σ_j |c_j| · sin(φ_j + h_j·X + k_j·Y + l_j·Z)
```

with √2 factor for proper normalization.

---

## Documentation Progress

### MIGRATION_REVIEW_STATUS.md
- [x] Created — Full module-by-module review with C↔Python mapping
- [x] Test results documented (119 passing, 1 warning)
- [x] Issues cataloged (P1: 2D ops LFS, P2: sorting warning, P3: compiler warnings)
- [x] Benchmark results included

### MIGRATION_STATUS.md
- [x] Up to date — All phases complete (7.3 Documentation polish)
- [x] 119 tests passing across 9 test suites
- [x] Benchmark results included

### MIGRATION_PLAN.md
- [x] All 7 phases complete
- [ ] Could be updated to show completion status with checkboxes

### PYTHON_REFERENCE_MAP.md
- [x] Created — Maps every Python function/class to C equivalent
- [x] Includes data structure mapping, function mapping, complexity notes
- [x] Notes on what has no C equivalent (SymPy-only functions)

### VALIDATION_PLAN.md
- [x] Created — Validation strategy with numerical tolerances
- [x] 6 test cases defined (lattice round-trip, ops parsing, star gen, coeff solver, e2e, gyroid)
- [x] Tolerances defined (1e-10 for metric, 1e-8 for phase factors, 1e-5 for field values)

### CHANGELOG_AGENT.md
- [x] Created — Session-by-session change log
- [x] Sessions 4-5 documented with implementation details

### User Guide
- [x] examples/INTRODUCTION.md — User-facing guide with use cases, getting started, examples

---

## Pending Documentation Work

1. **Formal theory documentation** — Create LaTeX-based documentation similar to DROPS:
   - Chapter 1: Introduction to SAFB
   - Chapter 2: Crystallographic Background (space groups, stars, reciprocal lattice)
   - Chapter 3: Mathematical Formulation (Fourier series, phase constraints, normalization)
   - Chapter 4: Numerical Methods (FFT, quadrature, coefficient solver)
   - Chapter 5: Library Architecture
   - Chapter 6: User Manual
   - Chapter 7: Examples and Workflows

2. **Read relevant papers** for theory-aligned documentation:
   - Wang & Edwards (1993) — Ia-3d gyroid phase
   - Matsen (2002) — Standard SCFT model
   - International Tables for Crystallography — Space group symmetry

3. **Generate figures**:
   - Architecture diagram (Python → C mapping)
   - Star generation visualization
   - Field morphology examples (gyroid, bcc, etc.)
   - Benchmark performance chart

4. **Add glossary** for crystallographic terms (space group, star, multiplicity, Miller index, etc.)

---

## Summary

- **MIGRATION_REVIEW_STATUS.md**: ✅ Created — comprehensive module review
- **MIGRATION_STATUS.md**: ✅ Up to date
- **MIGRATION_PLAN.md**: ✅ All phases complete
- **PYTHON_REFERENCE_MAP.md**: ✅ Complete function mapping
- **VALIDATION_PLAN.md**: ✅ Validation strategy documented
- **CHANGELOG_AGENT.md**: ✅ Session log maintained
- **User Guide**: ✅ examples/INTRODUCTION.md exists
- **Theory-aligned documentation**: ⏳ Pending (papers not yet read, LaTeX not yet generated)
