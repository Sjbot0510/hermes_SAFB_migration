# VALIDATION_PLAN.md

## Validation Strategy

Every C module must be validated against the Python reference by running the same input
and comparing numerical outputs within defined tolerances.

## Numerical Tolerances

| Type | Tolerance |
|---|---|
| Integer indices (h,k,l) | Exact match |
| Matrix elements (R) | Exact match (integers) |
| Metric inverse G⁻¹ | 1e-10 |
| q² values | 1e-10 |
| Phase factors | 1e-8 |
| Complex coefficients | 1e-8 magnitude, 1e-6 phase |
| Field values | 1e-5 relative |
| Star multiplicity | Exact match |
| Family keys | Exact string match |

## Test Cases

### Test 1: Lattice Info Round-Trip
**Input**: LatticeInfo(a=1.0, b=1.0, c=1.0, alpha=90, beta=90, gamma=90)
**Expected**: `direct_basis_matrix()` → identity-like matrix with columns [1,0,0], [0,1,0], [0,0,1]
**Python**: `lattice = LatticeInfo(1,1,1,90,90,90)` → `direct_basis_from_lattice_info(lattice)`
**C**: `LatticeInfo lat = lattice_info_new(1.0, 1.0, 1.0, 90.0, 90.0, 90.0); double* A = direct_basis_matrix(&lat);`
**Command**: `make test_lattice`
**Tolerance**: 1e-10

### Test 2: Space Group Op Parsing
**Input**: `symmetry_operators/3/Cubic/I_m_-3_m.txt`
**Expected**: Correct count of ops (48 for Im-3m), correct R and t matrices
**Python**: `ops = read_spacegroup_ops_txt("symmetry_operators/3/Cubic/I_m_-3_m.txt")` → 48 ops
**C**: `SymOp* ops = read_spacegroup_ops("symmetry_operators/3/Cubic/I_m_-3_m.txt", &n);` → n=48
**Command**: `make test_ops_parse`
**Tolerance**: Exact for R, 1e-10 for t

### Test 3: Star Generation
**Input**: P1 space group, N=5 modes, cubic lattice
**Expected**: First 5 stars sorted by q²: {100}, {010}, {001}, {110}, {101}
**Python**: `build_basis(ops, "P1", lattice, 5)` → 5 Star objects
**C**: `SAFBBasis* basis = build_basis(ops, "P1", &lattice, 5);` → basis->modes_count=5
**Command**: `make test_star_gen`
**Tolerance**: Exact for hkl tuples, 1e-10 for q²

### Test 4: Coefficient Solver
**Input**: Star with vectors and relationships
**Expected**: Coefficients satisfy all phase constraints c_k = ph * c_j
**Python**: `solve_star_coeffs(star, rels, ref_real=1.0)`
**C**: `complex* coeffs = solve_star_coeffs(star, rels, ref_real, &ncoeffs);`
**Command**: `make test_coeff_solver`
**Tolerance**: 1e-8

### Test 5: End-to-End Field Generation (Simple)
**Input**: P1 space group, random init, resol=0.1, grid ~10×10×10
**Expected**: Field values in [0,1], total mass preserved
**Python**: `engine.init_random()` → `engine.Output_field("test.vts", ...)`
**C**: `engine_init_random(engine, RNG_UNIFORM)` → `engine_output_field(engine, "test.vts", ...)`
**Command**: `make test_e2e_simple`
**Tolerance**: 1e-5 on field values, exact for dimensions

### Test 6: Gyroid (Ia-3d)
**Input**: Ia-3d, n_keep=10, cubic a=b=c=1, alpha=beta=gamma=90
**Expected**: 3D gyroid-like field with characteristic bicontinuous structure
**Python**: Full engine run → `Ia3d_init.vts`
**C**: Same → `Ia3d_init.vts`
**Command**: `make test_gyroid`
**Tolerance**: 1e-5 on field values

## How to Run Tests

```bash
source /sandbox/setup_build.sh
cd /sandbox/hermes_SAFB_migration
make test
```

Individual test:
```bash
make test_lattice       # Test 1
make test_ops_parse     # Test 2
make test_star_gen      # Test 3
make test_coeff_solver  # Test 4
make test_e2e_simple    # Test 5
make test_gyroid        # Test 6
```

## Validation Workflow

For each translated function:

1. Run Python reference: `python -c "import module; result = module.func(input)"`
2. Run C equivalent: `./test_<module> --func <name> --input <data>`
3. Diff outputs: `python -c "import numpy as np; a=np.loadtxt('py_out.txt'); b=np.loadtxt('c_out.txt'); print('max_diff:', np.max(np.abs(a-b)))"`
4. Log result in CHANGELOG_AGENT.md

## If Validation Fails

Document in MIGRATION_STATUS.md:
- Exact Python command and output
- Exact C command and output
- Difference analysis
- Root cause hypothesis
- Next debugging step
