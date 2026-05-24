# PYTHON_REFERENCE_MAP.md

Maps Python source → C target for every function, class, and data structure.

## Data Structures (domain.py)

### LatticeInfo (dataclass, frozen=True)

| Python | C Target |
|---|---|
| `a: float` | `double a;` |
| `b: float` | `double b;` |
| `c: float` | `double c;` |
| `alpha: float` | `double alpha;` |
| `beta: float` | `double beta;` |
| `gamma: float` | `double gamma;` |
| `dim: int = 3` | `int dim;` |

Factory: `LatticeInfo.from_2d(a, b, gamma)` → `lattice_info_new_2d(a, b, gamma)`

### Star (dataclass)

| Python | C Target |
|---|---|
| `family_key: str` | `char family_key[32];` |
| `q2: float` | `double q2;` |
| `star_close: bool` | `int star_close;` |
| `reason: str` | `const char* reason;` (literal) |
| `star_vectors: List[HKL]` | `int (*star_vectors)[3]; int star_vectors_count;` |
| `multiplicity: int` | `int multiplicity;` |
| `rels: Dict[Tuple, List]` | `/* TODO: graph adjacency representation */` |

### SAFBBasis (dataclass)

| Python | C Target |
|---|---|
| `space_group: str` | `char space_group[32];` |
| `centrosymmetric_group: bool` | `int centrosymmetric_group;` |
| `has_inversion_at_origin: bool` | `int has_inversion_at_origin;` |
| `inversion_info: List[Dict]` | `/* TODO: inversion info struct array */` |
| `additional_info: Optional[str]` | `const char* additional_info;` |
| `lattice: LatticeInfo` | `LatticeInfo lattice;` |
| `modes: List[Star]` | `Star* modes; int modes_count;` |

### InitializationResult (dataclass)

| Python | C Target |
|---|---|
| `space_group: str` | `char space_group[32];` |
| `lattice: LatticeInfo` | `LatticeInfo lattice;` |
| `SAFB: SAFBBasis` | `SAFBBasis* safb;` |
| `coeff: Dict[str, Dict[Tuple, complex]]` | `/* TODO: hash table of family_key → hkl→complex */` |
| `amplitudes: Dict[str, float]` | `/* TODO: hash table of family_key → float */` |

### ScatteringProfile

| Python | C Target |
|---|---|
| `hkl: np.ndarray` (N×3 int) | `int32_t* hkl; int num_peaks;` |
| `q: np.ndarray` (N float) | `double* q;` |
| `intensity: np.ndarray` (N float) | `double* intensity;` |

## symmetry_ops.py Functions

| Python Function | C Target | Complexity |
|---|---|---|
| `read_final_values(file_path)` | `parse_simulation_file()` | Low — file parser |
| `read_spacegroup_ops_txt(path)` | `read_spacegroup_ops()` | Medium — file parse + Fraction |
| `unique_rotations(ops)` | `unique_rotations()` | Low — array dedup |
| `star_from_hkl(hkl, rotations, Ginv)` | `generate_star()` | Low — matrix mult |
| `get_family_key_lexicographical(star_vectors)` | `canonical_family_key()` | Low — sorting |
| `frac_mod1(q)` | `frac_mod1()` | Low — math |
| `is_zero_mod1_vec(t)` | `is_zero_mod1_vec()` | Low — loop + frac_mod1 |
| `equal_int_mat(A, B)` | `equal_int_mat()` | Low — memcmp |
| `star_is_closed(star)` | `star_is_closed()` | Low — set lookup |
| `point_group_has_neg_identity(rotations)` | `has_neg_identity()` | Low — matrix compare |
| `find_inversion_ops(ops)` | `find_inversion_ops()` | Low — loop + negate check |
| `frac_to_float(q)` | `frac_to_float()` | Low — div |
| `phase_factor(hkl, t_frac)` | `phase_factor()` | Medium — dot + exp |
| `_metric_inverse(a,b,c,alpha,beta,gamma)` | `metric_inverse()` | Medium — matrix math |
| `_q2_metric(h,k,l,Ginv)` | `q2_metric()` | Low — quadratic form |
| `relationships_in_star(ops, star)` | `compute_star_relationships()` | High — graph + BFS |
| `family_planes_info(N, Ginv, Rs, ops, dim)` | `generate_family_planes()` | High — iterative search |
| `solve_star_coeffs(star, rels, ref_real)` | `solve_star_coeffs()` | High — BFS on graph |
| `direct_basis_from_lattice_info(info)` | `direct_basis_matrix()` | Medium — trig math |
| `lattice_params_from_basis(A)` | `lattice_params_from_basis()` | Medium — norm + acos |
| `float_to_miller_int(hB, ...)` | `float_to_miller_int()` | Medium — Fraction conversion |
| `transform_miller_between_lattices(...)` | `transform_miller_indices()` | Medium — matrix transform |

## space_group_plane_family.py Functions

| Python Function | C Target | Complexity |
|---|---|---|
| `build_basis(ops, sg_symbol, lattice, N)` | `build_basis()` | Medium — orchestrator |

## engine.py Class

### Engine Context
| Python | C Target | Status |
|---|---|---|
| `SpaceGroupInitializationEngine.__init__()` | `engine_create()` | ✅ Complete |
| (none — destructor) | `engine_free()` | ✅ Complete |

### Engine Methods
| Python | C Target | Status |
|---|---|---|
| `build_basis(N)` | `engine_build_basis()` | ✅ Complete |
| `init_random(distribution, dist_params, rng)` | `engine_random_init()` | ✅ Complete |
| `init_manual(amplitudes)` | `engine_manual_init()` | ✅ Complete |
| `init_from_file(filename, read_lattice_info, P)` | `engine_file_init()` | ✅ Complete |
| `transform_lattice_coordinate(lattice_A, P, result)` | `engine_transform_miller()` | ✅ Complete |
| `Output_field(filename, field_name, apply_tile, tile, result, resol, transform_coord)` | `engine_output_field()` | ✅ Complete |

### Pipeline Convenience
| Python | C Target | Status |
|---|---|---|
| (none — convenience wrapper) | `engine_full_pipeline()` | ✅ Complete |

## initializers.py Classes

| Python | C Target |
|---|---|
| `BaseInitializer._build_result()` | `build_result_from_coeffs()` |
| `RandomInitializer.initialize()` | `random_initialization()` |
| `ManualInitializer.initialize()` | `manual_initialization()` |
| `FileInitializer.initialize()` | `file_initialization()` |
| `read_scattering_data(filepath)` | `read_scattering_profile()` |

## field.py Functions

| Python Function | C Target | Complexity |
|---|---|---|
| `write_lattice_field_to_vts(...)` | `write_vts_field()` | Medium — VTK XML writer |
| `build_field(...)` | `generate_field()` | High — FFT + normalization |

## Analytic.py Functions

| Python Function | C Target | Complexity |
|---|---|---|
| `extract_basis(expr)` | *(no C equivalent needed)* — SymPy-only |
| `calculate_square_norm(basis_expr)` | `calculate_square_norm()` | Medium — numerical integration |
| `derive_analytical_star_function(...)` | *(no C equivalent needed)* — SymPy-only |

## Notes

- `extract_basis()`, `derive_analytical_star_function()` use SymPy symbolic math — **no C equivalent needed**.
- The VTK `.vts` writer can be a simple XML text writer in C.
- NumPy arrays → dynamically allocated 1D/2D arrays in C.
- Python `dict` → simple array + linear search for small tables, or a hash map for larger ones.
- `Fraction` → custom struct with integer numerator/denominator.
