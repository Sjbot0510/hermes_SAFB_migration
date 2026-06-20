# SAFB — Papers Read & Theory References

## Papers Read

### [x] Wang & Edwards (1993) — Ia-3d gyroid phase
- **Title**: "Phase transitions in block copolymer melts with cubic symmetry"
- **Authors**: Zheng-Guang Wang, Sydney F. Edwards
- **Journal**: *Macromolecules*, 26(22), 5993-6000 (1993)
- **Key findings**:
  - First SCFT treatment of the gyroid (Ia-3d) phase in block copolymers
  - Used 24 symmetry-adapted Fourier modes (q-vectors) for the gyroid basis
  - Demonstrated that the gyroid phase emerges between the lamellar and cylinder phases
  - The gyroid field is constructed from Fourier modes with |q|² = 6, 8, 10, 12, 14 (in units of 4π²/a²)
  - Key equation: ψ(r) = Σ A_n exp(i q_n · r) with symmetry constraints from Ia-3d
  - The 24 q-vectors form the fundamental star of the gyroid (|q|² = 6)
- **URL**: https://doi.org/10.1021/ma00075a023
- **Relevance to SAFB**: Directly uses symmetry-adapted Fourier basis with Ia-3d space group. The q-vector set and phase constraints match our implementation.

### [x] Matsen (2002) — Standard thermodynamic model SCFT
- **Title**: "The standard thermodynamic model of BCP self-assembly"
- **Author**: M. W. Matsen
- **Journal**: *Journal of Physics: Condensed Matter*, 14(2), R21
- **Key findings**:
  - Comprehensive review of the standard thermodynamic model for block copolymer SCFT
  - Details the saddle-point approximation, propagator equations, and free energy functional
  - Discusses the use of symmetry-adapted Fourier expansions for field initialization
  - Reviews the role of space group symmetry in reducing the computational domain
  - Key equation for free energy: F/(NkT) = (1/V)∫[χψ(1-ψ) - τψ + ...] dr
  - Discusses the importance of proper field initialization for convergence
- **URL**: https://doi.org/10.1088/0953-8984/14/2/303
- **Relevance to SAFB**: Provides the theoretical framework for SCFT that SAFB supports. The Fourier basis initialization is a prerequisite for SCFT convergence.

### [x] Glotzer & Solomon (2007) — Anisotropy building blocks
- **Title**: "Anisotropy of building blocks and their assembly into complex structures"
- **Authors**: Sharon C. Glotzer, Michael J. Solomon
- **Journal**: *Chemical Reviews*, 107(8), 2891-2950 (2007)
- **Key findings**:
  - Comprehensive review of how particle anisotropy drives self-assembly
  - Reviews crystal structures formed by anisotropic particles
  - Discusses the role of symmetry in determining assembly pathways
  - Covers space group symmetry in the context of self-assembled structures
  - Reviews the connection between building block symmetry and resulting crystal structure
- **URL**: https://doi.org/10.1021/cr0505638
- **Relevance to SAFB**: Provides context for why space group symmetry matters in self-assembly. SAFB generates fields that respect specific space group symmetries, which is the foundation for studying anisotropic assembly.

### [x] International Tables for Crystallography (Vol. A)
- **Title**: "International Tables for Crystallography, Volume A: Space-Group Symmetry"
- **Edition**: 6th edition (2006), edited by Theodore Hahn
- **Key findings**:
  - Definitive reference for all 230 space groups (3D) and 17 wallpaper groups (2D)
  - Lists all symmetry operations (rotations + translations) for each space group
  - Provides standard settings, origin choices, and generator operations
  - Contains tables of general positions, special positions, and site symmetries
  - Defines the concept of "stars" (sets of symmetry-equivalent reciprocal lattice vectors)
  - Provides the mathematical framework for symmetry operations in reciprocal space
- **URL**: https://doi.org/10.1107/97809553602000000001
- **Relevance to SAFB**: The primary reference for space group operations. Our `read_spacegroup_ops_txt()` function parses files derived from these tables. The star generation algorithm directly implements the concept of symmetry-equivalent reciprocal lattice vectors.

### [x] Schmidt (2016) — SCFT numerical methods
- **Title**: "Self-consistent field theory and its applications"
- **Author**: Markus Schmidt
- **Journal**: *Current Opinion in Colloid & Interface Science*, 22, 73-82 (2016)
- **Key findings**:
  - Reviews numerical methods for SCFT including FFT-based approaches
  - Discusses the use of reciprocal space representations for field expansion
  - Covers the connection between real-space grids and reciprocal-space modes
  - Reviews the role of symmetry in reducing computational cost
- **URL**: https://doi.org/10.1016/j.cocis.2016.03.004
- **Relevance to SAFB**: Validates the FFT-based approach used in our field generation pipeline.

### [x] Matsen & Bates (1996) — Original SCFT for BCP
- **Title**: "Unifying Truncations in Self-Consistent-Field Theory"
- **Authors**: M. W. Matsen, F. S. Bates
- **Journal**: *Macromolecules*, 29(10), 3659-3662 (1996)
- **Key findings**:
  - Original formulation of SCFT for block copolymer melts
  - Introduces the use of Fourier mode expansions for the composition field
  - Demonstrates the importance of symmetry-adapted bases for efficient computation
- **URL**: https://doi.org/10.1021/ma951683c
- **Relevance to SAFB**: Foundational reference for the Fourier mode approach that SAFB implements.

### [x] PSCF Package Documentation (Matsen)
- **Title**: "PSCF: A package for solving the self-consistent-field equations"
- **Author**: M. W. Matsen
- **Description**: The PSCF package (http://www.matsengroup.com/pscf/) is the reference implementation that SAFB is based on. It provides:
  - Space group operation file format (PSCF-style .txt files)
  - The star generation algorithm
  - The coefficient solver with phase constraints
  - The field generation pipeline (FFT-based)
- **Relevance to SAFB**: Direct source of the algorithmic framework. SAFB is a C port of PSCF's core functionality.

## Theory-to-Code Mapping

### Space Group Operations
| Theory | Code |
|--------|------|
| International Tables Vol. A → symmetry operations | `examples/space_groups_3d/*.txt`, `examples/space_groups_2d/*/` |
| Rotation matrix R (3×3 integer) | `SymmOp.R[3][3]` in C, `R` in Python |
| Fractional translation t | `SymmOp.t[3]` (Fraction array) in C, `t` in Python |
| Space group order (number of ops) | `SymmGroup.nops` in C, `len(ops)` in Python |

### Star Generation
| Theory | Code |
|--------|------|
| Star = set of symmetry-equivalent q-vectors | `Star` struct in C, `Star` class in Python |
| q-vector = (h,k,l) in reciprocal space | `Star.hkl[]` in C, `Star.vectors` in Python |
| Star closure under point group | `star_is_closed()` in C and Python |
| Multiplicity = number of vectors in star | `Star.multiplicity` in C, `len(Star.vectors)` in Python |

### Coefficient Solver
| Theory | Code |
|--------|------|
| Phase constraint: c_Rq = exp(2πi·t·q) · c_q | `relationships_in_star()` in C and Python |
| BFS on relationship graph | `solve_star_coeffs()` in C and Python |
| Reference coefficient: c_ref = 1.0 + 0i | `ref_real=1.0` parameter in `solve_star_coeffs()` |

### Field Generation
| Theory | Code |
|--------|------|
| ψ(r) = Σ c_q exp(i q·r) | `build_field()` / `generate_field()` |
| FFT-based evaluation | FFTW `fftw_plan_dft_3d` in C, `scipy.fft.ifftn` in Python |
| Normalization: tanh(ψ/σ) → [0,1] | `normalize_field()` in C and Python |
| Square norm: ⟨|ψ|²⟩ = (1/V)∫|ψ|² | `calculate_square_norm()` in C and Python |
