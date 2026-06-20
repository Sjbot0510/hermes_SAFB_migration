# Reference Materials for SAFB Documentation

## Primary References

### 1. Foundations of Crystallography with Computer Applications (3rd Edition, 2025)
**Authors:** Maureen M. Julian, Carla Slebodnick, Francis T. Julian  
**Publisher:** CRC Press (Taylor & Francis)  
**Pages:** 525  
**Path:** `/sandbox/hermes_SAFB_migration/ref/textbook_crystollography.pdf`

### 2. Symmetry-Adapted Fourier Series for the Wallpaper Groups (2012)
**Author:** Bart Verberck  
**Journal:** Symmetry, 4, 379-426  
**DOI:** 10.3390/sym4030379  
**Open Access:** Yes (CC BY 3.0)  
**Pages:** 48  
**Path:** `/sandbox/hermes_SAFB_migration/ref/Symmetry_Adapted_Fourier_Series_for_the_Wallpaper_Groups.pdf`

### 3. PSCF++ User Manual (2025)
**Author:** David M. Ceperley  
**URL:** https://dmorse.github.io/pscfpp-man/  
**Status:** Primary reference -- fetch from URL when needed.

---

## Key Theory from Crystallography Textbook

### Chapter 1: Lattices
- **2D Lattices:** 5 Bravais lattices (oblique, rectangular primitive, rectangular centered, square, hexagonal)
- **Basis Vectors:** Defined as a1, a2 in 2D (or a, b, c in 3D). Can always choose a1 parallel to x-axis.
- **Unit Cell:** Parallelogram defined by basis vectors. Area = sqrt(det(G)).
- **Metric Matrix G:**
  - 2D: G = [[a.a, a.b], [b.a, b.b]]
  - 3D: G = [[a.a, a.b, a.c], [b.a, b.b, b.c], [c.a, c.b, c.c]]
  - Volume = sqrt(det(G))
- **Transformation Matrix P:** (a2, b2, c2) = (a1, b1, c1) * P
  - Transformed coordinates: X2 = P^-1 * X1
  - Transformed G matrix: G2 = P^T * G1 * P
  - Area/volume ratio: A2/A1 = |det(P)|

### Chapter 3: Point Groups (2D)
- **2D Crystallographic Point Groups (10):** 1, 2, m, 2mm, 4, 4mm, 3, 3m, 6, 6mm
- **Crystallographic Restriction Theorem:** Only n-fold rotations with n = 1, 2, 3, 4, 6 are compatible with translational periodicity.
- **2D Crystal Systems (5):**
  - Oblique: no symmetry constraints on a, b, gamma (gamma != 90 deg)
  - Rectangular: a != b, gamma = 90 deg
  - Square: a = b, gamma = 90 deg
  - Hexagonal: a = b, gamma = 120 deg

### Chapter 4: Space Groups (2D -- Wallpaper Groups)
- **2D Space Groups (17):** p1, p2, pm, pg, cm, p2mm, p2mg, p2gg, c2mm, p3, p3m1, p31m, p4, p4mm, p4gm, p6, p6mm
- **5 Bravais Lattices in 2D:**
  - Oblique (p1, p2)
  - Rectangular primitive (pm, pg, p2mm, p2mg, p2gg)
  - Rectangular centered (cm, c2mm)
  - Square (p4, p4mm, p4gm)
  - Hexagonal (p3, p3m1, p31m, p6, p6mm)
- **Asymmetric Unit:** Smallest region that generates the entire pattern under symmetry operations
- **Space Group Hasse Diagram:** Shows subgroup relationships between wallpaper groups

### Chapter 5: Reciprocal Lattice
- **Reciprocal Basis Vectors (2D):**
  - b1 = 2pi/s * [(a2.a2)*a1 - (a1.a2)*a2]
  - b2 = 2pi/s * [(a1.a1)*a2 - (a1.a2)*a1]
  - where s = (a1.a1)(a2.a2) - (a1.a2)^2
  - Property: al . bm = 2pi*delta_lm
- **5 Bravais Lattices -- Basis and Reciprocal Vectors:**

| Lattice Type | a1 | a2 | b1 | b2 |
|---|---|---|---|---|
| Oblique | a(1,0) | (ax, ay) | (2pi/a)(1, -ax/ay) | (2pi/ay)(0,1) |
| Rectangular | a(1,0) | b(0,1) | (2pi/a)(1,0) | (2pi/b)(0,1) |
| Centered Rect. | (a/2,-b/2) | (a/2,b/2) | (2pi)(1/a,1/b) | (2pi)(1/a,-1/b) |
| Square | a(1,0) | a(0,1) | (2pi/a)(1,0) | (2pi/a)(0,1) |
| Hexagonal | a(1,0) | a(1/2,sqrt(3)/2) | (2pi/a)(1,-1/sqrt(3)) | (2pi/a)(0,2/sqrt(3)) |

- **d-spacing:** d_hk = 1/|hk1*b1 + hk2*b2|
- **Miller Indices:** (hk) in 2D -- reciprocals of fractional intercepts

### Chapter 7: Electron Density Maps / Fourier Series
- **Fourier Series in 2D:** rho(x,y) = Sum_hk F(hk) * exp(-2pi*i*(hx + ky))
- **Structure Factor:** F(hk) = Sum_j f_j * exp(-2pi*i*(hx_j + ky_j))
- **Reality Criterion:** For real rho(x,y), F(hk) = F*(-hk)
- **Friedel's Law:** |F(hk)| = |F(-h,-k)| for non-absorbing crystals

---

## Key Theory from Verberck (2012) -- Wallpaper Groups & Symmetry-Adapted Fourier Series

### Fundamental Fourier Framework
- **General Fourier Series:** f(r) = Sum_k c_k * exp(i*k*r)
  - k = k1*b1 + k2*b2 where b1, b2 are reciprocal basis vectors
  - c_k = (1/Omega) * integral_Omega f(r) * exp(-i*k*r) dr, where Omega is unit cell area

### Symmetry Operations on Fourier Coefficients

#### Rotation Axes
- **2-fold rotation (p2, rectangular/oblique):** c(k1,k2) = c(-k1,-k2)
- **3-fold rotation (p3, hexagonal):** c(k1,k2) = c(-k1+k2,-k1) = c(-k2,k1-k2) -- 3-cycle
- **4-fold rotation (p4, square):** c(k1,k2) = c(-k2,k1) = c(-k1,-k2) = c(k2,-k1) -- 4-cycle
- **6-fold rotation (p6, hexagonal):** c(k1,k2) = c(k2,-k1+k2) = c(-k1+k2,-k1) = c(-k1,-k2) = c(-k2,k1-k2) = c(k1-k2,k1) -- 6-cycle

#### Reflection Axes
- **Reflection about x-axis (pm, rectangular):** c(k1,k2) = c(k1,-k2)
- **Reflection about x-axis (hexagonal):** c(k1,k2) = c(k1,k1-k2)
- **Reflection about y-axis (square, p4mm):** c(k1,k2) = c(-k1,k2)
- **Reflection about diagonal (cm, centered rect.):** c(k1,k2) = c(k2,k1)

#### Glide Reflections
- **Glide along x-axis (pg, rectangular):** c(k1,k2) = (-1)^k1 * c(k1,-k2)
  - Implies: c(k1,0) = 0 for k1 odd
  - Phase factor: h(k1,k2) = (-1)^k1
- **Glide along y-axis (p2mg):** c(k1,k2) = (-1)^k2 * c(k1,-k2)
  - Phase factor: h(k1,k2) = (-1)^k2
- **Glide along diagonal (p2gg):** c(k1,k2) = (-1)^(k1+k2) * c(k1,-k2)
  - Phase factor: h(k1,k2) = (-1)^(k1+k2)

#### Centering (cm, c2mm -- centered rectangular)
- **Centering transformation:** r' = r + (a/2, b/2)
- **Index relation:** q1 = k1+k2, q2 = k1-k2 -> k1 = (q1+q2)/2, k2 = (q1-q2)/2
- **Centering condition:** (-1)^(q1+q2) * c(q1,q2) = c(q1,q2)
  - Implies: c(q1,q2) = 0 for q1+q2 odd
- **cm:** c(q1,q2) = c(q1,-q2) = (-1)^(q1+q2) * c(q1,q2)
- **c2mm:** c(q1,q2) = c(-q1,-q2) = c(q1,-q2) = c(-q1,q2) = (-1)^(q1+q2) * c(q1,q2)

### Minimal Domains for Independent Coefficients
- **D1:** Single point (k1,k2) -- cycle of 1 (identity)
- **D2:** Pairs {(k,-k), (-k,k)} -- cycle of 2
- **D3:** Triplets for 3-fold symmetry (hexagonal)
- **D4:** Quadruplets for 4-fold symmetry (square)
- **D6:** Sextuplets for 6-fold symmetry (hexagonal)
- **D12:** 12-cycles for p6mm (hexagonal with full symmetry)

### Key Explicit Fourier Expansions

#### p6 (hexagonal, 6-fold)
```
f(r) = c(0,0) + Sum_{(k1,k2) in D6} c(k1,k2) * [exp(i*[k1*b1+k2*b2]*r) + exp(i*[k2*b1+(-k1+k2)*b2]*r) + ... (6 terms)]
```
Domain D6: 0 < k1, 0 <= k2 < k1. Reduces coefficients by 6x.

#### p4 (square, 4-fold)
```
f(r) = c(0,0) + Sum_{(k1,k2) in D4} c(k1,k2) * [exp(i*[k1*b1+k2*b2]*r) + exp(i*[-k2*b1+k1*b2]*r) + exp(i*[-k1*b1-k2*b2]*r) + exp(i*[k2*b1-k1*b2]*r)]
```
Domain D4: k1 > 0, 0 <= k2 < k1. Reduces coefficients by 4x.

#### p2 (oblique/rectangular, 2-fold)
```
f(r) = c(0,0) + Sum_{(k1,k2) in D2} c(k1,k2) * [exp(i*[k1*b1+k2*b2]*r) + exp(i*[-k1*b1-k2*b2]*r)]
```
Domain D2: k1 > 0, or k1=0 and k2 > 0. Reduces by 2x.

#### pm (rectangular, mirror)
```
f(r) = c(0,0) + Sum_k1 [c(k1,0)*exp(i*k1*b1*r)] + 2*Sum_{k1,k2>0} Re[c(k1,k2)*exp(i*[k1*b1+k2*b2]*r)]
```
Domain: k2 >= 0. Reduces by 2x.

#### pg (rectangular, glide)
```
f(x,y) = c(0,0) + Sum_{k1 even} c(k1,0)*[cos(2pi*k1*x/a) + i*sin(2pi*k1*x/a)]
  + 2*Sum_{k2>0} { Sum_{k1 even} c(k1,k2)[cos*cos + i*sin*cos]
                 + Sum_{k1 odd}  c(k1,k2)[-sin*sin + i*cos*sin] }
```
Vanishing: c(k1,0) = 0 for k1 odd.

### Reality Criteria Summary (from Verberck Table 23)
- **Centrosymmetric groups** (have 2-fold at origin): all c_k must be real
- **Non-centrosymmetric groups:** c_k = c*(-k) relates coefficients within Dmin

Centrosymmetric wallpaper groups: p1, p2, pm, p2mm, p2mg, p2gg, c2mm, p4, p4mm, p4gm, p6, p6mm  
Non-centrosymmetric: pg, cm, p3, p3m1, p31m

### Fourier Coefficient Maps (Diffraction Patterns)
- Each wallpaper group has a characteristic "map" of equivalent coefficients in (k1,k2) space
- Maps show 100 for equivalent, 200 for doubly-equivalent, 1 for independent
- These maps are complementary to International Tables for Crystallography
- Maps can be used for identifying wallpaper groups from experimental data

---

## PSCF++ Manual Key Topics (fetch from URL when needed)

The PSCF++ manual covers:
1. **SCFT Theory:** Self-Consistent Field Theory for polymer systems
2. **Space Group Symmetry:** Implementation of space group operations in 3D
3. **Symmetry-Adapted Fourier Bases:** Construction of symmetry-adapted basis functions
4. **Periodic Functions:** Handling of periodic boundary conditions and Fourier representations
5. **Numerical Methods:** Real-space grid methods, FFT-based operations

**Key URLs:**
- Main: https://dmorse.github.io/pscfpp-man/
- SCFT Theory: /scft_theory.html
- Space Group Symmetry: /prdc_symmetry.html
- Symmetry-Adapted Fourier Bases: /prdc_basis.html
- Periodic Functions: /prdc_fourier.html
- Citation: /about_citation.html

---

## Mapping to SAFB C Code

### Theory to Code Mapping
| Theory Concept | SAFB C Module/Function | Description |
|---|---|---|
| Reciprocal lattice vectors | `reciprocal_lattice.c` | Computes b1, b2 from a1, a2 |
| Metric matrix G | `metric_matrix.c` | G_ij = ai.aj |
| 2D Bravais lattices | `lattice_type.c` | 5 lattice types: oblique, rectangular, centered, square, hexagonal |
| Wallpaper groups | `wallpaper_group.c` | 17 groups with symmetry operations |
| Symmetry-adapted bases | `sym_basis.c` | Constructs symmetry-adapted Fourier basis |
| Fourier coefficients | `fourier_coeffs.c` | Computes and manages Fourier coefficients |
| Space group operations | `space_group_ops.c` | Rotation, reflection, glide operations |
| Real-space grid | `real_space_grid.c` | Grid-based representation |
| FFT operations | `fft_operations.c` | Forward/inverse FFT |

### Key Mathematical Relationships
1. **Direct to Reciprocal:** b1 = 2pi*(a2 x n)/Omega, b2 = 2pi*(n x a1)/Omega (2D: cross products with normal vector n)
2. **Fourier Transform:** f(r) = Sum_k c_k * exp(i*k*r), c_k = (1/Omega)*integral f(r)*exp(-i*k*r) dr
3. **Symmetry constraint:** For operation R, c_k = chi(R) * c_{R^-1*k} where chi is the character
4. **Minimal domain:** Only independent coefficients stored; others generated by symmetry

---

## How to Use These References

1. **Primary reference:** PSCF++ manual (fetch from URL when writing specific sections)
2. **Crystallography textbook:** Use for foundational theory (lattices, metric matrices, reciprocal lattices, Miller indices, point groups, space groups)
3. **Verberck (2012):** Use for symmetry-adapted Fourier series -- explicit formulas, coefficient maps, minimal domains, reality criteria
4. **Large PDFs:** Only load when a specific section is needed; use the summary above for most tasks
