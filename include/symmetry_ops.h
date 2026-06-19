/**
 * symmetry_ops.h — Symmetry operation utilities for SAFB
 *
 * Translated from: Sg_init/symmetry_ops.py
 *
 * Contains:
 *   - Fraction arithmetic (rational numbers for symmetry translations)
 *   - Space group operation parsing (read_spacegroup_ops_txt)
 *   - Lattice basis helpers (direct_basis_from_lattice_info, lattice_params_from_basis)
 *   - Star generation and analysis
 *   - Phase constraint solver
 */

#ifndef SYMMETRY_OPS_H
#define SYMMETRY_OPS_H

#include <stdint.h>
#include <stddef.h>
#include "domain.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ========================================================================
 * Bounds constants
 * ======================================================================== */

#define MAX_ADJ 200          /* max adjacency entries per star vector */
#define MAX_WARNINGS 64      /* max warnings in solve_star_coeffs */

/* ========================================================================
 * Fraction — rational number for symmetry translations
 * ======================================================================== */

typedef struct {
    int num;  /* numerator */
    int den;  /* denominator (always > 0) */
} Fraction;

Fraction frac_new(int num, int den);
double frac_to_float(Fraction q);
Fraction frac_mod1(Fraction q);        /* reduce mod 1 into [0,1) */
int is_zero_mod1_vec(const Fraction t[3]);    /* all elements mod 1 == 0? */
int equal_int_mat(const int A[3][3], const int B[3][3]);

/* GCD helper */
int frac_gcd(int a, int b);

/* ========================================================================
 * Space group operation — one (R, t) pair
 * ======================================================================== */

#define MAX_OPS 1024

typedef struct {
    int R[3][3];              /* 3×3 integer rotation matrix */
    Fraction t[3];            /* 3-element fractional translation */
} SymmOp;

/* ========================================================================
 * Parsed space group / wallpaper group
 * ======================================================================== */

typedef struct {
    int dim;                  /* 2 or 3 */
    int count;                /* number of symmetry operations */
    SymmOp ops[MAX_OPS];
} SymmGroup;

SymmGroup read_spacegroup_ops_txt(const char *path);

/* ========================================================================
 * Unique rotations — extract unique rotation matrices from a SymmGroup
 * ======================================================================== */

int unique_rotations(const SymmGroup *sg, int result[/* MAX_OPS */][3][3]);

/* ========================================================================
 * Metric — inverse metric tensor from lattice parameters
 * ======================================================================== */

void metric_inverse(double a, double b, double c,
                    double alpha, double beta, double gamma,
                    double Ginv[3][3]);

double q2_metric(int h, int k, int l, const double Ginv[3][3]);

/* ========================================================================
 * Direct lattice basis from LatticeInfo
 * Returns 3×3 matrix A whose columns are the direct lattice vectors
 * ======================================================================== */

void direct_basis_from_lattice_info(const LatticeInfo *info, double A[3][3]);

/* ========================================================================
 * Lattice parameters from direct basis matrix
 * ======================================================================== */

void lattice_params_from_basis(const double A[3][3], double *a, double *b, double *c,
                                double *alpha, double *beta, double *gamma);

/* ========================================================================
 * Star generation — apply all rotations to an HKL vector
 * ======================================================================== */

int star_from_hkl(const int hkl[3],
                  const int rotations[/* n_rotations */][3][3],
                  int n_rotations,
                  const double Ginv[3][3],
                  int result[/* MAX_STAR_VECTORS */][3]);

/* ========================================================================
 * Family key — canonical representative of a star (lexicographic)
 * ======================================================================== */

void get_family_key_lexicographical(const int star_vectors[/* N */][3],
                                     int N,
                                     int key[3]);

/* ========================================================================
 * Star closure — does -G exist in the star?
 * ======================================================================== */

int star_is_closed(const int star[/* N */][3], int N);

/* ========================================================================
 * Point group check — is -I in the rotation set?
 * ======================================================================== */

int point_group_has_neg_identity(const int rotations[/* N */][3][3], int N);

/* ========================================================================
 * Inversion ops — find operations with R = -I
 * ======================================================================== */

int find_inversion_ops(const SymmGroup *sg,
                        Fraction t_result[/* MAX_OPS */][3],
                        int at_origin[/* MAX_OPS */]);

/* ========================================================================
 * Phase factor — exp(2πi * h·t) for a symmetry operation
 * ======================================================================== */

double phase_factor_real(const int hkl[3], const Fraction t[3]);
double phase_factor_imag(const int hkl[3], const Fraction t[3]);

/* ========================================================================
 * Star relationships — compute phase constraints within a star
 *
 * Returns relationship graph: for each Gk, which Gj maps to it with what phase.
 * Written into a flat structure:
 *   For each star vector j, rels[j] = list of (source_index, phase_real, phase_imag, op_index)
 *
 * Returns (rels, unary_contra_count, binary_contra_count)
 * ======================================================================== */

#define MAX_STAR_RELS_PER_VECTOR 100

typedef struct {
    int source;
    double phase_real;
    double phase_imag;
    int op_index;
} StarRel;

typedef struct {
    StarRel rels[MAX_STAR_VECTORS][MAX_STAR_RELS_PER_VECTOR];
    int rels_count[MAX_STAR_VECTORS];
    int unary_contra_count;
    int binary_contra_count;
} StarRelationships;

StarRelationships relationships_in_star(const SymmGroup *sg,
                                         const int star[/* N */][3],
                                         int N);

/* ========================================================================
 * Solve star coefficients — BFS over phase constraint graph
 *
 * Given a star and its relationship graph, compute complex coefficients
 * for each wave vector (up to a global phase reference).
 *
 * Returns: assigns coeffs[0..N-1] to complex values.
 * Warnings stored as strings.
 * ======================================================================== */

int solve_star_coeffs(const int star[/* N */][3], int N,
                       const StarRelationships *rels,
                       double ref_real,
                       double coeffs_real[/* N */],
                       double coeffs_imag[/* N */],
                       char warnings[/* MAX_WARNINGS */][256],
                       int max_warnings);

/* ========================================================================
 * Family planes — generate valid Miller index modes for a given space group
 *
 * Iteratively expands search range until N valid modes are found.
 * ======================================================================== */

/* Heap-allocated family planes info (avoids ~114MB stack allocation) */
typedef struct {
    int family_keys[MAX_MODES][3];   /* canonical HKL for each mode */
    int stars[MAX_MODES][MAX_STAR_VECTORS][3];  /* all vectors in each star */
    int star_counts[MAX_MODES];        /* how many vectors in each star */
    int count;                          /* number of valid modes found */
    StarRelationships* all_rels;        /* heap-allocated array */
} FamilyPlanesInfo;

FamilyPlanesInfo* family_planes_info(int N,
                                     const double Ginv[3][3],
                                     const SymmGroup *sg,
                                     int dim);
void family_planes_info_free(FamilyPlanesInfo* info);

/* ========================================================================
 * Read final simulation values
 * ======================================================================== */

typedef struct {
    int final_Nx, final_Ny, final_Nz;
    double final_Lx, final_Ly, final_Lz;
    double Nkai;
    double f;
    double free_energy;
} SimulationValues;

SimulationValues read_final_values(const char *file_path);

/* ========================================================================
 * Float to Miller integer conversion
 * ======================================================================== */

void float_to_miller_int(const double hB[3], int max_den, double tol,
                          int primitive,
                          int result[3]);

/* ========================================================================
 * Transform Miller indices between lattices
 * ======================================================================== */

void transform_miller_between_lattices(const LatticeInfo *lattice_A_info,
                                        const double P[3][3],
                                        const int existing_keys[][3],
                                        int num_keys,
                                        const double existing_coeffs_real[][MAX_STAR_VECTORS],
                                        const double existing_coeffs_imag[][MAX_STAR_VECTORS],
                                        const int existing_star_counts[],
                                        int (*trans_keys)[3],
                                        double trans_coeffs_real[MAX_STAR_VECTORS][MAX_STAR_VECTORS],
                                        double trans_coeffs_imag[MAX_STAR_VECTORS][MAX_STAR_VECTORS],
                                        int *trans_counts,
                                        double *new_lattice_a, double *new_lattice_b, double *new_lattice_c,
                                        double *new_lattice_alpha, double *new_lattice_beta, double *new_lattice_gamma);

#ifdef __cplusplus
}
#endif

#endif /* SYMMETRY_OPS_H */
