/**
 * domain.h — Data structures for SAFB (Symmetry-Adapted Fourier Basis)
 *
 * Translated from: Sg_init/domain.py
 *
 * Core data types: LatticeInfo, Star, SAFBBasis, InitializationResult,
 * ScatteringProfile.
 */

#ifndef DOMAIN_H
#define DOMAIN_H

#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ========================================================================
 * LatticeInfo — unit cell parameters
 * ======================================================================== */

typedef struct {
    double a, b, c;
    double alpha, beta, gamma;  /* degrees */
    int    dim;                  /* 2 or 3 */
} LatticeInfo;

LatticeInfo lattice_info_new(double a, double b, double c,
                             double alpha, double beta, double gamma, int dim);
LatticeInfo lattice_info_new_2d(double a, double b, double gamma);

/* ========================================================================
 * HKL — 3D reciprocal lattice vector
 * ======================================================================== */

typedef int HKL[3];

/* ========================================================================
 * Star — symmetry-equivalent set of reciprocal lattice vectors
 * ======================================================================== */

#define MAX_STAR_VECTORS 100

typedef struct {
    char family_key[32];
    double q2;
    int star_close;
    const char* reason;           /* literal string, not allocated */
    HKL star_vectors[MAX_STAR_VECTORS];
    int star_vectors_count;
    int multiplicity;
    /* rels: star → list of (source, phase, op_index)
     * Stored as flat array; use star_rels_lookup() to query.
     * Exact structure TBD during symmetry_ops translation.
     */
    void* rels;                   /* opaque — defined in symmetry_ops.h */
} Star;

/* ========================================================================
 * SAFBBasis — collection of valid Fourier modes for a space group
 * ======================================================================== */

#define MAX_MODES 500

typedef struct {
    char space_group[32];
    int  centrosymmetric_group;
    int  has_inversion_at_origin;
    const char* additional_info;   /* literal string */
    LatticeInfo lattice;
    Star modes[MAX_MODES];
    int modes_count;
} SAFBBasis;

/* ========================================================================
 * InitializationResult — output of amplitude assignment
 * ======================================================================== */

typedef struct {
    char space_group[32];
    LatticeInfo lattice;
    SAFBBasis* safb;
    void* coeff;                    /* opaque: hash table family_key → {hkl → complex} */
    void* amplitudes;               /* opaque: hash table family_key → float */
    int discretization_Na;
    int discretization_Nb;
    int discretization_Nc;
    /* trans_coeff and trans_lattice_params for lattice transforms */
    void* trans_coeff;
    double trans_lattice_a, trans_lattice_b, trans_lattice_c;
    double trans_lattice_alpha, trans_lattice_beta, trans_lattice_gamma;
} InitializationResult;

/* ========================================================================
 * ScatteringProfile — parsed experimental scattering data
 * ======================================================================== */

typedef struct {
    int32_t* hkl;       /* N×3 array of Miller indices */
    double* q;          /* N scattering vector magnitudes */
    double* intensity;  /* N scattering intensities */
    int num_peaks;
} ScatteringProfile;

ScatteringProfile scattering_profile_new(int N);
void scattering_profile_free(ScatteringProfile* sp);

#ifdef __cplusplus
}
#endif

#endif /* DOMAIN_H */
