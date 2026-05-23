/**
 * initializers.h — Amplitude assignment initializers for SAFB
 *
 * Translated from: Sg_init/initializers.py
 *
 * Contains:
 *   - ScatteringProfile parsing (read_scattering_data)
 *   - build_result_from_coeffs — solve phase constraints and build InitializationResult
 *   - random_initialization — random amplitude assignment
 *   - manual_initialization — manual amplitude assignment
 *   - file_initialization — amplitude from scattering data file
 */

#ifndef INITIALIZERS_H
#define INITIALIZERS_H

#include <stdint.h>
#include <stddef.h>
#include "domain.h"
#include "symmetry_ops.h"
#include "basis.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ========================================================================
 * Scattering data file parser
 * ======================================================================== */

ScatteringProfile read_scattering_data(const char *filepath);

/* ========================================================================
 * Coefficient storage — family_key → {hkl → complex}
 * ======================================================================== */

#define MAX_AMPLITUDE_KEYS 256
#define MAX_VEC_PER_STAR   64

typedef struct {
    int hkl[3];
    double real;
    double imag;
} ComplexCoeff;

typedef struct {
    char family_key[32];
    int count;
    ComplexCoeff coeffs[MAX_VEC_PER_STAR];
} FamilyCoeffs;

/* ========================================================================
 * Build InitializationResult from amplitudes — solve phase constraints
 * ======================================================================== */

/* Solve phase constraints for all modes and return coeffs + warnings */
int build_result_from_coeffs(
    const SAFBBasis *basis,
    const double amplitudes[],       /* amplitude_keys[i] → amplitudes[i] */
    const char *amplitude_keys[],    /* family_key strings */
    int n_amplitudes,
    FamilyCoeffs out_coeffs[],       /* output: family_key → {hkl → complex} */
    char warnings[][256],            /* warning strings */
    int max_warnings
);

/* ========================================================================
 * Random initialization
 * ======================================================================== */

typedef enum {
    DIST_UNIFORM = 0,
    DIST_NORMAL,
    DIST_LOGNORMAL,
    DIST_EXPONENTIAL
} DistributionType;

typedef struct {
    double loc;
    double scale;
} DistParams;

/* RNG state — minimal xorshift64 for portability */
typedef struct {
    uint64_t state[2];
} SAFBRNG;

void srand_init(SAFBRNG *rng, uint64_t seed);
double srand_next(SAFBRNG *rng);
double srand_normal(SAFBRNG *rng, double mean, double stddev);
double srand_uniform(SAFBRNG *rng);
double srand_exponential(SAFBRNG *rng, double lambda);
double srand_lognormal(SAFBRNG *rng, double mean_log, double sigma_log);

int random_initialization(
    const SAFBBasis *basis,
    int n_modes,
    DistributionType dist,
    DistParams params,
    SAFBRNG *rng,
    FamilyCoeffs out_coeffs[],
    char warnings[][256],
    int max_warnings
);

/* ========================================================================
 * Manual initialization
 * ======================================================================== */

int manual_initialization(
    const SAFBBasis *basis,
    const char *const amplitude_keys[],
    const double amplitudes[],
    int n_amplitudes,
    FamilyCoeffs out_coeffs[],
    char warnings[][256],
    int max_warnings
);

/* ========================================================================
 * File initialization — read amplitudes from scattering file
 * ======================================================================== */

int file_initialization(
    const SAFBBasis *basis,
    const char *scattering_file,
    const LatticeInfo *read_lattice_info,
    const double P[3][3],  /* change-of-basis matrix (NULL if same lattice) */
    int n_modes,
    FamilyCoeffs out_coeffs[],
    char warnings[][256],
    int max_warnings
);

/* ========================================================================
 * Build full InitializationResult (convenience wrapper)
 * ======================================================================== */

typedef struct {
    char space_group[32];
    LatticeInfo lattice;
    SAFBBasis safb;
    FamilyCoeffs coeffs[MAX_AMPLITUDE_KEYS];
    int n_coeffs;
    double amplitudes[MAX_AMPLITUDE_KEYS];
    char amplitude_keys[MAX_AMPLITUDE_KEYS][32];
    int n_amplitudes;
    int discretization_Na;
    int discretization_Nb;
    int discretization_Nc;
    int trans_valid;
    double trans_lattice_a, trans_lattice_b, trans_lattice_c;
    double trans_lattice_alpha, trans_lattice_beta, trans_lattice_gamma;
    double trans_coeffs_real[MAX_AMPLITUDE_KEYS][MAX_VEC_PER_STAR];
    double trans_coeffs_imag[MAX_AMPLITUDE_KEYS][MAX_VEC_PER_STAR];
    int trans_counts[MAX_AMPLITUDE_KEYS];
    int trans_keys[MAX_AMPLITUDE_KEYS][3];
} FullInitializationResult;

int build_initialization_result(
    const SAFBBasis *basis,
    const FamilyCoeffs coeffs[],
    int n_coeffs,
    const double amplitudes[],
    const char *amplitude_keys[],
    int n_amplitudes,
    FullInitializationResult *result
);

#ifdef __cplusplus
}
#endif

#endif /* INITIALIZERS_H */
