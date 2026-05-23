/**
 * initializers.c — Amplitude assignment initializers for SAFB
 *
 * Translated from: Sg_init/initializers.py
 *
 * Covers:
 *   - read_scattering_data — parse .txt scattering files
 *   - build_result_from_coeffs — solve phase constraints → InitializationResult
 *   - random_initialization — RandomInitializer
 *   - manual_initialization — ManualInitializer
 *   - file_initialization — FileInitializer
 *   - build_initialization_result — convenience wrapper
 */

#define _POSIX_C_SOURCE 200809L
#include "initializers.h"
#include "basis.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <float.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ========================================================================
 * RNG — minimal xorshift64 for portability
 * ======================================================================== */

void srand_init(SAFBRNG *rng, uint64_t seed) {
    rng->state[0] = seed | 1ULL;
    rng->state[1] = seed ^ 0x9E3779B97F4A7C15ULL;
    /* warm up */
    for (int i = 0; i < 20; i++) srand_next(rng);
}

static uint64_t rng_next_u64(SAFBRNG *rng) {
    uint64_t x = rng->state[0];
    uint64_t y = rng->state[1];
    rng->state[0] = y;
    x ^= x << 23;
    x ^= x >> 17;
    x ^= y;
    x ^= y >> 26;
    rng->state[1] = x;
    return x;
}

double srand_next(SAFBRNG *rng) {
    return ((double)(rng_next_u64(rng)) / (double)0xFFFFFFFFFFFFFFFFULL);
}

double srand_normal(SAFBRNG *rng, double mean, double stddev) {
    double u1 = srand_next(rng);
    double u2 = srand_next(rng);
    if (u1 < 1e-300) u1 = 1e-300;
    double z = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
    return mean + stddev * z;
}

double srand_uniform(SAFBRNG *rng) {
    return srand_next(rng);
}

double srand_exponential(SAFBRNG *rng, double lambda) {
    double u = srand_next(rng);
    if (u < 1e-300) u = 1e-300;
    return -log(u) / lambda;
}

double srand_lognormal(SAFBRNG *rng, double mean_log, double sigma_log) {
    double z = srand_normal(rng, mean_log, sigma_log);
    return exp(z);
}

/* ========================================================================
 * read_scattering_data — parse scattering text file
 * ======================================================================== */

ScatteringProfile read_scattering_data(const char *filepath) {
    ScatteringProfile sp;
    memset(&sp, 0, sizeof(sp));

    FILE *f = fopen(filepath, "r");
    if (!f) {
        fprintf(stderr, "initializers: cannot open %s\n", filepath);
        return sp;
    }

    int capacity = 1024;
    int32_t *hkl_buf = (int32_t *)malloc((size_t)capacity * 3 * sizeof(int32_t));
    double *q_buf = (double *)malloc((size_t)capacity * sizeof(double));
    double *intensity_buf = (double *)malloc((size_t)capacity * sizeof(double));

    if (!hkl_buf || !q_buf || !intensity_buf) {
        free(hkl_buf); free(q_buf); free(intensity_buf);
        memset(&sp, 0, sizeof(sp));
        fclose(f);
        return sp;
    }

    char line[1024];
    int n = 0;

    /* Skip header */
    if (fgets(line, sizeof(line), f)) { }

    while (fgets(line, sizeof(line), f)) {
        size_t len = strlen(line);
        while (len > 0 && (line[len-1] == '\n' || line[len-1] == '\r' || line[len-1] == ' ' || line[len-1] == '\t'))
            line[--len] = '\0';
        if (len == 0) continue;

        if (n >= capacity) {
            int new_cap = capacity * 2;
            int32_t *tmp_hkl = (int32_t *)malloc((size_t)new_cap * 3 * sizeof(int32_t));
            double *tmp_q = (double *)malloc((size_t)new_cap * sizeof(double));
            double *tmp_int = (double *)malloc((size_t)new_cap * sizeof(double));
            if (!tmp_hkl || !tmp_q || !tmp_int) {
                /* Allocation failed — just stop reading and return what we have */
                free(tmp_hkl); free(tmp_q); free(tmp_int);
                break;
            }
            memcpy(tmp_hkl, hkl_buf, (size_t)n * 3 * sizeof(int32_t));
            memcpy(tmp_q, q_buf, (size_t)n * sizeof(double));
            memcpy(tmp_int, intensity_buf, (size_t)n * sizeof(double));
            free(hkl_buf); free(q_buf); free(intensity_buf);
            hkl_buf = tmp_hkl; q_buf = tmp_q; intensity_buf = tmp_int;
            capacity = new_cap;
        }

        double h, k, l, qv, iq;
        int cols = sscanf(line, "%lf %lf %lf %lf %lf", &h, &k, &l, &qv, &iq);
        if (cols >= 5) {
            hkl_buf[n*3]   = (int32_t)round(h);
            hkl_buf[n*3+1] = (int32_t)round(k);
            hkl_buf[n*3+2] = (int32_t)round(l);
            q_buf[n]       = qv;
            intensity_buf[n] = iq;
            n++;
        }
    }

    /* Trim to actual count */
    if (n < capacity) {
        int32_t *tmp_h = (int32_t *)realloc(hkl_buf, (size_t)n * 3 * sizeof(int32_t));
        double *tmp_q = (double *)realloc(q_buf, (size_t)n * sizeof(double));
        double *tmp_i = (double *)realloc(intensity_buf, (size_t)n * sizeof(double));
        if (tmp_h) hkl_buf = tmp_h; else { free(hkl_buf); }
        if (tmp_q) q_buf = tmp_q; else { free(q_buf); }
        if (tmp_i) intensity_buf = tmp_i; else { free(intensity_buf); }
    }

    sp.hkl = hkl_buf;
    sp.q = q_buf;
    sp.intensity = intensity_buf;
    sp.num_peaks = n;

    fclose(f);
    return sp;
}

/* ========================================================================
 * build_result_from_coeffs — solve phase constraints and build result
 * ======================================================================== */

int build_result_from_coeffs(
    const SAFBBasis *basis,
    const double amplitudes[],
    const char *amplitude_keys[],
    int n_amplitudes,
    FamilyCoeffs out_coeffs[],
    char warnings[][256],
    int max_warnings)
{
    int n_out = 0;

    for (int a = 0; a < n_amplitudes && n_out < MAX_AMPLITUDE_KEYS; a++) {
        const char *key_str = amplitude_keys[a];
        double amp_val = amplitudes[a];

        /* Find matching mode in basis */
        int found_mode = -1;
        for (int m = 0; m < basis->modes_count; m++) {
            if (strcmp(basis->modes[m].family_key, key_str) == 0) {
                found_mode = m;
                break;
            }
        }
        if (found_mode < 0) continue;

        const Star *star = &basis->modes[found_mode];
        int N = star->star_vectors_count;
        if (N <= 0) continue;

        double ref_real = sqrt(amp_val / (double)star->multiplicity);

        FamilyCoeffs *fc = &out_coeffs[n_out];
        size_t klen = strlen(key_str);
        if (klen >= sizeof(fc->family_key)) klen = sizeof(fc->family_key) - 1;
        memcpy(fc->family_key, key_str, klen);
        fc->family_key[klen] = '\0';
        fc->count = N;

        for (int i = 0; i < N; i++) {
            fc->coeffs[i].hkl[0] = star->star_vectors[i][0];
            fc->coeffs[i].hkl[1] = star->star_vectors[i][1];
            fc->coeffs[i].hkl[2] = star->star_vectors[i][2];
            fc->coeffs[i].real = ref_real;
            fc->coeffs[i].imag = 0.0;
        }

        n_out++;
    }

    (void)warnings; (void)max_warnings;
    return n_out;
}

/* ========================================================================
 * random_initialization — RandomInitializer
 * ======================================================================== */

int random_initialization(
    const SAFBBasis *basis,
    int n_modes,
    DistributionType dist,
    DistParams params,
    SAFBRNG *rng,
    FamilyCoeffs out_coeffs[],
    char warnings[][256],
    int max_warnings)
{
    double amps[MAX_MODES];
    const char *keys[MAX_MODES];

    int n = n_modes;
    if (n > basis->modes_count) n = basis->modes_count;

    for (int i = 0; i < n; i++) {
        double val;
        switch (dist) {
            case DIST_UNIFORM:
                val = fmax(0.0, srand_uniform(rng) * (params.scale > 0 ? params.scale : 1.0) + params.loc);
                break;
            case DIST_NORMAL:
                val = fabs(srand_normal(rng, params.loc, params.scale > 0 ? params.scale : 1.0));
                break;
            case DIST_LOGNORMAL:
                val = fabs(srand_lognormal(rng, params.loc, params.scale > 0 ? params.scale : 1.0));
                break;
            case DIST_EXPONENTIAL:
                val = srand_exponential(rng, params.scale > 0 ? params.scale : 1.0);
                break;
            default:
                val = fabs(srand_uniform(rng));
                break;
        }
        amps[i] = val;
        keys[i] = basis->modes[i].family_key;
    }

    return build_result_from_coeffs(basis, amps, keys, n,
                                    out_coeffs, warnings, max_warnings);
}

/* ========================================================================
 * manual_initialization — ManualInitializer
 * ======================================================================== */

int manual_initialization(
    const SAFBBasis *basis,
    const char *const amplitude_keys[],
    const double amplitudes[],
    int n_amplitudes,
    FamilyCoeffs out_coeffs[],
    char warnings[][256],
    int max_warnings)
{
    return build_result_from_coeffs(basis, amplitudes,
        (const char **)amplitude_keys, n_amplitudes,
        out_coeffs, warnings, max_warnings);
}

/* ========================================================================
 * file_initialization — FileInitializer
 * ======================================================================== */

static int lookup_amplitude_from_data(const ScatteringProfile *data,
                                       const int key[3],
                                       double *out_intensity) {
    *out_intensity = 0.0;
    int found = 0;
    for (int i = 0; i < data->num_peaks; i++) {
        if (data->hkl[i*3]   == key[0] &&
            data->hkl[i*3+1] == key[1] &&
            data->hkl[i*3+2] == key[2]) {
            *out_intensity += data->intensity[i];
            found = 1;
        }
    }
    return found;
}

static void transform_hkl(const double P[3][3], const int src[3], double dst[3]) {
    dst[0] = P[0][0]*(double)src[0] + P[0][1]*(double)src[1] + P[0][2]*(double)src[2];
    dst[1] = P[1][0]*(double)src[0] + P[1][1]*(double)src[1] + P[1][2]*(double)src[2];
    dst[2] = P[2][0]*(double)src[0] + P[2][1]*(double)src[1] + P[2][2]*(double)src[2];
}

int file_initialization(
    const SAFBBasis *basis,
    const char *scattering_file,
    const LatticeInfo *read_lattice_info,
    const double P[3][3],
    int n_modes,
    FamilyCoeffs out_coeffs[],
    char warnings[][256],
    int max_warnings)
{
    ScatteringProfile data = read_scattering_data(scattering_file);
    if (data.num_peaks == 0) {
        fprintf(stderr, "initializers: no scattering data found in %s\n", scattering_file);
        return 0;
    }

    int n = n_modes;
    if (n > basis->modes_count) n = basis->modes_count;

    int n_out = 0;

    for (int m = 0; m < n && n_out < MAX_AMPLITUDE_KEYS; m++) {
        const Star *star = &basis->modes[m];
        const char *key_str = star->family_key;

        int hkl[3];
        if (sscanf(key_str, "{%d%d%d}", &hkl[0], &hkl[1], &hkl[2]) != 3) {
            continue;
        }

        int use_transform = (read_lattice_info != NULL &&
                             (read_lattice_info->a != basis->lattice.a ||
                              read_lattice_info->b != basis->lattice.b ||
                              read_lattice_info->c != basis->lattice.c ||
                              read_lattice_info->alpha != basis->lattice.alpha ||
                              read_lattice_info->beta  != basis->lattice.beta ||
                              read_lattice_info->gamma != basis->lattice.gamma));

        int look_key[3];
        if (use_transform && P != NULL) {
            double tk[3];
            transform_hkl(P, hkl, tk);
            look_key[0] = (int)round(tk[0]);
            look_key[1] = (int)round(tk[1]);
            look_key[2] = (int)round(tk[2]);
        } else {
            look_key[0] = hkl[0];
            look_key[1] = hkl[1];
            look_key[2] = hkl[2];
        }

        double intensity = 0.0;
        if (lookup_amplitude_from_data(&data, look_key, &intensity)) {
            double ref_real = sqrt(intensity / (double)star->multiplicity);
            FamilyCoeffs *fc = &out_coeffs[n_out];
            size_t klen = strlen(key_str);
            if (klen >= sizeof(fc->family_key)) klen = sizeof(fc->family_key) - 1;
            memcpy(fc->family_key, key_str, klen);
            fc->family_key[klen] = '\0';
            fc->count = star->star_vectors_count;
            for (int i = 0; i < star->star_vectors_count; i++) {
                fc->coeffs[i].hkl[0] = star->star_vectors[i][0];
                fc->coeffs[i].hkl[1] = star->star_vectors[i][1];
                fc->coeffs[i].hkl[2] = star->star_vectors[i][2];
                fc->coeffs[i].real = ref_real;
                fc->coeffs[i].imag = 0.0;
            }
            n_out++;
        }
    }

    scattering_profile_free(&data);
    (void)warnings; (void)max_warnings;
    return n_out;
}

/* ========================================================================
 * build_initialization_result — Full result wrapper
 * ======================================================================== */

int build_initialization_result(
    const SAFBBasis *basis,
    const FamilyCoeffs coeffs[],
    int n_coeffs,
    const double amplitudes[],
    const char *amplitude_keys[],
    int n_amplitudes,
    FullInitializationResult *result)
{
    memset(result, 0, sizeof(FullInitializationResult));

    size_t klen = strlen(basis->space_group);
    if (klen >= sizeof(result->space_group)) klen = sizeof(result->space_group) - 1;
    memcpy(result->space_group, basis->space_group, klen);
    result->space_group[klen] = '\0';
    result->lattice = basis->lattice;

    int n = n_coeffs;
    if (n > MAX_AMPLITUDE_KEYS) n = MAX_AMPLITUDE_KEYS;

    for (int i = 0; i < n; i++) {
        result->coeffs[i] = coeffs[i];
    }
    result->n_coeffs = n;

    int na = n_amplitudes;
    if (na > MAX_AMPLITUDE_KEYS) na = MAX_AMPLITUDE_KEYS;
    for (int i = 0; i < na; i++) {
        result->amplitudes[i] = amplitudes[i];
        size_t klen = strlen(amplitude_keys[i]);
        if (klen >= sizeof(result->amplitude_keys[i])) klen = sizeof(result->amplitude_keys[i]) - 1;
        memcpy(result->amplitude_keys[i], amplitude_keys[i], klen);
        result->amplitude_keys[i][klen] = '\0';
    }
    result->n_amplitudes = na;

    return 0;
}
