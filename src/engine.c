/**
 * engine.c — High-level SAFB initialization engine
 *
 * Translated from: Sg_init/engine.py
 *
 * Procedural implementation of SpaceGroupInitializationEngine methods.
 */

#include "engine.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

/* ========================================================================
 * Engine context management
 * ======================================================================== */

int engine_create(EngineContext *ctx,
                  const char *sg_symbol,
                  const char *symops_path,
                  const LatticeInfo *lattice,
                  int n_keep)
{
    if (!ctx || !sg_symbol || !symops_path || !lattice) {
        return -1;
    }

    memset(ctx, 0, sizeof(EngineContext));

    /* Copy space group symbol */
    strncpy(ctx->space_group, sg_symbol, sizeof(ctx->space_group) - 1);
    ctx->space_group[sizeof(ctx->space_group) - 1] = '\0';

    /* Copy lattice */
    ctx->lattice.a = lattice->a;
    ctx->lattice.b = lattice->b;
    ctx->lattice.c = lattice->c;
    ctx->lattice.alpha = lattice->alpha;
    ctx->lattice.beta = lattice->beta;
    ctx->lattice.gamma = lattice->gamma;
    ctx->lattice.dim = lattice->dim;

    ctx->n_keep = n_keep;

    /* Read symmetry operations */
    ctx->ops = read_spacegroup_ops_txt(symops_path);
    if (ctx->ops.count <= 0) {
        fprintf(stderr, "engine_create: failed to read ops from %s\n", symops_path);
        return -1;
    }

    /* Build basis */
    int ret = engine_build_basis(ctx, n_keep);
    if (ret != 0) {
        fprintf(stderr, "engine_create: failed to build basis\n");
        return -1;
    }

    return 0;
}

void engine_free(EngineContext *ctx)
{
    if (!ctx) return;
    /* basis contains embedded arrays, no free needed */
    memset(ctx, 0, sizeof(EngineContext));
}

int engine_build_basis(EngineContext *ctx, int N)
{
    if (!ctx) return -1;

    SAFBBasis result;
    int ret = basis_build(&ctx->ops, ctx->space_group, &ctx->lattice, N, &result);
    if (ret != 0) return -1;

    /* Copy result into context */
    memcpy(ctx->basis.space_group, result.space_group, sizeof(ctx->basis.space_group));
    ctx->basis.centrosymmetric_group = result.centrosymmetric_group;
    ctx->basis.has_inversion_at_origin = result.has_inversion_at_origin;
    ctx->basis.additional_info = result.additional_info; /* literal string, no copy */
    ctx->basis.lattice = result.lattice;
    ctx->basis.modes_count = result.modes_count;

    for (int i = 0; i < result.modes_count && i < MAX_MODES; i++) {
        ctx->basis.modes[i].q2 = result.modes[i].q2;
        ctx->basis.modes[i].star_close = result.modes[i].star_close;
        ctx->basis.modes[i].reason = result.modes[i].reason;
        ctx->basis.modes[i].rels = result.modes[i].rels;
        memcpy(ctx->basis.modes[i].family_key, result.modes[i].family_key,
               sizeof(ctx->basis.modes[i].family_key));
        ctx->basis.modes[i].multiplicity = result.modes[i].multiplicity;
        ctx->basis.modes[i].star_vectors_count = result.modes[i].star_vectors_count;
        memcpy(ctx->basis.modes[i].star_vectors, result.modes[i].star_vectors,
               sizeof(ctx->basis.modes[i].star_vectors));
    }

    return 0;
}

/* ========================================================================
 * Random initialization
 * ======================================================================== */

int engine_random_init(
    const EngineContext *ctx,
    FullInitializationResult *result)
{
    if (!ctx || !result) return -1;

    SAFBRNG rng;
    srand_init(&rng, 0); /* default seed 0, same as Python's default */

    FamilyCoeffs out_coeffs[MAX_MODES];
    char warnings[MAX_WARNINGS][256];

    int ret = random_initialization(
        &ctx->basis,
        ctx->n_keep,
        DIST_UNIFORM,
        (DistParams){0.0, 1.0},
        &rng,
        out_coeffs,
        warnings,
        MAX_WARNINGS
    );
    if (ret == 0) return -1;  // 0 = no coeffs matched = failure

    /* Build result from coeffs — construct proper pointer array for amplitude_keys */
    char raw_keys[MAX_MODES][32];
    const char *amplitude_keys[MAX_MODES];
    double amplitudes[MAX_MODES];
    memset(raw_keys, 0, sizeof(raw_keys));
    memset(amplitudes, 0, sizeof(amplitudes));
    int valid = ctx->basis.modes_count;
    if (valid > MAX_MODES) valid = MAX_MODES;
    for (int i = 0; i < valid; i++) {
        if (out_coeffs[i].count > 0) {
            strncpy(raw_keys[i], out_coeffs[i].family_key, 31);
            raw_keys[i][31] = '\0';
            amplitude_keys[i] = raw_keys[i];
            amplitudes[i] = out_coeffs[i].coeffs[0].real;
        }
    }

    int ret2 = build_initialization_result(
        &ctx->basis,
        out_coeffs,
        valid,
        amplitudes,
        amplitude_keys,
        valid,
        result
    );
    return ret2;
}

/* ========================================================================
 * Manual initialization
 * ======================================================================== */

int engine_manual_init(
    const EngineContext *ctx,
    const char *const amplitude_keys[],
    const double amplitudes[],
    int n_amplitudes,
    FullInitializationResult *result)
{
    if (!ctx || !amplitude_keys || !amplitudes || !result || n_amplitudes <= 0) return -1;

    FamilyCoeffs out_coeffs[MAX_MODES];
    char warnings[MAX_WARNINGS][256];

    int ret = manual_initialization(
        &ctx->basis,
        amplitude_keys,
        amplitudes,
        n_amplitudes,
        out_coeffs,
        warnings,
        MAX_WARNINGS
    );
    if (ret == 0) return -1;  // no matching keys found

    return build_initialization_result(
        &ctx->basis,
        out_coeffs,
        n_amplitudes,
        amplitudes,
        amplitude_keys,
        n_amplitudes,
        result
    );
}

/* ========================================================================
 * File initialization — replicate engine.init_from_file()
 * ======================================================================== */

int engine_file_init(
    const EngineContext *ctx,
    const char *scattering_file,
    const LatticeInfo *read_lattice_info,
    const double P[3][3],
    FullInitializationResult *result)
{
    if (!ctx || !scattering_file || !read_lattice_info || !result) return -1;

    FamilyCoeffs out_coeffs[MAX_MODES];
    char warnings[MAX_WARNINGS][256];

    int ret = file_initialization(
        &ctx->basis,
        scattering_file,
        read_lattice_info,
        P,
        ctx->n_keep,
        out_coeffs,
        warnings,
        MAX_WARNINGS
    );
    if (ret == 0) return -1;  // no matching peaks found

    /* For file init, amplitude_keys come from the file; use basis modes as fallback */
    char raw_keys[MAX_MODES][32];
    const char *amplitude_keys[MAX_MODES];
    double amplitudes[MAX_MODES];
    memset(raw_keys, 0, sizeof(raw_keys));
    memset(amplitudes, 0, sizeof(amplitudes));
    int na = ctx->basis.modes_count;
    if (na > MAX_MODES) na = MAX_MODES;
    for (int i = 0; i < na; i++) {
        strncpy(raw_keys[i], ctx->basis.modes[i].family_key, 31);
        raw_keys[i][31] = '\0';
        amplitude_keys[i] = raw_keys[i];
        amplitudes[i] = out_coeffs[i].coeffs[0].real;
    }

    return build_initialization_result(
        &ctx->basis,
        out_coeffs,
        na,
        amplitudes,
        amplitude_keys,
        na,
        result
    );
}

/* ========================================================================
 * Transform Miller indices between lattices
 * ======================================================================== */

int engine_transform_miller(
    const LatticeInfo *lattice_A,
    const double P[3][3],
    FullInitializationResult *result)
{
    if (!lattice_A || !P || !result) return -1;

    /* Collect existing keys and coeffs from the result */
    int num_keys = result->n_coeffs;
    if (num_keys <= 0) return -1;

    int existing_keys[MAX_AMPLITUDE_KEYS][3];
    double existing_coeffs_real[MAX_AMPLITUDE_KEYS][MAX_STAR_VECTORS];
    double existing_coeffs_imag[MAX_AMPLITUDE_KEYS][MAX_STAR_VECTORS];
    int existing_star_counts[MAX_AMPLITUDE_KEYS];

    memset(existing_coeffs_real, 0, sizeof(existing_coeffs_real));
    memset(existing_coeffs_imag, 0, sizeof(existing_coeffs_imag));

    for (int i = 0; i < num_keys && i < MAX_AMPLITUDE_KEYS; i++) {
        existing_keys[i][0] = result->trans_keys[i][0];
        existing_keys[i][1] = result->trans_keys[i][1];
        existing_keys[i][2] = result->trans_keys[i][2];
        int cnt = result->trans_counts[i];
        existing_star_counts[i] = cnt;
        for (int j = 0; j < MAX_STAR_VECTORS; j++) {
            existing_coeffs_real[i][j] = (j < cnt) ? result->trans_coeffs_real[i][j] : 0.0;
            existing_coeffs_imag[i][j] = (j < cnt) ? result->trans_coeffs_imag[i][j] : 0.0;
        }
    }

    /* Transform uses MAX_STAR_VECTORS (100) buffers — copy back after */
    double tmp_coeffs_real[MAX_AMPLITUDE_KEYS][MAX_STAR_VECTORS];
    double tmp_coeffs_imag[MAX_AMPLITUDE_KEYS][MAX_STAR_VECTORS];
    int tmp_counts[MAX_AMPLITUDE_KEYS];
    int tmp_keys[MAX_AMPLITUDE_KEYS][3];

    memset(tmp_coeffs_real, 0, sizeof(tmp_coeffs_real));
    memset(tmp_coeffs_imag, 0, sizeof(tmp_coeffs_imag));
    memset(tmp_counts, 0, sizeof(tmp_counts));
    memset(tmp_keys, 0, sizeof(tmp_keys));

    /* Transform — outputs go into temp buffers */
    transform_miller_between_lattices(
        (LatticeInfo*)lattice_A,
        P,
        existing_keys,
        num_keys,
        existing_coeffs_real,
        existing_coeffs_imag,
        existing_star_counts,
        tmp_keys,
        tmp_coeffs_real,
        tmp_coeffs_imag,
        tmp_counts,
        &result->trans_lattice_a,
        &result->trans_lattice_b,
        &result->trans_lattice_c,
        &result->trans_lattice_alpha,
        &result->trans_lattice_beta,
        &result->trans_lattice_gamma
    );

    /* Copy results back to result struct (which uses MAX_VEC_PER_STAR) */
    for (int i = 0; i < MAX_AMPLITUDE_KEYS; i++) {
        result->trans_keys[i][0] = tmp_keys[i][0];
        result->trans_keys[i][1] = tmp_keys[i][1];
        result->trans_keys[i][2] = tmp_keys[i][2];
        int cnt = tmp_counts[i];
        result->trans_counts[i] = cnt;
        int copy = cnt < MAX_VEC_PER_STAR ? cnt : MAX_VEC_PER_STAR;
        for (int j = 0; j < copy; j++) {
            result->trans_coeffs_real[i][j] = tmp_coeffs_real[i][j];
            result->trans_coeffs_imag[i][j] = tmp_coeffs_imag[i][j];
        }
        /* Zero out remaining */
        for (int j = copy; j < MAX_VEC_PER_STAR; j++) {
            result->trans_coeffs_real[i][j] = 0.0;
            result->trans_coeffs_imag[i][j] = 0.0;
        }
    }

    result->trans_valid = 1;
    return 0;
}

/* ========================================================================
 * Output field — write VTK file
 * ======================================================================== */

int engine_output_field(
    const char *filename,
    const char *field_name,
    int apply_tile,
    int tile_x, int tile_y, int tile_z,
    FullInitializationResult *result,
    double resol,
    int transform_coord)
{
    if (!filename || !field_name || !result) return -1;

    return build_field(
        filename,
        field_name,
        apply_tile,
        result,
        resol,
        transform_coord,
        tile_x, tile_y, tile_z
    );
}

/* ========================================================================
 * Full pipeline — convenience: create engine → random init → write field
 * ======================================================================== */

int engine_full_pipeline(
    const char *sg_symbol,
    const char *symops_path,
    const LatticeInfo *lattice,
    int n_keep,
    DistributionType dist,
    DistParams dist_params,
    uint64_t rng_seed,
    const char *output_vts,
    const char *field_name,
    double resol,
    int apply_tile,
    int tile_x, int tile_y, int tile_z)
{
    /* Suppress unused parameter warnings when dist/dist_params not used */
    (void)dist;
    (void)dist_params;
    (void)rng_seed;

    EngineContext ctx;
    int ret = engine_create(&ctx, sg_symbol, symops_path, lattice, n_keep);
    if (ret != 0) return -1;

    /* Random init */
    FullInitializationResult result;
    memset(&result, 0, sizeof(result));
    ret = engine_random_init(&ctx, &result);
    if (ret != 0) {
        engine_free(&ctx);
        return -1;
    }

    /* Output field */
    if (output_vts) {
        ret = engine_output_field(
            output_vts, field_name, apply_tile,
            tile_x, tile_y, tile_z, &result, resol, 0
        );
    }

    engine_free(&ctx);
    return ret;
}
