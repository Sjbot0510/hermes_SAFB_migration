/**
 * engine.h — High-level SAFB initialization engine
 *
 * Translated from: Sg_init/engine.py
 *
 * Procedural equivalents of SpaceGroupInitializationEngine methods:
 *   - engine_build_basis: build SAFB basis from ops, lattice, mode count
 *   - engine_random_init: random amplitude initialization
 *   - engine_manual_init: manual amplitude initialization
 *   - engine_file_init: file-based amplitude initialization
 *   - engine_build_field: generate field via iFFT + VTK export
 *   - engine_transform_miller: transform Miller indices between lattices
 *
 * Usage: chain the functions together to replicate the engine workflow.
 */

#ifndef ENGINE_H
#define ENGINE_H

#include <stdint.h>
#include <stddef.h>
#include "domain.h"
#include "initializers.h"
#include "symmetry_ops.h"
#include "basis.h"
#include "field.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ========================================================================
 * Engine context — holds reusable state (ops, lattice, etc.)
 * ======================================================================== */

typedef struct {
    SymmGroup ops;
    char space_group[32];
    LatticeInfo lattice;
    int n_keep;
    SAFBBasis basis;  /* pre-built basis */
} EngineContext;

/* Initialize engine context from a space group ops file and lattice */
int engine_create(EngineContext *ctx,
                  const char *sg_symbol,
                  const char *symops_path,
                  const LatticeInfo *lattice,
                  int n_keep);

/* Free engine context (currently no dynamically allocated members in basis) */
void engine_free(EngineContext *ctx);

/* ========================================================================
 * Build basis — replicate engine.build_basis()
 * ======================================================================== */

int engine_build_basis(EngineContext *ctx, int N);

/* ========================================================================
 * Random initialization — replicate engine.init_random()
 * ======================================================================== */

int engine_random_init(
    const EngineContext *ctx,
    FullInitializationResult *result
);

/* ========================================================================
 * Manual initialization — replicate engine.init_manual()
 *
 * amplitudes: array of amplitude values
 * amplitude_keys: parallel array of family_key strings
 * n_amplitudes: number of entries
 * ======================================================================== */

int engine_manual_init(
    const EngineContext *ctx,
    const char *const amplitude_keys[],
    const double amplitudes[],
    int n_amplitudes,
    FullInitializationResult *result
);

/* ========================================================================
 * File initialization — replicate engine.init_from_file()
 * ======================================================================== */

int engine_file_init(
    const EngineContext *ctx,
    const char *scattering_file,
    const LatticeInfo *read_lattice_info,
    const double P[3][3],  /* change-of-basis matrix (NULL = identity) */
    FullInitializationResult *result
);

/* ========================================================================
 * Transform Miller indices between lattices — replicate transform_lattice_coordinate()
 * ======================================================================== */

int engine_transform_miller(
    const LatticeInfo *lattice_A,
    const double P[3][3],  /* change-of-basis matrix (B = A @ P) */
    FullInitializationResult *result
);

/* ========================================================================
 * Output field — replicate Output_field() / build_field()
 * ======================================================================== */

int engine_output_field(
    const char *filename,
    const char *field_name,
    int apply_tile,
    int tile_x, int tile_y, int tile_z,
    FullInitializationResult *result,
    double resol,
    int transform_coord
);

/* ========================================================================
 * Full pipeline convenience — one call for the common random init workflow
 *
 * Returns 0 on success, -1 on error.
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
    int tile_x, int tile_y, int tile_z
);

#ifdef __cplusplus
}
#endif

#endif /* ENGINE_H */
