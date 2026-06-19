/**
 * demo_random_init.c — Random initialization via one-liner pipeline
 *
 * This is the fastest way to get a symmetric field:
 *   engine_full_pipeline() does everything in one call:
 *     1. Loads the space group from file
 *     2. Builds the SAFB basis (top N modes by q²)
 *     3. Assigns random amplitudes (deterministic via seed)
 *     4. Generates the real-space field via iFFT
 *     5. Writes VTK output
 *
 * Perfect for:
 *   - Quick prototyping
 *   - SCFT simulation initialization
 *   - Exploring different space groups
 *
 * Build:
 *   source /sandbox/setup_build.sh && gcc -o demo_random_init demo_random_init.c \
 *     ../src/domain.c ../src/symmetry_ops.c ../src/basis.c ../src/initializers.c \
 *     ../src/field.c ../src/engine.c -std=c11 -Wall -O2 -I../include \
 *     -lm -lfftw3 -L/sandbox/miniforge3/envs/build/lib
 *
 * Run:
 *   ./demo_random_init [space_group_symbol] [ops_file] [lattice_dim] [N] [output.vts]
 *
 * Example:
 *   ./demo_random_init Pm-3m ../examples/space_groups_3d/Cubic/P_m_-3_m.txt 4.0 5 field.vts
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../include/engine.h"

/* --- Space group files bundled for convenience --- */
#define P1_PATH     "../examples/space_groups_3d/Cubic/P_1.txt"
#define PM3M_PATH   "../examples/space_groups_3d/Cubic/P_m_-3_m.txt"
#define IA3D_PATH   "../examples/space_groups_3d/Cubic/I_a_-3_d.txt"
#define P4MMC_PATH  "../examples/space_groups_3d/Tetragonal/P_42_m_m_c.txt"

int main(int argc, char *argv[])
{
    printf("=== SAFB Demo: Random Init (One-Liner Pipeline) ===\n\n");

    /* --- Parse command-line arguments --- */
    const char *sg_symbol   = "Ia-3d";          /* default space group */
    const char *ops_path    = IA3D_PATH;        /* default ops file */
    double      lattice_dim = 4.0;              /* default lattice constant */
    int         n_modes     = 5;                /* default: top 5 modes */
    const char *output_file = "random_field.vts";
    int         rng_seed    = 42;               /* deterministic by default */

    if (argc >= 6) {
        sg_symbol     = argv[1];
        ops_path      = argv[2];
        lattice_dim   = atof(argv[3]);
        n_modes       = atoi(argv[4]);
        output_file   = argv[5];
        if (argc >= 7) rng_seed = atoi(argv[6]);
    }

    printf("Space group: %s\n", sg_symbol);
    printf("Ops file:    %s\n", ops_path);
    printf("Lattice:     cubic, a = %.1f\n", lattice_dim);
    printf("Modes (N):   %d\n", n_modes);
    printf("Output:      %s\n", output_file);
    printf("RNG seed:    %d\n\n", rng_seed);

    /* --- The one-liner! --- */
    /*
     * Distributions:
     *   DIST_UNIFORM  — uniform [0, 1]
     *   DIST_NORMAL   — normal (mean=loc, stddev=scale)
     *   DIST_EXPONENTIAL — exponential (rate=1/loc)
     *   DIST_LOGNORMAL — log-normal (mean/sigma of log)
     */

    int ret = engine_full_pipeline(
        sg_symbol,                     /* space group */
        ops_path,                      /* ops file */
        &(LatticeInfo){lattice_dim, lattice_dim, lattice_dim,
                       90.0, 90.0, 90.0, 3},  /* cubic lattice */
        n_modes,                       /* top N modes */
        DIST_NORMAL,                   /* random distribution */
        (DistParams){1.0, 0.5},        /* mean=1.0, stddev=0.5 */
        (uint64_t)rng_seed,            /* deterministic seed */
        output_file,                   /* output VTK file */
        "psi",                         /* field name */
        1.0,                           /* resolution */
        0, 1, 1, 1                     /* no tiling */
    );

    if (ret == 0) {
        printf("\n✅ Success! Field written to: %s\n", output_file);
        printf("Open in ParaView: File → Open → %s → Apply\n", output_file);
    } else {
        fprintf(stderr, "\n❌ Error: pipeline failed (ret=%d)\n", ret);
        return 1;
    }

    /* --- Bonus: Try a different space group --- */
    printf("\n--- Trying Pm-3m (simple cubic) ---\n");
    ret = engine_full_pipeline(
        "Pm-3m",
        PM3M_PATH,
        &(LatticeInfo){lattice_dim, lattice_dim, lattice_dim,
                       90.0, 90.0, 90.0, 3},
        3,                           /* fewer modes for this group */
        DIST_UNIFORM,
        (DistParams){0.0, 1.0},
        123,
        "pm3m_field.vts",
        "density",
        1.0, 0, 1, 1, 1            /* apply_tile, tile_x, tile_y, tile_z */
    );
    if (ret == 0) {
        printf("✅ Pm-3m field written to: pm3m_field.vts\n");
    } else {
        fprintf(stderr, "⚠ Pm-3m pipeline failed (ret=%d)\n", ret);
    }

    printf("\nDone. ✅\n");
    return 0;
}
