/**
 * demo_manual_init.c — Manual amplitude initialization for Ia-3d (gyroid)
 *
 * This demo shows how to:
 *   1. Load the Ia-3d space group (96 symmetry operations)
 *   2. Build the SAFB basis (5 modes: {222}, {400}, {420}, {440}, {531})
 *   3. Manually set equal amplitudes for all modes
 *   4. Generate the real-space field via iFFT
 *   5. Write the result as a VTK file (visualize in ParaView)
 *
 * Build:
 *   source /sandbox/setup_build.sh && gcc -o demo_manual_init demo_manual_init.c \
 *     ../src/domain.c ../src/symmetry_ops.c ../src/basis.c ../src/initializers.c \
 *     ../src/field.c ../src/engine.c -std=c11 -Wall -O2 -I../include \
 *     -lm -lfftw3 -L/sandbox/miniforge3/envs/build/lib
 *
 * Run:
 *   ./demo_manual_init
 *
 * Output:
 *   gyroid.vts — load in ParaView to see the gyroid bicontinuous structure
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../include/engine.h"
#include "../include/domain.h"

#define OPS_PATH "../examples/space_groups_3d/Cubic/I_a_-3_d.txt"

int main(void)
{
    printf("=== SAFB Demo: Manual Init for Ia-3d (Gyroid) ===\n\n");

    /* --- Step 1: Define cubic lattice (a = 4.0 Å) --- */
    LatticeInfo lattice = lattice_info_new(4.0, 4.0, 4.0,
                                            90.0, 90.0, 90.0, 3);
    printf("Lattice: cubic, a = %.1f\n", lattice.a);

    /* --- Step 2: Load space group and build basis --- */
    EngineContext ctx;
    int ret = engine_create(&ctx, "Ia-3d", OPS_PATH, &lattice, 5);
    if (ret != 0) {
        fprintf(stderr, "ERROR: failed to create engine\n");
        return 1;
    }

    printf("Space group: %s (%d ops, %d modes)\n",
           ctx.space_group, ctx.ops.count, ctx.basis.modes_count);
    printf("Centrosymmetric: %s\n",
           ctx.basis.centrosymmetric_group ? "yes" : "no");

    /* Print each mode */
    for (int i = 0; i < ctx.basis.modes_count; i++) {
        Star *m = &ctx.basis.modes[i];
        printf("  Mode %d: %-6s  q²=%.4f  multiplicity=%d  closed=%s\n",
               i, m->family_key, m->q2, m->multiplicity,
               m->star_close ? "yes" : "no");
    }
    printf("\n");

    /* --- Step 3: Manual initialization --- */
    /*
     * Set amplitude = 1.0 for every mode.
     * The family keys must match what the basis built.
     */
    const char *keys[] = {"{222}", "{400}", "{420}", "{440}", "{531}"};
    double amps[]    = {1.0, 1.0, 1.0, 1.0, 1.0};

    FullInitializationResult result;
    memset(&result, 0, sizeof(result));

    ret = engine_manual_init(&ctx,
                             (const char *const *)keys,
                             amps,
                             5,
                             &result);
    if (ret != 0) {
        fprintf(stderr, "ERROR: manual initialization failed\n");
        engine_free(&ctx);
        return 1;
    }

    printf("Amplitude keys and values:\n");
    for (int i = 0; i < 5; i++) {
        printf("  %-6s → %.2f\n", keys[i], amps[i]);
    }
    printf("\n");

    /* --- Step 4: Generate field (resolution = 1.0) --- */
    /* resol = 1.0 means grid points ≈ lattice_dim / resol = 4/1.0 = 4 per axis */
    ret = engine_output_field(
        "gyroid.vts",    /* output filename */
        "psi",           /* scalar field name */
        0,               /* no tiling */
        1, 1, 1,
        &result,         /* initialization result */
        1.0,             /* resolution */
        0                /* no lattice transform */
    );
    if (ret != 0) {
        fprintf(stderr, "ERROR: field generation failed\n");
        engine_free(&ctx);
        return 1;
    }

    printf("Field written to: gyroid.vts\n");
    printf("(Load in ParaView: File → Open → gyroid.vts → Apply)\n\n");

    /* --- Cleanup --- */
    engine_free(&ctx);

    printf("Done. ✅\n");
    return 0;
}
