/**
 * demo_2d.c — 2D square lattice (p4 symmetry)
 *
 * This demo shows:
 *   1. 2D lattice setup (dim=2 instead of 3)
 *   2. Loading a 2D wallpaper group
 *   3. Building a basis and generating a 2D field
 *
 * Use case: quasi-2D systems like confined block-copolymer films,
 * surface-confined polymer assemblies, or 2D photonic crystals.
 *
 * Build:
 *   source /sandbox/setup_build.sh && gcc -o demo_2d demo_2d.c \
 *     ../src/domain.c ../src/symmetry_ops.c ../src/basis.c ../src/initializers.c \
 *     ../src/field.c ../src/engine.c -std=c11 -Wall -O2 -I../include \
 *     -lm -lfftw3 -L/sandbox/miniforge3/envs/build/lib
 *
 * Run:
 *   ./demo_2d
 *
 * Output:
 *   p4_field.vts — 2D square-symmetric field (open in ParaView as 2D slice)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../include/engine.h"
#include "../include/domain.h"

#define P4_PATH "../examples/space_groups_2d/square/p_4.txt"

int main(void)
{
    printf("=== SAFB Demo: 2D Square Lattice (p4 symmetry) ===\n\n");

    /* --- Step 1: Define a 2D hexagonal (square) lattice --- */
    /*
     * For 2D, only a, b, and gamma matter.
     * alpha and beta are ignored (kept at 90 for API compatibility).
     */
    double a = 4.0;    /* unit cell length in x */
    double b = 4.0;    /* unit cell length in y */
    double gamma = 90.0;  /* angle between a and b (degrees) */

    LatticeInfo lattice = lattice_info_new_2d(a, b, gamma);
    printf("Lattice: 2D square, a = %.1f, b = %.1f, γ = %.1f°\n",
           lattice.a, lattice.b, lattice.gamma);
    printf("Lattice dim: %d\n\n", lattice.dim);

    /* --- Step 2: Load 2D space group (p4 = 4-fold rotation) --- */
    EngineContext ctx;
    int ret = engine_create(&ctx, "p4", P4_PATH, &lattice, 3);
    if (ret != 0) {
        fprintf(stderr, "ERROR: failed to create engine\n");
        return 1;
    }

    printf("Wallpaper group: %s (%d ops, %d modes)\n",
           ctx.space_group, ctx.ops.count, ctx.basis.modes_count);
    printf("Centrosymmetric: %s\n\n",
           ctx.basis.centrosymmetric_group ? "yes" : "no");

    /* Print modes */
    for (int i = 0; i < ctx.basis.modes_count; i++) {
        Star *m = &ctx.basis.modes[i];
        printf("  Mode %d: %-6s  q²=%.4f  multiplicity=%d\n",
               i, m->family_key, m->q2, m->multiplicity);
    }
    printf("\n");

    /* --- Step 3: Manual initialization --- */
    /* p4 generates modes like {100}, {110}, {200} */
    const char *keys[] = {"{100}", "{110}", "{200}"};
    double amps[]    = {1.0, 0.5, 0.3};  /* decreasing amplitude with q² */

    FullInitializationResult result;
    memset(&result, 0, sizeof(result));

    ret = engine_manual_init(&ctx,
                             (const char *const *)keys,
                             amps,
                             3,
                             &result);
    if (ret != 0) {
        fprintf(stderr, "ERROR: manual initialization failed\n");
        engine_free(&ctx);
        return 1;
    }

    printf("Amplitude keys:\n");
    for (int i = 0; i < 3; i++) {
        printf("  %-6s → %.1f\n", keys[i], amps[i]);
    }
    printf("\n");

    /* --- Step 4: Generate 2D field --- */
    ret = engine_output_field(
        "p4_field.vts",      /* output */
        "density",           /* field name */
        0,                   /* no tiling */
        1, 1, 1,
        &result,
        1.0,                 /* resolution */
        0                    /* no transform */
    );
    if (ret != 0) {
        fprintf(stderr, "ERROR: field generation failed\n");
        engine_free(&ctx);
        return 1;
    }

    printf("2D field written to: p4_field.vts\n");
    printf("(Open in ParaView and use the 'Slice' filter to view as 2D)\n\n");

    /* --- Cleanup --- */
    engine_free(&ctx);

    printf("Done. ✅\n");
    return 0;
}
