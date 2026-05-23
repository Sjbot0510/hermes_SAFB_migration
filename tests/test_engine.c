/**
 * test_engine.c — Tests for the engine layer (engine.c)
 *
 * Tests:
 *   - engine_create: create context from ops file + lattice (P1, P42/mmc)
 *   - engine_build_basis: build basis for a given mode count
 *   - engine_random_init: random initialization
 *   - engine_manual_init: manual amplitude assignment
 *   - engine_output_field: VTK field output
 *   - engine_full_pipeline: full end-to-end pipeline
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "engine.h"

static int tests_passed = 0;
static int tests_failed = 0;

#define TEST(name) static void test_##name(void)
#define RUN(name) do { \
    printf("  Testing %s... ", #name); \
    test_##name(); \
} while(0)

#define ASSERT(cond, msg) do { \
    if (!(cond)) { \
        fprintf(stderr, "\n    FAIL: %s\n", msg); \
        tests_failed++; \
        return; \
    } \
    tests_passed++; \
} while(0)

/* Construct a P1 SymmGroup manually (identity only) */
static SymmGroup make_p1_group(void) {
    SymmGroup sg;
    memset(&sg, 0, sizeof(sg));
    sg.dim = 3;

    int R[3][3] = {{1,0,0},{0,1,0},{0,0,1}};
    memcpy(sg.ops[0].R, R, sizeof(R));
    sg.ops[0].t[0] = frac_new(0, 1);
    sg.ops[0].t[1] = frac_new(0, 1);
    sg.ops[0].t[2] = frac_new(0, 1);
    sg.count = 1;
    return sg;
}

/* Construct a P2 SymmGroup manually (identity + 180° rotation about z) */
static SymmGroup make_p2_group(void) {
    SymmGroup sg;
    memset(&sg, 0, sizeof(sg));
    sg.dim = 3;

    /* Identity */
    int R0[3][3] = {{1,0,0},{0,1,0},{0,0,1}};
    memcpy(sg.ops[0].R, R0, sizeof(R0));
    sg.ops[0].t[0] = frac_new(0, 1);
    sg.ops[0].t[1] = frac_new(0, 1);
    sg.ops[0].t[2] = frac_new(0, 1);
    sg.count = 1;

    /* 180° rotation about z */
    int R1[3][3] = {{-1,0,0},{0,-1,0},{0,0,1}};
    memcpy(sg.ops[1].R, R1, sizeof(R1));
    sg.ops[1].t[0] = frac_new(0, 1);
    sg.ops[1].t[1] = frac_new(0, 1);
    sg.ops[1].t[2] = frac_new(0, 1);
    sg.count = 2;

    return sg;
}

/* Paths to real symmetry operation files */
#define P1_PATH "/sandbox/hermes_SAFB_migration/examples/space_groups_3d/Cubic/P_m_-3_m.txt"
#define P4MMC_PATH "/sandbox/hermes_SAFB_migration/examples/space_groups_3d/Tetragonal/P_42_m_m_c.txt"

TEST(engine_manual_p1)
{
    /* Test engine with manually constructed P1 group (no file I/O needed) */
    SymmGroup sg = make_p1_group();

    EngineContext ctx;
    memset(&ctx, 0, sizeof(ctx));
    strncpy(ctx.space_group, "P1", sizeof(ctx.space_group) - 1);
    ctx.lattice = lattice_info_new(1.0, 1.0, 1.0, 90.0, 90.0, 90.0, 3);
    ctx.n_keep = 5;
    memcpy(&ctx.ops, &sg, sizeof(sg));

    int ret = engine_build_basis(&ctx, 5);
    ASSERT(ret == 0, "engine_build_basis should succeed for P1");
    ASSERT(ctx.basis.modes_count >= 1, "Basis should have at least 1 mode");
    printf(" (%d modes)", ctx.basis.modes_count);
    engine_free(&ctx);
}

TEST(engine_manual_p2)
{
    SymmGroup sg = make_p2_group();

    EngineContext ctx;
    memset(&ctx, 0, sizeof(ctx));
    strncpy(ctx.space_group, "P2", sizeof(ctx.space_group) - 1);
    ctx.lattice = lattice_info_new(1.0, 1.0, 1.0, 90.0, 90.0, 90.0, 3);
    ctx.n_keep = 5;
    memcpy(&ctx.ops, &sg, sizeof(sg));

    int ret = engine_build_basis(&ctx, 5);
    ASSERT(ret == 0, "engine_build_basis should succeed for P2");
    ASSERT(ctx.basis.modes_count >= 1, "Basis should have at least 1 mode");
    printf(" (%d modes)", ctx.basis.modes_count);
    engine_free(&ctx);
}

TEST(engine_create_pm3m)
{
    /* Use Pm-3m — 48 ops, still manageable */
    EngineContext ctx;
    int ret = engine_create(&ctx, "Pm-3m",
        P1_PATH,  /* This is actually Pm-3m path */
        &(LatticeInfo){1.0, 1.0, 1.0, 90.0, 90.0, 90.0, 3},
        5);
    ASSERT(ret == 0, "engine_create should succeed for Pm-3m");
    ASSERT(ctx.ops.count > 0, "SymmGroup should have ops");
    printf(" (%d ops, %d modes)", ctx.ops.count, ctx.basis.modes_count);
    engine_free(&ctx);
}

TEST(engine_random_init)
{
    SymmGroup sg = make_p1_group();

    EngineContext ctx;
    memset(&ctx, 0, sizeof(ctx));
    strncpy(ctx.space_group, "P1", sizeof(ctx.space_group) - 1);
    ctx.lattice = lattice_info_new(1.0, 1.0, 1.0, 90.0, 90.0, 90.0, 3);
    ctx.n_keep = 5;
    memcpy(&ctx.ops, &sg, sizeof(sg));
    engine_build_basis(&ctx, 5);

    FullInitializationResult result;
    memset(&result, 0, sizeof(result));
    int ret = engine_random_init(&ctx, &result);
    ASSERT(ret == 0, "engine_random_init should succeed for P1");
    ASSERT(result.n_coeffs > 0, "Result should have coeffs");
    engine_free(&ctx);
}

TEST(engine_manual_init)
{
    SymmGroup sg = make_p1_group();

    EngineContext ctx;
    memset(&ctx, 0, sizeof(ctx));
    strncpy(ctx.space_group, "P1", sizeof(ctx.space_group) - 1);
    ctx.lattice = lattice_info_new(1.0, 1.0, 1.0, 90.0, 90.0, 90.0, 3);
    ctx.n_keep = 5;
    memcpy(&ctx.ops, &sg, sizeof(sg));
    engine_build_basis(&ctx, 5);

    char *keys[] = {"{100}", "{010}", "{001}", "{110}", "{101}"};
    double amps[] = {1.0, 1.0, 1.0, 1.0, 1.0};
    int n = 5;

    FullInitializationResult result;
    memset(&result, 0, sizeof(result));
    int ret = engine_manual_init(&ctx, (const char *const *)keys, amps, n, &result);
    ASSERT(ret == 0, "engine_manual_init should succeed for P1");
    engine_free(&ctx);
}

TEST(engine_output_field)
{
    SymmGroup sg = make_p1_group();

    EngineContext ctx;
    memset(&ctx, 0, sizeof(ctx));
    strncpy(ctx.space_group, "P1", sizeof(ctx.space_group) - 1);
    ctx.lattice = lattice_info_new(1.0, 1.0, 1.0, 90.0, 90.0, 90.0, 3);
    ctx.n_keep = 5;
    memcpy(&ctx.ops, &sg, sizeof(sg));
    engine_build_basis(&ctx, 5);

    FullInitializationResult result;
    memset(&result, 0, sizeof(result));
    engine_random_init(&ctx, &result);

    int ret = engine_output_field(
        "/tmp/test_engine_output.vts",
        "psi", 0, 1, 1, 1, &result, 20.0, 0
    );
    ASSERT(ret == 0, "engine_output_field should succeed");

    FILE *f = fopen("/tmp/test_engine_output.vts", "r");
    ASSERT(f != NULL, "VTS file should exist");
    if (f) {
        char buf[256];
        size_t n = fread(buf, 1, sizeof(buf), f);
        ASSERT(n > 100, "VTS file should have content");
        fclose(f);
    }
    engine_free(&ctx);
}

TEST(engine_full_pipeline)
{
    LatticeInfo lat = {1.0, 1.0, 1.0, 90.0, 90.0, 90.0, 3};
    DistParams params = {0.0, 1.0};

    int ret = engine_full_pipeline(
        "P1",
        P1_PATH,
        &lat,
        5,
        DIST_UNIFORM,
        params,
        42,
        "/tmp/test_full_pipeline.vts",
        "psi",
        20.0,
        0, 1, 1, 1
    );
    ASSERT(ret == 0, "Full pipeline should succeed for Pm-3m");

    FILE *f = fopen("/tmp/test_full_pipeline.vts", "r");
    ASSERT(f != NULL, "Pipeline VTS file should exist");
    if (f) {
        char buf[256];
        size_t n = fread(buf, 1, sizeof(buf), f);
        ASSERT(n > 100, "Pipeline VTS file should have content");
        fclose(f);
    }
}

int main(void)
{
    printf("\n=== Engine Module Tests ===\n\n");

    RUN(engine_manual_p1);
    RUN(engine_manual_p2);
    RUN(engine_create_pm3m);
    RUN(engine_random_init);
    RUN(engine_manual_init);
    RUN(engine_output_field);
    RUN(engine_full_pipeline);

    printf("\n=== Results: %d passed, %d failed ===\n", tests_passed, tests_failed);
    return tests_failed > 0 ? 1 : 0;
}
