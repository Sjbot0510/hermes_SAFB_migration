/**
 * benchmark.c — C SAFB benchmark suite
 *
 * Measures wall-clock time for each pipeline phase across space groups.
 * Also reports RSS memory usage.
 *
 * Usage: ./benchmark [space_group] [n_keep] [resolution] [n_iterations]
 *   Defaults: all space groups, n_keep=5, resol=0.5, 3 iterations
 */

#define _POSIX_C_SOURCE 200809L
#include "engine.h"
#include "analytic.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ========================================================================
 * Timing
 * ======================================================================== */

static double get_time_ms(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec * 1000.0 + ts.tv_nsec / 1e6;
}

static long get_rss_kb(void) {
    long rss = 0;
    FILE *f = fopen("/proc/self/status", "r");
    if (f) {
        char line[256];
        while (fgets(line, sizeof(line), f)) {
            if (strncmp(line, "VmRSS:", 6) == 0) {
                sscanf(line + 6, "%ld", &rss);
                break;
            }
        }
        fclose(f);
    }
    return rss;
}

/* ========================================================================
 * Test case definitions
 * ======================================================================== */

typedef struct {
    const char *sg_symbol;
    const char *symops_path;
    int dim;
    double a, b, c;
    double alpha, beta, gamma;
    int n_keep;
} TestCase;

static TestCase test_cases[] = {
    {"P4",        "examples/space_groups_2d/square/p_4.txt",        2, 4.0, 4.0, 1.0, 90.0, 90.0, 90.0, 5},
    {"Pm-3m",     "examples/space_groups_3d/Cubic/P_m_-3_m.txt",    3, 4.0, 4.0, 4.0, 90.0, 90.0, 90.0, 5},
    {"Ia-3d",     "examples/space_groups_3d/Cubic/I_a_-3_d.txt",    3, 4.0, 4.0, 4.0, 90.0, 90.0, 90.0, 5},
    {"Fm-3m",     "examples/space_groups_3d/Cubic/F_m_-3_m.txt",    3, 4.0, 4.0, 4.0, 90.0, 90.0, 90.0, 5},
    {"Ia-3d-10",  "examples/space_groups_3d/Cubic/I_a_-3_d.txt",    3, 4.0, 4.0, 4.0, 90.0, 90.0, 90.0, 10},
    {"Ia-3d-15",  "examples/space_groups_3d/Cubic/I_a_-3_d.txt",    3, 4.0, 4.0, 4.0, 90.0, 90.0, 90.0, 15},
};
static int n_test_cases = sizeof(test_cases) / sizeof(test_cases[0]);

/* ========================================================================
 * Benchmark one pipeline run
 * ======================================================================== */

typedef struct {
    double basis_ms;
    double init_ms;
    double field_ms;
    double analyt_ms;
    double total_ms;
    long rss_kb;
    int ok;
} BenchResult;

static int run_benchmark(const TestCase *tc,
                         double resol,
                         int n_iters,
                         BenchResult *out)
{
    memset(out, 0, sizeof(BenchResult));

    LatticeInfo lattice;
    if (tc->dim == 2) {
        lattice = lattice_info_new_2d(tc->a, tc->b, tc->gamma);
    } else {
        lattice = lattice_info_new(tc->a, tc->b, tc->c,
                                    tc->alpha, tc->beta, tc->gamma, tc->dim);
    }

    double total_basis = 0, total_init = 0, total_field = 0, total_analyt = 0;

    for (int iter = 0; iter < n_iters; iter++) {
        EngineContext ctx;
        memset(&ctx, 0, sizeof(ctx));

        double t0 = get_time_ms();
        int ret = engine_create(&ctx, tc->sg_symbol, tc->symops_path, &lattice, tc->n_keep);
        double basis_ms = get_time_ms() - t0;

        if (ret != 0) {
            fprintf(stderr, "  FAIL: engine_create for %s (iter %d)\n", tc->sg_symbol, iter);
            return -1;
        }

        /* Manual init */
        FullInitializationResult result;
        memset(&result, 0, sizeof(result));

        int n_modes = ctx.basis.modes_count;
        if (n_modes > MAX_MODES) n_modes = MAX_MODES;

        char mode_keys[MAX_MODES][32];
        double mode_amps[MAX_MODES];
        const char *keys[MAX_MODES];

        for (int i = 0; i < n_modes; i++) {
            strncpy(mode_keys[i], ctx.basis.modes[i].family_key, 31);
            mode_keys[i][31] = '\0';
            keys[i] = mode_keys[i];
            mode_amps[i] = 1.0;
        }

        t0 = get_time_ms();
        ret = engine_manual_init(&ctx, keys, mode_amps, n_modes, &result);
        double init_ms = get_time_ms() - t0;

        if (ret != 0) {
            fprintf(stderr, "  FAIL: manual_init for %s (iter %d)\n", tc->sg_symbol, iter);
            engine_free(&ctx);
            return -1;
        }

        /* Field generation */
        t0 = get_time_ms();
        char tmpfile[256];
        snprintf(tmpfile, sizeof(tmpfile), "/tmp/bench_tmp_%d_%d.vts", iter, (int)t0);

        ret = engine_output_field(
            tmpfile, "field", 0, 1, 1, 1,
            &result, resol, 0
        );
        double field_ms = get_time_ms() - t0;

        /* Analytical star function */
        t0 = get_time_ms();
        AnalyticBasis ab;
        memset(&ab, 0, sizeof(ab));
        ret = derive_analytical_star_function(&result, 32, &ab);
        double analyt_ms = get_time_ms() - t0;

        /* Clean up tmpfile */
        remove(tmpfile);

        engine_free(&ctx);

        total_basis += basis_ms;
        total_init += init_ms;
        total_field += field_ms;
        total_analyt += analyt_ms;
    }

    out->basis_ms = total_basis / n_iters;
    out->init_ms = total_init / n_iters;
    out->field_ms = total_field / n_iters;
    out->analyt_ms = total_analyt / n_iters;
    out->total_ms = out->basis_ms + out->init_ms + out->field_ms + out->analyt_ms;
    out->rss_kb = get_rss_kb();
    out->ok = 1;
    return 0;
}

/* ========================================================================
 * Main
 * ======================================================================== */

int main(int argc, char *argv[])
{
    printf("==========================================================\n");
    printf("  SAFB C Benchmark Suite\n");
    printf("  Times are per-run averages (3 iterations)\n");
    printf("==========================================================\n\n");

    int n_iters = 3;
    double resol = 0.5;

    if (argc > 1) {
        /* Specific space group requested */
        const char *target = argv[1];
        if (argc > 2) n_iters = atoi(argv[2]);
        if (n_iters < 1) n_iters = 1;
        if (argc > 3) resol = atof(argv[3]);
        if (resol <= 0) resol = 0.5;

        printf("Benchmarking %s (%d iterations, resol=%.2f)\n\n", target, n_iters, resol);

        for (int i = 0; i < n_test_cases; i++) {
            if (strcmp(test_cases[i].sg_symbol, target) == 0) {
                BenchResult res;
                if (run_benchmark(&test_cases[i], resol, n_iters, &res) == 0) {
                    printf("  %-10s  basis: %6.2fms  init: %6.2fms  "
                           "field: %6.2fms  analyt: %5.2fms  "
                           "total: %6.2fms  RSS: %ld KB\n",
                           test_cases[i].sg_symbol,
                           res.basis_ms, res.init_ms, res.field_ms,
                           res.analyt_ms, res.total_ms, res.rss_kb);
                }
                break;
            }
        }
        return 0;
    }

    /* Run all test cases */
    printf("Running all %d test cases (%d iterations each):\n\n", n_test_cases, n_iters);

    for (int i = 0; i < n_test_cases; i++) {
        BenchResult res;
        printf("Test %d/%d: %s (dim=%d, n_keep=%d)\n",
               i + 1, n_test_cases,
               test_cases[i].sg_symbol,
               test_cases[i].dim, test_cases[i].n_keep);

        if (run_benchmark(&test_cases[i], resol, n_iters, &res) == 0) {
            printf("  %-10s  basis: %6.2fms  init: %6.2fms  "
                   "field: %6.2fms  analyt: %5.2fms  "
                   "total: %6.2fms  RSS: %ld KB\n",
                   test_cases[i].sg_symbol,
                   res.basis_ms, res.init_ms, res.field_ms,
                   res.analyt_ms, res.total_ms, res.rss_kb);
        } else {
            printf("  %-10s  FAILED\n", test_cases[i].sg_symbol);
        }
    }

    printf("\n==========================================================\n");
    printf("  Benchmark complete.\n");
    printf("==========================================================\n");

    return 0;
}
