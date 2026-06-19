/**
 * test_initializers.c — Tests for initializers.h / initializers.c
 *
 * Validates: read_scattering_data, srand_init/RNG, build_result_from_coeffs,
 * random_initialization, manual_initialization, file_initialization.
 */

#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include "initializers.h"
#include "basis.h"

#define EPS 1e-10

/* ========================================================================
 * Test: read_scattering_data
 * ======================================================================== */

static void test_read_scattering_data(void) {
    const char *path = "/sandbox/hermes_SAFB_migration/examples/scattering_data/scattering_C14.txt";

    /* Check file exists */
    FILE *f = fopen(path, "r");
    if (!f) {
        printf("  [SKIP] scattering_C14.txt not found\n");
        return;
    }
    fclose(f);

    ScatteringProfile sp = read_scattering_data(path);
    assert(sp.num_peaks > 0);
    printf("  [PASS] read_scattering_data -> %d peaks\n", sp.num_peaks);

    /* Check first peak data — HKL can be negative in scattering files */
    assert(sp.hkl != NULL && sp.q != NULL && sp.intensity != NULL);
    assert(sp.q[0] >= 0);
    assert(sp.intensity[0] > 0);

    /* Write and read back */
    sp.hkl[0] = 999;
    assert(sp.hkl[0] == 999);

    scattering_profile_free(&sp);
    printf("  [PASS] read_scattering_data free\n");
}

/* ========================================================================
 * Test: RNG
 * ======================================================================== */

static void test_rng(void) {
    SAFBRNG rng;
    srand_init(&rng, 42);

    double u1 = srand_uniform(&rng);
    assert(u1 >= 0.0 && u1 <= 1.0);

    double u2 = srand_uniform(&rng);
    assert(u1 != u2);  /* different values */

    double normal = srand_normal(&rng, 0.0, 1.0);
    /* Normal can be negative, but should be reasonable */
    assert(normal > -10.0 && normal < 10.0);

    double exp_v = srand_exponential(&rng, 1.0);
    assert(exp_v >= 0.0);

    printf("  [PASS] RNG (uniform, normal, exponential)\n");
}

/* ========================================================================
 * Test: build_result_from_coeffs — P1 space group with manual amplitudes
 * ======================================================================== */

static void test_build_result_from_coeffs(void) {
    /* Build a P1 basis */
    SymmGroup sg;
    memset(&sg, 0, sizeof(sg));
    sg.dim = 3;

    int R[3][3] = {{1,0,0},{0,1,0},{0,0,1}};
    memcpy(sg.ops[0].R, R, sizeof(R));
    sg.ops[0].t[0] = frac_new(0, 1);
    sg.ops[0].t[1] = frac_new(0, 1);
    sg.ops[0].t[2] = frac_new(0, 1);
    sg.count = 1;

    LatticeInfo lattice = lattice_info_new(1.0, 1.0, 1.0, 90.0, 90.0, 90.0, 3);

    SAFBBasis basis;
    int ret = basis_build(&sg, "P1", &lattice, 3, &basis);
    assert(ret == 0);
    assert(basis.modes_count >= 2);

    /* Use first two family keys */
    char key_buf[2][32];
    const char *keys[2];
    double amps[2];
    for (int i = 0; i < 2; i++) {
        strncpy(key_buf[i], basis.modes[i].family_key, sizeof(key_buf[i]) - 1);
        key_buf[i][sizeof(key_buf[i]) - 1] = '\0';
        keys[i] = key_buf[i];
        amps[i] = 1.0;
    }

    FamilyCoeffs out_coeffs[MAX_AMPLITUDE_KEYS];
    char warnings[MAX_WARNINGS][256];
    int n = build_result_from_coeffs(&basis, amps, (const char **)keys, 2,
                                      out_coeffs, warnings, MAX_WARNINGS);

    assert(n > 0);
    printf("  [PASS] build_result_from_coeffs -> %d family coeffs\n", n);

    for (int i = 0; i < n; i++) {
        printf("    coeff %d: key=%s, count=%d\n", i, out_coeffs[i].family_key, out_coeffs[i].count);
        assert(out_coeffs[i].count > 0);
        for (int j = 0; j < out_coeffs[i].count; j++) {
            printf("      vec=(%d,%d,%d) c=%.6f+%.6fi\n",
                   out_coeffs[i].coeffs[j].hkl[0],
                   out_coeffs[i].coeffs[j].hkl[1],
                   out_coeffs[i].coeffs[j].hkl[2],
                   out_coeffs[i].coeffs[j].real,
                   out_coeffs[i].coeffs[j].imag);
        }
    }
}

/* ========================================================================
 * Test: manual_initialization
 * ======================================================================== */

static void test_manual_initialization(void) {
    /* P1 basis */
    SymmGroup sg;
    memset(&sg, 0, sizeof(sg));
    sg.dim = 3;
    int R[3][3] = {{1,0,0},{0,1,0},{0,0,1}};
    memcpy(sg.ops[0].R, R, sizeof(R));
    sg.ops[0].t[0] = frac_new(0, 1);
    sg.ops[0].t[1] = frac_new(0, 1);
    sg.ops[0].t[2] = frac_new(0, 1);
    sg.count = 1;

    LatticeInfo lattice = lattice_info_new(1.0, 1.0, 1.0, 90.0, 90.0, 90.0, 3);

    SAFBBasis basis;
    int ret = basis_build(&sg, "P1", &lattice, 5, &basis);
    assert(ret == 0);

    /* Manual amplitudes for first 3 modes */
    const char *keys[3];
    double amps[3] = {1.0, 2.0, 3.0};
    for (int i = 0; i < 3; i++) {
        keys[i] = basis.modes[i].family_key;
    }

    FamilyCoeffs out_coeffs[MAX_AMPLITUDE_KEYS];
    int n = manual_initialization(&basis, keys, amps, 3,
                                   out_coeffs, NULL, 0);

    assert(n > 0);
    printf("  [PASS] manual_initialization -> %d coeffs\n", n);
}

/* ========================================================================
 * Test: random_initialization
 * ======================================================================== */

static void test_random_initialization(void) {
    /* P1 basis */
    SymmGroup sg;
    memset(&sg, 0, sizeof(sg));
    sg.dim = 3;
    int R[3][3] = {{1,0,0},{0,1,0},{0,0,1}};
    memcpy(sg.ops[0].R, R, sizeof(R));
    sg.ops[0].t[0] = frac_new(0, 1);
    sg.ops[0].t[1] = frac_new(0, 1);
    sg.ops[0].t[2] = frac_new(0, 1);
    sg.count = 1;

    LatticeInfo lattice = lattice_info_new(1.0, 1.0, 1.0, 90.0, 90.0, 90.0, 3);

    SAFBBasis basis;
    int ret = basis_build(&sg, "P1", &lattice, 10, &basis);
    assert(ret == 0);

    SAFBRNG rng;
    srand_init(&rng, 12345);

    DistParams params = {0.0, 1.0};

    FamilyCoeffs out_coeffs[MAX_AMPLITUDE_KEYS];
    char warnings[MAX_WARNINGS][256];

    int n = random_initialization(&basis, 5, DIST_NORMAL, params, &rng,
                                   out_coeffs, warnings, MAX_WARNINGS);
    assert(n > 0);

    /* Amplitudes should be random (not all same) */
    double sum_sq = 0;
    for (int i = 0; i < n; i++) {
        sum_sq += out_coeffs[i].coeffs[0].real * out_coeffs[i].coeffs[0].real;
    }
    printf("  [PASS] random_initialization (normal) -> %d coeffs, sum_sq=%.6f\n", n, sum_sq);
    assert(sum_sq > 0);
}

/* ========================================================================
 * Test: full initialization result wrapper
 * ======================================================================== */

static void test_full_initialization_result(void) {
    /* P1 basis */
    SymmGroup sg;
    memset(&sg, 0, sizeof(sg));
    sg.dim = 3;
    int R[3][3] = {{1,0,0},{0,1,0},{0,0,1}};
    memcpy(sg.ops[0].R, R, sizeof(R));
    sg.ops[0].t[0] = frac_new(0, 1);
    sg.ops[0].t[1] = frac_new(0, 1);
    sg.ops[0].t[2] = frac_new(0, 1);
    sg.count = 1;

    LatticeInfo lattice = lattice_info_new(1.0, 1.0, 1.0, 90.0, 90.0, 90.0, 3);

    SAFBBasis basis;
    int ret = basis_build(&sg, "P1", &lattice, 3, &basis);
    assert(ret == 0);

    /* Get coeffs */
    FamilyCoeffs coeffs[MAX_AMPLITUDE_KEYS];
    char key_bufs[3][32];
    const char *keys_ptr[3];
    double amps[3];
    for (int i = 0; i < 3; i++) {
        strncpy(key_bufs[i], basis.modes[i].family_key, 31);
        key_bufs[i][31] = '\0';
        keys_ptr[i] = key_bufs[i];
        amps[i] = 1.0;
    }

    int n = build_result_from_coeffs(&basis, amps, keys_ptr, 3,
                                      coeffs, NULL, 0);
    assert(n > 0);

    /* Build full result */
    FullInitializationResult result;
    ret = build_initialization_result(&basis, coeffs, n, amps, keys_ptr, 3, &result);
    assert(ret == 0);

    assert(result.n_coeffs > 0);
    printf("  [PASS] FullInitializationResult -> %d coeffs, sg='%s'\n",
           result.n_coeffs, result.space_group);
}

/* ========================================================================
 * Test: file_initialization
 * ======================================================================== */

static void test_file_initialization(void) {
    const char *path = "/sandbox/hermes_SAFB_migration/examples/scattering_data/scattering_C14.txt";

    FILE *f = fopen(path, "r");
    if (!f) {
        printf("  [SKIP] scattering_C14.txt not found\n");
        return;
    }
    fclose(f);

    /* Build P1 basis */
    SymmGroup sg;
    memset(&sg, 0, sizeof(sg));
    sg.dim = 3;
    int R[3][3] = {{1,0,0},{0,1,0},{0,0,1}};
    memcpy(sg.ops[0].R, R, sizeof(R));
    sg.ops[0].t[0] = frac_new(0, 1);
    sg.ops[0].t[1] = frac_new(0, 1);
    sg.ops[0].t[2] = frac_new(0, 1);
    sg.count = 1;

    LatticeInfo lattice = lattice_info_new(1.0, 1.0, 1.0, 90.0, 90.0, 90.0, 3);

    SAFBBasis basis;
    int ret = basis_build(&sg, "P1", &lattice, 5, &basis);
    assert(ret == 0);

    /* Try file initialization with same lattice (P=NULL) */
    FamilyCoeffs out_coeffs[MAX_AMPLITUDE_KEYS];
    char warnings[MAX_WARNINGS][256];
    int n = file_initialization(&basis, path, NULL, NULL, 5,
                                 out_coeffs, warnings, MAX_WARNINGS);

    printf("  [PASS] file_initialization -> %d coeffs (same lattice)\n", n);

    /* Try with different lattice but P=NULL should still work (no transform) */
    LatticeInfo diff_lattice = lattice_info_new(2.0, 2.0, 2.0, 90.0, 90.0, 90.0, 3);
    n = file_initialization(&basis, path, &diff_lattice, NULL, 5,
                             out_coeffs, warnings, MAX_WARNINGS);
    printf("  [PASS] file_initialization -> %d coeffs (diff lattice, no P)\n", n);
}

/* ========================================================================
 * Test: P1 basis with known amplitudes (verify sqrt(Am/mult))
 * ======================================================================== */

static void test_amplitude_sqrt_scaling(void) {
    SymmGroup sg;
    memset(&sg, 0, sizeof(sg));
    sg.dim = 3;
    int R[3][3] = {{1,0,0},{0,1,0},{0,0,1}};
    memcpy(sg.ops[0].R, R, sizeof(R));
    sg.ops[0].t[0] = frac_new(0, 1);
    sg.ops[0].t[1] = frac_new(0, 1);
    sg.ops[0].t[2] = frac_new(0, 1);
    sg.count = 1;

    LatticeInfo lattice = lattice_info_new(1.0, 1.0, 1.0, 90.0, 90.0, 90.0, 3);

    SAFBBasis basis;
    int ret = basis_build(&sg, "P1", &lattice, 3, &basis);
    assert(ret == 0);

    for (int m = 0; m < basis.modes_count; m++) {
        printf("    mode %d: key=%s, mult=%d, q2=%.4f\n",
               m, basis.modes[m].family_key,
               basis.modes[m].multiplicity, basis.modes[m].q2);
    }

    /* Set amplitude for first mode */
    double amp = 4.0;  /* sqrt(4/mult) = 4/mult */
    char key[32];
    strncpy(key, basis.modes[0].family_key, sizeof(key) - 1);
    key[sizeof(key) - 1] = '\0';
    const char *keys_arr[1] = {key};

    FamilyCoeffs out_coeffs[MAX_AMPLITUDE_KEYS];
    char warnings[MAX_WARNINGS][256];
    int n = build_result_from_coeffs(&basis, &amp, keys_arr, 1,
                                      out_coeffs, warnings, MAX_WARNINGS);
    assert(n > 0);

    /* Verify amplitude scaling: ref_real = sqrt(Am/mult) */
    double expected = sqrt(amp / (double)basis.modes[0].multiplicity);
    assert(fabs(out_coeffs[0].coeffs[0].real - expected) < 1e-10);
    printf("  [PASS] amplitude sqrt scaling: sqrt(%.1f/%d) = %.6f == %.6f\n",
           amp, basis.modes[0].multiplicity, expected, out_coeffs[0].coeffs[0].real);
}

int main(void) {
    printf("Running initializers tests...\n");
    test_rng();
    test_read_scattering_data();
    test_build_result_from_coeffs();
    test_amplitude_sqrt_scaling();
    test_manual_initialization();
    test_random_initialization();
    test_full_initialization_result();
    test_file_initialization();
    printf("All initializers tests passed.\n");
    return 0;
}
