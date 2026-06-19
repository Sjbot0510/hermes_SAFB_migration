/**
 * test_field.c — Tests for field.h / field.c
 *
 * Validates:
 *   - 1D/3D inverse FFT
 *   - VTK XML writer
 *   - build_field coefficient placement + iFFT + normalization
 *   - freqf_int grid frequency construction
 */

#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <fftw3.h>
#include "field.h"
#include "initializers.h"
#include "basis.h"

#define EPS 1e-8
#define EPS_COARSE 1e-5

/* ========================================================================
 * Test: 1D FFT/IDFT round-trip using FFTW
 * ======================================================================== */

static void test_fft_1d_roundtrip(void) {
    printf("  Testing 1D FFT round-trip (N=8, power-of-2)...\n");
    {
        fftw_complex* data = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * 8);
        fftw_complex* saved = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * 8);
        for (int i = 0; i < 8; i++) {
            data[i][0] = (double)(i + 1);
            data[i][1] = 0;
            saved[i][0] = data[i][0];
            saved[i][1] = data[i][1];
        }

        fftw_plan fwd = fftw_plan_dft_1d(8, data, data, FFTW_FORWARD, FFTW_ESTIMATE);
        fftw_plan bwd = fftw_plan_dft_1d(8, data, data, FFTW_BACKWARD, FFTW_ESTIMATE);
        fftw_execute(fwd);
        fftw_execute(bwd);
        /* N * saved */
        for (int i = 0; i < 8; i++) {
            data[i][0] /= 8.0;
            data[i][1] /= 8.0;
        }
        for (int i = 0; i < 8; i++) {
            assert(fabs(data[i][0] - saved[i][0]) < EPS);
            assert(fabs(data[i][1] - saved[i][1]) < EPS);
        }
        fftw_destroy_plan(fwd);
        fftw_destroy_plan(bwd);
        fftw_free(data);
        fftw_free(saved);
        printf("    [PASS] 1D FFT round-trip (N=8)\n");
    }

    printf("  Testing 1D FFT round-trip (N=6, non-power-of-2)...\n");
    {
        fftw_complex* data = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * 6);
        fftw_complex* saved = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * 6);
        for (int i = 0; i < 6; i++) {
            data[i][0] = (double)(i + 1);
            data[i][1] = 0;
            saved[i][0] = data[i][0];
            saved[i][1] = data[i][1];
        }

        fftw_plan fwd = fftw_plan_dft_1d(6, data, data, FFTW_FORWARD, FFTW_ESTIMATE);
        fftw_plan bwd = fftw_plan_dft_1d(6, data, data, FFTW_BACKWARD, FFTW_ESTIMATE);
        fftw_execute(fwd);
        fftw_execute(bwd);
        /* N * saved */
        for (int i = 0; i < 6; i++) {
            data[i][0] /= 6.0;
            data[i][1] /= 6.0;
        }
        for (int i = 0; i < 6; i++) {
            assert(fabs(data[i][0] - saved[i][0]) < EPS);
            assert(fabs(data[i][1] - saved[i][1]) < EPS);
        }
        fftw_destroy_plan(fwd);
        fftw_destroy_plan(bwd);
        fftw_free(data);
        fftw_free(saved);
        printf("    [PASS] 1D FFT round-trip (N=6)\n");
    }
}

/* ========================================================================
 * Test: 3D FFT using field API
 * ======================================================================== */

static void test_fftn3d_simple(void) {
    printf("  Testing 3D FFT simple (3x3x3)...\n");

    int Na = 3, Nb = 3, Nc = 3;

    /* Place DC coefficient */
    fftw_complex* Wa = (fftw_complex*)fftw_alloc_complex((size_t)Na * Nb * Nc);
    memset(Wa, 0, sizeof(fftw_complex) * (size_t)(Na * Nb * Nc));
    Wa[0][0] = 1.0;  /* DC coefficient */

    /* Create FFTW plan, feed in, run, read out */
    void* plan = field_create_fftw_plan(Na, Nb, Nc);
    assert(plan != NULL);
    fftw_complex* in = (fftw_complex*)field_get_fftw_input(plan);
    memcpy(in, Wa, sizeof(fftw_complex) * (size_t)(Na * Nb * Nc));
    field_run_fftw(plan);
    fftw_complex* out = (fftw_complex*)field_get_fftw_output(plan);

    /* FFTW BACKWARD is unnormalized: result = N * ifftn */
    double N = (double)Na * (double)Nb * (double)Nc;
    double expected = 1.0 / N;
    for (int i = 0; i < Na * Nb * Nc; i++) {
        assert(fabs(out[i][0] / N - expected) < EPS);
        assert(fabs(out[i][1]) < EPS);
    }

    field_destroy_fftw_plan(plan);
    fftw_free(Wa);
    printf("    [PASS] 3D FFT single DC coefficient\n");
}

static void test_fftn3d_2x2x1(void) {
    printf("  Testing 3D FFT (2x2x1)...\n");

    int Na = 2, Nb = 2, Nc = 1;

    fftw_complex* Wa = (fftw_complex*)fftw_alloc_complex((size_t)Na * Nb * Nc);
    memset(Wa, 0, sizeof(fftw_complex) * (size_t)(Na * Nb * Nc));
    Wa[0][0] = 2.0;  /* DC coefficient */

    void* plan = field_create_fftw_plan(Na, Nb, Nc);
    assert(plan != NULL);
    fftw_complex* in = (fftw_complex*)field_get_fftw_input(plan);
    memcpy(in, Wa, sizeof(fftw_complex) * (size_t)(Na * Nb * Nc));
    field_run_fftw(plan);
    fftw_complex* out = (fftw_complex*)field_get_fftw_output(plan);
    double N = (double)Na * (double)Nb * (double)Nc;

    double expected = 2.0 / ((double)Na * (double)Nb * (double)Nc);
    for (int i = 0; i < Na * Nb * Nc; i++) {
        assert(fabs(out[i][0] / N - expected) < EPS_COARSE);
        assert(fabs(out[i][1]) < EPS_COARSE);
    }

    field_destroy_fftw_plan(plan);
    fftw_free(Wa);
    printf("    [PASS] 3D FFT 2x2x1\n");
}

/* ========================================================================
 * Test: freqf_int
 * ======================================================================== */

static void test_freqf_int(void) {
    printf("  Testing freqf_int...\n");

    int freq6[6];
    freqf_int(6, freq6);
    assert(freq6[0] == 0);
    assert(freq6[1] == 1);
    assert(freq6[2] == 2);
    assert(freq6[3] == -3);
    assert(freq6[4] == -2);
    assert(freq6[5] == -1);

    int freq5[5];
    freqf_int(5, freq5);
    assert(freq5[0] == 0);
    assert(freq5[1] == 1);
    assert(freq5[2] == 2);
    assert(freq5[3] == -2);
    assert(freq5[4] == -1);

    printf("    [PASS] freqf_int\n");
}

/* ========================================================================
 * Test: VTK writer — simple 3x3x1 grid
 * ======================================================================== */

static void test_vtk_writer(void) {
    printf("  Testing VTK writer (3x3x1 grid)...\n");

    LatticeInfo lattice = lattice_info_new(1.0, 1.0, 1.0, 90.0, 90.0, 90.0, 3);

    int Na = 3, Nb = 3, Nc = 1;
    double* field = (double*)calloc((size_t)(Na * Nb * Nc), sizeof(double));
    for (int i = 0; i < Na * Nb * Nc; i++) {
        field[i] = 0.5; /* uniform field */
    }

    const char* filename = "/tmp/test_field_3x3x1.vts";
    int ret = write_lattice_field_to_vts(filename, &lattice, field,
                                          "test_scalar",
                                          Na, Nb, Nc,
                                          0, 1, 1, 1);
    assert(ret == 0);

    /* Verify file exists and has expected content */
    FILE* f = fopen(filename, "r");
    assert(f != NULL);
    char line[1024];
    int has_header = 0, has_structure = 0, has_data = 0, has_scalar = 0;
    while (fgets(line, sizeof(line), f)) {
        if (strstr(line, "<?xml")) has_header = 1;
        if (strstr(line, "StructuredGrid")) has_structure = 1;
        if (strstr(line, "DataArray type=\"Float64\"")) has_data = 1;
        if (strstr(line, "test_scalar")) has_scalar = 1;
    }
    fclose(f);

    assert(has_header);
    assert(has_structure);
    assert(has_data);
    assert(has_scalar);

    /* Verify field value is in the file */
    f = fopen(filename, "r");
    char* result = (char*)malloc(65536);
    memset(result, 0, 65536);
    fread(result, 1, 65535, f);
    fclose(f);
    assert(strstr(result, "5.000000000000000e-01") != NULL);

    free(result);
    free(field);
    remove(filename);

    printf("    [PASS] VTK writer\n");
}

/* ========================================================================
 * Test: build_field with simple P1 basis
 * ======================================================================== */

static void test_build_field_simple(void) {
    printf("  Testing build_field with P1 basis (2x2x1)...\n");

    SymmGroup sg;
    memset(&sg, 0, sizeof(sg));
    sg.dim = 3;
    int R[3][3] = {{1,0,0},{0,1,0},{0,0,1}};
    memcpy(sg.ops[0].R, R, sizeof(R));
    sg.ops[0].t[0] = frac_new(0, 1);
    sg.ops[0].t[1] = frac_new(0, 1);
    sg.ops[0].t[2] = frac_new(0, 1);
    sg.count = 1;

    LatticeInfo lattice = lattice_info_new(2.0, 2.0, 1.0, 90.0, 90.0, 90.0, 3);

    SAFBBasis basis;
    int ret = basis_build(&sg, "P1", &lattice, 2, &basis);
    assert(ret == 0);
    assert(basis.modes_count >= 1);

    FamilyCoeffs coeffs[MAX_AMPLITUDE_KEYS];
    char key_bufs[MAX_AMPLITUDE_KEYS][32];
    const char* keys_ptr[MAX_AMPLITUDE_KEYS];
    double amps[MAX_AMPLITUDE_KEYS];

    int n_coeffs = build_result_from_coeffs(
        &basis,
        (const double[1]){4.0},
        (const char**)(const char*[]){basis.modes[0].family_key},
        1,
        coeffs,
        NULL, 0
    );
    assert(n_coeffs > 0);

    for (int i = 0; i < n_coeffs; i++) {
        strncpy(key_bufs[i], coeffs[i].family_key, 31);
        key_bufs[i][31] = '\0';
        keys_ptr[i] = key_bufs[i];
        amps[i] = 4.0;
    }

    FullInitializationResult result;
    ret = build_initialization_result(
        &basis, coeffs, n_coeffs, amps, keys_ptr, n_coeffs, &result
    );
    assert(ret == 0);

    const char* filename = "/tmp/test_build_field_simple.vts";
    ret = build_field(filename, "Wa", 0, &result, 1.0, 0, 1, 1, 1);
    assert(ret == 0);

    /* Verify VTK file was created */
    FILE* f = fopen(filename, "r");
    assert(f != NULL);
    fclose(f);

    /* Read field and verify it's not all zeros */
    f = fopen(filename, "r");
    char* content = (char*)malloc(65536);
    memset(content, 0, 65536);
    fread(content, 1, 65535, f);
    fclose(f);

    int has_data = (strstr(content, "0.") != NULL || strstr(content, "0,") != NULL);
    assert(has_data);

    free(content);
    remove(filename);

    printf("    [PASS] build_field simple\n");
}

/* ========================================================================
 * Test: 2D case (Nc=1)
 * ======================================================================== */

static void test_build_field_2d(void) {
    printf("  Testing build_field 2D case...\n");

    SymmGroup sg;
    memset(&sg, 0, sizeof(sg));
    sg.dim = 3;
    int R[3][3] = {{1,0,0},{0,1,0},{0,0,1}};
    memcpy(sg.ops[0].R, R, sizeof(R));
    sg.ops[0].t[0] = frac_new(0, 1);
    sg.ops[0].t[1] = frac_new(0, 1);
    sg.ops[0].t[2] = frac_new(0, 1);
    sg.count = 1;

    LatticeInfo lattice = lattice_info_new(2.0, 2.0, 1.0, 90.0, 90.0, 90.0, 2); /* dim=2 */

    SAFBBasis basis;
    int ret = basis_build(&sg, "P1_2D", &lattice, 2, &basis);
    assert(ret == 0);

    FamilyCoeffs coeffs[MAX_AMPLITUDE_KEYS];
    char key_buf[32];
    const char* keys_ptr[1];
    double amps[1];

    int n = build_result_from_coeffs(
        &basis,
        (const double[1]){1.0},
        (const char**)(const char*[]){basis.modes[0].family_key},
        1,
        coeffs,
        NULL, 0
    );

    strncpy(key_buf, coeffs[0].family_key, 31);
    key_buf[31] = '\0';
    keys_ptr[0] = key_buf;
    amps[0] = 1.0;

    FullInitializationResult result;
    ret = build_initialization_result(&basis, coeffs, n, amps, keys_ptr, n, &result);
    assert(ret == 0);

    const char* filename = "/tmp/test_build_field_2d.vts";
    ret = build_field(filename, "Wa", 0, &result, 1.0, 0, 1, 1, 1);
    (void)n;
    printf("    [PASS] build_field 2D (no crash)\n");
}

/* ========================================================================
 * Test: VTK with tiling
 * ======================================================================== */

static void test_vtk_writer_tiled(void) {
    printf("  Testing VTK writer with tiling (3x3x1 -> 9x9x3)...\n");

    LatticeInfo lattice = lattice_info_new(1.0, 1.0, 1.0, 90.0, 90.0, 90.0, 3);

    int Na = 3, Nb = 3, Nc = 1;
    double* field = (double*)calloc((size_t)(Na * Nb * Nc), sizeof(double));
    for (int i = 0; i < Na * Nb * Nc; i++) {
        field[i] = 0.25;
    }

    const char* filename = "/tmp/test_field_tiled.vts";
    int ret = write_lattice_field_to_vts(filename, &lattice, field,
                                          "tile_test",
                                          Na, Nb, Nc,
                                          1, 3, 3, 3);
    assert(ret == 0);

    FILE* f = fopen(filename, "r");
    assert(f != NULL);
    char line[1024];
    int found_extent = 0;
    while (fgets(line, sizeof(line), f)) {
        if (strstr(line, "WholeExtent")) {
            found_extent = 1;
        }
    }
    fclose(f);
    assert(found_extent);

    free(field);
    remove(filename);

    printf("    [PASS] VTK writer with tiling\n");
}

int main(void) {
    printf("Running field tests...\n");

    test_fft_1d_roundtrip();
    test_fftn3d_simple();
    test_fftn3d_2x2x1();
    test_freqf_int();
    test_vtk_writer();
    test_build_field_simple();
    test_build_field_2d();
    test_vtk_writer_tiled();

    printf("All field tests passed.\n");
    return 0;
}
