/**
 * field.c — Real-space field generation (iFFT via FFTW) + VTK XML export
 *
 * Translated from: Sg_init/field.py
 *
 * Contains:
 *   - FFTW 3D inverse FFT plans (full complex DFT)
 *   - write_lattice_field_to_vts — export scalar field to VTK .vts XML format
 *   - build_field — construct field from coefficients via 3D iFFT
 *   - Grid frequency construction (fftfreq equivalent)
 *   - Coefficient placement on reciprocal grid
 *
 * Dependencies: FFTW3 (installed via conda-forge/fftw)
 *               Link with: -lfftw3
 */

#define _POSIX_C_SOURCE 200809L
#include "field.h"
#include "initializers.h"
#include "basis.h"
#include <fftw3.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <float.h>
#include <ctype.h>
#include <sys/stat.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_TAU
#define M_TAU (2.0 * M_PI)
#endif

/* ========================================================================
 * FFTW context — stores plan and input/output arrays
 * ======================================================================== */

typedef struct {
    fftw_plan plan;
    fftw_complex* in;   /* complex input (frequency domain) */
    fftw_complex* out;  /* complex output (real space, take real part) */
    int Na, Nb, Nc;
} fftw_context_t;

/* ========================================================================
 * FFTW — 3D complex DFT
 *
 * Uses full complex-to-complex plans (fftw_plan_dft_3d) rather than
 * c2r because we work with general complex coefficients.
 * ======================================================================== */

void* field_create_fftw_plan(int Na, int Nb, int Nc) {
    size_t total = (size_t)Na * Nb * Nc;

    fftw_complex* in = (fftw_complex*)fftw_alloc_complex(total);
    fftw_complex* out = (fftw_complex*)fftw_alloc_complex(total);
    if (!in || !out) {
        fftw_free(in);
        fftw_free(out);
        return NULL;
    }

    /* Create plan: complex-to-complex inverse DFT */
    fftw_plan plan = fftw_plan_dft_3d(Na, Nb, Nc,
                                       in, out,
                                       FFTW_BACKWARD,
                                       FFTW_ESTIMATE);
    if (!plan) {
        fftw_free(in);
        fftw_free(out);
        return NULL;
    }

    fftw_context_t* ctx = (fftw_context_t*)malloc(sizeof(fftw_context_t));
    if (!ctx) {
        fftw_destroy_plan(plan);
        fftw_free(in);
        fftw_free(out);
        return NULL;
    }
    ctx->plan = plan;
    ctx->in = in;
    ctx->out = out;
    ctx->Na = Na;
    ctx->Nb = Nb;
    ctx->Nc = Nc;

    return (void*)ctx;
}

void field_run_fftw(void* plan) {
    if (!plan) return;
    fftw_context_t* ctx = (fftw_context_t*)plan;
    fftw_execute_dft(ctx->plan, ctx->in, ctx->out);
}

void field_destroy_fftw_plan(void* plan) {
    if (!plan) return;
    fftw_context_t* ctx = (fftw_context_t*)plan;
    fftw_destroy_plan(ctx->plan);
    fftw_free(ctx->in);
    fftw_free(ctx->out);
    free(ctx);
}

void* field_get_fftw_input(void* plan) {
    if (!plan) return NULL;
    return (void*)((fftw_context_t*)plan)->in;
}

void* field_get_fftw_output(void* plan) {
    if (!plan) return NULL;
    return (void*)((fftw_context_t*)plan)->out;
}

double* field_alloc_fftw_complex(int Na, int Nb, int Nc) {
    return (double*)fftw_alloc_complex((size_t)Na * Nb * Nc);
}

/* ========================================================================
 * Frequency grid construction — scipy.fft.fftfreq equivalent
 * ======================================================================== */

void freqf_int(int N, int* freq) {
    int threshold = N - N / 2;  /* ceil(N/2) */
    for (int i = 0; i < N; i++) {
        if (i < threshold) {
            freq[i] = i;
        } else {
            freq[i] = i - N;
        }
    }
}

/* ========================================================================
 * VTK XML Structured Grid Writer (.vts format)
 * ======================================================================== */

int write_lattice_field_to_vts(const char* filename,
                                const LatticeInfo* lattice_info,
                                const double* field,
                                const char* field_name,
                                int Na, int Nb, int Nc,
                                int apply_tile,
                                int tile_x, int tile_y, int tile_z) {
    int out_Na = apply_tile ? Na * tile_x : Na;
    int out_Nb = apply_tile ? Nb * tile_y : Nb;
    int out_Nc = apply_tile ? Nc * tile_z : Nc;

    /* Compute direct basis matrix */
    double A[3][3];
    direct_basis_from_lattice_info(lattice_info, A);

    /* Create output directory if needed */
    char dir_buf[1024];
    strncpy(dir_buf, filename, sizeof(dir_buf) - 1);
    dir_buf[sizeof(dir_buf) - 1] = '\0';
    char* last_slash = strrchr(dir_buf, '/');
    if (last_slash) {
        *last_slash = '\0';
        struct stat st;
        if (stat(dir_buf, &st) != 0) {
            /* mkdir not always portable, skip */
        }
    }

    FILE* f = fopen(filename, "w");
    if (!f) {
        fprintf(stderr, "field: cannot create %s\n", filename);
        return -1;
    }

    /* VTK XML Structured Grid header */
    fprintf(f, "<?xml version=\"1.0\"?>\n");
    fprintf(f, "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"\">\n");
    fprintf(f, "  <StructuredGrid WholeExtent=\"0 %d 0 %d 0 %d\">\n",
            out_Na - 1, out_Nb - 1, out_Nc - 1);

    /* Points — fractional coords mapped through direct basis */
    fprintf(f, "  <Points>\n");
    fprintf(f, "    <DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"0\">\n");

    int total_points = out_Na * out_Nb * out_Nc;
    fprintf(f, "      %d\n", total_points * 3);

    /* Write points: for each (i,j,k) compute lattice coords = A @ fractional */
    for (int k = 0; k < out_Nc; k++) {
        double w = (double)k / (double)out_Nc;
        for (int j = 0; j < out_Nb; j++) {
            double v = (double)j / (double)out_Nb;
            for (int i = 0; i < out_Na; i++) {
                double u = (double)i / (double)out_Na;
                double x = A[0][0] * u + A[0][1] * v + A[0][2] * w;
                double y = A[1][0] * u + A[1][1] * v + A[1][2] * w;
                double z = A[2][0] * u + A[2][1] * v + A[2][2] * w;
                fprintf(f, "%22.15e %22.15e %22.15e\n", x, y, z);
            }
        }
    }
    fprintf(f, "    </DataArray>\n");
    fprintf(f, "  </Points>\n");

    /* Scalar data — tile the field if requested */
    fprintf(f, "  <PointData Scalars=\"%s\">\n", field_name);
    fprintf(f, "    <DataArray type=\"Float64\" NumberOfComponents=\"1\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"0\">\n");
    fprintf(f, "      %d\n", total_points);

    for (int k = 0; k < out_Nc; k++) {
        int src_k = k % Nc;
        for (int j = 0; j < out_Nb; j++) {
            int src_j = j % Nb;
            for (int i = 0; i < out_Na; i++) {
                int src_i = i % Na;
                double val = field[src_k * Nb * Na + src_j * Na + src_i];
                fprintf(f, "      %22.15e\n", val);
            }
        }
    }
    fprintf(f, "    </DataArray>\n");
    fprintf(f, "  </PointData>\n");
    fprintf(f, "</StructuredGrid>\n");
    fprintf(f, "</VTKFile>\n");

    fclose(f);
    return 0;
}

/* ========================================================================
 * build_field — main entry point
 *
 * Translated from: Sg_init/field.py
 *
 * 1. Compute grid dimensions (Na, Nb, Nc) from lattice + resolution
 * 2. Build frequency arrays (H, K, L) using fftfreq
 * 3. Place coefficients from InitializationResult onto reciprocal grid
 * 4. Compute iFFT via FFTW
 * 5. Apply tanh normalization
 * 6. Write VTK XML output
 * ======================================================================== */

int build_field(const char* filename,
                const char* field_name,
                int apply_tile,
                FullInitializationResult* result,
                double resol,
                int transform_coord,
                int tile_x, int tile_y, int tile_z) {

    /* Resolve the lattice_info first */
    LatticeInfo lattice_buf;
    const LatticeInfo* lattice_info;

    if (transform_coord) {
        memset(&lattice_buf, 0, sizeof(lattice_buf));
        lattice_buf.a = result->trans_lattice_a;
        lattice_buf.b = result->trans_lattice_b;
        lattice_buf.c = result->trans_lattice_c;
        lattice_buf.alpha = result->trans_lattice_alpha;
        lattice_buf.beta = result->trans_lattice_beta;
        lattice_buf.gamma = result->trans_lattice_gamma;
        lattice_buf.dim = result->lattice.dim;
        lattice_info = &lattice_buf;
    } else {
        lattice_info = &result->lattice;
    }

    int dim = lattice_info->dim;

    /* Compute grid dimensions */
    int Na = (int)round(lattice_info->a / resol);
    Na += (Na % 2); /* ensure odd */

    int Nb = (int)round(lattice_info->b / resol);
    Nb += (Nb % 2);

    int Nc;
    if (dim == 2) {
        Nc = 1;
    } else {
        Nc = (int)round(lattice_info->c / resol);
        Nc += (Nc % 2);
    }

    /* Store discretization */
    if (!transform_coord) {
        result->discretization_Na = Na;
        result->discretization_Nb = Nb;
        result->discretization_Nc = Nc;
    }

    /* Build frequency arrays */
    int* H = (int*)malloc((size_t)Na * sizeof(int));
    int* K = (int*)malloc((size_t)Nb * sizeof(int));
    int* L = (int*)malloc((size_t)Nc * sizeof(int));
    if (!H || !K || !L) {
        free(H); free(K); free(L);
        fprintf(stderr, "field: allocation failed\n");
        return -1;
    }
    freqf_int(Na, H);
    freqf_int(Nb, K);
    if (Nc > 1) {
        freqf_int(Nc, L);
    } else {
        L[0] = 0;
    }

    /* Allocate and initialize reciprocal grid */
    fftw_complex* Wa = (fftw_complex*)fftw_alloc_complex((size_t)Na * Nb * Nc);
    if (!Wa) {
        free(H); free(K); free(L);
        fprintf(stderr, "field: allocation failed for Wa\n");
        return -1;
    }
    memset(Wa, 0, sizeof(fftw_complex) * (size_t)(Na * Nb * Nc));

    /* Place coefficients from InitializationResult onto reciprocal grid */
    for (int ci = 0; ci < result->n_coeffs && ci < MAX_AMPLITUDE_KEYS; ci++) {
        const FamilyCoeffs* fc = &result->coeffs[ci];
        for (int vi = 0; vi < fc->count && vi < MAX_VEC_PER_STAR; vi++) {
            const ComplexCoeff* cc = &fc->coeffs[vi];
            int hh = cc->hkl[0];
            int kk = cc->hkl[1];
            int ll = cc->hkl[2];

            /* Find matching position in reciprocal grid */
            for (int ii = 0; ii < Na; ii++) {
                if (H[ii] == hh) {
                    for (int jj = 0; jj < Nb; jj++) {
                        if (K[jj] == kk) {
                            for (int kk2 = 0; kk2 < Nc; kk2++) {
                                if (L[kk2] == ll) {
                                    Wa[ii * Nb * Nc + jj * Nc + kk2][0] = cc->real;
                                    Wa[ii * Nb * Nc + jj * Nc + kk2][1] = cc->imag;
                                    kk2 = Nc;  /* break inner */
                                    break;
                                }
                            }
                            jj = Nb;  /* break inner */
                        }
                    }
                    ii = Na;  /* break outer */
                }
            }
        }
    }

    /* Execute FFTW inverse plan */
    void* plan = field_create_fftw_plan(Na, Nb, Nc);
    if (!plan) {
        fftw_free(Wa);
        free(H); free(K); free(L);
        fprintf(stderr, "field: FFTW plan creation failed\n");
        return -1;
    }

    /* Copy input to plan's input buffer */
    fftw_complex* plan_in = (fftw_complex*)field_get_fftw_input(plan);
    memcpy(plan_in, Wa, sizeof(fftw_complex) * (size_t)(Na * Nb * Nc));

    /* Execute */
    field_run_fftw(plan);

    /* Get output (complex) and take real part */
    fftw_complex* plan_out = (fftw_complex*)field_get_fftw_output(plan);
    double* W_real = (double*)malloc((size_t)(Na * Nb * Nc) * sizeof(double));
    if (!W_real) {
        field_destroy_fftw_plan(plan);
        fftw_free(Wa);
        free(H); free(K); free(L);
        fprintf(stderr, "field: allocation failed for W_real\n");
        return -1;
    }

    /* FFTW BACKWARD gives unnormalized result: output = N * true_ifft(input) */
    /* So we need to divide by N = Na * Nb * Nc */
    double N = (double)Na * (double)Nb * (double)Nc;
    for (int i = 0; i < Na * Nb * Nc; i++) {
        W_real[i] = plan_out[i][0] / N;  /* real part / N */
    }

    /* Apply tanh normalization (same as Python) */
    double w_mean = 0.0;
    double w_sum_sq = 0.0;
    for (int i = 0; i < Na * Nb * Nc; i++) {
        w_mean += W_real[i];
        w_sum_sq += W_real[i] * W_real[i];
    }
    w_mean /= (Na * Nb * Nc);
    double w_std = sqrt(w_sum_sq / (Na * Nb * Nc) - w_mean * w_mean);

    double* W_normalized = (double*)malloc((size_t)(Na * Nb * Nc) * sizeof(double));
    if (!W_normalized) {
        field_destroy_fftw_plan(plan);
        fftw_free(Wa);
        free(H); free(K); free(L);
        free(W_real);
        fprintf(stderr, "field: allocation failed for W_normalized\n");
        return -1;
    }

    if (w_std < 1e-9) {
        for (int i = 0; i < Na * Nb * Nc; i++) {
            W_normalized[i] = w_mean;
        }
    } else {
        for (int i = 0; i < Na * Nb * Nc; i++) {
            W_normalized[i] = 0.5 * (1.0 + tanh((W_real[i] - w_mean) / w_std));
        }
    }

    /* Write VTK output */
    int ret = write_lattice_field_to_vts(filename, lattice_info, W_normalized,
                                         field_name, Na, Nb, Nc,
                                         apply_tile, tile_x, tile_y, tile_z);

    /* Cleanup */
    field_destroy_fftw_plan(plan);
    fftw_free(Wa);
    free(H); free(K); free(L);
    free(W_real);
    free(W_normalized);

    return ret;
}
