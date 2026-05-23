/**
 * field.c — Real-space field generation (iFFT) + VTK XML export
 *
 * Translated from: Sg_init/field.py
 *
 * Contains:
 *   - 3D inverse FFT (Cooley-Tukey radix-2 + DFT fallback)
 *   - Frequency grid construction (fftfreq equivalent)
 *   - VTK XML Structured Grid writer (.vts)
 *   - build_field — coefficient placement + iFFT + tanh normalization
 */

#define _POSIX_C_SOURCE 200809L
#include "field.h"
#include "initializers.h"
#include "basis.h"
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

void fftw3d_free(ComplexDouble* arr) {
    free(arr);
}

/* ========================================================================
 * Utility — check if N is a power of 2
 */

static int is_power_of_2(int n) {
    return n > 0 && (n & (n - 1)) == 0;
}

/* ========================================================================
 * 1D FFT — radix-2 Cooley-Tukey (for power-of-2 sizes)
 */

static void fft_radix2_recursive(ComplexDouble* data, int N, int direction) {
    if (N <= 1) return;

    if (N == 2) {
        double r0 = data[0].real, i0 = data[0].imag;
        double r1 = data[1].real, i1 = data[1].imag;
        data[0].real = r0 + r1;  data[0].imag = i0 + i1;
        data[1].real = r0 - r1;  data[1].imag = i0 - i1;
        return;
    }

    int half = N / 2;

    /* Split even/odd */
    ComplexDouble* even = (ComplexDouble*)malloc(sizeof(ComplexDouble) * (size_t)half);
    ComplexDouble* odd  = (ComplexDouble*)malloc(sizeof(ComplexDouble) * (size_t)half);
    if (!even || !odd) { free(even); free(odd); return; }

    for (int i = 0; i < half; i++) {
        even[i] = data[2 * i];
        odd[i]  = data[2 * i + 1];
    }

    fft_radix2_recursive(even, half, direction);
    fft_radix2_recursive(odd, half, direction);

    /* Butterfly */
    for (int k = 0; k < half; k++) {
        double angle = -direction * M_TAU * (double)k / (double)N;
        double cos_a = cos(angle);
        double sin_a = sin(angle);
        double t_r = cos_a * odd[k].real - sin_a * odd[k].imag;
        double t_i = cos_a * odd[k].imag + sin_a * odd[k].real;
        data[k].real     = even[k].real + t_r;
        data[k].imag     = even[k].imag + t_i;
        data[k + half].real = even[k].real - t_r;
        data[k + half].imag = even[k].imag - t_i;
    }

    free(even);
    free(odd);
}

/* ========================================================================
 * 1D DFT — O(N^2) fallback for non-power-of-2
 * ======================================================================== */

void fft_1d(ComplexDouble* data, int N, int direction) {
    if (N <= 0) return;

    if (is_power_of_2(N)) {
        fft_radix2_recursive(data, N, direction);
        return;
    }

    /* Fallback to DFT */
    ComplexDouble* temp = (ComplexDouble*)malloc(sizeof(ComplexDouble) * (size_t)N);
    if (!temp) return;
    memcpy(temp, data, sizeof(ComplexDouble) * (size_t)N);

    for (int k = 0; k < N; k++) {
        temp[k].real = 0;
        temp[k].imag = 0;
        for (int n = 0; n < N; n++) {
            double angle = -direction * M_TAU * (double)k * (double)n / (double)N;
            double cos_a = cos(angle);
            double sin_a = sin(angle);
            temp[k].real += cos_a * data[n].real - sin_a * data[n].imag;
            temp[k].imag += cos_a * data[n].imag + sin_a * data[n].real;
        }
    }
    memcpy(data, temp, sizeof(ComplexDouble) * (size_t)N);
    free(temp);
}

/* ========================================================================
 * 3D inverse FFT — apply 1D IFFT along each axis
 * ======================================================================== */

void fftn3d_idft(ComplexDouble* arr, int nx, int ny, int nz) {
    int total = nx * ny * nz;
    (void)total;

    /* We apply 1D FFT along the k-axis (fastest), then j, then i */
    /* Step 1: FFT along k-axis (index k varies fastest, i,j fixed) */
    {
        ComplexDouble* buf = (ComplexDouble*)malloc(sizeof(ComplexDouble) * (size_t)nz);
        if (!buf) return;
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                for (int k = 0; k < nz; k++) {
                    buf[k] = arr[i * ny * nz + j * nz + k];
                }
                fft_1d(buf, nz, -1); /* -1 for inverse (standard FFT has +1) */
                for (int k = 0; k < nz; k++) {
                    arr[i * ny * nz + j * nz + k] = buf[k];
                }
            }
        }
        free(buf);
    }

    /* Step 2: FFT along j-axis (i, k fixed) */
    {
        ComplexDouble* buf = (ComplexDouble*)malloc(sizeof(ComplexDouble) * (size_t)ny);
        if (!buf) return;
        for (int i = 0; i < nx; i++) {
            for (int k = 0; k < nz; k++) {
                for (int j = 0; j < ny; j++) {
                    buf[j] = arr[i * ny * nz + j * nz + k];
                }
                fft_1d(buf, ny, -1);
                for (int j = 0; j < ny; j++) {
                    arr[i * ny * nz + j * nz + k] = buf[j];
                }
            }
        }
        free(buf);
    }

    /* Step 3: FFT along i-axis (j, k fixed) */
    {
        ComplexDouble* buf = (ComplexDouble*)malloc(sizeof(ComplexDouble) * (size_t)nx);
        if (!buf) return;
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                for (int i = 0; i < nx; i++) {
                    buf[i] = arr[i * ny * nz + j * nz + k];
                }
                fft_1d(buf, nx, -1);
                for (int i = 0; i < nx; i++) {
                    arr[i * ny * nz + j * nz + k] = buf[i];
                }
            }
        }
        free(buf);
    }

    /* Scale by 1/(Nx * Ny * Nz) — this is the normalization used by scipy.fft.ifftn */
    {
        double scale = 1.0 / ((double)nx * (double)ny * (double)nz);
        for (int i = 0; i < total; i++) {
            arr[i].real *= scale;
            arr[i].imag *= scale;
        }
    }
}

/* ========================================================================
 * VTK XML Structured Grid Writer
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
        /* Try to create directory (ignore errors if it exists) */
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
 * Helper — map integer frequency index to FFT frequency (fftfreq equivalent)
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
 * build_field — main entry point
 *
 * Translated from: Sg_init/field.py
 *
 * 1. Compute grid dimensions (Na, Nb, Nc) from lattice + resolution
 * 2. Build frequency arrays (H, K, L) using fftfreq
 * 3. Place coefficients from InitializationResult onto reciprocal grid
 * 4. Compute iFFT via FFTW3D
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

    /* Resolve the lattice_info first — avoid dangling pointer to local tmp */
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

    /* Allocate reciprocal grid */
    ComplexDouble* Wa = fftw3d_alloc(Na, Nb, Nc);
    if (!Wa) {
        free(H); free(K); free(L);
        fprintf(stderr, "field: allocation failed for Wa\n");
        return -1;
    }
    /* Initialize to zero */
    memset(Wa, 0, sizeof(ComplexDouble) * (size_t)(Na * Nb * Nc));

    /* Build a lookup: family_key -> array of (hkl, real, imag) from coeffs */
    /* coeffs[] is already in FamilyCoeffs format */
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
                                    ComplexDouble* target = fftw3d_at(Wa, Na, Nb, Nc, ii, jj, kk2);
                                    target->real = cc->real;
                                    target->imag = cc->imag;
                                    /* Only match first L index (l is unique) */
                                    kk2 = Nc;
                                    break;
                                }
                            }
                            jj = Nb; /* break inner */
                        }
                    }
                    ii = Na; /* break outer */
                }
            }
        }
    }

    /* Compute iFFT */
    fftn3d_idft(Wa, Na, Nb, Nc);

    /* Compute real field — multiply by Na*Nb*Nc per Python's ifftn scaling */
    /* Python: W_real = np.real(ifftn(Wa)) * Na * Nb * Nc */
    /* But our FFT already divides by N, so we multiply back */
    double scale = (double)Na * (double)Nb * (double)Nc;

    /* Allocate real field */
    double* W_real = (double*)malloc((size_t)(Na * Nb * Nc) * sizeof(double));
    if (!W_real) {
        fftw3d_free(Wa);
        free(H); free(K); free(L);
        fprintf(stderr, "field: allocation failed for W_real\n");
        return -1;
    }
    for (int i = 0; i < Na * Nb * Nc; i++) {
        W_real[i] = Wa[i].real * scale;
    }
    fftw3d_free(Wa);

    /* tanh normalization */
    double w_mean = 0;
    for (int i = 0; i < Na * Nb * Nc; i++) {
        w_mean += W_real[i];
    }
    w_mean /= (double)(Na * Nb * Nc);

    double w_std_sq = 0;
    for (int i = 0; i < Na * Nb * Nc; i++) {
        double d = W_real[i] - w_mean;
        w_std_sq += d * d;
    }
    w_std_sq /= (double)(Na * Nb * Nc);
    double w_std = sqrt(w_std_sq);

    double* W_normalized = (double*)malloc((size_t)(Na * Nb * Nc) * sizeof(double));
    if (!W_normalized) {
        free(W_real);
        free(H); free(K); free(L);
        fprintf(stderr, "field: allocation failed for W_normalized\n");
        return -1;
    }

    if (w_std < 1e-9) {
        for (int i = 0; i < Na * Nb * Nc; i++) {
            W_normalized[i] = w_mean;
        }
    } else {
        for (int i = 0; i < Na * Nb * Nc; i++) {
            double x = 0.5 * (1.0 + tanh((W_real[i] - w_mean) / w_std));
            W_normalized[i] = x;
        }
    }

    free(W_real);
    free(H);
    free(K);
    free(L);

    /* Write VTK output */
    int ret = write_lattice_field_to_vts(
        filename,
        lattice_info,
        W_normalized,
        field_name,
        Na, Nb, Nc,
        apply_tile,
        tile_x, tile_y, tile_z
    );

    free(W_normalized);
    return ret;
}
