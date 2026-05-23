/**
 * field.h — Real-space field generation (iFFT) + VTK XML export
 *
 * Translated from: Sg_init/field.py
 *
 * Contains:
 *   - write_lattice_field_to_vts — export scalar field to VTK .vts XML format
 *   - build_field — construct field from coefficients via 3D iFFT
 *   - 3D inverse FFT (Cooley-Tukey, radix-2 or DFT fallback)
 *   - Grid frequency construction (fftfreq equivalent)
 *   - Coefficient placement on reciprocal grid
 */

#ifndef FIELD_H
#define FIELD_H

#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>
#include "domain.h"
#include "initializers.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ========================================================================
 * FFTW-free 3D inverse FFT
 *
 * For real-space grids where N is odd/even, we use a simple Cooley-Tukey
 * approach with radix-2 DFT fallback for non-power-of-2 sizes.
 *
 * Data is stored as complex double (real, imag) pairs in row-major order.
 * ======================================================================== */

/* Complex double array */
typedef struct {
    double real;
    double imag;
} ComplexDouble;

/* Allocate a 3D complex array of shape [Nx][Ny][Nz] */
static inline ComplexDouble* fftw3d_alloc(int nx, int ny, int nz) {
    return (ComplexDouble*)calloc((size_t)(nx * ny * nz), sizeof(ComplexDouble));
}
void fftw3d_free(ComplexDouble* arr);

/* Access element at (i, j, k) */
static inline ComplexDouble* fftw3d_at(ComplexDouble* arr, int nx, int ny, int nz,
                                        int i, int j, int k) {
    return &arr[i * ny * nz + j * nz + k];
}

/* 1D DFT — O(N^2) fallback for non-power-of-2 sizes */
void dft_1d(ComplexDouble* data, int N, int direction); /* direction: +1=DFT, -1=IDFT */

/* 1D FFT — O(N log N) for power-of-2; falls back to DFT */
void fft_1d(ComplexDouble* data, int N, int direction);

/* 3D inverse FFT along the k-axis (fastest-varying index) */
void fftn3d_idft(ComplexDouble* arr, int nx, int ny, int nz);

/* ========================================================================
 * Frequency grid construction — scipy.fft.fftfreq equivalent
 * ======================================================================== */

/* Build frequency array equivalent to numpy fftfreq(N, 1.0/N) */
void freqf_int(int N, int* freq);

/* ========================================================================
 * VTK XML Structured Grid Writer (.vts format)
 *
 * Outputs an ASCII VTK XML file compatible with ParaView.
 * ======================================================================== */

/*
 * Write a 3D scalar field to VTK XML Structured Grid format.
 *
 * Parameters:
 *   filename       — output file path
 *   lattice_info   — unit cell parameters
 *   field            — scalar field data (real double[3D], Nc x Nb x Na)
 *   field_name      — name of the scalar array
 *   Na, Nb, Nc     — grid dimensions (fastest to slowest: i, j, k)
 *   apply_tile     — whether to tile the field before writing
 *   tile_x, tile_y, tile_z — tiling factors (used only if apply_tile)
 *
 * Returns: 0 on success, -1 on error.
 */
int write_lattice_field_to_vts(const char* filename,
                                const LatticeInfo* lattice_info,
                                const double* field,
                                const char* field_name,
                                int Na, int Nb, int Nc,
                                int apply_tile,
                                int tile_x, int tile_y, int tile_z);

/* ========================================================================
 * Build field — main entry point, translated from field.py
 *
 * Parameters:
 *   filename       — output VTK file
 *   field_name     — name for the scalar field
 *   apply_tile     — whether to tile the output field
 *   result         — InitializationResult (contains coefficients)
 *   resol          — grid resolution (lattice_dim / resol = N points)
 *   transform_coord — if true, use trans_coeff and trans_lattice_params
 *   tile_x, tile_y, tile_z — tiling factors
 *
 * Returns: 0 on success, -1 on error.
 */
int build_field(const char* filename,
                const char* field_name,
                int apply_tile,
                FullInitializationResult* result,
                double resol,
                int transform_coord,
                int tile_x, int tile_y, int tile_z);

#ifdef __cplusplus
}
#endif

#endif /* FIELD_H */
