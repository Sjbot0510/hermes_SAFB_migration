/**
 * field.h — Real-space field generation (iFFT via FFTW) + VTK XML export
 *
 * Translated from: Sg_init/field.py
 *
 * Contains:
 *   - write_lattice_field_to_vts — export scalar field to VTK .vts XML format
 *   - build_field — construct field from coefficients via 3D iFFT
 *   - Grid frequency construction (fftfreq equivalent)
 *   - Coefficient placement on reciprocal grid
 *   - FFTW 3D complex DFT plans
 *
 * Dependencies: FFTW3 (installed via conda-forge/fftw)
 *               Link with: -lfftw3
 */

#ifndef FIELD_H
#define FIELD_H

#include <stdint.h>
#include <stddef.h>
#include "domain.h"
#include "initializers.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ========================================================================
 * FFTW — 3D complex DFT
 *
 * Uses FFTW's plan-based API for optimal performance on any grid size.
 * FFTW handles power-of-2, prime, and composite sizes efficiently.
 * ======================================================================== */

/* Create an inverse FFTW plan for a 3D complex array [Na][Nb][Nc] */
void* field_create_fftw_plan(int Na, int Nb, int Nc);

/* Execute the inverse FFT plan */
void field_run_fftw(void* plan);

/* Destroy the plan */
void field_destroy_fftw_plan(void* plan);

/* Allocate a contiguous complex double array for FFTW */
double* field_alloc_fftw_complex(int Na, int Nb, int Nc);

/* Get the input (complex) buffer from a plan for pre-fill */
void* field_get_fftw_input(void* plan);

/* Get the output (real) buffer from a plan after execution */
void* field_get_fftw_output(void* plan);

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
