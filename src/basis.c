/**
 * basis.c — SAFB basis construction
 *
 * Translated from: Sg_init/space_group_plane_family.py
 *
 * Core function: build_basis() — assembles SAFBBasis from symmetry ops,
 * lattice parameters, and desired mode count.
 */

#include "basis.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#define TAU (2.0 * M_PI)

/* ========================================================================
 * Helper: compute |q|^2 for an HKL vector given Ginv
 * ======================================================================== */

static double q2_for_hkl(const int hkl[3], const double Ginv[3][3]) {
    double v[3] = {(double)hkl[0], (double)hkl[1], (double)hkl[2]};
    double rv[3];
    rv[0] = Ginv[0][0]*v[0] + Ginv[0][1]*v[1] + Ginv[0][2]*v[2];
    rv[1] = Ginv[1][0]*v[0] + Ginv[1][1]*v[1] + Ginv[1][2]*v[2];
    rv[2] = Ginv[2][0]*v[0] + Ginv[2][1]*v[1] + Ginv[2][2]*v[2];
    return v[0]*rv[0] + v[1]*rv[1] + v[2]*rv[2];
}

/* ========================================================================
 * Helper: format a family key string like "{110}"
 * ======================================================================== */

static void format_family_key(const int key[3], char *buf, size_t buf_size) {
    snprintf(buf, buf_size, "{%d%d%d}", key[0], key[1], key[2]);
}

/* ========================================================================
 * Helper: check if a star (set of vectors) is closed under inversion
 * ======================================================================== */

static int is_star_closed(const int star[/* N */][3], int N) {
    for (int i = 0; i < N; i++) {
        int neg[3] = {-star[i][0], -star[i][1], -star[i][2]};
        int found = 0;
        for (int j = 0; j < N; j++) {
            if (star[j][0] == neg[0] && star[j][1] == neg[1] && star[j][2] == neg[2]) {
                found = 1;
                break;
            }
        }
        if (!found) return 0;
    }
    return 1;
}

/* ========================================================================
 * Helper: copy star vectors from FamilyPlanesInfo to Star array
 * ======================================================================== */

static void copy_star_to_mode(FamilyPlanesInfo* info, int mode_idx,
                              Star *star, const char *family_key) {
    star->family_key[0] = '\0';
    if (family_key) {
        size_t len = strlen(family_key);
        if (len >= sizeof(star->family_key)) len = sizeof(star->family_key) - 1;
        memcpy(star->family_key, family_key, len);
        star->family_key[len] = '\0';
    }
    star->q2 = 0.0;
    star->star_close = 0;
    star->reason = NULL;
    star->multiplicity = 0;
    star->star_vectors_count = info->star_counts[mode_idx];
    star->rels = NULL; /* opaque placeholder for now */

    for (int i = 0; i < star->star_vectors_count; i++) {
        star->star_vectors[i][0] = info->stars[mode_idx][i][0];
        star->star_vectors[i][1] = info->stars[mode_idx][i][1];
        star->star_vectors[i][2] = info->stars[mode_idx][i][2];
    }
}

/* ========================================================================
 * basis_build — main entry point, translated from space_group_plane_family.py
 * ======================================================================== */

int basis_build(const SymmGroup *ops,
                const char *sg_symbol,
                const LatticeInfo *lattice,
                int N,
                SAFBBasis *result) {
    int dim = lattice->dim;

    /* Step 1: extract unique rotation matrices */
    int unique_R[MAX_OPS][3][3];
    int n_unique = unique_rotations(ops, unique_R);

    /* Step 2: find inversion operations */
    Fraction inv_t[MAX_OPS][3];
    int inv_at_origin[MAX_OPS];
    int inv_count = find_inversion_ops(ops, inv_t, inv_at_origin);

    int centrosymmetric_gp = (inv_count > 0);

    /* Step 3: determine inversion info */
    int has_inversion_at_origin = 0;
    const char *additional_info = NULL;

    if (centrosymmetric_gp) {
        for (int i = 0; i < inv_count; i++) {
            if (inv_at_origin[i]) {
                has_inversion_at_origin = 1;
                break;
            }
        }

        if (has_inversion_at_origin) {
            additional_info = "Inversion at origin present and star is closed → you can form real SABF with coefficients of wavevector G = -G with definite parity (even/odd) directly.";
        } else {
            additional_info = "Group is centrosymmetric but inversion is about r0 = t/2 \xe2\x89\xa0 0; we only know that all the star is closed \xe2\x86\x92 you can form real SABF";
        }
    }

    /* Step 4: check point group for -I */
    int star_close_gp = point_group_has_neg_identity(unique_R, n_unique);
    const char *reason_gp = NULL;
    if (star_close_gp) {
        reason_gp = "Point group contains -I, it is a centrosymmetry group so R*G = -G for all G.";
    }

    /* Step 5: compute metric inverse */
    double Ginv[3][3];
    metric_inverse(lattice->a, lattice->b, lattice->c,
                   lattice->alpha, lattice->beta, lattice->gamma,
                   Ginv);
    fprintf(stderr, "DEBUG basis_build: metric_inverse done\n");

    /* Step 6: generate family planes (heap-allocated to avoid ~114MB stack) */
    fprintf(stderr, "DEBUG basis_build: calling family_planes_info(N=%d)\n", N);
    FamilyPlanesInfo* fam_info = family_planes_info(
        N, Ginv, ops, dim);

    if (!fam_info) return -1;

    /* Step 7: assemble modes */
    result->space_group[0] = '\0';
    if (sg_symbol) {
        strncpy(result->space_group, sg_symbol, sizeof(result->space_group) - 1);
        result->space_group[sizeof(result->space_group) - 1] = '\0';
    }
    result->centrosymmetric_group = centrosymmetric_gp;
    result->has_inversion_at_origin = has_inversion_at_origin;
    result->additional_info = additional_info;
    result->lattice = *lattice;
    result->modes_count = fam_info->count;

    for (int m = 0; m < fam_info->count; m++) {
        Star *mode = &result->modes[m];

        /* Format family key */
        char family_key_buf[32];
        format_family_key(fam_info->family_keys[m], family_key_buf, sizeof(family_key_buf));

        /* Copy star vectors */
        copy_star_to_mode(fam_info, m, mode, family_key_buf);

        /* Compute q2 */
        mode->q2 = q2_for_hkl(fam_info->family_keys[m], Ginv);

        /* Determine star closure for THIS specific star */
        if (star_close_gp) {
            mode->star_close = 1;
            mode->reason = reason_gp;
        } else {
            int closed = is_star_closed(
                fam_info->stars[m],
                fam_info->star_counts[m]);
            mode->star_close = closed;
            if (closed) {
                mode->reason = "For each G in the star, -G is also present.";
            } else {
                mode->reason = "No rotation maps the seed to its negative; -G is not in the star.";
            }
        }

        mode->multiplicity = mode->star_vectors_count;
    }

    family_planes_info_free(fam_info);

    return 0;
}
