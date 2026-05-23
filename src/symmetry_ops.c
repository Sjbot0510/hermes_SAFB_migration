/**
 * symmetry_ops.c — Symmetry operation utilities for SAFB
 *
 * Translated from: Sg_init/symmetry_ops.py
 *
 * Covers:
 *   - Fraction arithmetic for rational translation components
 *   - read_spacegroup_ops_txt — parse PSCF-style symmetry operation files
 *   - Fraction math: frac_mod1, is_zero_mod1_vec, frac_to_float
 *   - Lattice basis: direct_basis_from_lattice_info, lattice_params_from_basis
 *   - Star generation: star_from_hkl, get_family_key_lexicographical
 *   - Star analysis: star_is_closed, point_group_has_neg_identity, find_inversion_ops
 *   - Relationships: relationships_in_star, solve_star_coeffs
 *   - Family planes: family_planes_info
 *   - Misc: read_final_values, float_to_miller_int
 */

#define _POSIX_C_SOURCE 200809L
#include "symmetry_ops.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <float.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define TAU (2.0 * M_PI)

/* ========================================================================
 * Fraction arithmetic
 * ======================================================================== */

int frac_gcd(int a, int b) {
    a = (a < 0) ? -a : a;
    b = (b < 0) ? -b : b;
    while (b) {
        int t = b;
        b = a % b;
        a = t;
    }
    return a;
}

static void frac_reduce(Fraction *f) {
    if (f->den < 0) {
        f->num = -f->num;
        f->den = -f->den;
    }
    if (f->num == 0) {
        f->den = 1;
        return;
    }
    int g = frac_gcd(f->num, f->den);
    f->num /= g;
    f->den /= g;
}

Fraction frac_new(int num, int den) {
    Fraction f = {num, den};
    frac_reduce(&f);
    return f;
}

double frac_to_float(Fraction q) {
    return (double)q.num / (double)q.den;
}

Fraction frac_mod1(Fraction q) {
    /* reduce a rational mod 1 into [0,1)
     * e.g. 5/4 mod 1 = 1/4; -1/3 mod 1 = 2/3
     * x mod 1 = x - floor(x) */
    if (q.den <= 0) return frac_new(0, 1);
    int n = q.num % q.den;
    if (n < 0) n += q.den;
    return frac_new(n, q.den);
}

int is_zero_mod1_vec(const Fraction t[3]) {
    for (int i = 0; i < 3; i++) {
        Fraction m = frac_mod1(t[i]);
        if (m.num != 0) return 0;
    }
    return 1;
}

int equal_int_mat(const int A[3][3], const int B[3][3]) {
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            if (A[i][j] != B[i][j]) return 0;
    return 1;
}

/* ========================================================================
 * Helper: parse a token that may be a fraction like "1/4" or "-1/4"
 * ======================================================================== */
static Fraction parse_frac_token(const char *s) {
    const char *slash = strchr(s, '/');
    if (slash) {
        int a, b;
        /* parse numerator (may be negative) and denominator */
        char num_str[64], den_str[64];
        int num_len = (int)(slash - s);
        if (num_len <= 0) num_len = 1;
        if (num_len >= 64) num_len = 63;
        strncpy(num_str, s, num_len);
        num_str[num_len] = '\0';
        strncpy(den_str, slash + 1, 63);
        den_str[63] = '\0';
        a = atoi(num_str);
        b = atoi(den_str);
        if (b == 0) b = 1;
        return frac_new(a, b);
    }
    return frac_new(atoi(s), 1);
}

/* ========================================================================
 * read_spacegroup_ops_txt
 * ======================================================================== */

SymmGroup read_spacegroup_ops_txt(const char *path) {
    SymmGroup sg;
    memset(&sg, 0, sizeof(sg));
    sg.dim = 3;

    FILE *f = fopen(path, "r");
    if (!f) {
        fprintf(stderr, "symmetry_ops: cannot open %s\n", path);
        return sg;
    }

    /* Read all lines into buffer */
    char *all_lines[8192];
    int num_lines = 0;
    char buf[4096];
    while (num_lines < 8192 && fgets(buf, sizeof(buf), f)) {
        /* strip trailing whitespace */
        size_t len = strlen(buf);
        while (len > 0 && (buf[len-1] == '\n' || buf[len-1] == '\r' || buf[len-1] == ' ' || buf[len-1] == '\t'))
            buf[--len] = '\0';
        if (len == 0) continue;
        all_lines[num_lines] = strdup(buf);
        num_lines++;
    }
    fclose(f);

    /* Detect dimension */
    for (int i = 0; i < num_lines; i++) {
        if (strncmp(all_lines[i], "dim", 3) == 0) {
            int d = 3;
            if (sscanf(all_lines[i], "dim %d", &d) == 1) {
                if (d == 2) sg.dim = 2;
            }
            break;
        }
    }

    /* Tokenize all lines into data_rows */
    char *tokens[8192][6];
    int token_counts[8192];
    (void)token_counts;
    int num_data_rows = 0;

    for (int i = 0; i < num_lines; i++) {
        const char *ln = all_lines[i];
        if (strncmp(ln, "dim", 3) == 0) continue;
        if (strncmp(ln, "size", 4) == 0) continue;
        if (ln[0] == '[') continue;

        char tmp[256];
        strncpy(tmp, ln, 255);
        tmp[255] = '\0';
        int ncols = 0;
        char *tok = strtok(tmp, " \t");
        while (tok && ncols < 6) {
            tokens[num_data_rows][ncols] = tok;
            ncols++;
            tok = strtok(NULL, " \t");
        }
        if (ncols >= 2) num_data_rows++;
    }

    int nR = (sg.dim == 3) ? 3 : 2;
    int stride = nR + 1; /* nR rows for R, 1 row for t */

    for (int i = 0; i + stride - 1 < num_data_rows; i += stride) {
        if (sg.count >= MAX_OPS) break;

        /* Parse rotation matrix */
        int R[3][3];
        memset(R, 0, sizeof(R));
        for (int r = 0; r < nR; r++) {
            for (int c = 0; c < sg.dim; c++) {
                R[r][c] = atoi(tokens[i + r][c]);
            }
        }
        /* Fill in the 3rd row/col for 2D embedding */
        if (sg.dim == 2) {
            R[0][2] = 0;
            R[1][2] = 0;
            R[2][0] = 0;
            R[2][1] = 0;
            R[2][2] = 1;
        }

        /* Parse translation vector */
        Fraction t[3];
        for (int c = 0; c < sg.dim; c++) {
            t[c] = parse_frac_token(tokens[i + nR][c]);
        }
        t[2] = frac_new(0, 1); /* z = 0 for 2D */

        /* Store */
        SymmOp *op = &sg.ops[sg.count];
        memcpy(op->R, R, sizeof(R));
        op->t[0] = t[0];
        op->t[1] = t[1];
        op->t[2] = t[2];
        sg.count++;
    }

    /* Free line buffers */
    for (int i = 0; i < num_lines; i++) free(all_lines[i]);

    return sg;
}

/* ========================================================================
 * unique_rotations — extract unique R matrices from SymmGroup
 * ======================================================================== */

int unique_rotations(const SymmGroup *sg, int result[/* MAX_OPS */][3][3]) {
    int count = 0;
    for (int i = 0; i < sg->count && count < MAX_OPS; i++) {
        int found = 0;
        for (int j = 0; j < count; j++) {
            if (equal_int_mat(sg->ops[i].R, result[j])) {
                found = 1;
                break;
            }
        }
        if (!found) {
            memcpy(result[count], sg->ops[i].R, sizeof(int[3][3]));
            count++;
        }
    }
    return count;
}

/* ========================================================================
 * metric_inverse — compute inverse of the metric tensor G
 * ======================================================================== */

void metric_inverse(double a, double b, double c,
                    double alpha, double beta, double gamma,
                    double Ginv[3][3]) {
    double ca = cos(alpha * M_PI / 180.0);
    double cb = cos(beta  * M_PI / 180.0);
    double cg = cos(gamma * M_PI / 180.0);

    /* Direct metric tensor G (covariant) */
    double G[3][3];
    G[0][0] = a * a;          G[0][1] = a * b * cg;   G[0][2] = a * c * cb;
    G[1][0] = a * b * cg;     G[1][1] = b * b;        G[1][2] = b * c * ca;
    G[2][0] = a * c * cb;     G[2][1] = b * c * ca;   G[2][2] = c * c;

    /* Invert 3x3 matrix manually */
    double det = G[0][0] * (G[1][1]*G[2][2] - G[1][2]*G[2][1])
               - G[0][1] * (G[1][0]*G[2][2] - G[1][2]*G[2][0])
               + G[0][2] * (G[1][0]*G[2][1] - G[1][1]*G[2][0]);

    if (fabs(det) < 1e-30) {
        memset(Ginv, 0, sizeof(double) * 9);
        return;
    }

    double inv_det = 1.0 / det;
    Ginv[0][0] = (G[1][1]*G[2][2] - G[1][2]*G[2][1]) * inv_det;
    Ginv[0][1] = (G[0][2]*G[2][1] - G[0][1]*G[2][2]) * inv_det;
    Ginv[0][2] = (G[0][1]*G[1][2] - G[0][2]*G[1][1]) * inv_det;
    Ginv[1][0] = (G[1][2]*G[2][0] - G[1][0]*G[2][2]) * inv_det;
    Ginv[1][1] = (G[0][0]*G[2][2] - G[0][2]*G[2][0]) * inv_det;
    Ginv[1][2] = (G[0][2]*G[1][0] - G[0][0]*G[1][2]) * inv_det;
    Ginv[2][0] = (G[1][0]*G[2][1] - G[1][1]*G[2][0]) * inv_det;
    Ginv[2][1] = (G[0][1]*G[2][0] - G[0][0]*G[2][1]) * inv_det;
    Ginv[2][2] = (G[0][0]*G[1][1] - G[0][1]*G[1][0]) * inv_det;
}

double q2_metric(int h, int k, int l, const double Ginv[3][3]) {
    double v[3] = {(double)h, (double)k, (double)l};
    double rv[3];
    rv[0] = Ginv[0][0]*v[0] + Ginv[0][1]*v[1] + Ginv[0][2]*v[2];
    rv[1] = Ginv[1][0]*v[0] + Ginv[1][1]*v[1] + Ginv[1][2]*v[2];
    rv[2] = Ginv[2][0]*v[0] + Ginv[2][1]*v[1] + Ginv[2][2]*v[2];
    return v[0]*rv[0] + v[1]*rv[1] + v[2]*rv[2];
}

/* ========================================================================
 * direct_basis_from_lattice_info — build 3x3 direct basis matrix A
 * ======================================================================== */

void direct_basis_from_lattice_info(const LatticeInfo *info, double A[3][3]) {
    double a = info->a, b = info->b, c = info->c;
    double alpha = info->alpha * M_PI / 180.0;
    double beta  = info->beta  * M_PI / 180.0;
    double gamma = info->gamma * M_PI / 180.0;

    /* Direct lattice vectors in Cartesian frame */
    double a_vec[3] = {a, 0.0, 0.0};
    double b_vec[3] = {b * cos(gamma), b * sin(gamma), 0.0};

    double c_x = c * cos(beta);
    double sin_g = sin(gamma);
    double c_y;
    if (fabs(sin_g) < 1e-15) {
        c_y = 0.0;
    } else {
        c_y = c * (cos(alpha) - cos(beta) * cos(gamma)) / sin_g;
    }
    double c_z_sq = 1.0 - cos(beta)*cos(beta) - (cos(alpha) - cos(beta)*cos(gamma)) *
                    (cos(alpha) - cos(beta)*cos(gamma)) / (sin_g * sin_g);
    double c_z = (c_z_sq > 0.0) ? c * sqrt(c_z_sq) : 0.0;
    double c_vec[3] = {c_x, c_y, c_z};

    /* Columns are the basis vectors */
    A[0][0] = a_vec[0]; A[1][0] = a_vec[1]; A[2][0] = a_vec[2];
    A[0][1] = b_vec[0]; A[1][1] = b_vec[1]; A[2][1] = b_vec[2];
    A[0][2] = c_vec[0]; A[1][2] = c_vec[1]; A[2][2] = c_vec[2];
}

/* ========================================================================
 * lattice_params_from_basis — extract parameters from direct basis matrix A
 * ======================================================================== */

static double calc_angle(double u[3], double v[3]) {
    double dot = u[0]*v[0] + u[1]*v[1] + u[2]*v[2];
    double nu = sqrt(u[0]*u[0] + u[1]*u[1] + u[2]*u[2]);
    double nv = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    if (nu < 1e-30 || nv < 1e-30) return 90.0;
    double cuv = dot / (nu * nv);
    if (cuv > 1.0) cuv = 1.0;
    if (cuv < -1.0) cuv = -1.0;
    return atan2(sqrt(1.0 - cuv*cuv), cuv) * 180.0 / M_PI;
}

void lattice_params_from_basis(const double A[3][3], double *a, double *b, double *c,
                                double *alpha, double *beta, double *gamma) {
    double a_vec[3] = {A[0][0], A[1][0], A[2][0]};
    double b_vec[3] = {A[0][1], A[1][1], A[2][1]};
    double c_vec[3] = {A[0][2], A[1][2], A[2][2]};

    double norm_a = sqrt(a_vec[0]*a_vec[0] + a_vec[1]*a_vec[1] + a_vec[2]*a_vec[2]);
    double norm_b = sqrt(b_vec[0]*b_vec[0] + b_vec[1]*b_vec[1] + b_vec[2]*b_vec[2]);
    double norm_c = sqrt(c_vec[0]*c_vec[0] + c_vec[1]*c_vec[1] + c_vec[2]*c_vec[2]);

    *a = norm_a;
    *b = norm_b;
    *c = norm_c;

    *alpha = calc_angle(b_vec, c_vec);
    *beta  = calc_angle(a_vec, c_vec);
    *gamma = calc_angle(a_vec, b_vec);
}

/* ========================================================================
 * star_from_hkl — apply all rotations to an HKL vector
 * Returns count of unique vectors in the star
 * ======================================================================== */

int star_from_hkl(const int hkl[3],
                  const int rotations[/* n_rotations */][3][3],
                  int n_rotations,
                  const double Ginv[3][3],
                  int result[/* MAX_STAR_VECTORS */][3]) {
    int count = 0;

    for (int r = 0; r < n_rotations; r++) {
        /* v = R^T @ hkl */
        int v[3] = {0, 0, 0};
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                v[i] += rotations[r][j][i] * hkl[j];
            }
        }

        /* Check for duplicates */
        int dup = 0;
        for (int k = 0; k < count; k++) {
            if (result[k][0] == v[0] && result[k][1] == v[1] && result[k][2] == v[2]) {
                dup = 1;
                break;
            }
        }
        if (!dup) {
            result[count][0] = v[0];
            result[count][1] = v[1];
            result[count][2] = v[2];
            count++;
            if (count >= MAX_STAR_VECTORS) break;
        }
    }

    /* Sort by q^2, then lexicographically */
    for (int i = 0; i < count - 1; i++) {
        for (int j = i + 1; j < count; j++) {
            double q2_i = q2_metric(result[i][0], result[i][1], result[i][2], Ginv);
            double q2_j = q2_metric(result[j][0], result[j][1], result[j][2], Ginv);
            int swap = 0;
            if (q2_i > q2_j) {
                swap = 1;
            } else if (fabs(q2_i - q2_j) < 1e-12) {
                if (result[i][0] > result[j][0]) swap = 1;
                else if (result[i][0] == result[j][0]) {
                    if (result[i][1] > result[j][1]) swap = 1;
                    else if (result[i][1] == result[j][1]) {
                        if (result[i][2] > result[j][2]) swap = 1;
                    }
                }
            }
            if (swap) {
                int tmp[3] = {result[i][0], result[i][1], result[i][2]};
                result[i][0] = result[j][0]; result[i][1] = result[j][1]; result[i][2] = result[j][2];
                result[j][0] = tmp[0]; result[j][1] = tmp[1]; result[j][2] = tmp[2];
            }
        }
    }

    return count;
}

/* ========================================================================
 * get_family_key_lexicographical — canonical representative of a star
 * ======================================================================== */

void get_family_key_lexicographical(const int star_vectors[/* N */][3],
                                     int N,
                                     int key[3]) {
    if (N <= 0) {
        key[0] = key[1] = key[2] = 0;
        return;
    }
    int best = 0;
    for (int i = 1; i < N; i++) {
        int pos_best = (star_vectors[best][0] > 0) + (star_vectors[best][1] > 0) + (star_vectors[best][2] > 0);
        int pos_i   = (star_vectors[i][0] > 0) + (star_vectors[i][1] > 0) + (star_vectors[i][2] > 0);
        if (pos_i > pos_best) {
            best = i;
        } else if (pos_i == pos_best) {
            if (star_vectors[i][0] > star_vectors[best][0] ||
               (star_vectors[i][0] == star_vectors[best][0] &&
                (star_vectors[i][1] > star_vectors[best][1] ||
                (star_vectors[i][1] == star_vectors[best][1] &&
                 star_vectors[i][2] > star_vectors[best][2])))) {
                best = i;
            }
        }
    }
    key[0] = star_vectors[best][0];
    key[1] = star_vectors[best][1];
    key[2] = star_vectors[best][2];
}

/* ========================================================================
 * star_is_closed — does -G exist in the star?
 * ======================================================================== */

int star_is_closed(const int star[/* N */][3], int N) {
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
 * point_group_has_neg_identity — is -I in the rotation set?
 * ======================================================================== */

int point_group_has_neg_identity(const int rotations[/* N */][3][3], int N) {
    int negI[3][3];
    memset(negI, 0, sizeof(negI));
    negI[0][0] = -1; negI[1][1] = -1; negI[2][2] = -1;

    for (int i = 0; i < N; i++) {
        if (equal_int_mat(rotations[i], negI)) return 1;
    }
    return 0;
}

/* ========================================================================
 * find_inversion_ops — find operations with R = -I
 * ======================================================================== */

int find_inversion_ops(const SymmGroup *sg,
                        Fraction t_result[/* MAX_OPS */][3],
                        int at_origin[/* MAX_OPS */]) {
    int negI[3][3];
    memset(negI, 0, sizeof(negI));
    negI[0][0] = -1; negI[1][1] = -1; negI[2][2] = -1;

    int count = 0;
    for (int i = 0; i < sg->count && count < MAX_OPS; i++) {
        if (equal_int_mat(sg->ops[i].R, negI)) {
            t_result[count][0] = sg->ops[i].t[0];
            t_result[count][1] = sg->ops[i].t[1];
            t_result[count][2] = sg->ops[i].t[2];
            at_origin[count] = is_zero_mod1_vec(sg->ops[i].t);
            count++;
        }
    }
    return count;
}

/* ========================================================================
 * phase_factor — exp(2πi * h·t)
 * ======================================================================== */

double phase_factor_real(const int hkl[3], const Fraction t[3]) {
    double dot = 0.0;
    for (int i = 0; i < 3; i++) {
        dot += (double)hkl[i] * frac_to_float(t[i]);
    }
    double d = dot - floor(dot);
    return cos(TAU * d);
}

double phase_factor_imag(const int hkl[3], const Fraction t[3]) {
    double dot = 0.0;
    for (int i = 0; i < 3; i++) {
        dot += (double)hkl[i] * frac_to_float(t[i]);
    }
    double d = dot - floor(dot);
    return sin(TAU * d);
}

/* ========================================================================
 * relationships_in_star — compute phase constraints within a star
 * ======================================================================== */

StarRelationships relationships_in_star(const SymmGroup *sg,
                                         const int star[/* N */][3],
                                         int N) {
    StarRelationships sr;
    memset(&sr, 0, sizeof(sr));

    for (int op_i = 0; op_i < sg->count; op_i++) {
        const int (*R)[3] = sg->ops[op_i].R;
        const Fraction *t = sg->ops[op_i].t;

        for (int gj = 0; gj < N; gj++) {
            const int *Gj = star[gj];

            /* Gk = R^T @ Gj */
            int Gk[3] = {0, 0, 0};
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    Gk[i] += R[j][i] * Gj[j];
                }
            }

            /* Check if Gk is in the star */
            int gk_idx = -1;
            for (int k = 0; k < N; k++) {
                if (star[k][0] == Gk[0] && star[k][1] == Gk[1] && star[k][2] == Gk[2]) {
                    gk_idx = k;
                    break;
                }
            }
            if (gk_idx < 0) continue;

            /* Compute phase */
            double pr = phase_factor_real(Gj, t);
            double pi = phase_factor_imag(Gj, t);

            /* Unary contradiction check */
            if (gk_idx == gj) {
                double tol = 1e-14;
                if (fabs(pr - 1.0) > tol || fabs(pi) > tol) {
                    sr.unary_contra_count++;
                }
            }

            /* Store relationship */
            int count = sr.rels_count[gk_idx];
            if (count < MAX_STAR_RELS_PER_VECTOR) {
                sr.rels[gk_idx][count].source = gj;
                sr.rels[gk_idx][count].phase_real = pr;
                sr.rels[gk_idx][count].phase_imag = pi;
                sr.rels[gk_idx][count].op_index = op_i;
                sr.rels_count[gk_idx] = count + 1;
            }
        }
    }

    /* Binary contradictions */
    for (int gk_idx = 0; gk_idx < N; gk_idx++) {
        int src_indices[MAX_STAR_RELS_PER_VECTOR];
        double ph_real[MAX_STAR_RELS_PER_VECTOR];
        double ph_imag[MAX_STAR_RELS_PER_VECTOR];
        int src_count = sr.rels_count[gk_idx];

        for (int g = 0; g < src_count; g++) {
            src_indices[g] = sr.rels[gk_idx][g].source;
            ph_real[g] = sr.rels[gk_idx][g].phase_real;
            ph_imag[g] = sr.rels[gk_idx][g].phase_imag;
        }

        for (int s = 0; s < src_count; s++) {
            for (int e = s + 1; e < src_count; e++) {
                if (src_indices[e] == src_indices[s]) {
                    double diff = sqrt((ph_real[e]-ph_real[s])*(ph_real[e]-ph_real[s]) +
                                       (ph_imag[e]-ph_imag[s])*(ph_imag[e]-ph_imag[s]));
                    if (diff > 1e-14) {
                        sr.binary_contra_count++;
                    }
                }
            }
        }
    }

    return sr;
}

/* ========================================================================
 * solve_star_coeffs — BFS over phase constraint graph
 * ======================================================================== */

int solve_star_coeffs(const int star[/* N */][3], int N,
                       const StarRelationships *rels,
                       double ref_real,
                       double coeffs_real[/* N */],
                       double coeffs_imag[/* N */],
                       char warnings[/* MAX_WARNINGS */][256],
                       int max_warnings) {
    int assigned[MAX_STAR_VECTORS];
    memset(assigned, 0, sizeof(int) * N);

    int adj_dst[MAX_STAR_VECTORS][MAX_ADJ];
    double adj_pr[MAX_STAR_VECTORS][MAX_ADJ];
    double adj_pi[MAX_STAR_VECTORS][MAX_ADJ];
    int adj_count[MAX_STAR_VECTORS];
    memset(adj_count, 0, sizeof(adj_count));

    for (int k = 0; k < N; k++) {
        for (int ri = 0; ri < rels->rels_count[k]; ri++) {
            int source = rels->rels[k][ri].source;
            double pr = rels->rels[k][ri].phase_real;
            double pi = rels->rels[k][ri].phase_imag;

            if (adj_count[source] < MAX_ADJ) {
                adj_dst[source][adj_count[source]] = k;
                adj_pr[source][adj_count[source]] = pr;
                adj_pi[source][adj_count[source]] = pi;
                adj_count[source]++;
            }
        }
    }

    /* Pick reference: lexicographically largest tuple */
    int ref = 0;
    for (int i = 1; i < N; i++) {
        if (star[i][0] > star[ref][0] ||
           (star[i][0] == star[ref][0] && star[i][1] > star[ref][1]) ||
           (star[i][0] == star[ref][0] && star[i][1] == star[ref][1] && star[i][2] > star[ref][2])) {
            ref = i;
        }
    }

    double tol = 1e-14;
    int warn_count = 0;

    int queue[MAX_STAR_VECTORS];
    int qhead = 0, qtail = 0;

    coeffs_real[ref] = ref_real;
    coeffs_imag[ref] = 0.0;
    assigned[ref] = 1;
    queue[qtail++] = ref;

    while (qhead < qtail) {
        int u = queue[qhead++];
        double cu_r = coeffs_real[u];
        double cu_i = coeffs_imag[u];

        for (int a = 0; a < adj_count[u]; a++) {
            int v = adj_dst[u][a];
            double ph_r = adj_pr[u][a];
            double ph_i = adj_pi[u][a];

            double cv_r = ph_r * cu_r - ph_i * cu_i;
            double cv_i = ph_r * cu_i + ph_i * cu_r;

            if (!assigned[v]) {
                coeffs_real[v] = cv_r;
                coeffs_imag[v] = cv_i;
                assigned[v] = 1;
                queue[qtail++] = v;
            } else {
                double diff = sqrt((coeffs_real[v] - cv_r)*(coeffs_real[v] - cv_r) +
                                   (coeffs_imag[v] - cv_i)*(coeffs_imag[v] - cv_i));
                if (diff > tol && warn_count < max_warnings) {
                    snprintf(warnings[warn_count], 256,
                             "Inconsistent phase on %d->%d: have (%.6f,%.6f), propose (%.6f,%.6f)",
                             u, v, coeffs_real[v], coeffs_imag[v], cv_r, cv_i);
                    warn_count++;
                }
            }
        }
    }

    int assigned_count = 0;
    for (int i = 0; i < N; i++) {
        if (assigned[i]) assigned_count++;
    }
    if (assigned_count != N && warn_count < max_warnings) {
        snprintf(warnings[warn_count], 256,
                 "Star graph not fully connected; unassigned nodes present");
        warn_count++;
        for (int i = 0; i < N; i++) {
            if (!assigned[i]) {
                coeffs_real[i] = ref_real;
                coeffs_imag[i] = 0.0;
            }
        }
    }

    return warn_count;
}

/* ========================================================================
 * read_final_values — parse simulation output files
 * ======================================================================== */

SimulationValues read_final_values(const char *file_path) {
    SimulationValues sv;
    memset(&sv, 0, sizeof(sv));
    sv.final_Nx = sv.final_Ny = sv.final_Nz = 0;

    FILE *f = fopen(file_path, "r");
    if (!f) {
        fprintf(stderr, "symmetry_ops: cannot open %s\n", file_path);
        return sv;
    }

    char line[512];
    while (fgets(line, sizeof(line), f)) {
        size_t len = strlen(line);
        while (len > 0 && (line[len-1] == '\n' || line[len-1] == '\r' || line[len-1] == ' '))
            line[--len] = '\0';

        if (strncmp(line, "final Nx", 8) == 0)
            sscanf(line, "final Nx %d", &sv.final_Nx);
        else if (strncmp(line, "final Ny", 8) == 0)
            sscanf(line, "final Ny %d", &sv.final_Ny);
        else if (strncmp(line, "final Nz", 8) == 0)
            sscanf(line, "final Nz %d", &sv.final_Nz);
        else if (strncmp(line, "final Lx", 8) == 0)
            sscanf(line, "final Lx %lf", &sv.final_Lx);
        else if (strncmp(line, "final Ly", 8) == 0)
            sscanf(line, "final Ly %lf", &sv.final_Ly);
        else if (strncmp(line, "final Lz", 8) == 0)
            sscanf(line, "final Lz %lf", &sv.final_Lz);
        else if (strncmp(line, "Nkai", 4) == 0)
            sscanf(line, "Nkai %lf", &sv.Nkai);
        else if (strncmp(line, "f", 1) == 0) {
            if (line[1] == '\0' || line[1] == ' ' || line[1] == '=')
                sscanf(line, "f %lf", &sv.f);
        }
        else if (strncmp(line, "final free energy", 17) == 0)
            sscanf(line, "final free energy %lf", &sv.free_energy);
    }

    fclose(f);
    return sv;
}

/* ========================================================================
 * float_to_miller_int — convert float Miller vector to integer triple
 * ======================================================================== */

void float_to_miller_int(const double hB[3], int max_den, double tol,
                          int primitive,
                          int result[3]) {
    int all_zero = 1;
    for (int i = 0; i < 3; i++) {
        if (fabs(hB[i]) > tol) { all_zero = 0; break; }
    }
    if (all_zero) { result[0] = 0; result[1] = 0; result[2] = 0; return; }

    int fracs[3];
    int denoms[3];
    for (int i = 0; i < 3; i++) {
        if (fabs(hB[i]) < tol) {
            fracs[i] = 0; denoms[i] = 1;
        } else {
            double x = hB[i];
            if (x < 0) {
                fracs[i] = -(int)round(fabs(x) * max_den);
            } else {
                fracs[i] = (int)round(x * max_den);
            }
            denoms[i] = max_den;
        }
    }

    int lcm_val = 1;
    for (int i = 0; i < 3; i++) {
        lcm_val = lcm_val * denoms[i] / frac_gcd(lcm_val, denoms[i]);
    }

    for (int i = 0; i < 3; i++) {
        result[i] = fracs[i] * (lcm_val / denoms[i]);
    }

    if (primitive) {
        int g = 0;
        for (int i = 0; i < 3; i++) g = frac_gcd(g, result[i]);
        if (g > 1) {
            for (int i = 0; i < 3; i++) result[i] /= g;
        }
    }
}

/* ========================================================================
 * family_planes_info — generate valid Miller index modes
 *
 * This is a high-complexity function. For initial testing, we provide
 * a simplified version that works for small N values.
 * ======================================================================== */

FamilyPlanesInfo* family_planes_info(int N,
                                     const double Ginv[3][3],
                                     const SymmGroup *sg,
                                     int dim) {
    /* Allocate on heap to avoid ~114MB stack usage */
    FamilyPlanesInfo* info = (FamilyPlanesInfo*)malloc(sizeof(FamilyPlanesInfo));
    if (!info) return NULL;
    memset(info, 0, sizeof(FamilyPlanesInfo));

    info->all_rels = (StarRelationships*)calloc(MAX_MODES, sizeof(StarRelationships));
    if (!info->all_rels) {
        free(info);
        return NULL;
    }

    /* Build unique rotation set */
    int unique_R[MAX_OPS][3][3];
    int n_unique = unique_rotations(sg, unique_R);

    int current_max = 5;
    int step = 5;
    int hard_limit = 50;
    int seen_count = 0;

    /* Seen keys — track which family keys we've already accepted */
    int seen_keys[MAX_MODES][3];
    memset(seen_keys, 0, sizeof(seen_keys));

    /* Grid arrays (heap-allocated) */
    int grid_size_max = 2000000;
    int *grid_H = (int *)malloc(sizeof(int) * grid_size_max);
    int *grid_K = (int *)malloc(sizeof(int) * grid_size_max);
    int *grid_L = (int *)malloc(sizeof(int) * grid_size_max);
    double *grid_q2 = (double *)malloc(sizeof(double) * grid_size_max);

    while (seen_count < N) {
        if (current_max > hard_limit) {
            fprintf(stderr, "Warning: hard limit at max_index=%d without finding %d modes\n",
                    hard_limit, N);
            break;
        }

        int r = current_max;
        int gs = 0;

        for (int hh = -r; hh <= r && gs < grid_size_max; hh++) {
            for (int kk = -r; kk <= r && gs < grid_size_max; kk++) {
                if (dim == 2) {
                    int ll = 0;
                    if (hh == 0 && kk == 0 && ll == 0) continue;
                    grid_H[gs] = hh;
                    grid_K[gs] = kk;
                    grid_L[gs] = ll;
                    double v[3] = {(double)hh, (double)kk, (double)ll};
                    grid_q2[gs] = v[0]*(Ginv[0][0]*v[0]+Ginv[0][1]*v[1]+Ginv[0][2]*v[2]) +
                                  v[1]*(Ginv[1][0]*v[0]+Ginv[1][1]*v[1]+Ginv[1][2]*v[2]) +
                                  v[2]*(Ginv[2][0]*v[0]+Ginv[2][1]*v[1]+Ginv[2][2]*v[2]);
                    gs++;
                } else {
                    for (int ll = -r; ll <= r && gs < grid_size_max; ll++) {
                        if (hh == 0 && kk == 0 && ll == 0) continue;
                        grid_H[gs] = hh;
                        grid_K[gs] = kk;
                        grid_L[gs] = ll;
                        double v[3] = {(double)hh, (double)kk, (double)ll};
                        grid_q2[gs] = v[0]*(Ginv[0][0]*v[0]+Ginv[0][1]*v[1]+Ginv[0][2]*v[2]) +
                                      v[1]*(Ginv[1][0]*v[0]+Ginv[1][1]*v[1]+Ginv[1][2]*v[2]) +
                                      v[2]*(Ginv[2][0]*v[0]+Ginv[2][1]*v[1]+Ginv[2][2]*v[2]);
                        gs++;
                    }
                }
            }
        }

        /* Sort by q2, then lexicographic */
        for (int i = 0; i < gs - 1; i++) {
            for (int j = i + 1; j < gs; j++) {
                int swap = 0;
                if (grid_q2[i] > grid_q2[j]) swap = 1;
                else if (fabs(grid_q2[i] - grid_q2[j]) < 1e-12) {
                    if (grid_H[i] > grid_H[j]) swap = 1;
                    else if (grid_H[i] == grid_H[j]) {
                        if (grid_K[i] > grid_K[j]) swap = 1;
                        else if (grid_K[i] == grid_K[j]) {
                            if (grid_L[i] > grid_L[j]) swap = 1;
                        }
                    }
                }
                if (swap) {
                    int th = grid_H[i]; grid_H[i] = grid_H[j]; grid_H[j] = th;
                    int tk = grid_K[i]; grid_K[i] = grid_K[j]; grid_K[j] = tk;
                    int tl = grid_L[i]; grid_L[i] = grid_L[j]; grid_L[j] = tl;
                    double tq = grid_q2[i]; grid_q2[i] = grid_q2[j]; grid_q2[j] = tq;
                }
            }
        }

        for (int g = 0; g < gs; g++) {
            if (seen_count >= N) break;

            int hkl[3] = {grid_H[g], grid_K[g], grid_L[g]};

            int star_buf[MAX_STAR_VECTORS][3];
            int star_count = star_from_hkl(hkl, unique_R, n_unique, Ginv, star_buf);

            int key[3];
            get_family_key_lexicographical(star_buf, star_count, key);

            /* Check if already seen */
            int dup = 0;
            for (int s = 0; s < seen_count; s++) {
                if (seen_keys[s][0] == key[0] && seen_keys[s][1] == key[1] && seen_keys[s][2] == key[2]) {
                    dup = 1;
                    break;
                }
            }
            if (dup) continue;

            /* Compute relationships */
            StarRelationships rels = relationships_in_star(sg, star_buf, star_count);

            /* Skip if there are contradictions */
            if (rels.unary_contra_count > 0 || rels.binary_contra_count > 0) {
                continue;
            }

            /* Store */
            int m = seen_count;
            memcpy(info->family_keys[m], key, 3 * sizeof(int));
            info->star_counts[m] = star_count;
            memcpy(info->stars[m], star_buf, star_count * 3 * sizeof(int));
            info->all_rels[m] = rels;
            seen_keys[m][0] = key[0];
            seen_keys[m][1] = key[1];
            seen_keys[m][2] = key[2];
            seen_count++;
        }

        current_max += step;
    }

    info->count = seen_count;

    free(grid_H);
    free(grid_K);
    free(grid_L);
    free(grid_q2);

    return info;
}

/* ========================================================================
 * family_planes_info_free — release heap allocation
 * ======================================================================== */

void family_planes_info_free(FamilyPlanesInfo* info) {
    if (!info) return;
    free(info->all_rels);
    free(info);
}

/* ========================================================================
 * transform_miller_between_lattices
 * ======================================================================== */

void transform_miller_between_lattices(const LatticeInfo *lattice_A_info,
                                        const double P[3][3],
                                        const int existing_keys[/* K */][3],
                                        int num_keys,
                                        const double existing_coeffs_real[/* K */][MAX_STAR_VECTORS],
                                        const double existing_coeffs_imag[/* K */][MAX_STAR_VECTORS],
                                        const int existing_star_counts[/* K */],
                                        int (*trans_keys)[3],
                                        double trans_coeffs_real[MAX_STAR_VECTORS][MAX_STAR_VECTORS],
                                        double trans_coeffs_imag[MAX_STAR_VECTORS][MAX_STAR_VECTORS],
                                        int *trans_counts,
                                        double *new_lattice_a, double *new_lattice_b, double *new_lattice_c,
                                        double *new_lattice_alpha, double *new_lattice_beta, double *new_lattice_gamma) {
    double A[3][3];
    direct_basis_from_lattice_info(lattice_A_info, A);

    double B[3][3];
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
            B[i][j] = A[i][0]*P[0][j] + A[i][1]*P[1][j] + A[i][2]*P[2][j];
        }

    double PT[3][3];
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            PT[i][j] = P[j][i];

    int out_count = 0;
    for (int f = 0; f < num_keys && out_count < MAX_STAR_VECTORS; f++) {
        double hA[3] = {
            (double)existing_keys[f][0],
            (double)existing_keys[f][1],
            (double)existing_keys[f][2]
        };
        double hB[3] = {
            PT[0][0]*hA[0] + PT[0][1]*hA[1] + PT[0][2]*hA[2],
            PT[1][0]*hA[0] + PT[1][1]*hA[1] + PT[1][2]*hA[2],
            PT[2][0]*hA[0] + PT[2][1]*hA[1] + PT[2][2]*hA[2]
        };

        float_to_miller_int(hB, 8, 1e-8, 0, trans_keys[out_count]);
        trans_counts[out_count] = 1;

        for (int v = 0; v < existing_star_counts[f]; v++) {
            trans_coeffs_real[out_count][v] = existing_coeffs_real[f][v];
            trans_coeffs_imag[out_count][v] = existing_coeffs_imag[f][v];
        }
        out_count++;
    }

    lattice_params_from_basis(B, new_lattice_a, new_lattice_b, new_lattice_c,
                             new_lattice_alpha, new_lattice_beta, new_lattice_gamma);
}
