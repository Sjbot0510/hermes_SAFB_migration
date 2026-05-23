/**
 * test_symmetry_ops.c — Tests for symmetry_ops.h / symmetry_ops.c
 *
 * Validates: fraction math, lattice basis, star generation, phase factors,
 * space group parsing, and star relationships.
 */

#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include "symmetry_ops.h"

#define EPS 1e-10
#define EPS_TOL 1e-14

static void test_frac_new(void) {
    Fraction f = frac_new(1, 4);
    assert(f.num == 1 && f.den == 4);
    assert(fabs(frac_to_float(f) - 0.25) < EPS);

    Fraction g = frac_new(2, 8);  /* should reduce to 1/4 */
    assert(g.num == 1 && g.den == 4);

    Fraction h = frac_new(-1, 3);
    assert(h.num == -1 && h.den == 3);
    assert(fabs(frac_to_float(h) - (-1.0/3.0)) < EPS);

    printf("  [PASS] frac_new reduction\n");
}

static void test_frac_mod1(void) {
    Fraction f = frac_mod1(frac_new(5, 4));
    assert(f.num == 1 && f.den == 4);  /* 5/4 mod 1 = 1/4 */

    Fraction g = frac_mod1(frac_new(-1, 3));
    assert(g.num == 2 && g.den == 3);  /* -1/3 mod 1 = 2/3 */

    Fraction h = frac_mod1(frac_new(3, 1));
    assert(h.num == 0 && h.den == 1);  /* 3 mod 1 = 0 */

    printf("  [PASS] frac_mod1\n");
}

static void test_is_zero_mod1_vec(void) {
    Fraction t1[3] = {frac_new(0,1), frac_new(0,1), frac_new(0,1)};
    assert(is_zero_mod1_vec(t1) == 1);

    Fraction t2[3] = {frac_new(1,2), frac_new(0,1), frac_new(0,1)};
    assert(is_zero_mod1_vec(t2) == 0);

    Fraction t3[3] = {frac_new(2,1), frac_new(0,1), frac_new(0,1)};
    assert(is_zero_mod1_vec(t3) == 1);  /* 2 mod 1 = 0 */

    printf("  [PASS] is_zero_mod1_vec\n");
}

static void test_equal_int_mat(void) {
    int A[3][3] = {{1,0,0},{0,1,0},{0,0,1}};
    int B[3][3] = {{1,0,0},{0,1,0},{0,0,1}};
    int C[3][3] = {{-1,0,0},{0,-1,0},{0,0,-1}};
    assert(equal_int_mat(A, B) == 1);
    assert(equal_int_mat(A, C) == 0);
    printf("  [PASS] equal_int_mat\n");
}

static void test_metric_inverse(void) {
    /* Cubic lattice: a=b=c=1, alpha=beta=gamma=90 */
    double Ginv[3][3];
    metric_inverse(1.0, 1.0, 1.0, 90.0, 90.0, 90.0, Ginv);

    /* Should be identity */
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            double expected = (i == j) ? 1.0 : 0.0;
            assert(fabs(Ginv[i][j] - expected) < EPS);
        }
    }
    printf("  [PASS] metric_inverse (cubic)\n");
}

static void test_q2_metric(void) {
    double Ginv[3][3];
    metric_inverse(1.0, 1.0, 1.0, 90.0, 90.0, 90.0, Ginv);

    assert(fabs(q2_metric(1,0,0,Ginv) - 1.0) < EPS);
    assert(fabs(q2_metric(1,1,0,Ginv) - 2.0) < EPS);
    assert(fabs(q2_metric(1,1,1,Ginv) - 3.0) < EPS);
    printf("  [PASS] q2_metric (cubic)\n");
}

static void test_direct_basis_from_lattice_info(void) {
    LatticeInfo info = lattice_info_new(1.0, 1.0, 1.0, 90.0, 90.0, 90.0, 3);
    double A[3][3];
    direct_basis_from_lattice_info(&info, A);

    /* Cubic: A = identity */
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
            double expected = (i == j) ? 1.0 : 0.0;
            assert(fabs(A[i][j] - expected) < EPS);
        }
    printf("  [PASS] direct_basis_from_lattice_info (cubic)\n");
}

static void test_lattice_params_from_basis(void) {
    /* b_vec = (0.5, sqrt(3)/2, 0) => |b_vec| = 1, gamma = 60° */
    double A[3][3] = {{1.0, 0.5, 0.0},
                       {0.0, 0.8660254037844386, 0.0},
                       {0.0, 0.0, 1.0}};

    double a, b, c, alpha, beta, gamma;
    lattice_params_from_basis(A, &a, &b, &c, &alpha, &beta, &gamma);

    assert(fabs(a - 1.0) < EPS);
    assert(fabs(b - 1.0) < EPS);
    assert(fabs(c - 1.0) < EPS);
    assert(fabs(gamma - 60.0) < EPS);  /* cos(gamma) = 0.5 => gamma = 60° */
    printf("  [PASS] lattice_params_from_basis\n");
}

static void test_star_from_hkl(void) {
    double Ginv[3][3];
    metric_inverse(1.0, 1.0, 1.0, 90.0, 90.0, 90.0, Ginv);

    /* Identity rotations only */
    int R[1][3][3];
    memset(R, 0, sizeof(R));
    R[0][0][0] = 1; R[0][1][1] = 1; R[0][2][2] = 1;

    int hkl[3] = {1, 0, 0};
    int result[MAX_STAR_VECTORS][3];
    int count = star_from_hkl(hkl, R, 1, Ginv, result);
    assert(count == 1);
    assert(result[0][0] == 1 && result[0][1] == 0 && result[0][2] == 0);

    /* 4-fold rotation about z */
    int R4[4][3][3];
    memset(R4, 0, sizeof(R4));
    R4[0][0][0]=1; R4[0][1][1]=1; R4[0][2][2]=1;  /* 0°   */
    R4[1][0][1]=1; R4[1][1][0]=-1; R4[1][2][2]=1; /* 90°  */
    R4[2][0][0]=-1; R4[2][1][1]=-1; R4[2][2][2]=1;/* 180° */
    R4[3][0][1]=-1; R4[3][1][0]=1; R4[3][2][2]=1;  /* 270° */

    count = star_from_hkl(hkl, R4, 4, Ginv, result);
    assert(count == 4);  /* (1,0,0), (0,1,0), (-1,0,0), (0,-1,0) */
    printf("  [PASS] star_from_hkl (4-fold)\n");
}

static void test_family_key_lexicographical(void) {
    int vecs[4][3] = {{-1,0,0},{0,1,0},{1,0,0},{0,-1,0}};
    int key[3];
    get_family_key_lexicographical(vecs, 4, key);
    assert(key[0] == 1 && key[1] == 0 && key[2] == 0);  /* (1,0,0) has most positives */
    printf("  [PASS] get_family_key_lexicographical\n");
}

static void test_star_is_closed(void) {
    int star[4][3] = {{1,0,0},{-1,0,0},{0,1,0},{0,-1,0}};
    assert(star_is_closed(star, 4) == 1);

    int star2[2][3] = {{1,0,0},{0,1,0}};
    assert(star_is_closed(star2, 2) == 0);
    printf("  [PASS] star_is_closed\n");
}

static void test_point_group_has_neg_identity(void) {
    int Rs[2][3][3];
    memset(Rs, 0, sizeof(Rs));
    Rs[0][0][0]=-1; Rs[0][1][1]=-1; Rs[0][2][2]=-1; /* -I */
    Rs[1][0][0]=1; Rs[1][1][1]=1; Rs[1][2][2]=1;    /* I */
    assert(point_group_has_neg_identity(Rs, 2) == 1);

    int Rs2[1][3][3];
    memset(Rs2, 0, sizeof(Rs2));
    Rs2[0][0][0]=1; Rs2[0][1][1]=1; Rs2[0][2][2]=1;
    assert(point_group_has_neg_identity(Rs2, 1) == 0);
    printf("  [PASS] point_group_has_neg_identity\n");
}

static void test_phase_factor(void) {
    Fraction t[3] = {frac_new(0,1), frac_new(0,1), frac_new(0,1)};
    int hkl[3] = {1, 0, 0};
    assert(fabs(phase_factor_real(hkl, t) - 1.0) < EPS);
    assert(fabs(phase_factor_imag(hkl, t)) < EPS);

    /* t = (1/4, 0, 0), h = (1,0,0) => exp(2πi * 1/4) = i */
    Fraction t2[3] = {frac_new(1,4), frac_new(0,1), frac_new(0,1)};
    assert(fabs(phase_factor_real(hkl, t2)) < EPS);
    assert(fabs(phase_factor_imag(hkl, t2) - 1.0) < EPS);
    printf("  [PASS] phase_factor\n");
}

static void test_find_inversion_ops(void) {
    SymmGroup sg;
    memset(&sg, 0, sizeof(sg));
    sg.dim = 3;

    /* Identity op */
    int R1[3][3];
    memset(R1, 0, sizeof(R1));
    R1[0][0]=1; R1[1][1]=1; R1[2][2]=1;
    sg.ops[0].R[0][0]=R1[0][0]; sg.ops[0].R[0][1]=R1[0][1]; sg.ops[0].R[0][2]=R1[0][2];
    sg.ops[0].R[1][0]=R1[1][0]; sg.ops[0].R[1][1]=R1[1][1]; sg.ops[0].R[1][2]=R1[1][2];
    sg.ops[0].R[2][0]=R1[2][0]; sg.ops[0].R[2][1]=R1[2][1]; sg.ops[0].R[2][2]=R1[2][2];
    sg.ops[0].t[0]=frac_new(0,1); sg.ops[0].t[1]=frac_new(0,1); sg.ops[0].t[2]=frac_new(0,1);
    sg.count = 1;

    /* Inversion op: R = -I, t = (0,0,0) */
    int R2[3][3];
    memset(R2, 0, sizeof(R2));
    R2[0][0]=-1; R2[1][1]=-1; R2[2][2]=-1;
    sg.ops[1].R[0][0]=R2[0][0]; sg.ops[1].R[0][1]=R2[0][1]; sg.ops[1].R[0][2]=R2[0][2];
    sg.ops[1].R[1][0]=R2[1][0]; sg.ops[1].R[1][1]=R2[1][1]; sg.ops[1].R[1][2]=R2[1][2];
    sg.ops[1].R[2][0]=R2[2][0]; sg.ops[1].R[2][1]=R2[2][1]; sg.ops[1].R[2][2]=R2[2][2];
    sg.ops[1].t[0]=frac_new(0,1); sg.ops[1].t[1]=frac_new(0,1); sg.ops[1].t[2]=frac_new(0,1);
    sg.count = 2;

    Fraction t_result[MAX_OPS][3];
    int at_origin[MAX_OPS];
    int inv_count = find_inversion_ops(&sg, t_result, at_origin);
    assert(inv_count == 1);
    assert(at_origin[0] == 1);
    printf("  [PASS] find_inversion_ops\n");
}

static void test_lattice_roundtrip(void) {
    /* Start with non-orthogonal lattice */
    LatticeInfo info;
    info = lattice_info_new(2.0, 3.0, 4.0, 90.0, 90.0, 60.0, 3);

    double A[3][3];
    direct_basis_from_lattice_info(&info, A);

    double a, b, c, alpha, beta, gamma;
    lattice_params_from_basis(A, &a, &b, &c, &alpha, &beta, &gamma);

    assert(fabs(a - 2.0) < EPS);
    assert(fabs(b - 3.0) < EPS);
    assert(fabs(c - 4.0) < EPS);
    assert(fabs(gamma - 60.0) < EPS);
    printf("  [PASS] lattice roundtrip (a=2,b=3,c=4,gamma=60°)\n");
}

static void test_unique_rotations(void) {
    SymmGroup sg;
    memset(&sg, 0, sizeof(sg));
    sg.dim = 3;

    /* Two identical rotations, one different */
    int R1[3][3]; memset(R1, 0, sizeof(R1));
    R1[0][0]=1; R1[1][1]=1; R1[2][2]=1;

    int R2[3][3]; memset(R2, 0, sizeof(R2));
    R2[0][0]=-1; R2[1][1]=-1; R2[2][2]=-1;

    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
            sg.ops[0].R[i][j] = R1[i][j];
            sg.ops[1].R[i][j] = R1[i][j];
            sg.ops[2].R[i][j] = R2[i][j];
            sg.ops[i].t[0]=frac_new(0,1); sg.ops[i].t[1]=frac_new(0,1); sg.ops[i].t[2]=frac_new(0,1);
        }
    sg.count = 3;

    int result[MAX_OPS][3][3];
    int n = unique_rotations(&sg, result);
    assert(n == 2);
    printf("  [PASS] unique_rotations\n");
}

int main(void) {
    printf("Running symmetry_ops tests...\n");
    test_frac_new();
    test_frac_mod1();
    test_is_zero_mod1_vec();
    test_equal_int_mat();
    test_metric_inverse();
    test_q2_metric();
    test_direct_basis_from_lattice_info();
    test_lattice_params_from_basis();
    test_star_from_hkl();
    test_family_key_lexicographical();
    test_star_is_closed();
    test_point_group_has_neg_identity();
    test_phase_factor();
    test_find_inversion_ops();
    test_lattice_roundtrip();
    test_unique_rotations();
    printf("All symmetry_ops tests passed.\n");
    return 0;
}
