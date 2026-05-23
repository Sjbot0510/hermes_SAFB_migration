/**
 * test_domain.c — Tests for domain.h / domain.c
 *
 * Validates: LatticeInfo construction, 2D factory,
 * ScatteringProfile allocation/free.
 */

#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "domain.h"

#define EPS 1e-10

static void test_lattice_info_3d(void)
{
    LatticeInfo l = lattice_info_new(1.0, 1.0, 1.0, 90.0, 90.0, 90.0, 3);
    assert(fabs(l.a - 1.0) < EPS);
    assert(fabs(l.b - 1.0) < EPS);
    assert(fabs(l.c - 1.0) < EPS);
    assert(fabs(l.alpha - 90.0) < EPS);
    assert(fabs(l.beta - 90.0) < EPS);
    assert(fabs(l.gamma - 90.0) < EPS);
    assert(l.dim == 3);
    printf("  [PASS] lattice_info_new (cubic)\n");
}

static void test_lattice_info_2d(void)
{
    LatticeInfo l = lattice_info_new_2d(1.0, 1.0, 120.0);
    assert(fabs(l.a - 1.0) < EPS);
    assert(fabs(l.b - 1.0) < EPS);
    assert(fabs(l.c - 1.0) < EPS);  /* dummy */
    assert(fabs(l.alpha - 90.0) < EPS);
    assert(fabs(l.beta - 90.0) < EPS);
    assert(fabs(l.gamma - 120.0) < EPS);
    assert(l.dim == 2);
    printf("  [PASS] lattice_info_new_2d (hexagonal)\n");
}

static void test_scattering_profile(void)
{
    int N = 5;
    ScatteringProfile sp = scattering_profile_new(N);
    assert(sp.num_peaks == N);
    assert(sp.hkl != NULL);
    assert(sp.q != NULL);
    assert(sp.intensity != NULL);

    /* Write and read back a value to verify memory is usable */
    sp.hkl[0] = 1;
    sp.hkl[1] = 2;
    sp.hkl[2] = 3;
    sp.q[0] = 1.234;
    sp.intensity[0] = 56.789;
    assert(sp.hkl[0] == 1);
    assert(fabs(sp.q[0] - 1.234) < EPS);
    assert(fabs(sp.intensity[0] - 56.789) < EPS);

    scattering_profile_free(&sp);
    assert(sp.hkl == NULL);
    assert(sp.q == NULL);
    assert(sp.intensity == NULL);
    printf("  [PASS] ScatteringProfile alloc/free\n");
}

int main(void)
{
    printf("Running domain tests...\n");
    test_lattice_info_3d();
    test_lattice_info_2d();
    test_scattering_profile();
    printf("All domain tests passed.\n");
    return 0;
}
