
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <fftw3.h>
#include "field.h"

int main() {
    printf("Debug test: 3x3x3 DC coefficient\n");

    int Na = 3, Nb = 3, Nc = 3;

    fftw_complex* Wa = (fftw_complex*)fftw_alloc_complex((size_t)Na * Nb * Nc);
    memset(Wa, 0, sizeof(fftw_complex) * (size_t)(Na * Nb * Nc));
    Wa[0][0] = 1.0;
    Wa[0][1] = 0.0;

    printf("Wa[0] before plan: real=%.6f imag=%.6f\n", Wa[0][0], Wa[0][1]);

    void* plan = field_create_fftw_plan(Na, Nb, Nc);
    assert(plan != NULL);
    
    fftw_complex* in = (fftw_complex*)field_get_fftw_input(plan);
    memcpy(in, Wa, sizeof(fftw_complex) * (size_t)(Na * Nb * Nc));
    
    printf("in[0] after memcpy: real=%.6f imag=%.6f\n", in[0][0], in[0][1]);
    
    field_run_fftw(plan);
    
    double* W_real = (double*)field_get_fftw_output(plan);
    double N = (double)Na * (double)Nb * (double)Nc;
    
    printf("\nOutput (first 5 elements):\n");
    for (int i = 0; i < 5; i++) {
        printf("  [%d] raw=%.6f / N=%.6f\n", i, W_real[i], W_real[i] / N);
    }
    
    double expected = 1.0 / N;
    printf("\nexpected = %.6f\n", expected);
    printf("W_real[0] raw = %.6f\n", W_real[0]);
    printf("W_real[0] / N = %.6f\n", W_real[0] / N);
    printf("diff = %.15e\n", fabs(W_real[0] / N - expected));
    
    field_destroy_fftw_plan(plan);
    fftw_free(Wa);
    
    if (fabs(W_real[0] / N - expected) < 1e-8) {
        printf("PASS!\n");
    } else {
        printf("FAIL!\n");
    }
    return 0;
}
