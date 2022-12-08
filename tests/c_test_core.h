#ifndef C_TEST_CORE_H_DEFINED
#define C_TEST_CORE_H_DEFINED

#include <stdbool.h>
#include <complex.h>

#define INDEX(i, j, lda) ((lda) * (j) + (i))

bool is_mtx_equal(int m, int n, const double *x, const double *y, double tol);
bool is_cmplx_mtx_equal(int m, int n, const double complex *x, 
    const double complex *y, double tol);

#endif