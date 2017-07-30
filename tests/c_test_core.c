// c_test_core.c

#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "c_test_core.h"

bool isSeeded = false;

void make_rand_mtx(int m, int n, double *x) {
    int i, j;
    double denom;
    if (!isSeeded) {
        srand(time(NULL));
        isSeeded = true;
    }
    denom = (double)RAND_MAX;
    for (j = 0; j < n; ++j)
        for (i = 0; i < m; ++i)
            x[INDEX(i,j,m)] = rand() / denom;
}


bool is_dbl_mtx_equal(int m, int n, const double *x, const double *y, 
                      double tol)
{
    int i, j;
    for (j = 0; j < n; ++j) {
        for (i = 0; i < m; ++i) {
            if (abs(x[INDEX(i,j,m)] - y[INDEX(i,j,m)]) > tol) {
                return false;
            }
        }
    }
    return true;
}

bool is_cmplx_mtx_equal(int m, int n, const double complex *x, 
                        const double complex *y, double tol)
{
    int i, j;
    for (j = 0; j < n; ++j) {
        for (i = 0; i < m; ++i) {
            if (cabs(x[INDEX(i,j,m)] - y[INDEX(i,j,m)]) > tol) return false;
        }
    }
    return true;
}
