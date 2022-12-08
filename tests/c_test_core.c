#include "c_test_core.h"
#include <math.h>

bool is_mtx_equal(int m, int n, const double *x, const double *y, double tol)
{
    int i, j;
    for (j = 0; j < n; ++j)
    {
        for (i = 0; i < m; ++i)
        {
            if (fabs(x[INDEX(i,j,m)] - y[INDEX(i,j,m)]) > tol) return false;
        }
    }
    return true;
}

bool is_cmplx_mtx_equal(int m, int n, const double complex *x, 
    const double complex *y, double tol)
{
    int i, j;
    double xr, xi, yr, yi;
    for (j = 0; j < n; ++j)
    {
        for (i = 0; i < m; ++i)
        {
            xr = creal(x[INDEX(i,j,m)]);
            xi = cimag(x[INDEX(i,j,m)]);
            yr = creal(y[INDEX(i,j,m)]);
            yi = cimag(y[INDEX(i,j,m)]);
            if (fabs(xr - yr) > tol || fabs(xi - yi) > tol) return false;
        }
    }
    return true;
}