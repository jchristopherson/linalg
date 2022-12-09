#include "c_test_core.h"
#include <math.h>
#include <stdlib.h>

double drand(double low, double high)
{
    return ((double)rand() * (high - low)) / (double)RAND_MAX + low;
}

double complex dcrand(double low, double high)
{
    double a, b;
    a = drand(low, high);
    b = drand(low, high);
    return a + b * I;
}




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

void mtx_mult(int m, int n, int p, const double *x, const double *y,
    double *z)
{
    int i, j, k;
    double temp;
    for (j = 0; j < n; ++j)
    {
        for (i = 0; i < m; ++i)
        {
            temp = 0.0;
            for (k = 0; k < p; ++k) temp += x[INDEX(i,k,m)] * y[INDEX(k,j,p)];
            z[INDEX(i,j,m)] = temp;
        }
    }
}

void cmplx_mtx_mult(int m, int n, int p, const double complex *x,
    const double complex *y, double complex *z)
{
    int i, j, k;
    double complex temp;
    for (j = 0; j < n; ++j)
    {
        for (i = 0; i < n; ++i)
        {
            temp = 0.0 + I * 0.0;
            for (k = 0; k < p; ++k) temp += x[INDEX(i,k,m)] * y[INDEX(k,j,p)];
            z[INDEX(i,j,m)] = temp;
        }
    }
}


void promote_diagonal(int n, const double *x, double *z)
{
    int i, j;
    for (j = 0; j < n; ++j)
    {
        for (i = 0; i < n; ++i)
        {
            z[INDEX(i,j,n)] = i == j ? x[j] : 0.0;
        }
    }
}

void cmplx_promote_diagonal(int n, const double complex *x, double complex *z)
{
    int i, j;
    for (j = 0; j < n; ++j)
    {
        for (i = 0; i < n; ++i)
        {
            z[INDEX(i,j,n)] = i == j ? x[j] : 0.0 + I * 0.0;
        }
    }
}


void transpose(int m, int n, const double *x, double *z)
{
    int i, j;
    for (j = 0; j < n; ++j)
    {
        for (i = 0; i < m; ++i)
        {
            z[INDEX(j,i,n)] = x[INDEX(i,j,m)];
        }
    }
}

void cmplx_transpose(int m, int n, const double complex *x, double complex *z)
{
    int i, j;
    for (j = 0; j < n; ++j)
    {
        for (i = 0; i < m; ++i)
        {
            z[INDEX(j,i,n)] = x[INDEX(i,j,m)];
        }
    }
}

void conj_transpose(int m, int n, const double complex *x, double complex *z)
{
    int i, j;
    for (j = 0; j < n; ++j)
    {
        for (i = 0; i < m; ++i)
        {
            z[INDEX(j,i,n)] = conj(x[INDEX(i,j,m)]);
        }
    }
}

void create_matrix(int m, int n, double *x)
{
    int i;
    for (i = 0; i < m * n; ++i)
    {
        x[i] = drand(-1.0, 1.0);
    }
}

void cmplx_create_matrix(int m, int n, double complex *x)
{
    int i;
    for (i = 0; i < m * n; ++i)
    {
        x[i] = dcrand(-1.0, 1.0);
    }
}



void create_symmetric_matrix(int n, double *x)
{
    int i, j;
    double *temp1, *temp2;
    temp1 = (double*)malloc((size_t)(n * n * sizeof(double)));
    if (!temp1) return;
    temp2 = (double*)malloc((size_t)(n * n * sizeof(double)));
    if (!temp2) return;
    create_matrix(n, n, temp1);
    transpose(n, n, temp1, temp2);
    mtx_mult(n, n, n, temp1, temp2, x);
    free(temp1);
    free(temp2);
}

void cmplx_create_symmetric_matrix(int n, double complex *x)
{
    int i, j;
    double complex *temp1, *temp2;
    temp1 = (double complex*)malloc((size_t)(n * n * sizeof(double complex)));
    if (!temp1) return;
    temp2 = (double complex*)malloc((size_t)(n * n * sizeof(double complex)));
    if (!temp2) return;
    cmplx_create_matrix(n, n, temp1);
    cmplx_transpose(n, n, temp1, temp2);
    cmplx_mtx_mult(n, n, n, temp1, temp2, x);
    free(temp1);
    free(temp2);
}



void identity_matrix(int n, double *x)
{
    int i, j;
    for (j = 0; j < n; ++j)
    {
        for (i = 0; i < n; ++i)
        {
            x[INDEX(i,j,n)] = i == j ? 1.0 : 0.0;
        }
    }
}

void cmplx_identity_matrix(int n, double complex *x)
{
    int i, j;
    for (j = 0; j < n; ++j)
    {
        for (i = 0; i < n; ++i)
        {
            x[INDEX(i,j,n)] = (i == j ? 1.0 : 0.0) + I * 0.0;
        }
    }
}



void copy_matrix(int m, int n, const double *src, double *dst)
{
    int i, j;
    for (j = 0; j < n; ++j)
    {
        for (i = 0; i < m; ++i)
        {
            dst[INDEX(i,j,m)] = src[INDEX(i,j,m)];
        }
    }
}

void cmplx_copy_matrix(int m, int n, const double complex *src, 
    double complex *dst)
{
    int i, j;
    for (j = 0; j < n; ++j)
    {
        for (i = 0; i < m; ++i)
        {
            dst[INDEX(i,j,m)] = src[INDEX(i,j,m)];
        }
    }
}




void create_array(int n, double *x)
{
    int i;
    for (i = 0; i < n; ++i) x[i] = drand(-1.0, 1.0);
}

void cmplx_create_array(int n, double complex *x)
{
    int i;
    for (i = 0; i < n; ++i) x[i] = dcrand(-1.0, 1.0);
}




void copy_array(int n, const double *src, double *dst)
{
    int i;
    for (i = 0; i < n; ++i) dst[i] = src[i];
}

void cmplx_copy_array(int n, const double complex *src, double complex *dst)
{
    int i;
    for (i = 0; i < n; ++i) dst[i] = src[i];
}


void add_matrix(int m, int n, double alpha, const double *x, double beta,
    const double *y, double *z)
{
    int i, j;
    for (j = 0; j < n; ++j)
    {
        for (i = 0; i < m; ++i)
        {
            z[INDEX(i,j,m)] = alpha * x[INDEX(i,j,m)] + beta * y[INDEX(i,j,m)];
        }
    }
}

void cmplx_add_matrix(int m, int n, double complex alpha, 
    const double complex *x, double complex beta, const double complex *y,
    double complex *z)
{
    int i, j;
    for (j = 0; j < n; ++j)
    {
        for (i = 0; i < m; ++i)
        {
            z[INDEX(i,j,m)] = alpha * x[INDEX(i,j,m)] + beta * y[INDEX(i,j,m)];
        }
    }
}

void scale_matrix(int m, int n, double x, double *y)
{
    int i;
    for (i = 0; i < m * n; ++i) y[i] *= x;
}

void cmplx_scale_matrix(int m, int n, double complex x, double complex *y)
{
    int i;
    for (i = 0; i < m * n; ++i) y[i] *= x;
}



void rank1_update(int m, int n, const double *x, const double *y, double *z)
{
    int i, j;
    for (j = 0; j < n; ++j)
    {
        for (i = 0; i < m; ++i)
        {
            z[INDEX(i,j,m)] = x[i] * y[j];
        }
    }
}

void cmplx_rank1_update(int m, int n, const double complex *x, 
    const double complex *y, double complex *z)
{
    int i, j;
    for (j = 0; j < n; ++j)
    {
        for (i = 0; i < m; ++i)
        {
            z[INDEX(i,j,m)] = x[i] * conj(y[j]);
        }
    }
}



double product(int n, const double *x)
{
    double ans = 1.0;
    int i;
    for (i = 0; i < n; ++i) ans *= x[i];
    return ans;
}

double complex cmplx_product(int n, const double complex * x)
{
    double complex ans = 1.0;
    int i;
    for (i = 0; i < n; ++i) ans *= x[i];
    return ans;
}





double sum(int n, const double *x)
{
    double ans = 0.0;
    int i;
    for (i = 0; i < n; ++i) ans += x[i];
    return ans;
}

double complex cmplx_sum(int n, const double complex *x)
{
    double complex ans = 0.0 + 0.0 * I;
    int i;
    for (i = 0; i < n; ++i) ans += x[i];
    return ans;
}
