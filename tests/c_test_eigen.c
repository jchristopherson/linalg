// c_test_eigen.c

#include "c_test_core.h"


bool test_eigen_symm() {
    // Local Variables
    const int n = 100;
    const double tol = 1.0e-8;

    double a[n*n], a1[n*n], vecs[n*n], x[n*n], y[n*n], vals[n];
    int i;
    bool rst;

    // Initialization
    rst = true;
    make_rand_mtx(n, n, a1);
    mtx_mult_(false, true, n, n, n, 1.0, a1, n, a1, n, 0.0, a);
    for (i = 0; i < n * n; ++i) vecs[i] = a[i];

    // Compute the eigenvalues and vectors of A
    eigen_symm_(n, true, vecs, vals, NULL);

    // Compute VECS * VALS where VALS is a diagonal matrix.
    diag_mtx_rmult_(n, n, n, 1.0, vecs, vals, 0.0, x);

    // Compute A * VECS, and then test
    mtx_mult_(false, false, n, n, n, 1.0, a, n, vecs, n, 0.0, y);
    if (!is_dbl_mtx_equal(n, n, x, y, tol)) {
        rst = false;
        printf("Test Failed: Symmetric Eigen\n");
    }

    // End
    return rst;
}



bool test_eigen_asymm() {
    // Local Variables
    const int n = 100;
    const double tol = 1.0e-8;

    double a[n*n];
    double complex vecs[n*n], x[n*n], y[n*n], vals[n], ac[n*n];
    bool rst;
    int i;

    // Initialization
    rst = true;
    make_rand_mtx(n, n, a);
    for (i = 0; i < n * n; ++i) ac[i] = a[i] + 0.0 * _Complex_I;

    // Compute the eigenvalues and vectors
    eigen_asymm_(n, a, vals, vecs, NULL);

    // Compute VECS * VALS, where VALS is a diagonal matrix
    diag_cmtx_rmult_(n, n, n, 1.0, vecs, vals, 0.0, x);

    // Compute A * VECS, and then test
    cmtx_mult_(false, false, n, n, n, 1.0, ac, n, vecs, n, 0.0, y);
    if (!is_cmplx_mtx_equal(n, n, x, y, tol)) {
        rst = false;
        printf("Test Failed: Asymmetric Eigen\n");
    }

    // End
    return rst;
}



bool test_eigen_gen() {
    // Local Variables
    const int n = 100;
    const double tol = 1.0e-8;

    double a[n*n], b[n*n], beta[n];
    double complex vals[n], vecs[n*n], x[n*n], y[n*n], ac[n*n], bc[n*n], t[n*n];
    int i;
    bool rst;

    // Initialization
    rst = true;
    make_rand_mtx(n, n, a);
    make_rand_mtx(n, n, b);
    for (i = 0; i < n * n; ++i) {
        ac[i] = a[i] + 0.0 * _Complex_I;
        bc[i] = b[i] + 0.0 * _Complex_I;
    }

    // Compute the eigenvalues and vectors
    eigen_gen_(n, a, b, vals, beta, vecs, NULL);
    for (i = 0; i < n; ++i) vals[i] /= beta[i];

    // Compute X = A * VECS
    cmtx_mult_(false, false, n, n, n, 1.0, ac, n, vecs, n, 0.0, x);

    // And Y = B * (VECS * VALS)
    diag_cmtx_rmult_(n, n, n, 1.0, vecs, vals, 0.0, t);
    cmtx_mult_(false, false, n, n, n, 1.0, bc, n, t, n, 0.0, y);
    if (!is_cmplx_mtx_equal(n, n, x, y, tol)) {
        rst = false;
        printf("Test Failed: Geenralized Eigen\n");
    }

    // End
    return rst;
}

