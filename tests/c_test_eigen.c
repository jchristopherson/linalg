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

    double a[n*n], a1[n*n];
    double complex vecs[n*n], x[n*n], y[n*n], vals[n];
    bool rst;
    int i;

    // Initialization
    rst = true;
    make_rand_mtx(n, n, a);
    for (i = 0; i < n * n; ++i) a1[i] = a[i];

    // Compute the eigenvalues and vectors
    eigen_asymm_(n, a, vals, vecs, NULL);

    // Compute VECS * VALS, where VALS is a diagonal matrix
    diag_cmtx_rmult_(n, n, n, 1.0, vecs, vals, 0.0, x);

    // Compute A * VECS, and then test

    // End
    return rst;
}