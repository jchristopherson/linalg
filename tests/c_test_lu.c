// c_test_lu.c

#include "c_test_core.h"

bool test_lu_factor() {
    // Local Variables
    const int n = 100;
    const int nn = 10000;
    const double tol = 1.0e-8;

    double a[nn], a1[nn], l[nn], u[nn], p[nn], b[nn];
    int i, ipvt[n];
    bool rst;

    // Initialization
    rst = true;
    make_rand_mtx(n, n, a);
    for (i = 0; i < nn; ++i) a1[i] = a[i];

    // Compute the factorization
    lu_factor_(n, n, a1, n, ipvt, NULL);

    // Extract L, U, and P to determine if P * A = L * U
    for (i = 0; i < nn; ++i) l[i] = a1[i];
    form_lu_(n, l, ipvt, u, p);

    // Compute A1 = P * A
    mtx_mult_(false, false, n, n, n, 1.0, p, n, a, n, 0.0, a1);

    // Compute B = L * U
    mtx_mult_(false, false, n, n, n, 1.0, l, n, u, n, 0.0, b);

    // Test
    if (!is_dbl_mtx_equal(n, n, a1, b, tol)) {
        rst = false;
        printf("Test Failed: LU Factorization Test\n");
    }

    // End
    return rst;
}


bool test_lu_solve() {
    // Local Variables
    const int n = 100;
    const int nrhs = 20;
    const int nn = 10000;
    const int nnrhs = 2000;
    const double tol = 1.0e-8;

    double a[nn], a1[nn], b[nnrhs], x[nnrhs], b1[nnrhs];
    int i, ipvt[n];
    bool rst;

    // Initialization
    rst = true;
    make_rand_mtx(n, n, a);
    make_rand_mtx(n, nrhs, b);
    for (i = 0; i < nn; ++i) a1[i] = a[i];
    for (i = 0; i < nnrhs; ++i) x[i] = b[i];

    // Compute the factorization
    lu_factor_(n, n, a1, n, ipvt, NULL);

    // Solve for X
    solve_lu_(n, nrhs, a1, ipvt, x);

    // Test by determining if A * X = B
    mtx_mult_(false, false, n, nrhs, n, 1.0, a, n, x, n, 0.0, b1);
    if (!is_dbl_mtx_equal(n, nrhs, b, b1, tol)) {
        rst = false;
        printf("Test Failed: LU Factorization & Solution Test\n");
    }

    // End
    return rst;
}
