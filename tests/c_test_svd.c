// c_test_svd.c

#include "c_test_core.h"

bool test_svd() {
    // Local Variables
    const int m = 100;
    const int n = 100;
    const double tol = 1.0e-8;

    double a[m*n], a1[m*n], u1[m*m], vt1[n*n], s1[n], svt[m*n];
    bool rst;
    int i, mn;

    // Initialization
    rst = true;
    mn = MIN(m, n);
    make_rand_mtx(m, n, a);
    for (i = 0; i < m * n; ++i) a1[i] = a[i];

    // Compute the SVD of A
    svd(m, n, a1, s1, u1, vt1, NULL);

    // Ensure A = U * S * V**T
    diag_mtx_mult(m, n, n, 1.0, s1, vt1, 0.0, svt);
    mtx_mult(false, false, m, n, m, 1.0, u1, m, svt, m, 0.0, a1);

    // Test
    if (!is_dbl_mtx_equal(m, n, a1, a, tol)) {
        rst = false;
        printf("Test Failed: Singular Value Decomposition\n");
    }

    // End
    return rst;
}



bool test_svd_od() {
    // Local Variables
    const int m = 200;
    const int n = 100;
    const double tol = 1.0e-8;

    double a[m*n], a1[m*n], u1[m*m], vt1[n*n], s1[n], svt[m*n];
    bool rst;
    int i, mn;

    // Initialization
    rst = true;
    mn = MIN(m, n);
    make_rand_mtx(m, n, a);
    for (i = 0; i < m * n; ++i) a1[i] = a[i];

    // Compute the SVD of A
    svd(m, n, a1, s1, u1, vt1, NULL);

    // Ensure A = U * S * V**T
    diag_mtx_mult(m, n, n, 1.0, s1, vt1, 0.0, svt);
    mtx_mult(false, false, m, n, m, 1.0, u1, m, svt, m, 0.0, a1);

    // Test
    if (!is_dbl_mtx_equal(m, n, a1, a, tol)) {
        rst = false;
        printf("Test Failed: Overdetermined Singular Value Decomposition\n");
    }

    // End
    return rst;
}



bool test_svd_ud() {
    // Local Variables
    const int m = 100;
    const int n = 200;
    const double tol = 1.0e-8;

    double a[m*n], a1[m*n], u1[m*m], vt1[n*n], s1[n], svt[m*n];
    bool rst;
    int i, mn;

    // Initialization
    rst = true;
    mn = MIN(m, n);
    make_rand_mtx(m, n, a);
    for (i = 0; i < m * n; ++i) a1[i] = a[i];

    // Compute the SVD of A
    svd(m, n, a1, s1, u1, vt1, NULL);

    // Ensure A = U * S * V**T
    diag_mtx_mult(m, n, n, 1.0, s1, vt1, 0.0, svt);
    mtx_mult(false, false, m, n, m, 1.0, u1, m, svt, m, 0.0, a1);

    // Test
    if (!is_dbl_mtx_equal(m, n, a1, a, tol)) {
        rst = false;
        printf("Test Failed: Underdetermined Singular Value Decomposition\n");
    }

    // End
    return rst;
}
