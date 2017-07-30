// c_test_inv.c

#include "c_test_core.h"


bool test_pinv() {
    // Local Variables
    const int m = 100;
    const int n = 100;
    const int nrhs = 20;
    const double tol = 1.0e-8;

    double a[m*n], a1[m * n], ainv[n *m], b[m*nrhs], x[n*nrhs], ans[m * nrhs];
    int i;
    bool rst;

    // Initialization
    rst = true;
    make_rand_mtx(m, n, a);
    make_rand_mtx(m, nrhs, b);
    for (i = 0; i < m * n; ++i) a1[i] = a[i];

    // Compute the inverse
    mtx_pinverse_(m, n, a1, ainv, NULL);

    // Compute X = inv(A) * B
    mtx_mult_(false, false, n, nrhs, m, 1.0, ainv, n, b, m, 0.0, x);

    // Test A *X = B
    mtx_mult_(false, false, m, nrhs, n, 1.0, a, m, x, n, 0.0, ans);
    if (!is_dbl_mtx_equal(m, nrhs, b, ans, tol)) {
        rst = false;
        printf("Test Failed: Pseudo-Inverse\n");
    }

    // End
    return rst;
}



bool test_pinv_od() {
    // Local Variables
    const int m = 200;
    const int n = 100;
    const int nrhs = 20;
    const double tol = 1.0e-8;

    double a[m*n], a1[m * n], ainv[n *m], b[m*nrhs], x[n*nrhs], ans[m * nrhs];
    int i;
    bool rst;

    // Initialization
    rst = true;
    make_rand_mtx(m, n, a);
    make_rand_mtx(m, nrhs, b);
    for (i = 0; i < m * n; ++i) a1[i] = a[i];

    // Compute the inverse
    mtx_pinverse_(m, n, a1, ainv, NULL);

    // Compute X = inv(A) * B
    mtx_mult_(false, false, n, nrhs, m, 1.0, ainv, n, b, m, 0.0, x);

    // Test A *X = B
    mtx_mult_(false, false, m, nrhs, n, 1.0, a, m, x, n, 0.0, ans);
    if (!is_dbl_mtx_equal(m, nrhs, b, ans, tol)) {
        rst = false;
        printf("Test Failed: Overdetermined Pseudo-Inverse\n");
    }

    // End
    return rst;
}



bool test_inv() {
    // Local Variables
    const int m = 100;
    const int nrhs = 20;
    const double tol = 1.0e-8;

    double a[m*m], a1[m * m], ainv[m *m];
    int i;
    bool rst;

    // Initialization
    rst = true;
    make_rand_mtx(m, m, a);
    for (i = 0; i < m * m; ++i) a1[i] = a[i];

    // Compute the inverse
    mtx_inverse_(m, a1, NULL);

    // Compute the psuedo-inverse
    mtx_pinverse_(m, m, a, ainv, NULL);

    // Test 
    if (!is_dbl_mtx_equal(m, m, a1, ainv, tol)) {
        rst = false;
        printf("Test Failed: Matrix Inverse\n");
    }

    // End
    return rst;
}
