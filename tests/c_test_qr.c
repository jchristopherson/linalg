// c_test_qr.c

#include "c_test_core.h"

bool test_qr_factor() {
    // Local Variables
    const int m = 100;
    const int n = 100;
    const double tol = 1.0e-8;

    double *a, *r1, *r2, *q1, *q2, *p2, *tau1, *tau2, *a1, *a2;
    int i, mn, pvt2[n];
    bool rst;

    // Initialization
    rst = true;
    mn = MIN(m, n);
    a = (double*)malloc((size_t)(m * n * sizeof(double)));
    r1 = (double*)malloc((size_t)(m * n * sizeof(double)));
    r2 = (double*)malloc((size_t)(m * n * sizeof(double)));
    q1 = (double*)malloc((size_t)(m * m * sizeof(double)));
    q2 = (double*)malloc((size_t)(m * m * sizeof(double)));
    p2 = (double*)malloc((size_t)(n * n * sizeof(double)));
    tau1 = (double*)malloc((size_t)(mn * sizeof(double)));
    tau2 = (double*)malloc((size_t)(mn * sizeof(double)));
    a1 = (double*)malloc((size_t)(m * n * sizeof(double)));
    a2 = (double*)malloc((size_t)(m * n * sizeof(double)));
    make_rand_mtx(m, n, a);
    for (i = 0; i < m * n; ++i) r1[i] = r2[i] = a[i];

    // Compute the factorization of A
    qr_factor_(m, n, r1, mn, tau1, NULL);

    // Extract Q and R, and then check that Q * R = A
    form_qr_(m, n, r1, mn, tau1, q1, NULL);
    mtx_mult_(false, false, m, n, m, 1.0, q1, m, r1, m, 0.0, a1);
    if (!is_dbl_mtx_equal(m, n, a, a1, tol)) {
        rst = false;
        printf("Test Failed: QR Factorization Test 1\n");
    }

    // Compute the QR factorization of A with pivoting
    for (i = 0; i < n; ++i) pvt2[i] = 0;
    qr_factor_pivot_(m, n, r2, mn, tau2, pvt2, NULL);

    // Extract Q and R, and then check that Q * R = A * P
    form_qr_pivot_(m, n, r2, mn, tau2, pvt2, q2, p2, NULL);
    mtx_mult_(false, false, m, n, n, 1.0, a, m, p2, n, 0.0, a1);
    mtx_mult_(false, false, m, n, m, 1.0, q2, m, r2, m, 0.0, a2);
    if (!is_dbl_mtx_equal(m, n, a1, a2, tol)) {
        rst = false;
        printf("Test Failed: QR Factorization Test 2\n");
    }

    // End
    free(a);
    free(r1);
    free(r2);
    free(q1);
    free(q2);
    free(p2);
    free(tau1);
    free(tau2);
    free(a1);
    free(a2);
    return rst;
}


bool test_qr_factor_od() {
    // Local Variables
    const int m = 200;
    const int n = 100;
    const double tol = 1.0e-8;

    double *a, *r1, *r2, *q1, *q2, *p2, *tau1, *tau2, *a1, *a2;
    int i, mn, pvt2[n];
    bool rst;

    // Initialization
    rst = true;
    mn = MIN(m, n);
    a = (double*)malloc((size_t)(m * n * sizeof(double)));
    r1 = (double*)malloc((size_t)(m * n * sizeof(double)));
    r2 = (double*)malloc((size_t)(m * n * sizeof(double)));
    q1 = (double*)malloc((size_t)(m * m * sizeof(double)));
    q2 = (double*)malloc((size_t)(m * m * sizeof(double)));
    p2 = (double*)malloc((size_t)(n * n * sizeof(double)));
    tau1 = (double*)malloc((size_t)(mn * sizeof(double)));
    tau2 = (double*)malloc((size_t)(mn * sizeof(double)));
    a1 = (double*)malloc((size_t)(m * n * sizeof(double)));
    a2 = (double*)malloc((size_t)(m * n * sizeof(double)));
    make_rand_mtx(m, n, a);
    for (i = 0; i < m * n; ++i) r1[i] = r2[i] = a[i];

    // Compute the factorization of A
    qr_factor_(m, n, r1, mn, tau1, NULL);

    // Extract Q and R, and then check that Q * R = A
    form_qr_(m, n, r1, mn, tau1, q1, NULL);
    mtx_mult_(false, false, m, n, m, 1.0, q1, m, r1, m, 0.0, a1);
    if (!is_dbl_mtx_equal(m, n, a, a1, tol)) {
        rst = false;
        printf("Test Failed: Overdetermined QR Factorization Test 1\n");
    }

    // Compute the QR factorization of A with pivoting
    for (i = 0; i < n; ++i) pvt2[i] = 0;
    qr_factor_pivot_(m, n, r2, mn, tau2, pvt2, NULL);

    // Extract Q and R, and then check that Q * R = A * P
    form_qr_pivot_(m, n, r2, mn, tau2, pvt2, q2, p2, NULL);
    mtx_mult_(false, false, m, n, n, 1.0, a, m, p2, n, 0.0, a1);
    mtx_mult_(false, false, m, n, m, 1.0, q2, m, r2, m, 0.0, a2);
    if (!is_dbl_mtx_equal(m, n, a1, a2, tol)) {
        rst = false;
        printf("Test Failed: Overdetermined QR Factorization Test 2\n");
    }

    // End
    free(a);
    free(r1);
    free(r2);
    free(q1);
    free(q2);
    free(p2);
    free(tau1);
    free(tau2);
    free(a1);
    free(a2);
    return rst;
}



bool test_qr_factor_ud() {
    // Local Variables
    const int m = 100;
    const int n = 200;
    const double tol = 1.0e-8;

    double *a, *r1, *r2, *q1, *q2, *p2, *tau1, *tau2, *a1, *a2;
    int i, mn, pvt2[n];
    bool rst;

    // Initialization
    rst = true;
    mn = MIN(m, n);
    a = (double*)malloc((size_t)(m * n * sizeof(double)));
    r1 = (double*)malloc((size_t)(m * n * sizeof(double)));
    r2 = (double*)malloc((size_t)(m * n * sizeof(double)));
    q1 = (double*)malloc((size_t)(m * m * sizeof(double)));
    q2 = (double*)malloc((size_t)(m * m * sizeof(double)));
    p2 = (double*)malloc((size_t)(n * n * sizeof(double)));
    tau1 = (double*)malloc((size_t)(mn * sizeof(double)));
    tau2 = (double*)malloc((size_t)(mn * sizeof(double)));
    a1 = (double*)malloc((size_t)(m * n * sizeof(double)));
    a2 = (double*)malloc((size_t)(m * n * sizeof(double)));
    make_rand_mtx(m, n, a);
    for (i = 0; i < m * n; ++i) r1[i] = r2[i] = a[i];

    // Compute the factorization of A
    qr_factor_(m, n, r1, mn, tau1, NULL);

    // Extract Q and R, and then check that Q * R = A
    form_qr_(m, n, r1, mn, tau1, q1, NULL);
    mtx_mult_(false, false, m, n, m, 1.0, q1, m, r1, m, 0.0, a1);
    if (!is_dbl_mtx_equal(m, n, a, a1, tol)) {
        rst = false;
        printf("Test Failed: Underdetermined QR Factorization Test 1\n");
    }

    // Compute the QR factorization of A with pivoting
    for (i = 0; i < n; ++i) pvt2[i] = 0;
    qr_factor_pivot_(m, n, r2, mn, tau2, pvt2, NULL);

    // Extract Q and R, and then check that Q * R = A * P
    form_qr_pivot_(m, n, r2, mn, tau2, pvt2, q2, p2, NULL);
    mtx_mult_(false, false, m, n, n, 1.0, a, m, p2, n, 0.0, a1);
    mtx_mult_(false, false, m, n, m, 1.0, q2, m, r2, m, 0.0, a2);
    if (!is_dbl_mtx_equal(m, n, a1, a2, tol)) {
        rst = false;
        printf("Test Failed: Underdetermined QR Factorization Test 2\n");
    }

    // End
    free(a);
    free(r1);
    free(r2);
    free(q1);
    free(q2);
    free(p2);
    free(tau1);
    free(tau2);
    free(a1);
    free(a2);
    return rst;
}
