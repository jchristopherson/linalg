// c_test_qr.c

#include "c_test_core.h"

bool test_qr_factor() {
    // Local Variables
    const int m = 100;
    const int n = 100;
    const double tol = 1.0e-8;

    double a[m*n], r1[m*n], r2[m*n], q1[m*m], q2[m*m], p2[n*n], a1[m*n], 
        a2[m*n], tau1[MIN(m,n)], tau2[MIN(m,n)];
    int i, mn, pvt2[n];
    bool rst;

    // Initialization
    rst = true;
    mn = MIN(m, n);
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
    return rst;
}


bool test_qr_factor_od() {
    // Local Variables
    const int m = 200;
    const int n = 100;
    const double tol = 1.0e-8;

    double a[m*n], r1[m*n], r2[m*n], q1[m*m], q2[m*m], p2[n*n], a1[m*n], 
        a2[m*n], tau1[MIN(m,n)], tau2[MIN(m,n)];
    int i, mn, pvt2[n];
    bool rst;

    // Initialization
    rst = true;
    mn = MIN(m, n);
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
    return rst;
}



bool test_qr_factor_ud() {
    // Local Variables
    const int m = 100;
    const int n = 200;
    const double tol = 1.0e-8;

    double a[m*n], r1[m*n], r2[m*n], q1[m*m], q2[m*m], p2[n*n], a1[m*n], 
        a2[m*n], tau1[MIN(m,n)], tau2[MIN(m,n)];
    int i, mn, pvt2[n];
    bool rst;

    // Initialization
    rst = true;
    mn = MIN(m, n);
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
    return rst;
}



bool test_qr_mult() {
    // Local Variables
    const int m = 100;
    const int n = 100;
    const double tol = 1e-8;

    double a[m*n], r[m*n], c1[m*n], c2[m*n], ans[m*n], q[m*m], tau[MIN(m,n)];
    bool rst;
    int i, mn;

    // Initialization
    rst = true;
    mn = MIN(m, n);
    make_rand_mtx(m, n, a);
    make_rand_mtx(m, n, c1);
    for (i = 0; i < m * n; ++i) c2[i] = c1[i];

    // Generate the QR factorization of A
    qr_factor_(m, n, a, mn, tau, NULL);
    for (i = 0; i < m * n; ++i) r[i] = a[i];
    form_qr_(m, n, a, mn, tau, q, NULL);

    // Compute C = Q * C
    mult_qr_(false, m, n, r, mn, tau, c1, NULL);

    // Compute ANS = Q * C
    mtx_mult_(false, false, m, n, m, 1.0, q, m, c2, m, 0.0, ans);

    // Test
    if (!is_dbl_mtx_equal(m, n, c1, ans, tol)) {
        rst = false;
        printf("Test Failed: QR Multiplication Test 1\n");
    }

    // Compute C = Q**T * C
    for (i = 0; i < m * n; ++i) c1[i] = c2[i];
    mult_qr_(true, m, n, r, mn, tau, c1, NULL);

    // Compute ANS = Q**T * C
    mtx_mult_(true, false, m, n, m, 1.0, q, m, c2, m, 0.0, ans);

    // Test
    if (!is_dbl_mtx_equal(m, n, c1, ans, tol)) {
        rst = false;
        printf("Test Failed: QR Multiplication Test 2\n");
    }

    // End
    return rst;
}



bool test_qr_mult_od() {
    // Local Variables
    const int m = 200;
    const int n = 100;
    const double tol = 1e-8;

    double a[m*n], r[m*n], c1[m*n], c2[m*n], ans[m*n], q[m*m], tau[MIN(m,n)];
    bool rst;
    int i, mn;

    // Initialization
    rst = true;
    mn = MIN(m, n);
    make_rand_mtx(m, n, a);
    make_rand_mtx(m, n, c1);
    for (i = 0; i < m * n; ++i) c2[i] = c1[i];

    // Generate the QR factorization of A
    qr_factor_(m, n, a, mn, tau, NULL);
    for (i = 0; i < m * n; ++i) r[i] = a[i];
    form_qr_(m, n, a, mn, tau, q, NULL);

    // Compute C = Q * C
    mult_qr_(false, m, n, r, mn, tau, c1, NULL);

    // Compute ANS = Q * C
    mtx_mult_(false, false, m, n, m, 1.0, q, m, c2, m, 0.0, ans);

    // Test
    if (!is_dbl_mtx_equal(m, n, c1, ans, tol)) {
        rst = false;
        printf("Test Failed: Overdetermined QR Multiplication Test 1\n");
    }

    // Compute C = Q**T * C
    for (i = 0; i < m * n; ++i) c1[i] = c2[i];
    mult_qr_(true, m, n, r, mn, tau, c1, NULL);

    // Compute ANS = Q**T * C
    mtx_mult_(true, false, m, n, m, 1.0, q, m, c2, m, 0.0, ans);

    // Test
    if (!is_dbl_mtx_equal(m, n, c1, ans, tol)) {
        rst = false;
        printf("Test Failed: Overdetermined QR Multiplication Test 2\n");
    }

    // End
    return rst;
}



bool test_qr_mult_ud() {
    // Local Variables
    const int m = 100;
    const int n = 200;
    const double tol = 1e-8;

    double a[m*n], r[m*n], c1[m*n], c2[m*n], ans[m*n], q[m*m], tau[MIN(m,n)];
    bool rst;
    int i, mn;

    // Initialization
    rst = true;
    mn = MIN(m, n);
    make_rand_mtx(m, n, a);
    make_rand_mtx(m, n, c1);
    for (i = 0; i < m * n; ++i) c2[i] = c1[i];

    // Generate the QR factorization of A
    qr_factor_(m, n, a, mn, tau, NULL);
    for (i = 0; i < m * n; ++i) r[i] = a[i];
    form_qr_(m, n, a, mn, tau, q, NULL);

    // Compute C = Q * C
    mult_qr_(false, m, n, r, mn, tau, c1, NULL);

    // Compute ANS = Q * C
    mtx_mult_(false, false, m, n, m, 1.0, q, m, c2, m, 0.0, ans);

    // Test
    if (!is_dbl_mtx_equal(m, n, c1, ans, tol)) {
        rst = false;
        printf("Test Failed: Underdetermined QR Multiplication Test 1\n");
    }

    // Compute C = Q**T * C
    for (i = 0; i < m * n; ++i) c1[i] = c2[i];
    mult_qr_(true, m, n, r, mn, tau, c1, NULL);

    // Compute ANS = Q**T * C
    mtx_mult_(true, false, m, n, m, 1.0, q, m, c2, m, 0.0, ans);

    // Test
    if (!is_dbl_mtx_equal(m, n, c1, ans, tol)) {
        rst = false;
        printf("Test Failed: Underdetermined QR Multiplication Test 2\n");
    }

    // End
    return rst;
}

