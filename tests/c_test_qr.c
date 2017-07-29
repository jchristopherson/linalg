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



bool test_qr_solve_no_pivot() {
    // Local Variables
    const int m = 100;
    const int n = 80;
    const int nrhs = 20;
    const double tol = 1.0e-8;

    double a[m*n], a1[m*n], b[m*nrhs], b1[m*nrhs], ans1[m*nrhs], x1[n*nrhs],
        tau[MIN(m,n)], b2a[m], b2[m], ans2[m], x2[n];
    int i, j, mn;
    bool rst;

    // Initialization
    rst = true;
    mn = MIN(m, n);
    make_rand_mtx(m, n, a);
    make_rand_mtx(m, nrhs, b);
    make_rand_mtx(m, 1, b2a);
    for (i = 0; i < m * n; ++i) a1[i] = a[i];
    for (i = 0; i < m * nrhs; ++i) b1[i] = b[i];
    for (i = 0; i < m; ++i) b2[i] = b2a[i];

    // Compute the QR factorization of A
    qr_factor_(m, n, a1, mn, tau, NULL);

    // Solve the system of equations
    solve_qr_(m, n, nrhs, a1, mn, tau, b1, NULL);
    for (j = 0; j < nrhs; ++j)
        for (i = 0; i < n; ++i)
            x1[INDEX(i,j,n)] = b1[INDEX(i,j,m)];
    
    // Test
    mtx_mult_(false, false, m, nrhs, n, 1.0, a, m, x1, n, 0.0, ans1);
    if (!is_dbl_mtx_equal(m, nrhs, ans1, b, tol)) {
        rst = false;
        printf("Test Failed: QR Solution Test 1, No Pivoting\n");
    }

    // Solve the system of equations
    solve_qr_(m, n, 1, a1, mn, tau, b2, NULL);
    for (i = 0; i < n; ++i) x2[i] = b2[i];

    // Test 2
    mtx_mult_(false, false, m, 1, n, 1.0, a, m, x2, n, 0.0, ans2);
    if (!is_dbl_mtx_equal(m, 1, ans2, b2a, tol)) {
        rst = false;
        printf("Test Failed: QR Solution Test 2, No Pivoting\n");
    }

    // End
    return rst;
}



bool test_qr_solve_pivot() {
    // Local Variables
    const int m = 100;
    const int n = 80;
    const int nrhs = 20;
    const double tol = 1.0e-8;

    double a[m*n], a1[m*n], b[m*nrhs], b1[m*nrhs], ans1[m*nrhs], x1[n*nrhs],
        tau[MIN(m,n)], b2a[m], b2[m], ans2[m], x2[n];
    int i, j, mn, pvt[n];
    bool rst;

    // Initialization
    rst = true;
    mn = MIN(m, n);
    make_rand_mtx(m, n, a);
    make_rand_mtx(m, nrhs, b);
    make_rand_mtx(m, 1, b2a);
    for (i = 0; i < m * n; ++i) a1[i] = a[i];
    for (i = 0; i < m * nrhs; ++i) b1[i] = b[i];
    for (i = 0; i < m; ++i) b2[i] = b2a[i];
    for (i = 0; i < n; ++i) pvt[i] = 0;

    // Compute the QR factorization of A
    qr_factor_pivot_(m, n, a1, mn, tau, pvt, NULL);

    // Solve the system of equations
    solve_qr_pivot_(m, n, nrhs, a, mn, tau, pvt, m, b1, NULL);
    for (j = 0; j < nrhs; ++j)
        for (i = 0; i < n; ++i)
            x1[INDEX(i,j,n)] = b1[INDEX(i,j,m)];
    
    // Test
    mtx_mult_(false, false, m, nrhs, n, 1.0, a, m, x1, n, 0.0, ans1);
    if (!is_dbl_mtx_equal(m, nrhs, ans1, b, tol)) {
        rst = false;
        printf("Test Failed: QR Solution Test 1, With Pivoting\n");
    }

    // Solve the system of equations
    solve_qr_pivot_(m, n, 1, a1, mn, tau, pvt, m, b2, NULL);
    for (i = 0; i < n; ++i) x2[i] = b2[i];

    // Test 2
    mtx_mult_(false, false, m, 1, n, 1.0, a, m, x2, n, 0.0, ans2);
    if (!is_dbl_mtx_equal(m, 1, ans2, b2a, tol)) {
        rst = false;
        printf("Test Failed: QR Solution Test 2, With Pivoting\n");
    }

    // End
    return rst;
}
