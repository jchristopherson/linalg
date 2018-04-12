// c_test_qr.c

#include "c_test_core.h"

bool test_qr_factor() {
    // Local Variables
    const int m = 100;
    const int n = 100;
    const double tol = 1.0e-8;

    double a[m*n], r1[m*n], r2[m*n], q1[m*m], q2[m*m], p2[n*n], a1[m*n], 
        a2[m*n], tau1[MIN(m,n)], tau2[MIN(m,n)];
    int i, pvt2[n];
    bool rst;

    // Initialization
    rst = true;
    make_rand_mtx(m, n, a);
    for (i = 0; i < m * n; ++i) r1[i] = r2[i] = a[i];

    // Compute the factorization of A
    qr_factor(m, n, r1, tau1, NULL);

    // Extract Q and R, and then check that Q * R = A
    form_qr(m, n, r1, tau1, q1, NULL);
    mtx_mult(false, false, m, n, m, 1.0, q1, m, r1, m, 0.0, a1);
    if (!is_dbl_mtx_equal(m, n, a, a1, tol)) {
        rst = false;
        printf("Test Failed: QR Factorization Test 1\n");
    }

    // Compute the QR factorization of A with pivoting
    for (i = 0; i < n; ++i) pvt2[i] = 0;
    qr_factor_pivot(m, n, r2, tau2, pvt2, NULL);

    // Extract Q and R, and then check that Q * R = A * P
    form_qr_pivot(m, n, r2, tau2, pvt2, q2, p2, NULL);
    mtx_mult(false, false, m, n, n, 1.0, a, m, p2, n, 0.0, a1);
    mtx_mult(false, false, m, n, m, 1.0, q2, m, r2, m, 0.0, a2);
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
    int i, pvt2[n];
    bool rst;

    // Initialization
    rst = true;
    make_rand_mtx(m, n, a);
    for (i = 0; i < m * n; ++i) r1[i] = r2[i] = a[i];

    // Compute the factorization of A
    qr_factor(m, n, r1, tau1, NULL);

    // Extract Q and R, and then check that Q * R = A
    form_qr(m, n, r1, tau1, q1, NULL);
    mtx_mult(false, false, m, n, m, 1.0, q1, m, r1, m, 0.0, a1);
    if (!is_dbl_mtx_equal(m, n, a, a1, tol)) {
        rst = false;
        printf("Test Failed: Overdetermined QR Factorization Test 1\n");
    }

    // Compute the QR factorization of A with pivoting
    for (i = 0; i < n; ++i) pvt2[i] = 0;
    qr_factor_pivot(m, n, r2, tau2, pvt2, NULL);

    // Extract Q and R, and then check that Q * R = A * P
    form_qr_pivot(m, n, r2, tau2, pvt2, q2, p2, NULL);
    mtx_mult(false, false, m, n, n, 1.0, a, m, p2, n, 0.0, a1);
    mtx_mult(false, false, m, n, m, 1.0, q2, m, r2, m, 0.0, a2);
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
    int i, pvt2[n];
    bool rst;

    // Initialization
    rst = true;
    make_rand_mtx(m, n, a);
    for (i = 0; i < m * n; ++i) r1[i] = r2[i] = a[i];

    // Compute the factorization of A
    qr_factor(m, n, r1, tau1, NULL);

    // Extract Q and R, and then check that Q * R = A
    form_qr(m, n, r1, tau1, q1, NULL);
    mtx_mult(false, false, m, n, m, 1.0, q1, m, r1, m, 0.0, a1);
    if (!is_dbl_mtx_equal(m, n, a, a1, tol)) {
        rst = false;
        printf("Test Failed: Underdetermined QR Factorization Test 1\n");
    }

    // Compute the QR factorization of A with pivoting
    for (i = 0; i < n; ++i) pvt2[i] = 0;
    qr_factor_pivot(m, n, r2, tau2, pvt2, NULL);

    // Extract Q and R, and then check that Q * R = A * P
    form_qr_pivot(m, n, r2, tau2, pvt2, q2, p2, NULL);
    mtx_mult(false, false, m, n, n, 1.0, a, m, p2, n, 0.0, a1);
    mtx_mult(false, false, m, n, m, 1.0, q2, m, r2, m, 0.0, a2);
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
    int i;

    // Initialization
    rst = true;
    make_rand_mtx(m, n, a);
    make_rand_mtx(m, n, c1);
    for (i = 0; i < m * n; ++i) c2[i] = c1[i];

    // Generate the QR factorization of A
    qr_factor(m, n, a, tau, NULL);
    for (i = 0; i < m * n; ++i) r[i] = a[i];
    form_qr(m, n, a, tau, q, NULL);

    // Compute C = Q * C
    mult_qr(false, m, n, r, tau, c1, NULL);

    // Compute ANS = Q * C
    mtx_mult(false, false, m, n, m, 1.0, q, m, c2, m, 0.0, ans);

    // Test
    if (!is_dbl_mtx_equal(m, n, c1, ans, tol)) {
        rst = false;
        printf("Test Failed: QR Multiplication Test 1\n");
    }

    // Compute C = Q**T * C
    for (i = 0; i < m * n; ++i) c1[i] = c2[i];
    mult_qr(true, m, n, r, tau, c1, NULL);

    // Compute ANS = Q**T * C
    mtx_mult(true, false, m, n, m, 1.0, q, m, c2, m, 0.0, ans);

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
    int i;

    // Initialization
    rst = true;
    make_rand_mtx(m, n, a);
    make_rand_mtx(m, n, c1);
    for (i = 0; i < m * n; ++i) c2[i] = c1[i];

    // Generate the QR factorization of A
    qr_factor(m, n, a, tau, NULL);
    for (i = 0; i < m * n; ++i) r[i] = a[i];
    form_qr(m, n, a, tau, q, NULL);

    // Compute C = Q * C
    mult_qr(false, m, n, r, tau, c1, NULL);

    // Compute ANS = Q * C
    mtx_mult(false, false, m, n, m, 1.0, q, m, c2, m, 0.0, ans);

    // Test
    if (!is_dbl_mtx_equal(m, n, c1, ans, tol)) {
        rst = false;
        printf("Test Failed: Overdetermined QR Multiplication Test 1\n");
    }

    // Compute C = Q**T * C
    for (i = 0; i < m * n; ++i) c1[i] = c2[i];
    mult_qr(true, m, n, r, tau, c1, NULL);

    // Compute ANS = Q**T * C
    mtx_mult(true, false, m, n, m, 1.0, q, m, c2, m, 0.0, ans);

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
    int i;

    // Initialization
    rst = true;
    make_rand_mtx(m, n, a);
    make_rand_mtx(m, n, c1);
    for (i = 0; i < m * n; ++i) c2[i] = c1[i];

    // Generate the QR factorization of A
    qr_factor(m, n, a, tau, NULL);
    for (i = 0; i < m * n; ++i) r[i] = a[i];
    form_qr(m, n, a, tau, q, NULL);

    // Compute C = Q * C
    mult_qr(false, m, n, r, tau, c1, NULL);

    // Compute ANS = Q * C
    mtx_mult(false, false, m, n, m, 1.0, q, m, c2, m, 0.0, ans);

    // Test
    if (!is_dbl_mtx_equal(m, n, c1, ans, tol)) {
        rst = false;
        printf("Test Failed: Underdetermined QR Multiplication Test 1\n");
    }

    // Compute C = Q**T * C
    for (i = 0; i < m * n; ++i) c1[i] = c2[i];
    mult_qr(true, m, n, r, tau, c1, NULL);

    // Compute ANS = Q**T * C
    mtx_mult(true, false, m, n, m, 1.0, q, m, c2, m, 0.0, ans);

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
    int i, j;
    bool rst;

    // Initialization
    rst = true;
    make_rand_mtx(m, n, a);
    make_rand_mtx(m, nrhs, b);
    make_rand_mtx(m, 1, b2a);
    for (i = 0; i < m * n; ++i) a1[i] = a[i];
    for (i = 0; i < m * nrhs; ++i) b1[i] = b[i];
    for (i = 0; i < m; ++i) b2[i] = b2a[i];

    // Compute the QR factorization of A
    qr_factor(m, n, a1, tau, NULL);

    // Solve the system of equations
    solve_qr(m, n, nrhs, a1, tau, b1, NULL);
    for (j = 0; j < nrhs; ++j)
        for (i = 0; i < n; ++i)
            x1[INDEX(i,j,n)] = b1[INDEX(i,j,m)];
    
    // Test
    mtx_mult(false, false, m, nrhs, n, 1.0, a, m, x1, n, 0.0, ans1);
    if (!is_dbl_mtx_equal(m, nrhs, ans1, b, tol)) {
        rst = false;
        printf("Test Failed: QR Solution Test 1, No Pivoting\n");
    }

    // Solve the system of equations
    solve_qr(m, n, 1, a1, tau, b2, NULL);
    for (i = 0; i < n; ++i) x2[i] = b2[i];

    // Test 2
    mtx_mult(false, false, m, 1, n, 1.0, a, m, x2, n, 0.0, ans2);
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
    const int n = 100;
    const int nrhs = 20;
    const double tol = 1.0e-8;

    double a[m*n], a1[m*n], b[m*nrhs], b1[m*nrhs], ans1[m*nrhs], x1[n*nrhs],
        tau[MIN(m,n)];
    int i, j, pvt[n];
    bool rst;

    // Initialization
    rst = true;
    make_rand_mtx(m, n, a);
    make_rand_mtx(m, nrhs, b);
    for (i = 0; i < m * n; ++i) a1[i] = a[i];
    for (i = 0; i < m * nrhs; ++i) b1[i] = b[i];
    for (i = 0; i < n; ++i) pvt[i] = 0;

    // Compute the QR factorization of A
    qr_factor_pivot(m, n, a1, tau, pvt, NULL);

    // Solve the system of equations
    solve_qr_pivot(m, n, nrhs, a, tau, pvt, b1, NULL);
    for (j = 0; j < nrhs; ++j)
        for (i = 0; i < n; ++i)
            x1[INDEX(i,j,n)] = b1[INDEX(i,j,m)];
    
    // Test
    mtx_mult(false, false, m, nrhs, n, 1.0, a, m, x1, n, 0.0, ans1);
    if (!is_dbl_mtx_equal(m, nrhs, ans1, b, tol)) {
        rst = false;
        printf("Test Failed: QR Solution Test 1, With Pivoting\n");
    }

    // End
    return rst;
}



bool test_qr_solve_pivot_ud() {
    // Local Variables
    const int m = 100;
    const int n = 200;
    const int nrhs = 20;
    const double tol = 1.0e-8;

    double a[m*n], a1[m*n], b[m*nrhs], b1[n*nrhs], ans1[m*nrhs], x1[n*nrhs],
        tau[MIN(m,n)];
    int i, j, pvt[n];
    bool rst;

    // Initialization
    rst = true;
    make_rand_mtx(m, n, a);
    make_rand_mtx(m, nrhs, b);
    for (i = 0; i < m * n; ++i) a1[i] = a[i];
    for (j = 0; j < nrhs; ++j)
        for (i = 0; i < m; ++i)
            b1[INDEX(i,j,n)] = b[INDEX(i,j,m)];
    for (i = 0; i < n; ++i) pvt[i] = 0;

    // Compute the QR factorization of A
    qr_factor_pivot(m, n, a1, tau, pvt, NULL);

    // Solve the system of equations
    solve_qr_pivot(m, n, nrhs, a, tau, pvt, b1, NULL);
    for (j = 0; j < nrhs; ++j)
        for (i = 0; i < n; ++i)
            x1[INDEX(i,j,n)] = b1[INDEX(i,j,n)];
    
    // Test
    mtx_mult(false, false, m, nrhs, n, 1.0, a, m, x1, n, 0.0, ans1);
    if (!is_dbl_mtx_equal(m, nrhs, ans1, b, tol)) {
        rst = false;
        printf("Test Failed: Underdetermined QR Solution Test 1, With Pivoting\n");
    }

    // End
    return rst;
}


bool test_qr_update() {
    // Local Variables
    const int m = 200;
    const int n = 100;
    const double tol = 1.0e-8;

    double a[m*n], a1[m*n], r[m*n], q[m*m], u[m], v[n], tau[n];
    bool rst;
    int i;

    // Initialization
    rst = true;
    make_rand_mtx(m, n, a);
    make_rand_mtx(m, 1, u);
    make_rand_mtx(n, 1, v);
    for (i = 0; i < m * n; ++i) r[i] = a1[i] = a[i];

    // Compute the QR factorization of A
    qr_factor(m, n, r, tau, NULL);

    // Form Q & R
    form_qr(m, n, r, tau, q, NULL);

    // Compute the rank 1 update A1 = A + u * v**T
    rank1_update(m, n, 1.0, u, v, a1);

    // Use the QR update to update the original R & Q
    qr_rank1_update(m, n, q, r, u, v, NULL);

    // Test that A1 = Q * R
    mtx_mult(false, false, m, n, m, 1.0, q, m, r, m, 0.0, a);
    if (!is_dbl_mtx_equal(m, n, a, a1, tol)) {
        rst = false;
        printf("Test Failed: Rank 1 QR Update\n");
    }

    // End
    return rst;
}
