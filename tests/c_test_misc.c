// c_test_misc.c

#include "c_test_core.h"


bool test_diagonal_mtx_mult() {
    // Local Variables
    const int m = 30;
    const int n = 30;
    const int k = 30;
    const double tol = 1e-12;
    const double alpha = 0.5;
    const double beta = 0.25;

    int i, j;
    double c1[m*n], ans1[m*n], b1[k*n], d1[m*k], d1v[k];
    bool rst;

    // Initialization
    rst = true;
    make_rand_mtx(m, n, c1);
    make_rand_mtx(k, n, b1);
    make_rand_mtx(k, 1, d1v);

    for (j = 0; j < m*k; ++j) d1[j] = 0.0;
    for (j = 0; j < k; ++j) d1[INDEX(j,j,m)] = d1v[k];

    // Compute C1 = D1 * B1 + C1
    for (i = 0; i < m*n; ++i) ans1[i] = c1[i];
    diag_mtx_mult_(m, n, k, alpha, d1v, b1, beta, c1);
    mtx_mult_(false, false, m, n, k, alpha, d1, m, b1, k, beta, ans1);
    if (!is_dbl_mtx_equal(m, n, ans1, c1, tol)) {
        rst = false;
        printf("Test Failed: Diagonal Matrix Multiply Test 1\n");
    }

    // End
    return rst;
}


bool test_rank1_update() {
    // Local Variables
    const int m = 50;
    const int n = 20;
    const double alpha = 0.5;
    const double tol = 1e-12;

    double a[m*n], b[m*n], x[m], y[n];
    bool rst = true;
    int i, j;

    // Initialization
    make_rand_mtx(m, n, a);
    make_rand_mtx(m, 1, x);
    make_rand_mtx(n, 1, y);

    // Define the solution
    mtx_mult_(false, true, m, n, 1, alpha, x, m, y, n, 0.0, b);
    for (j = 0; j < n; ++j)
        for (i = 0; i < m; ++i)
            b[INDEX(i,j,m)] += a[INDEX(i,j,m)];
    
    // Test
    rank1_update_(m, n, alpha, x, y, a);
    if (!is_dbl_mtx_equal(m, n, a, b, tol)) {
        rst = false;
        printf("Test Failed: Rank 1 Update\n");
    }

    // End
    return rst;
}


bool test_rank() {
    // Local Variables
    const int m = 8;
    const int n = 6;
    double a[48] = {64.0, 9.0, 17.0, 40.0, 32.0, 41.0, 49.0, 8.0, 2.0, 55.0, 
        47.0, 26.0, 34.0, 23.0, 15.0, 58.0, 3.0, 54.0, 46.0, 27.0, 35.0, 22.0, 
        14.0, 59.0, 61.0, 12.0, 20.0, 37.0, 29.0, 44.0, 52.0, 5.0, 60.0, 13.0, 
        21.0, 36.0, 28.0, 45.0, 53.0, 4.0, 6.0, 51.0, 43.0, 30.0, 38.0, 19.0, 
        11.0, 62.0};
    bool rst = true;

    // The rank of A should be 3
    if (mtx_rank_(m, n, a, NULL) != 3) {
        rst = false;
        printf("Test Failed: Matrix Rank");
    }

    // End
    return rst;
}


bool test_tri_mtx_mult() {
    // Local Variables
    const int n = 100;
    const double alpha = 1.5;
    const double beta = -3.0;
    const double tol = 1.0e-12;

    bool check, rst;
    int i, j;
    double a[n*n], b[n*n], bans[n*n];

    // Initialization
    check = true;
    make_rand_mtx(n, n, a);
    for (j = 0; j < n; ++j)
        for (i = j + 1; i < n; ++i)
            a[INDEX(i,j,n)] = 0.0;
    
    // Test 1 (beta = 0)
    mtx_mult_(true, false, n, n, n, alpha, a, n, a, n, 0.0, bans);
    tri_mtx_mult_(true, n, alpha, a, 0.0, b, NULL);
    if (!is_dbl_mtx_equal(n, n, b, bans, tol)) {
        check = false;
        printf("Test Failed: Triangular Matrix Update - Test 1A\n");
    }
    rst = check;

    // Test 2 (beta /= 0)
    check = true;
    mtx_mult_(true, false, n, n, n,alpha, a, n, a, n, beta, bans);
    tri_mtx_mult_(true, n, alpha, a, beta, b, NULL);
    if (!is_dbl_mtx_equal(n, n, b, bans, tol)) {
        check = false;
        rst = false;
        printf("Test Failed: Triangular Matrix Update - Test 1B\n");
    }

    // End
    return rst;
}
