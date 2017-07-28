// c_test_misc.c

#include "c_test_core.h"


bool test_diagonal_mtx_mult() {
    // Local Variables
    const int m = 30;
    const int n = 30;
    const int k = 30;
    const int mn = 900;
    const int kn = 900;
    const int mk = 900;
    const double tol = 1e-8;
    const double alpha = 0.5;
    const double beta = 0.25;

    int i, j;
    double c1[mn], ans1[mn], b1[kn], b2[kn], ans2[kn], d1[mk], d1v[k];
    bool rst;

    // Initialization
    rst = true;
    make_rand_mtx(m, n, c1);
    make_rand_mtx(k, n, b1);
    make_rand_mtx(k, 1, d1v);

    for (j = 0; j < mk; ++j) d1[j] = 0.0;
    for (j = 0; j < k; ++j) d1[INDEX(j,j,m)] = d1v[k];

    // Compute C1 = D1 * B1 + C1
    diag_mtx_mult_(false, m, n, alpha, k, d1v, k, m, b1, beta, c1, NULL);
    mtx_mult_(false, false, m, n, k, alpha, d1, m, b1, k, beta, ans1);
    if (!is_dbl_mtx_equal(m, n, ans1, c1, tol)) {
        rst = false;
        printf("Test Failed: Diagonal Matrix Multiply Test 1\n");
    }


    // End
    return rst;
}