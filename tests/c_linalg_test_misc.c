#include "c_linalg_test.h"
#include "c_test_core.h"
#include "linalg.h"

bool test_diagonal_mtx_mult()
{
    const int m = 30;
    const int n = 30;
    const int k = 30;
    const double tol = 1.0e-8;
    const double alpha = 0.5;
    const double beta = 0.25;
    const int mn = m * n;
    const int kn = k * n;
    const int mk = m * k;
    const int kk = k * k;
    double c1[mn], ans1[mn], b1[kn], b2[kn], ans2[kn], d1[mk], d1v[k], D1[kk],
        temp1[kn];
    bool rst, check;
    int flag;

    // Initialization
    rst = true;
    create_matrix(m, n, c1);
    create_matrix(k, n, b1);
    create_array(k, d1v);
    copy_array(k, d1v, d1);
    promote_diagonal(k, d1v, D1);

    // Compute C1 = alpha * D1 * B1 + beta * C1
    mtx_mult(k, n, k, D1, b1, temp1);   // TEMP1 = D1 * B1
    add_matrix(m, n, alpha, temp1, beta, c1, ans1); // ANS1 = alpha * D1 * B1 + beta * C1

    flag = la_diag_mtx_mult(true, false, m, n, k, alpha, d1, b1, k, beta, c1, m);

    check = is_mtx_equal(m, n, c1, ans1, tol);
    if (!check) rst = false;
    

    // End
    return rst;
}




bool test_rank1_update()
{
    //
}