#include "c_linalg_test.h"
#include "c_test_core.h"
#include "linalg.h"
#include <math.h>

bool test_diagonal_mtx_mult()
{
    const int m = 30;
    const int n = 30;
    const int k = 30;
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
    if (flag != LA_NO_ERROR) rst = false;
    check = is_mtx_equal(m, n, c1, ans1, DBL_TOL);
    if (!check) rst = false;

    // End
    return rst;
}

bool test_cmplx_diagonal_mtx_mult()
{
    const int m = 30;
    const int n = 30;
    const int k = 30;
    const double complex alpha = 0.5 + 0.0 * I;
    const double complex beta = 0.25 + 0.0 * I;
    const int mn = m * n;
    const int kn = k * n;
    const int mk = m * k;
    const int kk = k * k;
    double complex c1[mn], ans1[mn], b1[kn], b2[kn], ans2[kn], d1[mk], d1v[k], 
        D1[kk], temp1[kn];
    bool rst, check;
    int flag;

    // Initialization
    rst = true;
    cmplx_create_matrix(m, n, c1);
    cmplx_create_matrix(k, n, b1);
    cmplx_create_array(k, d1v);
    cmplx_copy_array(k, d1v, d1);
    cmplx_promote_diagonal(k, d1v, D1);

    // Compute C1 = alpha * D1 * B1 + beta * C1
    cmplx_mtx_mult(k, n, k, D1, b1, temp1);   // TEMP1 = D1 * B1
    cmplx_add_matrix(m, n, alpha, temp1, beta, c1, ans1); // ANS1 = alpha * D1 * B1 + beta * C1

    flag = la_diag_mtx_mult_cmplx(true, LA_NO_OPERATION, m, n, k, alpha, d1, b1, k, beta, c1, m);
    if (flag != LA_NO_ERROR) rst = false;
    check = is_cmplx_mtx_equal(m, n, c1, ans1, DBL_TOL);
    if (!check) rst = false;
    
    // End
    return rst;
}



bool test_rank1_update()
{
    // Variables
    const int m = 50;
    const int n = 20;
    const double alpha = 0.5;
    const int mn = m * n;
    double a[mn], a1[mn], b[mn], x[m], y[n], temp1[mn];
    bool rst;
    int flag;

    // Initialization
    rst = true;
    create_matrix(m, n, a);
    copy_matrix(m, n, a, a1);
    create_array(m, x);
    create_array(n, y);

    // Compute the solution: A = alpha * x * y**T + A
    rank1_update(m, n, x, y, temp1);
    add_matrix(m, n, alpha, temp1, 1.0, a1, a1);    // A = A + alpha * TEMP1

    // Test
    flag = la_rank1_update(m, n, alpha, x, y, a, m);
    if (flag != LA_NO_ERROR) rst = false;
    if (!is_mtx_equal(m, n, a, a1, DBL_TOL)) rst = false;

    // End
    return rst;
}

bool test_cmplx_rank1_update()
{
    // Variables
    const int m = 50;
    const int n = 20;
    const double complex alpha = 0.5 + 0.0 * I;
    const double complex one = 1.0 + 0.0 * I;
    const int mn = m * n;
    double complex a[mn], a1[mn], b[mn], x[m], y[n], temp1[mn];
    bool rst;
    int flag;

    // Initialization
    rst = true;
    cmplx_create_matrix(m, n, a);
    cmplx_copy_matrix(m, n, a, a1);
    cmplx_create_array(m, x);
    cmplx_create_array(n, y);

    // Compute the solution: A = alpha * x * y**T + A
    cmplx_rank1_update(m, n, x, y, temp1);
    cmplx_add_matrix(m, n, alpha, temp1, one, a1, a1); // A = A + alpha * TEMP1

    // Test
    flag = la_rank1_update_cmplx(m, n, alpha, x, y, a, m);
    if (flag != LA_NO_ERROR) rst = false;
    if (!is_cmplx_mtx_equal(m, n, a, a1, DBL_TOL)) rst = false;

    // End
    return rst;
}




bool test_trace()
{
    // Variables
    const int n = 20;
    const int nn = n * n;
    double xd[n], x[nn], z, ans;
    bool rst;
    int flag;

    // Initialization
    rst = true;
    create_array(n, xd);
    promote_diagonal(n, xd, x);
    ans = sum(n, xd);

    // Test
    flag = la_trace(n, n, x, n, &z);
    if (flag != LA_NO_ERROR) rst = false;
    if (fabs(z - ans) > DBL_TOL) rst = false;

    // End
    return rst;
}

bool test_cmplx_trace()
{
    // Variables
    const int n = 20;
    const int nn = n * n;
    double complex xd[n], x[nn], z, ans;
    bool rst;
    int flag;

    // Initialization
    rst = true;
    cmplx_create_array(n, xd);
    cmplx_promote_diagonal(n, xd, x);
    ans = cmplx_sum(n, xd);

    // Test
    flag = la_trace_cmplx(n, n, x, n, &z);
    if (flag != LA_NO_ERROR) rst = false;
    if (fabs(creal(z) - creal(ans)) > DBL_TOL || 
        fabs(cimag(z) - cimag(ans)) > DBL_TOL) rst = false;

    // End
    return rst;
}


bool test_matrix_mulitply()
{
    // Variables
    const int m = 40;
    const int n = 20;
    const int p = 30;
    const int mn = m * n;
    const int mp = m * p;
    const int pn = p * n;
    const double alpha = 1.0;
    const double beta = 0.0;
    double x[mp], y[pn], z[mn], ans[mn];
    bool rst;
    int flag;

    // Initialization
    rst = true;
    create_matrix(m, p, x);
    create_matrix(p, n, y);
    mtx_mult(m, n, p, x, y, ans);

    // Test
    flag = la_mtx_mult(false, false, m, n, p, alpha, x, m, y, p, beta, z, m);
    if (flag != LA_NO_ERROR) rst = false;
    if (!is_mtx_equal(m, n, ans, z, DBL_TOL)) rst = false;

    // End
    return rst;
}

bool test_cmplx_matrix_mulitply()
{
    // Variables
    const int m = 40;
    const int n = 20;
    const int p = 30;
    const int mn = m * n;
    const int mp = m * p;
    const int pn = p * n;
    const double complex alpha = 1.0 + 0.0 * I;
    const double complex beta = 0.0 + 0.0 * I;
    double complex x[mp], y[pn], z[mn], ans[mn];
    bool rst;
    int flag;

    // Initialization
    rst = true;
    cmplx_create_matrix(m, p, x);
    cmplx_create_matrix(p, n, y);
    cmplx_mtx_mult(m, n, p, x, y, ans);

    // Test
    flag = la_mtx_mult_cmplx(LA_NO_OPERATION, LA_NO_OPERATION, m, n, p, alpha, 
        x, m, y, p, beta, z, m);
    if (flag != LA_NO_ERROR) rst = false;
    if (!is_cmplx_mtx_equal(m, n, ans, z, DBL_TOL)) rst = false;

    // End
    return rst;
}




bool test_triangular_matrix_multiply()
{
    // Variables
    const int n = 50;
    const int nn = n * n;
    const double alpha = 0.5;
    const double beta = -1.0;
    double x[nn], y[nn], y1[nn];
    bool rst;
    int flag;

    // Initialization
    rst = true;
    create_triangular_matrix(true, n, x);
    create_matrix(n, n, y);
    copy_matrix(n, n, y, y1);

    // Compute the solution (Y = alpha * X**T * X + beta * Y)
    flag = la_mtx_mult(true, false, n, n, n, alpha, x, n, x, n, beta, y1, n);
    if (flag != LA_NO_ERROR) rst = false;

    // Test
    flag = la_tri_mtx_mult(true, alpha, n, x, n, beta, y, n);
    if (flag != LA_NO_ERROR) rst = false;
    if (!is_mtx_equal(n, n, y, y1, DBL_TOL)) rst = false;

    // End
    return rst;
}

bool test_cmplx_triangular_matrix_multiply()
{
    // Variables
    const int n = 50;
    const int nn = n * n;
    const double complex alpha = 0.5 + 0.0 * I;
    const double complex beta = -1.0 + 0.0 * I;
    double complex x[nn], y[nn], y1[nn];
    bool rst;
    int flag;

    // Initialization
    rst = true;
    cmplx_create_triangular_matrix(true, n, x);
    cmplx_create_matrix(n, n, y);
    cmplx_copy_matrix(n, n, y, y1);

    // Compute the solution (Y = alpha * X**T * X + beta * Y)
    flag = la_mtx_mult_cmplx(LA_TRANSPOSE, LA_NO_OPERATION, n, n, n, alpha, x, 
        n, x, n, beta, y1, n);
    if (flag != LA_NO_ERROR) rst = false;

    // Test
    flag = la_tri_mtx_mult_cmplx(true, alpha, n, x, n, beta, y, n);
    if (flag != LA_NO_ERROR) rst = false;
    if (!is_cmplx_mtx_equal(n, n, y, y1, DBL_TOL)) rst = false;

    // End
    return rst;
}

