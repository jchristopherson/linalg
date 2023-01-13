#include "linalg.h"
#include "c_linalg_test.h"
#include "c_test_core.h"


bool test_lu()
{
    // Variables
    const int n = 50;
    const int nrhs = 20;
    const int nn = n * n;
    const int nnrhs = n * nrhs;
    double a[nn], a1[nn], b[nnrhs], b1[nnrhs], bref[nnrhs];
    int pvt[n];
    bool rst;
    int flag;

    // Initialization
    rst = true;
    create_matrix(n, n, a);
    copy_matrix(n, n, a, a1);
    create_matrix(n, nrhs, b);
    copy_matrix(n, nrhs, b, b1);

    // Factor A
    flag = la_lu_factor(n, n, a, n, pvt);
    if (flag != LA_NO_ERROR) rst = false;

    // Solve A * X = B - store results in B
    flag = la_solve_lu(n, nrhs, a, n, pvt, b, n);
    if (flag != LA_NO_ERROR) rst = false;

    // Compare by multiplying the solutions
    mtx_mult(n, nrhs, n, a1, b, bref);
    if (!is_mtx_equal(n, nrhs, b1, bref, DBL_TOL)) rst = false;

    // End
    return rst;
}

bool test_cmplx_lu()
{
    // Variables
    const int n = 50;
    const int nrhs = 20;
    const int nn = n * n;
    const int nnrhs = n * nrhs;
    const double complex zero = 0.0 + 0.0 * I;
    const double complex one = 1.0 + 0.0 * I;
    double complex a[nn], a1[nn], b[nnrhs], b1[nnrhs], bref[nnrhs];
    int pvt[n];
    bool rst;
    int flag;

    // Initialization
    rst = true;
    cmplx_create_matrix(n, n, a);
    cmplx_copy_matrix(n, n, a, a1);
    cmplx_create_matrix(n, nrhs, b);
    cmplx_copy_matrix(n, nrhs, b, b1);

    // Factor A
    flag = la_lu_factor_cmplx(n, n, a, n, pvt);
    if (flag != LA_NO_ERROR) rst = false;

    // Solve A * X = B - store results in B
    flag = la_solve_lu_cmplx(n, nrhs, a, n, pvt, b, n);
    if (flag != LA_NO_ERROR) rst = false;

    // Compare by multiplying the solutions
    flag = la_mtx_mult_cmplx(LA_NO_OPERATION, LA_NO_OPERATION, n, nrhs, n, one, 
        a1, n, b, n, zero, bref, n);
    if (!is_cmplx_mtx_equal(n, nrhs, b1, bref, DBL_TOL)) rst = false;

    // End
    return rst;
}




bool test_qr()
{
    // Variables
    const int m = 50;
    const int n = 50;
    const int nrhs = 20;
    const int mn = m * n;
    const int mnrhs = m * nrhs;
    const int minmn = MIN(m, n);
    const double zero = 0.0;
    const double one = 1.0;
    double a[mn], a1[mn], b[mnrhs], b1[mnrhs], tau[minmn], bref[mnrhs];
    bool rst;
    int flag;

    // Initialization
    rst = true;
    create_matrix(m, n, a);
    copy_matrix(m, n, a, a1);
    create_matrix(m, nrhs, b);
    copy_matrix(m, nrhs, b, b1);

    // Factor
    flag = la_qr_factor(m, n, a, m, tau);
    if (flag != LA_NO_ERROR) rst = false;

    // Solve
    flag = la_solve_qr(m, n, nrhs, a, m, tau, b, m);
    if (flag != LA_NO_ERROR) rst = false;

    // Test by ensuring A * X = B
    flag = la_mtx_mult(false, false, m, nrhs, n, one, a1, m, b, m, zero, bref, m);
    if (flag != LA_NO_ERROR) rst = false;
    if (!is_mtx_equal(m, nrhs, b1, bref, DBL_TOL)) rst = false;

    // End
    return rst;
}

bool test_cmplx_qr()
{
    // Variables
    const int m = 50;
    const int n = 50;
    const int nrhs = 20;
    const int mn = m * n;
    const int mnrhs = m * nrhs;
    const int minmn = MIN(m, n);
    const double complex zero = 0.0 + 0.0 * I;
    const double complex one = 1.0 + 0.0 * I;
    double complex a[mn], a1[mn], b[mnrhs], b1[mnrhs], tau[minmn], bref[mnrhs];
    bool rst;
    int flag;

    // Initialization
    rst = true;
    cmplx_create_matrix(m, n, a);
    cmplx_copy_matrix(m, n, a, a1);
    cmplx_create_matrix(m, nrhs, b);
    cmplx_copy_matrix(m, nrhs, b, b1);

    // Factor
    flag = la_qr_factor_cmplx(m, n, a, m, tau);
    if (flag != LA_NO_ERROR) rst = false;

    // Solve
    flag = la_solve_qr_cmplx(m, n, nrhs, a, m, tau, b, m);
    if (flag != LA_NO_ERROR) rst = false;

    // Test by ensuring A * X = B
    flag = la_mtx_mult_cmplx(LA_NO_OPERATION, LA_NO_OPERATION, m, nrhs, n, one, 
        a1, m, b, m, zero, bref, m);
    if (flag != LA_NO_ERROR) rst = false;
    if (!is_cmplx_mtx_equal(m, nrhs, b1, bref, DBL_TOL)) rst = false;

    // End
    return rst;
}




bool test_qr_pivot()
{
    // Variables
    const int m = 50;
    const int n = 50;
    const int nrhs = 20;
    const int mn = m * n;
    const int mnrhs = m * nrhs;
    const int minmn = MIN(m, n);
    const double zero = 0.0;
    const double one = 1.0;
    double a[mn], a1[mn], b[mnrhs], b1[mnrhs], tau[minmn], bref[mnrhs];
    bool rst;
    int flag, pvt[n];

    // Initialization
    rst = true;
    create_matrix(m, n, a);
    copy_matrix(m, n, a, a1);
    create_matrix(m, nrhs, b);
    copy_matrix(m, nrhs, b, b1);
    zero_int_array(n, pvt);

    // Factor
    flag = la_qr_factor_pvt(m, n, a, m, tau, pvt);
    if (flag != LA_NO_ERROR) rst = false;

    // Solve
    flag = la_solve_qr_pvt(m, n, nrhs, a, m, tau, pvt, b, m);
    if (flag != LA_NO_ERROR) rst = false;

    // Test by ensuring A * X = B
    flag = la_mtx_mult(false, false, m, nrhs, n, one, a1, m, b, m, zero, bref, m);
    if (flag != LA_NO_ERROR) rst = false;
    if (!is_mtx_equal(m, nrhs, b1, bref, DBL_TOL)) rst = false;

    // End
    return rst;
}

bool test_cmplx_qr_pivot()
{
    // Variables
    const int m = 50;
    const int n = 50;
    const int nrhs = 20;
    const int mn = m * n;
    const int mnrhs = m * nrhs;
    const int minmn = MIN(m, n);
    const double complex zero = 0.0;
    const double complex one = 1.0;
    double complex a[mn], a1[mn], b[mnrhs], b1[mnrhs], tau[minmn], bref[mnrhs];
    bool rst;
    int flag, pvt[n];

    // Initialization
    rst = true;
    cmplx_create_matrix(m, n, a);
    cmplx_copy_matrix(m, n, a, a1);
    cmplx_create_matrix(m, nrhs, b);
    cmplx_copy_matrix(m, nrhs, b, b1);
    zero_int_array(n, pvt);

    // Factor
    flag = la_qr_factor_cmplx_pvt(m, n, a, m, tau, pvt);
    if (flag != LA_NO_ERROR) rst = false;

    // Solve
    flag = la_solve_qr_cmplx_pvt(m, n, nrhs, a, m, tau, pvt, b, m);
    if (flag != LA_NO_ERROR) rst = false;

    // Test by ensuring A * X = B
    flag = la_mtx_mult_cmplx(LA_NO_OPERATION, LA_NO_OPERATION, m, nrhs, n, one, 
        a1, m, b, m, zero, bref, m);
    if (flag != LA_NO_ERROR) rst = false;
    if (!is_cmplx_mtx_equal(m, nrhs, b1, bref, DBL_TOL)) rst = false;

    // End
    return rst;
}




bool test_qr_rank1_update()
{
    // Variables
    const int m = 50;
    const int n = 30;
    const int mn = m * n;
    const int minmn = MIN(m, n);
    const int mm = m * m;
    const double zero = 0.0;
    const double one = 1.0;
    double a[mn], u[m], v[n], tau[minmn], q[mm], a1[mn], ans[mn];
    bool rst;
    int flag;

    // Initialization
    rst = true;
    create_matrix(m, n, a);
    copy_matrix(m, n, a, a1);
    create_array(m, u);
    create_array(n, v);

    // Compute the QR factorization of A and then explicity form Q & R
    flag = la_qr_factor(m, n, a, m, tau);
    if (flag != LA_NO_ERROR) rst = false;

    flag = la_form_qr(true, m, n, a, m, tau, q, m);    // R is stored in A
    if (flag != LA_NO_ERROR) rst = false;

    // Perform the rank 1 update to the original matrix A
    flag = la_rank1_update(m, n, one, u, v, a1, m);
    if (flag != LA_NO_ERROR) rst = false;

    // Perform the rank 1 update to the QR factored matrices
    flag = la_qr_rank1_update(m, n, q, m, a, m, u, v);
    if (flag != LA_NO_ERROR) rst = false;

    // See if A1 = Q1 * R1
    flag = la_mtx_mult(false, false, m, n, m, one, q, m, a, m, zero, ans, m);
    if (!is_mtx_equal(m, n, ans, a1, DBL_TOL)) rst = false;

    // End
    return rst;
}

bool test_cmplx_qr_rank1_update()
{
    // Variables
    const int m = 50;
    const int n = 30;
    const int mn = m * n;
    const int minmn = MIN(m, n);
    const int mm = m * m;
    const double complex zero = 0.0 + 0.0 * I;
    const double complex one = 1.0 + 0.0 * I;
    double complex a[mn], u[m], v[n], tau[minmn], q[mm], a1[mn], ans[mn];
    bool rst;
    int flag;

    // Initialization
    rst = true;
    cmplx_create_matrix(m, n, a);
    cmplx_copy_matrix(m, n, a, a1);
    cmplx_create_array(m, u);
    cmplx_create_array(n, v);

    // Compute the QR factorization of A and then explicity form Q & R
    // flag = la_qr_factor_cmplx(m, n, a, m, tau);
    // if (flag != LA_NO_ERROR) rst = false;

    // flag = la_form_qr_cmplx(true, m, n, a, m, tau, q, m);    // R is stored in A
    // if (flag != LA_NO_ERROR) rst = false;

    // // Perform the rank 1 update to the original matrix A
    // flag = la_rank1_update_cmplx(m, n, one, u, v, a1, m);
    // if (flag != LA_NO_ERROR) rst = false;

    // // Perform the rank 1 update to the QR factored matrices
    // flag = la_qr_rank1_update_cmplx(m, n, q, m, a, m, u, v);
    // if (flag != LA_NO_ERROR) rst = false;

    // // See if A1 = Q1 * R1
    // flag = la_mtx_mult_cmplx(false, false, m, n, m, one, q, m, a, m, zero, ans, m);
    // if (!is_cmplx_mtx_equal(m, n, ans, a1, DBL_TOL)) rst = false;

    // End
    return rst;
}





bool test_cholesky()
{
    // Variables
    const int n = 50;
    const int nrhs = 20;
    const int nn = n * n;
    const int nnrhs = n * nrhs;
    const double zero = 0.0;
    const double one = 1.0;
    double a[nn], a1[nn], b[nnrhs], b1[nnrhs], bref[nnrhs];
    bool rst;
    int flag;

    // Initialization
    rst = true;
    create_symmetric_matrix(n, a);
    copy_matrix(n, n, a, a1);
    create_matrix(n, nrhs, b);
    copy_matrix(n, nrhs, b, b1);

    // Factor A
    flag = la_cholesky_factor(true, n, a, n);
    if (flag != LA_NO_ERROR) rst = false;

    // Solve A * X = B for X
    flag = la_solve_cholesky(true, n, nrhs, a, n, b, n);
    if (flag != LA_NO_ERROR) rst = false;

    // Check A * X = B
    flag = la_mtx_mult(false, false, n, nrhs, n, one, a1, n, b, n, zero, bref, n);
    if (flag != LA_NO_ERROR) rst = false;

    // Test
    if (!is_mtx_equal(n, nrhs, b1, bref, DBL_TOL)) rst = false;

    // End
    return rst;
}

bool test_cmplx_cholesky()
{
    // Variables
    const int n = 50;
    const int nrhs = 20;
    const int nn = n * n;
    const int nnrhs = n * nrhs;
    const double complex zero = 0.0 + 0.0 * I;
    const double complex one = 1.0 + 0.0 * I;
    double complex a[nn], a1[nn], b[nnrhs], b1[nnrhs], bref[nnrhs];
    bool rst;
    int flag;

    // Initialization
    rst = true;
    cmplx_create_symmetric_matrix(n, a);
    cmplx_copy_matrix(n, n, a, a1);
    cmplx_create_matrix(n, nrhs, b);
    cmplx_copy_matrix(n, nrhs, b, b1);

    // Factor A
    flag = la_cholesky_factor_cmplx(true, n, a, n);
    if (flag != LA_NO_ERROR) rst = false;

    // Solve A * X = B for X
    flag = la_solve_cholesky_cmplx(true, n, nrhs, a, n, b, n);
    if (flag != LA_NO_ERROR) rst = false;

    // Check A * X = B
    flag = la_mtx_mult_cmplx(LA_NO_OPERATION, LA_NO_OPERATION, n, nrhs, n, one, 
        a1, n, b, n, zero, bref, n);
    if (flag != LA_NO_ERROR) rst = false;

    // Test
    if (!is_cmplx_mtx_equal(n, nrhs, b1, bref, DBL_TOL)) rst = false;

    // End
    return rst;
}




bool test_cholesky_rank1_update()
{
    // Variables
    const int n = 50;
    const int nn = n * n;
    const double zero = 0.0;
    const double one = 1.0;
    double a[nn], a1[nn], u[n], c[nn];
    bool rst;
    int flag;

    // Initialization
    rst = true;
    create_symmetric_matrix(n, a);
    copy_matrix(n, n, a, a1);
    create_array(n, u);

    // Compute the factorization
    flag = la_cholesky_factor(true, n, a, n);
    if (flag != LA_NO_ERROR) rst = false;

    // Perform the rank 1 update to matrix A
    flag = la_rank1_update(n, n, one, u, u, a1, n);
    if (flag != LA_NO_ERROR) rst = false;

    // Update the factorization
    flag = la_cholesky_rank1_update(n, a, n, u);
    if (flag != LA_NO_ERROR) rst = false;

    // Ensure R**T * R = A
    flag = la_mtx_mult(true, false, n, n, n, one, a, n, a, n, zero, c, n);
    if (flag != LA_NO_ERROR) rst = false;
    if (!is_mtx_equal(n, n, c, a1, DBL_TOL)) rst = false;

    // End
    return rst;
}

bool test_cmplx_cholesky_rank1_update()
{
    // Variables
    const int n = 50;
    const int nn = n * n;
    const double complex zero = 0.0 + 0.0 * I;
    const double complex one = 1.0 + 0.0 * I;
    double complex a[nn], a1[nn], u[n], c[nn];
    bool rst;
    int flag;

    // Initialization
    rst = true;
    cmplx_create_symmetric_matrix(n, a);
    cmplx_copy_matrix(n, n, a, a1);
    cmplx_create_array(n, u);

    // Compute the factorization
    flag = la_cholesky_factor_cmplx(true, n, a, n);
    if (flag != LA_NO_ERROR) rst = false;

    // Perform the rank 1 update to matrix A
    flag = la_rank1_update_cmplx(n, n, one, u, u, a1, n);
    if (flag != LA_NO_ERROR) rst = false;

    // Update the factorization
    flag = la_cholesky_rank1_update_cmplx(n, a, n, u);
    if (flag != LA_NO_ERROR) rst = false;

    // Ensure R**T * R = A
    flag = la_mtx_mult_cmplx(LA_TRANSPOSE, LA_NO_OPERATION, n, n, n, 
        one, a, n, a, n, zero, c, n);
    if (flag != LA_NO_ERROR) rst = false;
    if (!is_cmplx_mtx_equal(n, n, c, a1, DBL_TOL)) rst = false;

    // End
    return rst;
}




bool test_svd()
{
    // Variables
    const int m = 60;
    const int n = 50;
    const int mn = m * n;
    const int mm = m * m;
    const int nn = n * n;
    const int minmn = MIN(m, n);
    const double zero = 0.0;
    const double one = 1.0;
    double a[mn], a1[mn], u[mm], s[minmn], vt[nn], us[mn];
    bool rst;
    int flag;

    // Initialization
    rst = true;
    create_matrix(m, n, a);
    copy_matrix(m, n, a, a1);

    // Compute the SVD
    flag = la_svd(m, n, a, m, s, u, m, vt, n);
    if (flag != LA_NO_ERROR) rst = false;

    // Compute U * S * V**T
    flag = la_diag_mtx_mult(false, false, m, n, minmn, one, s, u, m, zero, us, m);
    if (flag != LA_NO_ERROR) rst = false;

    // Compute (U * S) * V**T
    flag = la_mtx_mult(false, false, m, n, n, one, us, m, vt, n, zero, a, m);
    if (flag != LA_NO_ERROR) rst = false;

    // Test
    if (!is_mtx_equal(m, n, a, a1, DBL_TOL)) rst = false;
    return rst;
}

bool test_cmplx_svd()
{
    // Variables
    const int m = 60;
    const int n = 50;
    const int mn = m * n;
    const int mm = m * m;
    const int nn = n * n;
    const int minmn = MIN(m, n);
    const double complex zero = 0.0 + 0.0 * I;
    const double complex one = 1.0 + 0.0 * I;
    double complex a[mn], a1[mn], u[mm], vt[nn], us[mn];
    double s[minmn];
    bool rst;
    int flag;

    // Initialization
    rst = true;
    cmplx_create_matrix(m, n, a);
    cmplx_copy_matrix(m, n, a, a1);

    // Compute the SVD
    flag = la_svd_cmplx(m, n, a, m, s, u, m, vt, n);
    if (flag != LA_NO_ERROR) rst = false;

    // Compute U * S * V**T
    flag = la_diag_mtx_mult_mixed(false, LA_NO_OPERATION, m, n, minmn, one, s, 
        u, m, zero, us, m);
    if (flag != LA_NO_ERROR) rst = false;

    // Compute (U * S) * V**T
    flag = la_mtx_mult_cmplx(LA_NO_OPERATION, LA_NO_OPERATION, m, n, n, one, 
        us, m, vt, n, zero, a, m);
    if (flag != LA_NO_ERROR) rst = false;

    // Test
    if (!is_cmplx_mtx_equal(m, n, a, a1, DBL_TOL)) rst = false;
    return rst;
}





bool test_inverse()
{
    // Variables
    const int n = 50;
    const int nn = n * n;
    double a[nn], a1[nn], ainv[nn];
    bool rst;
    int flag;

    // Initialization
    rst = true;
    create_matrix(n, n, a);
    copy_matrix(n, n, a, a1);

    // Compute the inverse (traditional)
    flag = la_inverse(n, a, n);
    if (flag != LA_NO_ERROR) rst = false;

    // Compute the inverse using the Moore-Penrose approach
    flag = la_pinverse(n, n, a1, n, ainv, n);
    if (flag != LA_NO_ERROR) rst = false;

    // Test
    if (!is_mtx_equal(n, n, a, ainv, DBL_TOL)) rst = false;
    return rst;
}

bool test_cmplx_inverse()
{
    // Variables
    const int n = 50;
    const int nn = n * n;
    double complex a[nn], a1[nn], ainv[nn];
    bool rst;
    int flag;

    // Initialization
    rst = true;
    cmplx_create_matrix(n, n, a);
    cmplx_copy_matrix(n, n, a, a1);

    // Compute the inverse (traditional)
    flag = la_inverse_cmplx(n, a, n);
    if (flag != LA_NO_ERROR) rst = false;

    // Compute the inverse using the Moore-Penrose approach
    flag = la_pinverse_cmplx(n, n, a1, n, ainv, n);
    if (flag != LA_NO_ERROR) rst = false;

    // Test
    if (!is_cmplx_mtx_equal(n, n, a, ainv, DBL_TOL)) rst = false;
    return rst;
}





bool test_lq()
{
    // Variables
    const int m = 50;
    const int n = 50;
    const int nrhs = 20;
    const int mn = m * n;
    const int mnrhs = m * nrhs;
    const int nnrhs = n * nrhs;
    const int minmn = MIN(m, n);
    const double zero = 0.0;
    const double one = 1.0;
    double a[mn], a1[mn], x[nnrhs], tau[minmn], b[mnrhs], bref[mnrhs];
    bool rst;
    int i, j, flag;

    // Initialization
    rst = true;
    create_matrix(m, n, a);
    copy_matrix(m, n, a, a1);
    create_matrix(m, nrhs, bref);
    copy_matrix(m, nrhs, bref, x);

    // Factor
    flag = la_lq_factor(m, n, a, m, tau);
    if (flag != LA_NO_ERROR) rst = false;

    // Solve
    flag = la_solve_lq(m, n, nrhs, a, m, tau, x, n);
    if (flag != LA_NO_ERROR) rst = false;

    // Test by ensuring A * X = B
    flag = la_mtx_mult(false, false, m, nrhs, n, one, a1, m, x, n, zero, b, m);
    if (flag != LA_NO_ERROR) rst = false;
    if (!is_mtx_equal(m, nrhs, b, bref, DBL_TOL)) rst = false;

    // End
    return rst;
}

bool test_cmplx_lq()
{
    // Variables
    const int m = 50;
    const int n = 50;
    const int nrhs = 20;
    const int mn = m * n;
    const int mnrhs = m * nrhs;
    const int nnrhs = n * nrhs;
    const int minmn = MIN(m, n);
    const double complex zero = 0.0 + 0.0 * I;
    const double complex one = 1.0 + 1.0 * I;
    double complex a[mn], a1[mn], x[nnrhs], tau[minmn], b[mnrhs], bref[mnrhs];
    bool rst;
    int i, j, flag;

    // Initialization
    rst = true;
    cmplx_create_matrix(m, n, a);
    cmplx_copy_matrix(m, n, a, a1);
    cmplx_create_matrix(m, nrhs, bref);
    cmplx_copy_matrix(m, nrhs, bref, x);

    // Factor
    flag = la_lq_factor_cmplx(m, n, a, m, tau);
    if (flag != LA_NO_ERROR) rst = false;

    // Solve
    flag = la_solve_lq_cmplx(m, n, nrhs, a, m, tau, x, n);
    if (flag != LA_NO_ERROR) rst = false;

    // Test by ensuring A * X = B
    flag = la_mtx_mult_cmplx(LA_NO_OPERATION, LA_NO_OPERATION, m, nrhs, n, 
        one, a1, m, x, n, zero, b, m);
    if (flag != LA_NO_ERROR) rst = false;
    if (!is_cmplx_mtx_equal(m, nrhs, b, bref, DBL_TOL)) rst = false;

    // End
    return rst;
}
