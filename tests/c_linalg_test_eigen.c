#include "linalg.h"
#include "c_linalg_test.h"
#include "c_test_core.h"


bool test_eigen_symm()
{
    // Variables
    const int n = 50;
    const int nn = n * n;
    const double tol = 1.0e-8;
    double zero = 0.0;
    double one = 1.0;
    double a[nn], a1[nn], vals[n], t[nn], t1[nn];
    bool rst;
    int flag;

    // Initialization
    rst = true;
    create_symmetric_matrix(n, a);
    copy_matrix(n, n, a, a1);

    // Compute the eigenvalues and eigenvectors
    flag = la_eigen_symm(true, n, a, n, vals);
    if (flag != LA_NO_ERROR) rst = false;

    // Compute A1 * VECS = VECS * VALS
    flag = la_mtx_mult(false, false, n, n, n, one, a1, n, a, n, zero, t, n);
    if (flag != LA_NO_ERROR) rst = false;

    flag = la_diag_mtx_mult(false, false, n, n, n, one, vals, a, n, zero, t1, n);
    if (flag != LA_NO_ERROR) rst = false;

    // Test
    if (!is_mtx_equal(n, n, t, t1, tol)) rst = false;
    return rst;
}





bool test_eigen_asymm()
{
    // Variables
    const int n = 50;
    const int nn = n * n;
    const double tol = 1.0e-8;
    const double zero = 0.0;
    const double one = 1.0;
    double a[nn];
    double complex a1[nn], vals[n], vecs[nn], t[nn], t1[nn];
    bool rst;
    int flag;

    // Initialization
    rst = true;
    create_matrix(n, n, a);
    to_complex(nn, a, a1);

    // Compute the eigenvalues and eigenvectors
    flag = la_eigen_asymm(true, n, a, n, vals, vecs, n);
    if (flag != LA_NO_ERROR) rst = false;

    // Compute A1 * VECS = VECS * VALS
    flag = la_mtx_mult_cmplx(LA_NO_OPERATION, LA_NO_OPERATION, n, n, n, one, 
        a1, n, vecs, n, zero, t, n);
    if (flag != LA_NO_ERROR) rst = false;

    flag = la_diag_mtx_mult_cmplx(false, LA_NO_OPERATION, n, n, n, one, vals, 
        vecs, n, zero, t1, n);
    if (flag != LA_NO_ERROR) rst = false;

    // Test
    if (!is_cmplx_mtx_equal(n, n, t, t1, tol)) rst = false;
    return rst;
}

bool test_cmplx_eigen_asymm()
{
    // Variables
    const int n = 50;
    const int nn = n * n;
    const double tol = 1.0e-8;
    const double complex zero = 0.0 + 0.0 * I;
    const double complex one = 1.0 + 0.0 * I;
    double complex a[nn], a1[nn], vals[n], vecs[nn], t[nn], t1[nn];
    bool rst;
    int flag;

    // Initialization
    rst = true;
    cmplx_create_matrix(n, n, a);
    cmplx_copy_matrix(n, n, a, a1);

    // Compute the eigenvalues and eigenvectors
    flag = la_eigen_cmplx(true, n, a, n, vals, vecs, n);
    if (flag != LA_NO_ERROR) rst = false;

    // Compute A1 * VECS = VECS * VALS
    flag = la_mtx_mult_cmplx(LA_NO_OPERATION, LA_NO_OPERATION, n, n, n, one, 
        a1, n, vecs, n, zero, t, n);
    if (flag != LA_NO_ERROR) rst = false;

    flag = la_diag_mtx_mult_cmplx(false, LA_NO_OPERATION, n, n, n, one, vals, 
        vecs, n, zero, t1, n);
    if (flag != LA_NO_ERROR) rst = false;

    // Test
    if (!is_cmplx_mtx_equal(n, n, t, t1, tol)) rst = false;
    return rst;
}




bool test_eigen_gen()
{
    // Variables
    const int n = 50;
    const int nn = n * n;
    const double tol = 1.0e-8;
    const double complex zero = 0.0 + 0.0 * I;
    const double complex one = 1.0 + 0.0 * I;
    double a[nn], b[nn], beta[n];
    double complex alpha[n], vecs[nn], vals[n], ac[nn], bc[nn], av[nn], bv[nn],
        bvn[nn];
    bool rst;
    int i, flag;

    // Initialization
    rst = true;
    create_matrix(n, n, a);
    create_matrix(n, n, b);
    to_complex(nn, a, ac);
    to_complex(nn, b, bc);

    // Compute the eigenvalues and eigenvectors
    flag = la_eigen_gen(true, n, a, n, b, n, alpha, beta, vecs, n);
    if (flag != LA_NO_ERROR) rst = false;

    // Compute alpha / beta - may over or underflow
    for (i = 0; i < n; ++i) vals[i] = alpha[i] / beta[i];

    // Compute A * VECS = B * VECS * VALS
    flag = la_mtx_mult_cmplx(LA_NO_OPERATION, LA_NO_OPERATION, n, n, n, one,
        ac, n, vecs, n, zero, av, n);
    if (flag != LA_NO_ERROR) rst = false;
    
    flag = la_mtx_mult_cmplx(LA_NO_OPERATION, LA_NO_OPERATION, n, n, n, one,
        bc, n, vecs, n, zero, bv, n);
    if (flag != LA_NO_ERROR) rst = false;
    flag = la_diag_mtx_mult_cmplx(false, LA_NO_OPERATION, n, n, n, one, vals,
        bv, n, zero, bvn, n);
    if (flag != LA_NO_ERROR) rst = false;

    // Test
    if (!is_cmplx_mtx_equal(n, n, av, bvn, tol)) rst = false;
    return rst;
}
