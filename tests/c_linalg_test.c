#include <stdbool.h>
#include "c_linalg_test.h"

int main()
{
    // Variables
    bool check;
    int flag = 0;

    // Tests
    // check = test_diagonal_mtx_mult();
    // if (!check) flag = 1;
    // check = test_cmplx_diagonal_mtx_mult();
    // if (!check) flag = -1;

    // check = test_rank1_update();
    // if (!check) flag = 2;
    // check = test_cmplx_rank1_update();
    // if (!check) flag = -2;

    // check = test_trace();
    // if (!check) flag = 3;
    // check = test_cmplx_trace();
    // if (!check) flag = -3;

    // check = test_triangular_matrix_multiply();
    // if (!check) flag = 4;
    // check = test_cmplx_triangular_matrix_multiply();
    // if (!check) flag = -4;

    // check = test_lu();
    // if (!check) flag = 5;
    // check = test_cmplx_lu();
    // if (!check) flag = -5;

    // check = test_qr();
    // if (!check) flag = 6;
    // check = test_cmplx_qr();
    // if (!check) flag = -6;

    // check = test_qr_pivot();
    // if (!check) flag = 7;
    // check = test_cmplx_qr_pivot();
    // if (!check) flag = -7;

    // check = test_qr_rank1_update();
    // if (!check) flag = 8;
    // check = test_cmplx_qr_rank1_update();
    // if (!check) flag = -8;

    // check = test_cholesky();
    // if (!check) flag = 9;
    // check = test_cmplx_cholesky();
    // if (!check) flag = -9;

    check = test_cholesky_rank1_update();
    if (!check) flag = 10;
    // check = test_cmplx_cholesky_rank1_update();
    // if (!check) flag = -10;
    
    check = test_svd();
    if (!check) flag = 11;
    check = test_cmplx_svd();
    if (!check) flag = -11;

    check = test_inverse();
    if (!check) flag = 12;
    check = test_cmplx_inverse();
    if (!check) flag = -12;

    check = test_eigen_symm();
    if (!check) flag = 13;

    check = test_eigen_asymm();
    if (!check) flag = 14;
    check = test_cmplx_eigen_asymm();
    if (!check) flag = -14;

    check = test_eigen_gen();
    if (!check) flag = 15;

    check = test_lq();
    if (!check) flag = 16;

    // check = test_cmplx_lq();
    // if (!check) flag = 17;

    // End
    return flag;
}
