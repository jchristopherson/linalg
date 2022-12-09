#ifndef C_LINALG_TEST_H_DEFINED
#define C_LINALG_TEST_H_DEFINED

#include <stdbool.h>
#include <complex.h>

// c_linalg_test_misc.c
bool test_diagonal_mtx_mult();
bool test_cmplx_diagonal_mtx_mult();
bool test_rank1_update();
bool test_cmplx_rank1_update();
bool test_trace();
bool test_cmplx_trace();
bool test_matrix_mulitply();
bool test_cmplx_matrix_mulitply();
bool test_triangular_matrix_multiply();
bool test_cmplx_triangular_matrix_multiply();

// c_linalg_test_factor.c
bool test_lu();
bool test_cmplx_lu();
bool test_qr();
bool test_cmplx_qr();
bool test_qr_pivot();
bool test_cmplx_qr_pivot();
bool test_qr_rank1_update();
bool test_cmplx_qr_rank1_update();
bool test_cholesky();

#endif
