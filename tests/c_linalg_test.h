#ifndef C_LINALG_TEST_H_DEFINED
#define C_LINALG_TEST_H_DEFINED

#include <stdbool.h>
#include <complex.h>

bool test_diagonal_mtx_mult();
bool test_cmplx_diagonal_mtx_mult();
bool test_rank1_update();
bool test_cmplx_rank1_update();
bool test_trace();
bool test_cmplx_trace();

#endif
