// c_test_core.h
#ifndef C_TEST_CORE_INCLUDED
#define C_TEST_CORE_INCLUDED

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include "linalg.h"

#define INDEX(i, j, m) ((j) * (m) + (i))
#define MIN(x, y) ((x) < (y) ? (x) : (y))

void make_rand_mtx(int m, int n, double *x);
bool is_dbl_mtx_equal(int m, int n, const double *x, const double *y, 
                      double tol);

// Misc. Routines
bool test_diagonal_mtx_mult();
bool test_rank1_update();
bool test_rank();
bool test_tri_mtx_mult();

// LU Factorization Routines
bool test_lu_factor();
bool test_lu_solve();

// QR Factorization Routines
bool test_qr_factor();
bool test_qr_factor_od();
bool test_qr_factor_ud();
bool test_qr_mult();
bool test_qr_mult_od();
bool test_qr_mult_ud();
bool test_qr_solve_no_pivot();
bool test_qr_solve_pivot();
bool test_qr_solve_pivot_ud();

#endif
