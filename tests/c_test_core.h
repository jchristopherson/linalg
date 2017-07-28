// c_test_core.h
#ifndef C_TEST_CORE_INCLUDED
#define C_TEST_CORE_INCLUDED

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include "linalg.h"

#define INDEX(i, j, m) ((j) * (m) + (i))

void make_rand_mtx(int m, int n, double *x);
bool is_dbl_mtx_equal(int m, int n, const double *x, const double *y, 
                      double tol);

// Misc. Routines
bool test_diagonal_mtx_mult();
bool test_rank1_update();
bool test_rank();

#endif
