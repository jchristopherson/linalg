// c_test_core.h
#ifndef C_TEST_CORE_INCLUDED
#define C_TEST_CORE_INCLUDED

#include <stdbool.h>

#define INDEX(i, j, m) ((j) * (m) + (i))

void make_rand_mtx(int m, int n, double *x);
bool is_dbl_mtx_equal(int m, int n, const double *x, const double *y, 
                      double tol);

#endif
