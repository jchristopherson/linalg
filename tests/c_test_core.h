#ifndef C_TEST_CORE_H_DEFINED
#define C_TEST_CORE_H_DEFINED

#include <stdbool.h>
#include <complex.h>

#define INDEX(i, j, lda) ((lda) * (j) + (i))

// Tests to see if two matrices are equal within the specified tolerance.
bool is_mtx_equal(int m, int n, const double *x, const double *y, double tol);
bool is_cmplx_mtx_equal(int m, int n, const double complex *x, 
    const double complex *y, double tol);

// Multiplies two matrices: Z = X * Y.
//
// x is M-by-P
// y is P-by-N
// z is M-by-N
void mtx_mult(int m, int n, int p, const double *x, const double *y,
    double *z);
void cmplx_mtx_mult(int m, int n, int p, const double complex *x,
    const double complex *y, double complex *z);

// Promotes an N-element vector to an N-by-N diagonal matrix.
void promote_diagonal(int n, const double *x, double *z);
void cmplx_promote_diagonal(int n, const double complex *x, double complex *z);

// Transposes a matrix
void transpose(int m, int n, const double *x, double *z);
void cmplx_transpose(int m, int n, const double complex *x, double complex *z);
void conj_transpose(int m, int n, const double complex *x, double complex *z);

// Create a random M-by-N matrix
void create_matrix(int m, int n, double *x);
void cmplx_create_matrix(int m, int n, double complex *x);

// Create an N-by-N symmetric matrix
void create_symmetric_matrix(int n, double *x);
void cmplx_create_symmetric_matrix(int n, double complex *x);

// Create an N-by-N identity matrix
void identity_matrix(int n, double *x);
void cmplx_identity_matrix(int n, double complex *x);

// Copy an M-by-N matrix
void copy_matrix(int m, int n, const double *src, double *dst);
void cmplx_copy_matrix(int m, int n, const double complex *src, 
    double complex *dst);

#endif