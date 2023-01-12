#ifndef C_TEST_CORE_H_DEFINED
#define C_TEST_CORE_H_DEFINED

#include <stdbool.h>
#include <complex.h>

#define INDEX(i, j, lda) ((lda) * (j) + (i))
#define MIN(a, b)((a) < (b) ? (a) : (b))
#define DBL_TOL 1.0e-6

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

// Create an N-element array.
void create_array(int n, double *x);
void cmplx_create_array(int n, double complex *x);

// Copy an N-element array.
void copy_array(int n, const double *src, double *dst);
void cmplx_copy_array(int n, const double complex *src, double complex *dst);

// Adds two M-by-N matrices: Z = alpha * A + beta * C
void add_matrix(int m, int n, double alpha, const double *x, double beta,
    const double *y, double *z);
void cmplx_add_matrix(int m, int n, double complex alpha, 
    const double complex *x, double complex beta, const double complex *y,
    double complex *z);

// Scales a matrix by a scalar.
void scale_matrix(int m, int n, double x, double *y);
void cmplx_scale_matrix(int m, int n, double complex x, double complex *y);

// Performs a rank 1 update X * Y**T where X is an M-element array and 
// Y is an N-element array.  The complex-valued case uses a conjugate
// transpose.
void rank1_update(int m, int n, const double *x, const double *y, double *z);
void cmplx_rank1_update(int m, int n, const double complex *x, 
    const double complex *y, double complex *z);

// Computes the product of the elements in an array.
double product(int n, const double *x);
double complex cmplx_product(int n, const double complex * x);

// Computes the sum of the elements in an array.
double sum(int n, const double *x);
double complex cmplx_sum(int n, const double complex *x);

// Creates a triangular matrix.
void create_triangular_matrix(bool upper, int n, double *x);
void cmplx_create_triangular_matrix(bool upper, int n, double complex *x);

// Zeros an array
void zero_int_array(int n, int *x);

// Convert from double to complex
void to_complex(int n, const double *src, double complex *dst);

#endif