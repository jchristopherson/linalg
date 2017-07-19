// linalg.h
#ifndef LINALG_H_DEFINED
#define LINALG_H_DEFINED

#include <stdbool.h>
#include <complex.h>
#include "ferror/include/ferror.h"

#ifdef __cplusplus
extern "C" {
#endif

/** @brief Performs the matrix operation: 
 *  C = alpha * op(A) * op(B) + beta * C.
 *
 * @param transa Set to true if op(A) == A**T; else, set to false if
 *  op(A) == A.
 * @param transa Set to true if op(B) == B**T; else, set to false if
 *  op(B) == B.
 * @param m The number of rows in matrix C, and the number of rows in
 *  matrix op(A).
 * @param n The number of columns in matrix C, and the number of columns
 *  in matrix op(B).
 * @param k The number of columns in matrix op(A), and the number of
 *  rows in the matrix op(B).
 * @param alpha The scalar multiplier to matrix A.
 * @param a The M-by-K matrix A.
 * @param lda The leading dimension of matrix A.  If @p transa is true, 
 *  this value must be at least MAX(1, K); else, if @p transa is false, this
 *  value must be at least MAX(1, M).
 * @param b The K-by-N matrix B.
 * @param ldb The leading dimension of matrix B.  If @p transb is true,
 *  this value must be at least MAX(1, N); else, if @p transb is false, this
 *  value must be at least MAX(1, K).
 * @param beta The scalar multiplier to matrix C.
 * @param c The M-by-N matrix C.
 */
void mtx_mult(bool transa, bool transb, int m, int n, int k, double alpha, 
              const double *a, int lda, const double *b, double beta, int ldb,
              double *c);

/** @brief Computes the matrix operation: C = alpha * A * op(B) + beta * C.
 *
 * @param trans Set to true if op(B) == B**T; else, set to false if
 *  op(B) == B.
 * @param m The number of rows in matrix C.
 * @param n The number of columns in matrix C.
 * @param alpha The scalar multiplier to matrix A.
 * @param na The length of @p a.
 * @param a A MIN(M,P)-element array containing the diagonal elements 
 *  of matrix A.
 * @param mb The number of rows in matrix B.
 * @param nb The number of columns in matrix B.
 * @param b The LDB-by-TDB matrix B where (LDB = leading dimension of B,
 *  and TDB = trailing dimension of B):
 *  - @p trans == true: LDB = N, TDB = P
 *  - @p trans == false: LDB = P, TDB = N
 * @param beta The scalar multiplier to matrix C.
 * @param c THe M-by-N matrix C.
 */
void diag_mtx_mult(bool trans, int m, int n, double alpha, int na, 
                   const double *a, int mb, int nb, const double *b, 
                   double beta, double *c);

/** @brief Computes the matrix operation: C = alpha * A * op(B) + beta * C,
 *  where A and C are complex-valued.
 *
 * @param trans Set to true if op(B) == B**T; else, set to false if
 *  op(B) == B.
 * @param m The number of rows in matrix C.
 * @param n The number of columns in matrix C.
 * @param alpha The scalar multiplier to matrix A.
 * @param na The length of @p a.
 * @param a A MIN(M,P)-element array containing the diagonal elements 
 *  of matrix A.
 * @param mb The number of rows in matrix B.
 * @param nb The number of columns in matrix B.
 * @param b The LDB-by-TDB matrix B where (LDB = leading dimension of B,
 *  and TDB = trailing dimension of B):
 *  - @p trans == true: LDB = N, TDB = P
 *  - @p trans == false: LDB = P, TDB = N
 * @param beta The scalar multiplier to matrix C.
 * @param c THe M-by-N matrix C.
 */
void diag_mtx_mult_cmplx(bool trans, int m, int n, double alpha, int na,
                         const double complex *a, int mb, int nb, 
                         const double *b, double beta, double complex *c);

/** @brief Performs the rank-1 update to matrix A such that:
 * A = alpha * X * Y**T + A, where A is an M-by-N matrix, alpha is a scalar,
 * X is an M-element array, and N is an N-element array.
 *
 * @param m The number of elements in @p x, and the number of rows in
 *  matrix @p a.
 * @param n The number of elements in @p y, and the number of columns in
 *  matrix @p a.
 * @param alpha The scalar multiplier.
 * @param x An M-element array.
 * @param y An N-element array.
 * @param a On input, the M-by-N matrix to update.  On output, the
 *  updated M-by-N matrix.
 */
void rank1_update(int m, int n, double alpha, const double *x, const double *y,
                  double *a);

/** @brief Computes the trace of a matrix (the sum of the main diagonal
 * elements).
 *
 * @param m The number of rows in the matrix.
 * @param n The number of columns in the matrix.
 * @param x The matrix on which to operate.
 *
 * @return The trace of @p x.
 */
double trace(int m, int n, const double *x);

/** @brief Computes the rank of a matrix.
 *
 * @param m The number of rows in the matrix.
 * @param n The number of columns in the matrix.
 * @param a On input, the M-by-N matrix of interest.  On output, the
 *  contents of the matrix are overwritten.
 * @param err A pointer to the C error handler object.  If no error
 *  handling is desired, simply pass NULL, and errors will be dealt with
 *  by the default internal error handler.  Possible errors that may be
 *  encountered are as follows.
 *  - LA_ARRAY_SIZE_ERROR: Occurs if any of the input arrays are not sized
 *      appropriately.
 *  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
 *      there is insufficient memory available.
 *  - LA_CONVERGENCE_ERROR: Occurs as a warning if the QR iteration process
 *      could not converge to a zero value.
 */
int mtx_rank(int m, int n, double *a, errorhandler err);

/** @brief Computes the determinant of a square matrix.
 *
 * @param n The dimension of the matrix.
 * @param a On input, the N-by-N matrix on which to operate.  On
 * output the contents are overwritten by the LU factorization of the
 * original matrix.
 * @param err A pointer to the C error handler object.  If no error
 *  handling is desired, simply pass NULL, and errors will be dealt with
 *  by the default internal error handler.  Possible errors that may be
 *  encountered are as follows.
 *  - LA_ARRAY_SIZE_ERROR: Occurs if the input matrix is not square.
 *  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
 *      there is insufficient memory available.
 */
double det(int n, double *a, errorhandler err);

/** @brief Swaps the contents of two arrays.
 *
 * @param n The number of elements either array.
 * @param x One of the N-element arrays.
 * @param y The other N-element array.
 */
void swap(int n, double *x, double *y);

/** brief Computes the triangular matrix operation: 
 * B = alpha * A**T * A + beta * B, or B = alpha * A * A**T + beta * B, 
 * where A is a triangular matrix.
 *
 * @param upper Set to true if matrix A is upper triangular, and 
 *  B = alpha * A**T * A + beta * B is to be calculated; else, set to false
 *  if A is lower triangular, and B = alpha * A * A**T + beta * B is to
 *  be computed.
 * @param n The size of the matrix.
 * @param alpha A scalar multiplier.
 * @param a The N-by-N triangular matrix.  Notice, if @p upper is true
 *  only the upper triangular portion of this matrix is referenced; else,
 *  if @p upper is false, only the lower triangular portion of this matrix
 *  is referenced.
 * @param beta A scalar multiplier.
 * @param b On input, the N-by-N matrix B.  On output, the N-by-N
 *  solution matrix.
 * @param err A pointer to the C error handler object.  If no error
 *  handling is desired, simply pass NULL, and errors will be dealt with
 *  by the default internal error handler.  Possible errors that may be
 *  encountered are as follows.
 *  - LA_ARRAY_SIZE_ERROR: Occurs if any of the input arrays are not sized
 *      appropriately.
 */
void tri_mtx_mult(bool upper, int n, double alpha, const double *a, double beta,
                  double *b, errorhandler err);

/** @brief Computes the LU factorization of an M-by-N matrix.
 *
 * @param m The number of rows in the matrix.
 * @param n The number of columns in the matrix.
 * @param a On input, the M-by-N matrix on which to operate.  On
 * output, the LU factored matrix in the form [L\\U] where the unit diagonal
 * elements of L are not stored.
 * @param ni The number of elements in the pivot array @p ipvt.  This
 *  value must be equal to MIN(M, N).
 * @param ipvt An MIN(M, N)-element array used to track row-pivot
 *  operations.  The array stored pivot information such that row I is
 *  interchanged with row IPVT(I).
 * @param err A pointer to the C error handler object.  If no error
 *  handling is desired, simply pass NULL, and errors will be dealt with
 *  by the default internal error handler.  Possible errors that may be
 *  encountered are as follows.
 *  - LA_ARRAY_SIZE_ERROR: Occurs if the pivot array is not sized 
 *      appropriately.
 *  - LA_SINGULAR_MATRIX_ERROR: Occurs as a warning if @p a is found to be
 *      singular.
 */
void lu_factor(int m, int n, double *a, int ni, int *ipvt, errorhandler err);

/** @brief Extracts the L, U, and P matrices from the output of the
 * @ref lu_factor routine.
 *
 * @param n The dimension of the original matrix.
 * @param lu On input, the N-by-N matrix as output by
 *  @ref lu_factor.  On output, the N-by-N lower triangular matrix L.
 * @param ipvt The N-element pivot array as output by
 *  @ref lu_factor.
 * @param u An N-by-N matrix where the U matrix will be written.
 * @param p An N-by-N matrix where the row permutation matrix will be
 *  written.
 *
 * @par Remarks
 * This routine allows extraction of the actual "L", "U", and "P" matrices
 * of the decomposition.  To use these matrices to solve the system A*X = B,
 * the following approach is used.
 *
 * 1. First, solve the linear system: L*Y = P*B for Y.
 * 2. Second, solve the linear system: U*X = Y for X.
 *
 * Notice, as both L and U are triangular in structure, the above equations
 * can be solved by forward and backward substitution.
 */
void form_lu(int n, double *lu, const int *ipvt, double *u, double *p);

/** @brief Computes the QR factorization of an M-by-N matrix without
 * pivoting.
 *
 * @param m The number of rows in the matrix.
 * @param n The number of columns in the matrix.
 * @param a On input, the M-by-N matrix to factor.  On output, the
 *  elements on and above the diagonal contain the MIN(M, N)-by-N upper
 *  trapezoidal matrix R (R is upper triangular if M >= N).  The elements
 *  below the diagonal, along with the array @p tau, represent the
 *  orthogonal matrix Q as a product of elementary reflectors.
 * @param nt The number of elements in the scalar factor array @p tau.
 *  This value must be equal to MIN(M, N).
 * @param tau A MIN(M, N)-element array used to store the scalar
 *  factors of the elementary reflectors.
 * @param err A pointer to the C error handler object.  If no error
 *  handling is desired, simply pass NULL, and errors will be dealt with
 *  by the default internal error handler.  Possible errors that may be
 *  encountered are as follows.
 *  - LA_ARRAY_SIZE_ERROR: Occurs if the scalar factor array is not sized 
 *      appropriately.
 *  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
 *      there is insufficient memory available.
 */
void qr_factor(int m, int n, double *a, int nt, double *tau, errorhandler err);

/** @brief Computes the QR factorization of an M-by-N matrix with column
 * pivoting such that A * P = Q * R.
 *
 * @param m The number of rows in the matrix.
 * @param n The number of columns in the matrix.
 * @param a On input, the M-by-N matrix to factor.  On output, the
 *  elements on and above the diagonal contain the MIN(M, N)-by-N upper
 *  trapezoidal matrix R (R is upper triangular if M >= N).  The elements
 *  below the diagonal, along with the array @p tau, represent the
 *  orthogonal matrix Q as a product of elementary reflectors.
 * @param nt The number of elements in the scalar factor array @p tau.
 *  This value must be equal to MIN(M, N).
 * @param tau A MIN(M, N)-element array used to store the scalar
 *  factors of the elementary reflectors.
 * @param jpvt On input, an N-element array that if JPVT(I) .ne. 0,
 *  the I-th column of A is permuted to the front of A * P; if JPVT(I) = 0,
 *  the I-th column of A is a free column.  On output, if JPVT(I) = K, then
 *  the I-th column of A * P was the K-th column of A.
 * @param err A pointer to the C error handler object.  If no error
 *  handling is desired, simply pass NULL, and errors will be dealt with
 *  by the default internal error handler.  Possible errors that may be
 *  encountered are as follows.
 *  - LA_ARRAY_SIZE_ERROR: Occurs if the scalar factor array is not sized 
 *      appropriately.
 *  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
 *      there is insufficient memory available.
 */
void qr_factor_pivot(int m, int n, double *a, int nt, double *tau, int *jpvt,
                     errorhandler err);

/** @brief Forms the full M-by-M orthogonal matrix Q from the elementary
 * reflectors returned by the base QR factorization algorithm.
 *
 * @param m The number of rows in the original matrix.
 * @param n The number of columns in the original matrix.
 * @param r On input, an M-by-N matrix where the elements below the
 *  diagonal contain the elementary reflectors generated from the QR
 *  factorization.  On and above the diagonal, the matrix contains the
 *  matrix R.  On output, the elements below the diagonal are zeroed such
 *  that the remaining matrix is simply the M-by-N matrix R.
 * @param nt The number of elements in the scalar factor array @p tau.
 *  This value must be equal to MIN(M, N).
 * @param tau A MIN(M, N)-element array containing the scalar factors of
 *  each elementary reflector defined in @p r.
 * @param q An M-by-M matrix where the full orthogonal matrix Q will be
 *  written.  In the event that M > N, Q may be supplied as M-by-N, and
 *  therefore only return the useful submatrix Q1 (Q = [Q1, Q2]) as the
 *  factorization can be written as Q * R = [Q1, Q2] * [R1; 0].
 * @param err A pointer to the C error handler object.  If no error
 *  handling is desired, simply pass NULL, and errors will be dealt with
 *  by the default internal error handler.  Possible errors that may be
 *  encountered are as follows.
 *  - LA_ARRAY_SIZE_ERROR: Occurs if the scalar factor array is not sized 
 *      appropriately.
 *  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
 *      there is insufficient memory available.
 */
void form_qr(int m, int n, double *r, int nt, const double *tau, double *q,
             errorhandler err);

/** @brief Forms the full M-by-M orthogonal matrix Q from the elementary
 * reflectors returned by the base QR factorization algorithm.
 *
 * @param m The number of rows in the original matrix.
 * @param n The number of columns in the original matrix.
 * @param r On input, an M-by-N matrix where the elements below the
 *  diagonal contain the elementary reflectors generated from the QR
 *  factorization.  On and above the diagonal, the matrix contains the
 *  matrix R.  On output, the elements below the diagonal are zeroed such
 *  that the remaining matrix is simply the M-by-N matrix R.
 * @param nt The number of elements in the scalar factor array @p tau.
 *  This value must be equal to MIN(M, N).
 * @param tau A MIN(M, N)-element array containing the scalar factors of
 *  each elementary reflector defined in @p r.
 * @param pvt An N-element column pivot array as returned by the QR
 *  factorization.
 * @param q An M-by-M matrix where the full orthogonal matrix Q will be
 *  written.  In the event that M > N, Q may be supplied as M-by-N, and
 *  therefore only return the useful submatrix Q1 (Q = [Q1, Q2]) as the
 *  factorization can be written as Q * R = [Q1, Q2] * [R1; 0].
 * @param p An N-by-N matrix where the pivot matrix will be written.
 * @param err A pointer to the C error handler object.  If no error
 *  handling is desired, simply pass NULL, and errors will be dealt with
 *  by the default internal error handler.  Possible errors that may be
 *  encountered are as follows.
 *  - LA_ARRAY_SIZE_ERROR: Occurs if the scalar factor array is not sized 
 *      appropriately.
 *  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
 *      there is insufficient memory available.
 */
void form_qr_pivot(int m, int n, double *r, int nt, const double *tau,
                   const int *pvt, double *q, double *p, errorhandler err);

/** @brief Multiplies a general matrix by the orthogonal matrix Q from a QR
 * factorization such that: C = op(Q) * C.
 *
 * @param trans Set to true to apply Q**T; else, set to false.
 * @param m The number of rows in the matrix @p c.
 * @param n The number of columns in the matrix @p c.
 * @param q On input, an M-by-M matrix containing the elementary
 *  reflectors output from the QR factorization.    Notice, the contents of 
 *  this matrix are restored on exit.
 *  that the remaining matrix is simply the M-by-N matrix R.
 * @param nt The number of elements in the scalar factor array @p tau.
 *  This value must be equal to MIN(M, N).
 * @param tau A MIN(M,N)-element array containing the scalar factors of 
 *  each elementary reflector defined in @p a.
 * @param c On input, the M-by-N matrix C.  On output, the product
 *  of the orthogonal matrix Q and the original matrix C.
 * @param err A pointer to the C error handler object.  If no error
 *  handling is desired, simply pass NULL, and errors will be dealt with
 *  by the default internal error handler.  Possible errors that may be
 *  encountered are as follows.
 *  - LA_ARRAY_SIZE_ERROR: Occurs if the scalar factor array is not sized 
 *      appropriately.
 *  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
 *      there is insufficient memory available.
 */
void mult_qr(bool trans, int m, int n, double *q, int nt, const double *tau,
             double *c, errorhandler err);

/** @brief Computes the rank 1 update to an M-by-N QR factored matrix A
 * (M >= N) where A = Q * R, and A1 = A + U * V**T such that A1 = Q1 * R1.
 *
 * @param m The number of rows in the original matrix.
 * @param n The number of columns in the original matrix.
 * @param q On input, the original M-by-M orthogonal matrix Q.  On
 *  output, the updated matrix Q1.
 * @param r On input, the M-by-N matrix R.  On output, the updated
 *  matrix R1.
 * @param u On input, the M-element U update vector.  On output,
 *  the original content of the array is overwritten.
 * @param v On input, the N-element V update vector.  On output,
 *  the original content of the array is overwritten.
 * @param err A pointer to the C error handler object.  If no error
 *  handling is desired, simply pass NULL, and errors will be dealt with
 *  by the default internal error handler.  Possible errors that may be
 *  encountered are as follows.
 *  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
 *      there is insufficient memory available.
 */
void qr_rank1_update(int m, int n, double *q, double *r, double *u, double *v, 
                     errorhandler err);

/** @brief Computes the Cholesky factorization of a symmetric, positive
 * definite matrix.
 *
 * @param n The dimension of the matrix.
 * @param a On input, the N-by-N matrix to factor.  On output, the
 *  factored matrix is returned in either the upper or lower triangular
 *  portion of the matrix, dependent upon the value of @p upper.
 * @param upper An optional input that, if specified, provides control
 *  over whether the factorization is computed as A = U**T * U (set to
 *  true), or as A = L * L**T (set to false).  The default value is true
 *  such that A = U**T * U.
 * @param err A pointer to the C error handler object.  If no error
 *  handling is desired, simply pass NULL, and errors will be dealt with
 *  by the default internal error handler.  Possible errors that may be
 *  encountered are as follows.
 *  - LA_MATRIX_FORMAT_ERROR: Occurs if @p a is not positive definite.
 */
void cholesky_factor(int n, double *a, bool upper, errorhandler err);

/** @brief Computes the rank 1 update to a Cholesky factored matrix (upper
 * triangular).
 *
 * @param n The dimension of the matrix.
 * @param r On input, the N-by-N upper triangular matrix R.  On
 *  output, the updated matrix R1.
 * @param u On input, the N-element update vector U.  On output,
 *  the rotation sines used to transform R to R1.
 * @param err A pointer to the C error handler object.  If no error
 *  handling is desired, simply pass NULL, and errors will be dealt with
 *  by the default internal error handler.  Possible errors that may be
 *  encountered are as follows.
 *  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
 *      there is insufficient memory available.
 */
void cholesky_rank1_update(int n, double *r, double *u, errorhandler err);

/** @brief Computes the rank 1 downdate to a Cholesky factored matrix (upper
 * triangular).
 *
 * @param n The dimension of the matrix.
 * @param r On input, the N-by-N upper triangular matrix R.  On
 *  output, the updated matrix R1.
 * @param u On input, the N-element update vector U.  On output,
 *  the rotation sines used to transform R to R1.
 * @param err A pointer to the C error handler object.  If no error
 *  handling is desired, simply pass NULL, and errors will be dealt with
 *  by the default internal error handler.  Possible errors that may be
 *  encountered are as follows.
 *  - LA_ARRAY_SIZE_ERROR: Occurs if any of the input array sizes are 
 *      incorrect.
 *  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
 *      there is insufficient memory available.
 *  - LA_MATRIX_FORMAT_ERROR: Occurs if the downdated matrix is not 
 *      positive definite.
 *  - LA_SINGULAR_MATRIX_ERROR: Occurs if @p r is singular.
 */
void cholesky_rank1_downdate(int n, double *r, double *u, errorhandler err);

/** @brief Factors an upper trapezoidal matrix by means of orthogonal
 * transformations such that A = R * Z = (R 0) * Z.  Z is an orthogonal
 * matrix of dimension N-by-N, and R is an M-by-M upper triangular
 * matrix.
 *
 * @param m The number of rows in the original matrix.
 * @param n The number of columns in the original matrix.
 * @param a On input, the M-by-N upper trapezoidal matrix to factor.
 *  On output, the leading M-by-M upper triangular part of the matrix
 *  contains the upper triangular matrix R, and elements N-L+1 to N of the
 *  first M rows of A, with the array @p tau, represent the orthogonal
 *  matrix Z as a product of M elementary reflectors.
 * @param tau An M-element array used to store the scalar
 *  factors of the elementary reflectors.
 * @param err A pointer to the C error handler object.  If no error
 *  handling is desired, simply pass NULL, and errors will be dealt with
 *  by the default internal error handler.  Possible errors that may be
 *  encountered are as follows.
 *  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
 *      there is insufficient memory available.
 */
void rz_factor(int m, int n, double *a, double *tau, errorhandler err);

/** @brief Multiplies a general matrix by the orthogonal matrix Z from an
 * RZ factorization such that: C = op(Z) * C.
 *
 * @param trans Set to true to apply Z**T; else, set to false.
 * @param m The number of rows in the matrix @p c.
 * @param n The number of columns in the matrix @p c.
 * @param l The number of columns in matrix @p a containing the
 *  meaningful part of the Householder vectors (M >= L >= 0).
 * @param a On input, the M-by-M matrix Z as output by @p rz_factor.
 *  The matrix is used as in-place storage during execution; however, the
 *  contents of the matrix are restored on exit.
 * @param tau An M-element array containing the scalar factors of the
 *  elementary reflectors found in @p a.
 * @param c On input, the M-by-N matrix C.  On output, the product
 *  of the orthogonal matrix Z and the original matrix C.
 * @param err A pointer to the C error handler object.  If no error
 *  handling is desired, simply pass NULL, and errors will be dealt with
 *  by the default internal error handler.  Possible errors that may be
 *  encountered are as follows.
 *  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
 *      there is insufficient memory available.
 */
void mult_rz(bool trans, int m, int n, int l, double *a, const double *tau, 
             double *c, errorhandler err);

/** @brief Computes the singular value decomposition of a matrix A.  The
 *  SVD is defined as: A = U * S * V**T, where U is an M-by-M orthogonal
 *  matrix, S is an M-by-N diagonal matrix, and V is an N-by-N orthogonal
 *  matrix.
 *
 * @param m The number of rows in the original matrix.
 * @param n The number of columns in the original matrix.
 * @param a On input, the M-by-N matrix to factor.  The matrix is
 *  overwritten on output.
 *  that the remaining matrix is simply the M-by-N matrix R.
 * @param ns The number of elements in the singular value array @p s.
 *  This value must be equal to MIN(M, N).
 * @param s A MIN(M, N)-element array containing the singular values
 *  of @p a sorted in descending order.
 * @param u An M-by-M matrix that on output contains the left singular
 *  vectors (matrix U in the decomposition: A = U * S * V**T)
 * @param vt An N-by-N matrix that on output contains the right
 *  singular vectors (matrix V**T in the decomposition: A = U * S * V**T).
 * @param err A pointer to the C error handler object.  If no error
 *  handling is desired, simply pass NULL, and errors will be dealt with
 *  by the default internal error handler.  Possible errors that may be
 *  encountered are as follows.
 *  - LA_ARRAY_SIZE_ERROR: Occurs if the singular value array is not sized 
 *      appropriately.
 *  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
 *      there is insufficient memory available.
 *  - LA_CONVERGENCE_ERROR: Occurs as a warning if the QR iteration process
 *      could not converge to a zero value.
 */
void svd(int m, int n, double *a, int ns, double *s, double *u, double *vt,
         errorhandler err);

/** @brief Solves one of the matrix equations: op(A) * X = alpha * B, where 
 * A is a triangular matrix.
 *
 * @param upper Set to true if A is an upper triangular matrix; else,
 *  set to false if A is a lower triangular matrix.
 * @param trans Set to true if op(A) = A**T; else, set to false if
 *  op(A) = A.
 * @param nounit Set to true if A is not a unit-diagonal matrix (ones on
 *  every diagonal element); else, set to false if A is a unit-diagonal
 *  matrix.
 * @param n The dimension of the triangular matrix @p a.
 * @param nrhs The number of right-hand-side vectors (number of columns
 *  in matrix @p b).
 * @param alpha The scalar multiplier to B.
 * @param a N-by-N triangular matrix on which to operate.
 * @param b On input, the N-by-NRHS right-hand-side.  On output, the
 *  N-by-NRHS solution.
 */
void solve_triangular_system(bool upper, bool trans, bool nounit, int n, 
                             int nrhs, double alpha, const double *a, 
                             double *b);

/** @brief Solves a system of LU-factored equations.
 *
 * @param n The dimension of the original matrix @p a.
 * @param nrhs The number of right-hand-side vectors (number of columns
 *  in matrix @p b).
 * @param a The N-by-N LU factored matrix as output by @ref lu_factor.
 * @param ipvt The N-element pivot array as output by @ref lu_factor.
 * @param b On input, the N-by-NRHS right-hand-side matrix.  On
 *  output, the N-by-NRHS solution matrix.
 */
void solve_lu(int n, int nrhs, const double *a, const int *ipvt, double *b);

/** @brief Solves a system of M QR-factored equations of N unknowns where
 * M >= N.
 *
 * @param m The number of rows in the original coefficient matrix A.
 * @param n The number of columns in the original coefficient matrix A.
 * @param nrhs The number of right-hand-side vectors (number of columns
 *  in matrix @p b).
 * @param a On input, the M-by-N QR factored matrix as returned by
 *  @ref qr_factor.  On output, the contents of this matrix are restored.
 *  Notice, M must be greater than or equal to N.
 * @param tau A MIN(M, N)-element array containing the scalar factors of
 *  the elementary reflectors as returned by @ref qr_factor.
 * @param b On input, the M-by-NRHS right-hand-side matrix.  On output,
 *  the first N columns are overwritten by the solution matrix X.
 * @param err A pointer to the C error handler object.  If no error
 *  handling is desired, simply pass NULL, and errors will be dealt with
 *  by the default internal error handler.  Possible errors that may be
 *  encountered are as follows.
 *  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
 *      there is insufficient memory available.
 */
void solve_qr(int m, int n, int nrhs, double *a, const double *tau, double *b,
              errorhandler err);

/** @brief Solves a system of M QR-factored equations of N unknowns where the
 * QR factorization made use of column pivoting.
 *
 * @param m The number of rows in the original coefficient matrix A.
 * @param n The number of columns in the original coefficient matrix A.
 * @param nrhs The number of right-hand-side vectors (number of columns
 *  in matrix @p b).
 * @param a On input, the M-by-N QR factored matrix as returned by
 *  @ref qr_factor.  On output, the contents of this matrix are altered.
 * @param nt The number of elements in the scalar factor array @p tau.
 *  This value must be equal to MIN(M, N).
 * @param tau A MIN(M, N)-element array containing the scalar factors of
 *  the elementary reflectors as returned by @ref qr_factor.
 * @param jpvt An N-element array, as output by @ref qr_factor, used to
 *  track the column pivots.
 * @param mb The number of rows in the matrix @p b.  This value must be
 *  equal to MAX(M, N).
 * @param b On input, the MAX(M, N)-by-NRHS matrix where the first M
 *  rows contain the right-hand-side matrix B.  On output, the first N rows
 *  are overwritten by the solution matrix X.
 * @param err A pointer to the C error handler object.  If no error
 *  handling is desired, simply pass NULL, and errors will be dealt with
 *  by the default internal error handler.  Possible errors that may be
 *  encountered are as follows.
 *  - LA_ARRAY_SIZE_ERROR: Occurs if any of the input arrays are not sized 
 *      appropriately.
 *  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
 *      there is insufficient memory available.
 */
void solve_qr_pivot(int m, int n, int nrhs, double *a, int nt, 
                    const double *tau, const int *jpvt, int mb, double *b,
                    errorhandler err);

/** @brief Solves a system of Cholesky factored equations.
 *
 * @param upper Set to true if the original matrix A was factored such
 *  that A = U**T * U; else, set to false if the factorization of A was
 *  A = L**T * L.
 * @param n The dimension of the original matrix @p a.
 * @param nrhs The number of right-hand-side vectors (number of columns
 *  in matrix @p b).
 * @param a The N-by-N Cholesky factored matrix.
 * @param b On input, the N-by-NRHS right-hand-side matrix B.  On
 *  output, the solution matrix X.
 */
void solve_cholesky(bool upper, int n, int nrhs, const double *a, double *b);

/** @brief Computes the inverse of a square matrix.
 *
 * @param n The dimension of the matrix.
 * @param a On input, the N-by-N matrix to invert.  On output, the
 *  inverted matrix.
 * @param err A pointer to the C error handler object.  If no error
 *  handling is desired, simply pass NULL, and errors will be dealt with
 *  by the default internal error handler.  Possible errors that may be
 *  encountered are as follows.
 *  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
 *      there is insufficient memory available.
 *  - LA_SINGULAR_MATRIX_ERROR: Occurs if the input matrix is singular.
 */
void mtx_inverse(int n, double *a, errorhandler err);

/** @brief Computes the Moore-Penrose pseudo-inverse of a M-by-N matrix
 * using the singular value decomposition of the matrix.
 *
 * @param m The number of rows in the matrix to invert.
 * @param n The number of columns in the matrix to invert.
 * @param a On input, the M-by-N matrix to invert.  The matrix is
 *  overwritten on output.
 * @param ainv The N-by-M matrix where the pseudo-inverse of @p a
 *  will be written.
 * @param err A pointer to the C error handler object.  If no error
 *  handling is desired, simply pass NULL, and errors will be dealt with
 *  by the default internal error handler.  Possible errors that may be
 *  encountered are as follows.
 *  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
 *      there is insufficient memory available.
 *  - LA_CONVERGENCE_ERROR: Occurs as a warning if the QR iteration process
 *      could not converge to a zero value.
 */
void mtx_pinverse(int m, int n, double *a, double *ainv, errorhandler err);

/** @brief Solves the overdetermined or underdetermined system (A*X = B) of
 * M equations of N unknowns using a complete orthogonal factorization of
 * matrix A.
 *
 * @param m The number of rows in the original coefficient matrix A.
 * @param n The number of columns in the original coefficient matrix A.
 * @param nrhs The number of right-hand-side vectors (number of columns
 *  in matrix @p b).
 * @param a On input, the M-by-N matrix A.  On output, the matrix
 *  is overwritten by the details of its complete orthogonal factorization.
 * @param mb The number of rows in the matrix @p b.  This value must be
 *  equal to MAX(M, N).
 * @param b If M >= N, the M-by-NRHS matrix B.  On output, the first
 *  N rows contain the N-by-NRHS solution matrix X.  If M < N, an
 *  N-by-NRHS matrix with the first M rows containing the matrix B.  On
 *  output, the N-by-NRHS solution matrix X.
 * @param err A pointer to the C error handler object.  If no error
 *  handling is desired, simply pass NULL, and errors will be dealt with
 *  by the default internal error handler.  Possible errors that may be
 *  encountered are as follows.
 *  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
 *      there is insufficient memory available.
 */
void solve_least_squares(int m, int n, int nrhs, double *a, int mb, double *b,
                         errorhandler err);

/** @brief Computes the eigenvalues, and optionally the eigenvectors of a
 * real, symmetric matrix.
 *
 * @param n The dimension of the matrix.
 * @param vecs Set to true to compute the eigenvectors as well as the
 *  eigenvalues; else, set to false to just compute the eigenvalues.
 * @param a On input, the N-by-N symmetric matrix on which to
 *  operate.  On output, and if @p vecs is set to true, the matrix will
 *  contain the eigenvectors (one per column) corresponding to each
 *  eigenvalue in @p vals.  If @p vecs is set to false, the lower triangular
 *  portion of the matrix is overwritten.
 * @param vals An N-element array that will contain the eigenvalues
 *  sorted into ascending order.
 * @param err A pointer to the C error handler object.  If no error
 *  handling is desired, simply pass NULL, and errors will be dealt with
 *  by the default internal error handler.  Possible errors that may be
 *  encountered are as follows.
 *  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
 *      there is insufficient memory available.
 *  - LA_CONVERGENCE_ERROR: Occurs if the algorithm failed to converge.
 */
void eigen_symm(int n, bool vecs, double *a, double *vals, errorhandler err);

/** @brief Computes the eigenvalues, and the right eigenvectors of a square 
 *  matrix.
 *
 * @param n The dimension of the matrix.
 * @param a On input, the N-by-N matrix on which to operate.  On
 *  output, the contents of this matrix are overwritten.
 * @param vals An N-element array containing the eigenvalues of the
 *  matrix on output.  The eigenvalues are not sorted.
 * @param vecs An N-by-N matrix containing the right eigenvectors 
 *  (one per column) on output.
 * @param err A pointer to the C error handler object.  If no error
 *  handling is desired, simply pass NULL, and errors will be dealt with
 *  by the default internal error handler.  Possible errors that may be
 *  encountered are as follows.
 *  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
 *      there is insufficient memory available.
 *  - LA_CONVERGENCE_ERROR: Occurs if the algorithm failed to converge.
 */
void eigen_asymm(int n, double *a, double complex *vals, double complex *vecs,
                 errorhandler err);

/** @brief Computes the eigenvalues, and optionally the right eigenvectors of
 * a square matrix assuming the structure of the eigenvalue problem is
 * A*X = lambda*B*X.
 *
 * @param a On input, the N-by-N matrix A.  On output, the contents
 *  of this matrix are overwritten.
 * @param b On input, the N-by-N matrix B.  On output, the contents
 *  of this matrix are overwritten.
 * @param alpha An N-element array that, on output, contains the 
 *  numerator of the eigenvalue ration ALPHA / BETA.  Computation of this
 *  ratio isn't necessarily as trivial as it seems as it is entirely 
 *  possible, and likely, that ALPHA / BETA can overflow or underflow.  With
 *  that said, the values in ALPHA will always be less than and usually 
 *  comparable with the NORM(A).
 * @param beta An N-element array that, on output, contains the
 *  denominator used to determine the eigenvalues as ALPHA / BETA.  The 
 *  values in this array will always be less than and usually comparable
 *  with the NORM(B).
 * @param vecs An N-by-N matrix containing the right eigenvectors 
 *  (one per column) on output.
 * @param err A pointer to the C error handler object.  If no error
 *  handling is desired, simply pass NULL, and errors will be dealt with
 *  by the default internal error handler.  Possible errors that may be
 *  encountered are as follows.
 *  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
 *      there is insufficient memory available.
 *  - LA_CONVERGENCE_ERROR: Occurs if the algorithm failed to converge.
 */
void eigen_gen(int n, double *a, double *b, double complex *alpha, double *beta,
               double complex *vecs, errorhandler err);

#ifdef __cplusplus
}
#endif
#endif // LINALG_H_DEFINED