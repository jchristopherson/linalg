#ifndef LINALG_H_
#define LINALG_H_

#include <stdbool.h>
#include <complex.h>

#define LA_NO_OPERATION             0
#define LA_TRANSPOSE                1
#define LA_HERMITIAN_TRANSPOSE      2
#define LA_NO_ERROR                 0
#define LA_INVALID_INPUT_ERROR      101
#define LA_ARRAY_SIZE_ERROR         102
#define LA_SINGULAR_MATRIX_ERROR    103
#define LA_MATRIX_FORMAT_ERROR      104
#define LA_OUT_OF_MEMORY_ERROR      105
#define LA_CONVERGENCE_ERROR        106
#define LA_INVALID_OPERATION_ERROR  107

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Performs the rank-1 update to matrix A such that:
 * A = alpha * X * Y**T + A, where A is an M-by-N matrix, alpha is a scalar,
 * X is an M-element array, and N is an N-element array.
 *
 * @param m The number of rows in the matrix.
 * @param n The number of columns in the matrix.
 * @param alpha The scalar multiplier.
 * @param x An M-element array.
 * @param y An N-element array.
 * @param a On input, the M-by-N matrix to update.  On output, the
 *  updated M-by-N matrix.
 * @param lda The leading dimension of matrix A.
 *
 * @return An error code.  The following codes are possible.
 *  - LA_NO_ERROR: No error occurred.  Successful operation.
 *  - LA_INVALID_INPUT_ERROR: Occurs if @p lda is not correct.
 */
int la_rank1_update(int m, int n, double alpha, const double *x, 
    const double *y, double *a, int lda);

/**
 * Performs the rank-1 update to matrix A such that:
 * A = alpha * X * Y**T + A, where A is an M-by-N matrix, alpha is a scalar,
 * X is an M-element array, and N is an N-element array.
 *
 * @param m The number of rows in the matrix.
 * @param n The number of columns in the matrix.
 * @param alpha The scalar multiplier.
 * @param x An M-element array.
 * @param y An N-element array.
 * @param a On input, the M-by-N matrix to update.  On output, the
 *  updated M-by-N matrix.
 * @param lda The leading dimension of matrix A.
 *
 * @return An error code.  The following codes are possible.
 *  - LA_NO_ERROR: No error occurred.  Successful operation.
 *  - LA_INVALID_INPUT_ERROR: Occurs if @p lda is not correct.
 */
int la_rank1_update_cmplx(int m, int n, double complex alpha,
    const double complex *x, const double complex *y, double complex *a,
    int lda);

/**
 * Computes the trace of a matrix (the sum of the main diagonal
 * elements).
 *
 * @param m The number of rows in the matrix.
 * @param n The number of columns in the matrix.
 * @param a The M-by-N matrix on which to operate.
 * @param lda The leading dimension of the matrix.
 * @param rst The results of the operation.
 *
 * @return An error code.  The following codes are possible.
 *  - LA_NO_ERROR: No error occurred.  Successful operation.
 *  - LA_INVALID_INPUT_ERROR: Occurs if @p lda is not correct.
 */
int la_trace(int m, int n, const double *a, int lda, double *rst);

/**
 * Computes the trace of a matrix (the sum of the main diagonal
 * elements).
 *
 * @param m The number of rows in the matrix.
 * @param n The number of columns in the matrix.
 * @param a The M-by-N matrix on which to operate.
 * @param lda The leading dimension of the matrix.
 * @param rst The results of the operation.
 *
 * @return An error code.  The following codes are possible.
 *  - LA_NO_ERROR: No error occurred.  Successful operation.
 *  - LA_INVALID_INPUT_ERROR: Occurs if @p lda is not correct.
 */
int la_trace_cmplx(int m, int n, const double complex *a, int lda,
    double complex *rst);

/**
 * Computes the matrix operation C = alpha * op(A) * op(B) + beta * C.
 * 
 * @param transa Set to true to compute op(A) as the transpose of A; else,
 *   set to false to compute op(A) as A.
 * @param transb Set to true to compute op(B) as the transpose of B; else,
 *   set to false to compute op(B) as B.
 * @param m The number of rows in @p c.
 * @param n The number of columns in @p c.
 * @param k The interior dimension of the product @p a and @p b.
 * @param alpha A scalar multiplier.
 * @param a If @p transa is true, this matrix must be @p k by @p m; else,
 *   if @p transa is false, this matrix must be @p m by @p k.
 * @param lda The leading dimension of matrix @p a.
 * @param b If @p transb is true, this matrix must be @p n by @p k; else,
 *   if @p transb is false, this matrix must be @p k by @p n.
 * @param ldb The leading dimension of matrix @p b.
 * @param beta A scalar multiplier.
 * @param c The @p m by @p n matrix C.
 * @param ldc The leading dimension of matrix @p c.
 * 
 * @return An error code.  The following codes are possible.
 *  - LA_NO_ERROR: No error occurred.  Successful operation.
 *  - LA_INVALID_INPUT_ERROR: Occurs if @p lda, @p ldb, or @p ldc are not
 *      correct.
 */
int la_mtx_mult(bool transa, bool transb, int m, int n, int k, double alpha,
    const double *a, int lda, const double *b, int ldb, double beta, 
    double *c, int ldc);

/**
 * Computes the matrix operation C = alpha * op(A) * op(B) + beta * C.
 *
 * @param opa Set to LA_TRANSPOSE to compute op(A) as a direct transpose of A,
 *  set to LA_HERMITIAN_TRANSPOSE to compute op(A) as the Hermitian transpose
 *  of A, otherwise, set to LA_NO_OPERATION to compute op(A) as A.
 * @param opb Set to TLA_RANSPOSE to compute op(B) as a direct transpose of B,
 *  set to LA_HERMITIAN_TRANSPOSE to compute op(B) as the Hermitian transpose
 *  of B, otherwise, set to LA_NO_OPERATION to compute op(B) as B.
 * @param mThe number of rows in @p c.
 * @param n The number of columns in @p c.
 * @param k The interior dimension of the product @p a and @p b.
 * @param alpha A scalar multiplier.
 * @param a If @p opa is LA_TRANSPOSE or LA_HERMITIAN_TRANSPOSE, this matrix 
 *  must be @p k by @p m; else, this matrix must be @p m by @p k.
 * @param lda The leading dimension of matrix @p a.
 * @param b If @p opb is LA_TRANSPOSE or LA_HERMITIAN_TRANSPOSE, this matrix 
 *  must be @p n by @p k; else, this matrix must be @p k by @p n.
 * @param ldb The leading dimension of matrix @p b.
 * @param beta A scalar multiplier.
 * @param c The @p m by @p n matrix C.
 * @param ldc The leading dimension of matrix @p c.
 * 
 * @return An error code.  The following codes are possible.
 *  - LA_NO_ERROR: No error occurred.  Successful operation.
 *  - LA_INVALID_INPUT_ERROR: Occurs if @p lda, @p ldb, or @p ldc are not
 *      correct.
 */
int la_mtx_mult_cmplx(int opa, int opb, int m, int n, int k, 
    double complex alpha, const double complex *a, int lda,
    const double complex *b, int ldb, double complex beta, double complex *c,
    int ldc);

/**
 * Computes the matrix operation: C = alpha * A * op(B) + beta * C,
 * or C = alpha * op(B) * A + beta * C.
 *
 * @param lside Set to true to apply matrix A from the left; else, set
 *  to false to apply matrix A from the left.
 * @param trans Set to true if op(B) == B**T; else, set to false if
 *  op(B) == B.
 * @param m The number of rows in the matrix C.
 * @param n The number of columns in the matrix C.
 * @param k The inner dimension of the matrix product A * op(B).
 * @param alpha A scalar multiplier.
 * @param a A P-element array containing the diagonal elements of matrix A
 *  where P = MIN(@p m, @p k) if @p lside is true; else, P = MIN(@p n, @p k)
 *  if @p lside is false.
 * @param b The LDB-by-TDB matrix B where (LDB = leading dimension of B,
 *  and TDB = trailing dimension of B):
 *  - @p lside == true & @p trans == true: LDB = @p n, TDB = @p k
 *  - @p lside == true & @p trans == false: LDB = @p k, TDB = @p n
 *  - @p lside == false & @p trans == true: LDB = @p k, TDB = @p m
 *  - @p lside == false & @p trans == false: LDB = @p m, TDB = @p k
 * @param ldb The leading dimension of matrix B.
 * @param beta A scalar multiplier.
 * @param c The @p m by @p n matrix C.
 * @param ldc The leading dimension of matrix C.
 *
 * @return An error code.  The following codes are possible.
 *  - LA_NO_ERROR: No error occurred.  Successful operation.
 *  - LA_INVALID_INPUT_ERROR: Occurs if @p ldb, or @p ldc are not
 *      correct.
 *  - LA_ARRAY_SIZE_ERROR: Occurs if any of the input array sizes are
 *      incorrect.
 */
int la_diag_mtx_mult(bool lside, bool transb, int m, int n, int k, 
    double alpha, const double *a, const double *b, int ldb, double beta, 
    double *c, int ldc);

/**
 * Computes the matrix operation: C = alpha * A * op(B) + beta * C,
 * or C = alpha * op(B) * A + beta * C.
 *
 * @param lside Set to true to apply matrix A from the left; else, set
 *  to false to apply matrix A from the left.
 * @param opb Set to LA_TRANSPOSE to compute op(B) as a direct transpose of B,
 *  set to LA_HERMITIAN_TRANSPOSE to compute op(B) as the Hermitian transpose
 *  of B, otherwise, set to LA_NO_OPERATION to compute op(B) as B.
 * @param m The number of rows in the matrix C.
 * @param n The number of columns in the matrix C.
 * @param k The inner dimension of the matrix product A * op(B).
 * @param alpha A scalar multiplier.
 * @param a A P-element array containing the diagonal elements of matrix A
 *  where P = MIN(@p m, @p k) if @p lside is true; else, P = MIN(@p n, @p k)
 *  if @p lside is false.
 * @param b The LDB-by-TDB matrix B where (LDB = leading dimension of B,
 *  and TDB = trailing dimension of B):
 *  - @p lside == true & @p trans == true: LDB = @p n, TDB = @p k
 *  - @p lside == true & @p trans == false: LDB = @p k, TDB = @p n
 *  - @p lside == false & @p trans == true: LDB = @p k, TDB = @p m
 *  - @p lside == false & @p trans == false: LDB = @p m, TDB = @p k
 * @param ldb The leading dimension of matrix B.
 * @param beta A scalar multiplier.
 * @param c The @p m by @p n matrix C.
 * @param ldc The leading dimension of matrix C.
 * 
 * @return An error code.  The following codes are possible.
 *  - LA_NO_ERROR: No error occurred.  Successful operation.
 *  - LA_INVALID_INPUT_ERROR: Occurs if @p ldb, or @p ldc are not
 *      correct.
 *  - LA_ARRAY_SIZE_ERROR: Occurs if any of the input array sizes are
 *      incorrect.
 */
int la_diag_mtx_mult_cmplx(bool lside, int opb, int m, int n, int k, 
    double complex alpha, const double complex *a, const double complex *b, 
    int ldb, double complex beta, double complex *c, int ldc);

/**
 * Computes the matrix operation: C = alpha * A * op(B) + beta * C,
 * or C = alpha * op(B) * A + beta * C.
 *
 * @param lside Set to true to apply matrix A from the left; else, set
 *  to false to apply matrix A from the left.
 * @param opb Set to LA_TRANSPOSE to compute op(B) as a direct transpose of B,
 *  set to LA_HERMITIAN_TRANSPOSE to compute op(B) as the Hermitian transpose
 *  of B, otherwise, set to LA_NO_OPERATION to compute op(B) as B.
 * @param m The number of rows in the matrix C.
 * @param n The number of columns in the matrix C.
 * @param k The inner dimension of the matrix product A * op(B).
 * @param alpha A scalar multiplier.
 * @param a A P-element array containing the diagonal elements of matrix A
 *  where P = MIN(@p m, @p k) if @p lside is true; else, P = MIN(@p n, @p k)
 *  if @p lside is false.
 * @param b The LDB-by-TDB matrix B where (LDB = leading dimension of B,
 *  and TDB = trailing dimension of B):
 *  - @p lside == true & @p trans == true: LDB = @p n, TDB = @p k
 *  - @p lside == true & @p trans == false: LDB = @p k, TDB = @p n
 *  - @p lside == false & @p trans == true: LDB = @p k, TDB = @p m
 *  - @p lside == false & @p trans == false: LDB = @p m, TDB = @p k
 * @param ldb The leading dimension of matrix B.
 * @param beta A scalar multiplier.
 * @param c The @p m by @p n matrix C.
 * @param ldc The leading dimension of matrix C.
 * 
 * @return An error code.  The following codes are possible.
 *  - LA_NO_ERROR: No error occurred.  Successful operation.
 *  - LA_INVALID_INPUT_ERROR: Occurs if @p ldb, or @p ldc are not
 *      correct.
 *  - LA_ARRAY_SIZE_ERROR: Occurs if any of the input array sizes are
 *      incorrect.
 */
int la_diag_mtx_mult_mixed(bool lside, int opb, int m, int n, int k, 
    double complex alpha, const double *a, const double complex *b, 
    int ldb, double complex beta, double complex *c, int ldc);

/**
 * Computes the rank of a matrix.
 *
 * @param m The number of rows in the matrix.
 * @param n The number of columns in the matrix.
 * @param a The M-by-N matrix.  The matrix is overwritten as part of this
 *  operation.
 * @param lda The leading dimension of matrix A.
 * @param rnk The rank of @p a.
 *
 * @return An error code.  The following codes are possible.
 *  - LA_NO_ERROR: No error occurred.  Successful operation.
 *  - LA_INVALID_INPUT_ERROR: Occurs if @p lda is not correct.
 *  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
 *      there is insufficient memory available.
 *  - LA_CONVERGENCE_ERROR: Occurs as a warning if the QR iteration process
 *      could not converge to a zero value.
 */
int la_rank(int m, int n, double *a, int lda, int *rnk);

/**
 * Computes the rank of a matrix.
 *
 * @param m The number of rows in the matrix.
 * @param n The number of columns in the matrix.
 * @param a The M-by-N matrix.  The matrix is overwritten as part of this
 *  operation.
 * @param lda The leading dimension of matrix A.
 * @param rnk The rank of @p a.
 *
 * @return An error code.  The following codes are possible.
 *  - LA_NO_ERROR: No error occurred.  Successful operation.
 *  - LA_INVALID_INPUT_ERROR: Occurs if @p lda is not correct.
 *  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
 *      there is insufficient memory available.
 *  - LA_CONVERGENCE_ERROR: Occurs as a warning if the QR iteration process
 *      could not converge to a zero value.
 */
int la_rank_cmplx(int m, int n, double complex *a, int lda, int *rnk);

/**
 * Computes the determinant of a square matrix.
 *
 * @param n The dimension of the matrix.
 * @param a The N-by-N matrix.  The matrix is overwritten on output.
 * @param lda The leading dimension of the matrix.
 * @param d The determinant of @p a.
 *
 * @return An error code.  The following codes are possible.
 *  - LA_NO_ERROR: No error occurred.  Successful operation.
 *  - LA_ARRAY_SIZE_ERROR: Occurs if any of the input arrays are not sized
 *      appropriately.
 *  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
 *      there is insufficient memory available.
 */
int la_det(int n, double *a, int lda, double *d);

/**
 * Computes the determinant of a square matrix.
 *
 * @param n The dimension of the matrix.
 * @param a The N-by-N matrix.  The matrix is overwritten on output.
 * @param lda The leading dimension of the matrix.
 * @param d The determinant of @p a.
 *
 * @return An error code.  The following codes are possible.
 *  - LA_NO_ERROR: No error occurred.  Successful operation.
 *  - LA_ARRAY_SIZE_ERROR: Occurs if any of the input arrays are not sized
 *      appropriately.
 *  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
 *      there is insufficient memory available.
 */
int la_det_cmplx(int n, double complex *a, int lda, double complex *d);

/**
 * Computes the triangular matrix operation:
 * B = alpha * A**T * A + beta * B, or B = alpha * A * A**T + beta * B,
 * where A is a triangular matrix.
 *
 * @param upper Set to true if matrix A is upper triangular, and
 *  B = alpha * A**T * A + beta * B is to be calculated; else, set to false
 *  if A is lower triangular, and B = alpha * A * A**T + beta * B is to
 *  be computed.
 * @param alpha A scalar multiplier.
 * @param n The dimension of the matrix.
 * @param a The @p n by @p n triangular matrix A.  Notice, if @p upper is
 *  true, only the upper triangular portion of this matrix is referenced;
 *  else, if @p upper is false, only the lower triangular portion of this
 *  matrix is referenced.
 * @param lda The leading dimension of matrix A.
 * @param beta A scalar multiplier.
 * @param b On input, the @p n by @p n matrix B.  On output, the @p n by
 *  @p n resulting matrix.
 * @param ldb The leading dimension of matrix B.
 *
 * @return An error code.  The following codes are possible.
 *  - LA_NO_ERROR: No error occurred.  Successful operation.
 *  - LA_INVALID_INPUT_ERROR: Occurs if @p lda or @p ldb are not correct.
 */
int la_tri_mtx_mult(bool upper, double alpha, int n, const double *a, int lda,
    double beta, double *b, int ldb);

/**
 * Computes the triangular matrix operation:
 * B = alpha * A**T * A + beta * B, or B = alpha * A * A**T + beta * B,
 * where A is a triangular matrix.
 *
 * @param upper Set to true if matrix A is upper triangular, and
 *  B = alpha * A**T * A + beta * B is to be calculated; else, set to false
 *  if A is lower triangular, and B = alpha * A * A**T + beta * B is to
 *  be computed.
 * @param alpha A scalar multiplier.
 * @param n The dimension of the matrix.
 * @param a The @p n by @p n triangular matrix A.  Notice, if @p upper is
 *  true, only the upper triangular portion of this matrix is referenced;
 *  else, if @p upper is false, only the lower triangular portion of this
 *  matrix is referenced.
 * @param lda The leading dimension of matrix A.
 * @param beta A scalar multiplier.
 * @param b On input, the @p n by @p n matrix B.  On output, the @p n by
 *  @p n resulting matrix.
 * @param ldb The leading dimension of matrix B.
 *
 * @return An error code.  The following codes are possible.
 *  - LA_NO_ERROR: No error occurred.  Successful operation.
 *  - LA_INVALID_INPUT_ERROR: Occurs if @p lda or @p ldb are not correct.
 */
int la_tri_mtx_mult_cmplx(bool upper, double complex alpha, int n, 
    const double complex *a, int lda, double complex beta, 
    double complex *b, int ldb);

/**
 * Computes the LU factorization of an M-by-N matrix.
 *
 * @param m The number of rows in the matrix.
 * @param n The number of columns in the matrix.
 * @param a On input, the M-by-N matrix on which to operate.  On
 * output, the LU factored matrix in the form [L\\U] where the unit diagonal
 * elements of L are not stored.
 * @param lda The leading dimension of matrix A.
 * @param ipvt An MIN(M, N)-element array used to track row-pivot
 *  operations.  The array stored pivot information such that row I is
 *  interchanged with row IPVT(I).
 *
 * @return An error code.  The following codes are possible.
 *  - LA_NO_ERROR: No error occurred.  Successful operation.
 *  - LA_INVALID_INPUT_ERROR: Occurs if @p lda is not correct.
 *  - LA_SINGULAR_MATRIX_ERROR: Occurs as a warning if @p a is found to be
 *      singular.
 */
int la_lu_factor(int m, int n, double *a, int lda, int *ipvt);

/**
 * Computes the LU factorization of an M-by-N matrix.
 *
 * @param m The number of rows in the matrix.
 * @param n The number of columns in the matrix.
 * @param a On input, the M-by-N matrix on which to operate.  On
 * output, the LU factored matrix in the form [L\\U] where the unit diagonal
 * elements of L are not stored.
 * @param lda The leading dimension of matrix A.
 * @param ipvt An MIN(M, N)-element array used to track row-pivot
 *  operations.  The array stored pivot information such that row I is
 *  interchanged with row IPVT(I).
 *
 * @return An error code.  The following codes are possible.
 *  - LA_NO_ERROR: No error occurred.  Successful operation.
 *  - LA_INVALID_INPUT_ERROR: Occurs if @p lda is not correct.
 *  - LA_SINGULAR_MATRIX_ERROR: Occurs as a warning if @p a is found to be
 *      singular.
 */
int la_lu_factor_cmplx(int m, int n, double complex *a, int lda, int *ipvt);

/**
 * Extracts the L, U, and P matrices from the LU factorization
 * output from la_lu_factor.
 *
 * @param n The dimension of the input matrix.
 * @param a On input, the N-by-N matrix as output by
 *  @ref la_lu_factor.  On output, the N-by-N lower triangular matrix L.
 * @param lda The leading dimension of @p a.
 * @param ipvt The N-element pivot array as output by
 *  @ref la_lu_factor.
 * @param u An N-by-N matrix where the U matrix will be written.
 * @param ldu The leading dimension of @p u.
 * @param p An N-by-N matrix where the row permutation matrix will be
 *  written.
 * @param ldp The leading dimension of @p p.
 *
 * @return An error code.  The following codes are possible.
 *  - LA_NO_ERROR: No error occurred.  Successful operation.
 *  - LA_INVALID_INPUT_ERROR: Occurs if @p lda, @p ldu, or @p ldp is not 
 *      correct.
 */
int la_form_lu(int n, double *a, int lda, int *ipvt, double *u, int ldu, 
    double *p, int ldp);

/**
 * Extracts the L, U, and P matrices from the LU factorization
 * output from la_lu_factor_cmplx.
 *
 * @param n The dimension of the input matrix.
 * @param a On input, the N-by-N matrix as output by
 *  @ref la_lu_factor_cmplx.  On output, the N-by-N lower triangular matrix L.
 * @param lda The leading dimension of @p a.
 * @param ipvt The N-element pivot array as output by
 *  @ref la_lu_factor_cmplx.
 * @param u An N-by-N matrix where the U matrix will be written.
 * @param ldu The leading dimension of @p u.
 * @param p An N-by-N matrix where the row permutation matrix will be
 *  written.
 * @param ldp The leading dimension of @p p.
 *
 * @return An error code.  The following codes are possible.
 *  - LA_NO_ERROR: No error occurred.  Successful operation.
 *  - LA_INVALID_INPUT_ERROR: Occurs if @p lda, @p ldu, or @p ldp is not 
 *      correct.
 */
int la_form_lu_cmplx(int n, double complex *a, int lda, int *ipvt, 
    double complex *u, int ldu, double *p, int ldp);

/**
 * Computes the QR factorization of an M-by-N matrix without
 * pivoting.
 *
 * @param m The number of rows in the matrix.
 * @param n The number of columns in the matrix.
 * @param a  On input, the M-by-N matrix to factor.  On output, the
 *  elements on and above the diagonal contain the MIN(M, N)-by-N upper
 *  trapezoidal matrix R (R is upper triangular if M >= N).  The elements
 *  below the diagonal, along with the array @p tau, represent the
 *  orthogonal matrix Q as a product of elementary reflectors.
 * @param lda The leading dimension of matrix A.
 * @param tau A MIN(M, N)-element array used to store the scalar
 *  factors of the elementary reflectors.
 *
 * @return An error code.  The following codes are possible.
 *  - LA_NO_ERROR: No error occurred.  Successful operation.
 *  - LA_INVALID_INPUT_ERROR: Occurs if @p lda is not correct.
 *  - LA_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory 
 *      available.
 */
int la_qr_factor(int m, int n, double *a, int lda, double *tau);

/**
 * Computes the QR factorization of an M-by-N matrix without
 * pivoting.
 *
 * @param m The number of rows in the matrix.
 * @param n The number of columns in the matrix.
 * @param a  On input, the M-by-N matrix to factor.  On output, the
 *  elements on and above the diagonal contain the MIN(M, N)-by-N upper
 *  trapezoidal matrix R (R is upper triangular if M >= N).  The elements
 *  below the diagonal, along with the array @p tau, represent the
 *  orthogonal matrix Q as a product of elementary reflectors.
 * @param lda The leading dimension of matrix A.
 * @param tau A MIN(M, N)-element array used to store the scalar
 *  factors of the elementary reflectors.
 *
 * @return An error code.  The following codes are possible.
 *  - LA_NO_ERROR: No error occurred.  Successful operation.
 *  - LA_INVALID_INPUT_ERROR: Occurs if @p lda is not correct.
 *  - LA_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory 
 *      available.
 */
int la_qr_factor_cmplx(int m, int n, double complex *a, int lda, 
    double complex *tau);

/**
 * Computes the QR factorization of an M-by-N matrix with column pivoting.
 *
 * @param m The number of rows in the matrix.
 * @param n The number of columns in the matrix.
 * @param a  On input, the M-by-N matrix to factor.  On output, the
 *  elements on and above the diagonal contain the MIN(M, N)-by-N upper
 *  trapezoidal matrix R (R is upper triangular if M >= N).  The elements
 *  below the diagonal, along with the array @p tau, represent the
 *  orthogonal matrix Q as a product of elementary reflectors.
 * @param lda The leading dimension of matrix A.
 * @param tau A MIN(M, N)-element array used to store the scalar
 *  factors of the elementary reflectors.
 * @param jpvt On input, an N-element array that if JPVT(I) .ne. 0,
 *  the I-th column of A is permuted to the front of A * P; if JPVT(I) = 0,
 *  the I-th column of A is a free column.  On output, if JPVT(I) = K, then
 *  the I-th column of A * P was the K-th column of A.
 *
 * @return An error code.  The following codes are possible.
 *  - LA_NO_ERROR: No error occurred.  Successful operation.
 *  - LA_INVALID_INPUT_ERROR: Occurs if @p lda is not correct.
 *  - LA_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory 
 *      available.
 */
int la_qr_factor_pvt(int m, int n, double *a, int lda, double *tau, int *jpvt);

/**
 * Computes the QR factorization of an M-by-N matrix with column pivoting.
 *
 * @param m The number of rows in the matrix.
 * @param n The number of columns in the matrix.
 * @param a  On input, the M-by-N matrix to factor.  On output, the
 *  elements on and above the diagonal contain the MIN(M, N)-by-N upper
 *  trapezoidal matrix R (R is upper triangular if M >= N).  The elements
 *  below the diagonal, along with the array @p tau, represent the
 *  orthogonal matrix Q as a product of elementary reflectors.
 * @param lda The leading dimension of matrix A.
 * @param tau A MIN(M, N)-element array used to store the scalar
 *  factors of the elementary reflectors.
 * @param jpvt On input, an N-element array that if JPVT(I) .ne. 0,
 *  the I-th column of A is permuted to the front of A * P; if JPVT(I) = 0,
 *  the I-th column of A is a free column.  On output, if JPVT(I) = K, then
 *  the I-th column of A * P was the K-th column of A.
 *
 * @return An error code.  The following codes are possible.
 *  - LA_NO_ERROR: No error occurred.  Successful operation.
 *  - LA_INVALID_INPUT_ERROR: Occurs if @p lda is not correct.
 *  - LA_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory 
 *      available.
 */
int la_qr_factor_cmplx_pvt(int m, int n, double complex *a, int lda,
    double complex *tau, int *jpvt);

/**
 * Forms the full M-by-M orthogonal matrix Q from the elementary
 * reflectors returned by the base QR factorization algorithm.
 *
 * @param fullq Set to true to always return the full Q matrix; else,
 *  set to false, and in the event that M > N, Q may be supplied as M-by-N,
 *  and therefore only return the useful submatrix Q1 (Q = [Q1, Q2]) as the
 *  factorization can be written as Q * R = [Q1, Q2] * [R1; 0].
 * @param m The number of rows in R.
 * @param n The number of columns in R.
 * @param r On input, the M-by-N factored matrix as returned by the
 *  QR factorization routine.  On output, the upper triangular matrix R.
 * @param ldr The leading dimension of matrix R.
 * @param tau A MIN(M, N)-element array containing the scalar factors of
 *  each elementary reflector defined in @p r.
 * @param q An M-by-M matrix where the full Q matrix will be written.
 *  In the event that @p fullq is set to false, and M > N, this matrix need
 *  only by M-by-N.
 * @param ldq The leading dimension of matrix Q.
 *
 * @return An error code.  The following codes are possible.
 *  - LA_NO_ERROR: No error occurred.  Successful operation.
 *  - LA_INVALID_INPUT_ERROR: Occurs if @p lda is not correct.
 *  - LA_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory 
 *      available.
 */
int la_form_qr(bool fullq, int m, int n, double *r, int ldr, const double *tau,
    double *q, int ldq);

/**
 * Forms the full M-by-M orthogonal matrix Q from the elementary
 * reflectors returned by the base QR factorization algorithm.
 *
 * @param fullq Set to true to always return the full Q matrix; else,
 *  set to false, and in the event that M > N, Q may be supplied as M-by-N,
 *  and therefore only return the useful submatrix Q1 (Q = [Q1, Q2]) as the
 *  factorization can be written as Q * R = [Q1, Q2] * [R1; 0].
 * @param m The number of rows in R.
 * @param n The number of columns in R.
 * @param r On input, the M-by-N factored matrix as returned by the
 *  QR factorization routine.  On output, the upper triangular matrix R.
 * @param ldr The leading dimension of matrix R.
 * @param tau A MIN(M, N)-element array containing the scalar factors of
 *  each elementary reflector defined in @p r.
 * @param q An M-by-M matrix where the full Q matrix will be written.
 *  In the event that @p fullq is set to false, and M > N, this matrix need
 *  only by M-by-N.
 * @param ldq The leading dimension of matrix Q.
 *
 * @return An error code.  The following codes are possible.
 *  - LA_NO_ERROR: No error occurred.  Successful operation.
 *  - LA_INVALID_INPUT_ERROR: Occurs if @p lda is not correct.
 *  - LA_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory 
 *      available.
 */
int la_form_qr_cmplx(bool fullq, int m, int n, double complex *r, int ldr,
    const double complex *tau, double complex *q, int ldq);

/**
 * Forms the full M-by-M orthogonal matrix Q from the elementary
 * reflectors returned by the base QR factorization algorithm.  This
 * routine also inflates the pivot array into an N-by-N matrix P such
 * that A * P = Q * R.
 *
 * @param fullq Set to true to always return the full Q matrix; else,
 *  set to false, and in the event that M > N, Q may be supplied as M-by-N,
 *  and therefore only return the useful submatrix Q1 (Q = [Q1, Q2]) as the
 *  factorization can be written as Q * R = [Q1, Q2] * [R1; 0].
 * @param m The number of rows in R.
 * @param n The number of columns in R.
 * @param r On input, the M-by-N factored matrix as returned by the
 *  QR factorization routine.  On output, the upper triangular matrix R.
 * @param ldr The leading dimension of matrix R.
 * @param tau A MIN(M, N)-element array containing the scalar factors of
 *  each elementary reflector defined in @p r.
 * @param pvt An N-element array containing the pivot information from
 *  the QR factorization.
 * @param q An M-by-M matrix where the full Q matrix will be written.
 *  In the event that @p fullq is set to false, and M > N, this matrix need
 *  only by M-by-N.
 * @param ldq The leading dimension of matrix Q.
 * @param p An N-by-N matrix where the pivot matrix P will be written.
 * @param ldp The leading dimension of matrix P.
 *
 * @return An error code.  The following codes are possible.
 *  - LA_NO_ERROR: No error occurred.  Successful operation.
 *  - LA_INVALID_INPUT_ERROR: Occurs if @p lda is not correct.
 *  - LA_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory 
 *      available.
 */
int la_form_qr_pvt(bool fullq, int m, int n, double *r, int ldr, 
    const double *tau, const int *pvt, double *q, int ldq, double *p, int ldp);

/**
 * Forms the full M-by-M orthogonal matrix Q from the elementary
 * reflectors returned by the base QR factorization algorithm.  This
 * routine also inflates the pivot array into an N-by-N matrix P such
 * that A * P = Q * R.
 *
 * @param fullq Set to true to always return the full Q matrix; else,
 *  set to false, and in the event that M > N, Q may be supplied as M-by-N,
 *  and therefore only return the useful submatrix Q1 (Q = [Q1, Q2]) as the
 *  factorization can be written as Q * R = [Q1, Q2] * [R1; 0].
 * @param m The number of rows in R.
 * @param n The number of columns in R.
 * @param r On input, the M-by-N factored matrix as returned by the
 *  QR factorization routine.  On output, the upper triangular matrix R.
 * @param ldr The leading dimension of matrix R.
 * @param tau A MIN(M, N)-element array containing the scalar factors of
 *  each elementary reflector defined in @p r.
 * @param pvt An N-element array containing the pivot information from
 *  the QR factorization.
 * @param q An M-by-M matrix where the full Q matrix will be written.
 *  In the event that @p fullq is set to false, and M > N, this matrix need
 *  only by M-by-N.
 * @param ldq The leading dimension of matrix Q.
 * @param p An N-by-N matrix where the pivot matrix P will be written.
 * @param ldp The leading dimension of matrix P.
 *
 * @return An error code.  The following codes are possible.
 *  - LA_NO_ERROR: No error occurred.  Successful operation.
 *  - LA_INVALID_INPUT_ERROR: Occurs if @p lda is not correct.
 *  - LA_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory 
 *      available.
 */
int la_form_qr_cmplx_pvt(bool fullq, int m, int n, double complex *r, int ldr,
    const double complex *tau, const int *pvt, double complex *q, int ldq,
    double complex *p, int ldp);

/**
 * Multiplies a general matrix by the orthogonal matrix Q from a QR
 * factorization such that: C = op(Q) * C, or C = C * op(Q).
 *
 * @param lside Set to true to apply Q or Q**T from the left; else, set
 *  to false to apply Q or Q**T from the right.
 * @param trans Set to true to apply Q**T; else, set to false.
 * @param m The number of rows in matrix C.
 * @param n The number of columns in matrix C.
 * @param k The number of elementary reflectors whose product defines 
 *  the matrix Q.
 * @param a On input, an LDA-by-K matrix containing the elementary
 *  reflectors output from the QR factorization.  If @p lside is set to
 *  true, LDA = M, and M >= K >= 0; else, if @p lside is set to false,
 *  LDA = N, and N >= K >= 0.  Notice, the contents of this matrix are
 *  restored on exit.
 * @param lda The leading dimension of matrix A.
 * @param tau A K-element array containing the scalar factors of each
 *  elementary reflector defined in @p a.
 * @param c On input, the M-by-N matrix C.  On output, the product
 *  of the orthogonal matrix Q and the original matrix C.
 * @param ldc THe leading dimension of matrix C.
 *
 * @return An error code.  The following codes are possible.
 *  - LA_NO_ERROR: No error occurred.  Successful operation.
 *  - LA_INVALID_INPUT_ERROR: Occurs if @p lda is not correct.
 *  - LA_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory 
 *      available.
 */
int la_mult_qr(bool lside, bool trans, int m, int n, int k, double *a, int lda,
    const double *tau, double *c, int ldc);

/**
 * Multiplies a general matrix by the orthogonal matrix Q from a QR
 * factorization such that: C = op(Q) * C, or C = C * op(Q).
 *
 * @param lside Set to true to apply Q or Q**H from the left; else, set
 *  to false to apply Q or Q**H from the right.
 * @param trans Set to true to apply Q**H; else, set to false.
 * @param m The number of rows in matrix C.
 * @param n The number of columns in matrix C.
 * @param k The number of elementary reflectors whose product defines 
 *  the matrix Q.
 * @param a On input, an LDA-by-K matrix containing the elementary
 *  reflectors output from the QR factorization.  If @p lside is set to
 *  true, LDA = M, and M >= K >= 0; else, if @p lside is set to false,
 *  LDA = N, and N >= K >= 0.  Notice, the contents of this matrix are
 *  restored on exit.
 * @param lda The leading dimension of matrix A.
 * @param tau A K-element array containing the scalar factors of each
 *  elementary reflector defined in @p a.
 * @param c On input, the M-by-N matrix C.  On output, the product
 *  of the orthogonal matrix Q and the original matrix C.
 * @param ldc THe leading dimension of matrix C.
 *
 * @return An error code.  The following codes are possible.
 *  - LA_NO_ERROR: No error occurred.  Successful operation.
 *  - LA_INVALID_INPUT_ERROR: Occurs if @p lda is not correct.
 *  - LA_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory 
 *      available.
 */
int la_mult_qr_cmplx(bool lside, bool trans, int m, int n, int k, 
    double complex *a, int lda, const double complex *tau, double complex *c,
    int ldc);

/**
 * Computes the rank 1 update to an M-by-N QR factored matrix A
 * (M >= N) where A = Q * R, and A1 = A + U * V**T such that A1 = Q1 * R1.
 *
 * @param m The number of rows in R.
 * @param n The number of columns in R.
 * @param q On input, the original M-by-K orthogonal matrix Q.  On
 *  output, the updated matrix Q1.
 * @param ldq The leading dimension of matrix Q.
 * @param r On input, the M-by-N matrix R.  On output, the updated
 *  matrix R1.
 * @param ldr The leading dimension of matrix R.
 * @param u On input, the M-element U update vector.  On output,
 *  the original content of the array is overwritten.
 * @param v On input, the N-element V update vector.  On output,
 *  the original content of the array is overwritten.
 *
 * @return An error code.  The following codes are possible.
 *  - LA_NO_ERROR: No error occurred.  Successful operation.
 *  - LA_INVALID_INPUT_ERROR: Occurs if @p ldq or @p ldr is not correct.
 *  - LA_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory 
 *      available.
 */
int la_qr_rank1_update(int m, int n, double *q, int ldq, double *r, int ldr,
    double *u, double *v);

/**
 * Computes the rank 1 update to an M-by-N QR factored matrix A
 * (M >= N) where A = Q * R, and A1 = A + U * V**T such that A1 = Q1 * R1.
 *
 * @param m The number of rows in R.
 * @param n The number of columns in R.
 * @param q On input, the original M-by-K orthogonal matrix Q.  On
 *  output, the updated matrix Q1.
 * @param ldq The leading dimension of matrix Q.
 * @param r On input, the M-by-N matrix R.  On output, the updated
 *  matrix R1.
 * @param ldr The leading dimension of matrix R.
 * @param u On input, the M-element U update vector.  On output,
 *  the original content of the array is overwritten.
 * @param v On input, the N-element V update vector.  On output,
 *  the original content of the array is overwritten.
 *
 * @return An error code.  The following codes are possible.
 *  - LA_NO_ERROR: No error occurred.  Successful operation.
 *  - LA_INVALID_INPUT_ERROR: Occurs if @p ldq or @p ldr is not correct.
 *  - LA_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory 
 *      available.
 */
int la_qr_rank1_update_cmplx(int m, int n, double complex *q, int ldq,
    double complex *r, int ldr, double complex *u, double complex *v);

/**
 * Computes the Cholesky factorization of a symmetric, positive
 * definite matrix.
 *
 * @param upper Set to true to compute the upper triangular factoriztion
 *  A = U**T * U; else, set to false to compute the lower triangular
 *  factorzation A = L * L**T.
 * @param n The dimension of matrix A.
 * @param a On input, the N-by-N matrix to factor.  On output, the
 *  factored matrix is returned in either the upper or lower triangular
 *  portion of the matrix, dependent upon the value of @p upper.
 * @param lda The leading dimension of matrix A.
 *
 * @return An error code.  The following codes are possible.
 *  - LA_NO_ERROR: No error occurred.  Successful operation.
 *  - LA_INVALID_INPUT_ERROR: Occurs if @p lda is not correct.
 *  - LA_MATRIX_FORMAT_ERROR: Occurs if @p a is not positive definite.
 */
int la_cholesky_factor(bool upper, int n, double *a, int lda);

/**
 * Computes the Cholesky factorization of a symmetric, positive
 * definite matrix.
 *
 * @param upper Set to true to compute the upper triangular factoriztion
 *  A = U**T * U; else, set to false to compute the lower triangular
 *  factorzation A = L * L**T.
 * @param n The dimension of matrix A.
 * @param a On input, the N-by-N matrix to factor.  On output, the
 *  factored matrix is returned in either the upper or lower triangular
 *  portion of the matrix, dependent upon the value of @p upper.
 * @param lda The leading dimension of matrix A.
 *
 * @return An error code.  The following codes are possible.
 *  - LA_NO_ERROR: No error occurred.  Successful operation.
 *  - LA_INVALID_INPUT_ERROR: Occurs if @p lda is not correct.
 *  - LA_MATRIX_FORMAT_ERROR: Occurs if @p a is not positive definite.
 */
int la_cholesky_factor_cmplx(bool upper, int n, double complex *a, int lda);

/**
 * Computes the rank 1 update to a Cholesky factored matrix (upper
 * triangular).
 *
 * @param n The dimension of the matrix.
 * @param r On input, the N-by-N upper triangular matrix R.  On
 *  output, the updated matrix R1.
 * @param ldr The leading dimension of matrix R.
 * @param u On input, the N-element update vector U.  On output,
 *  the rotation sines used to transform R to R1.
 *
 * @return An error code.  The following codes are possible.
 *  - LA_NO_ERROR: No error occurred.  Successful operation.
 *  - LA_INVALID_INPUT_ERROR: Occurs if @p ldr is not correct.
 *  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
 *      there is insufficient memory available.
 */
int la_cholesky_rank1_update(int n, double *r, int ldr, double *u);

/**
 * Computes the rank 1 update to a Cholesky factored matrix (upper
 * triangular).
 *
 * @param n The dimension of the matrix.
 * @param r On input, the N-by-N upper triangular matrix R.  On
 *  output, the updated matrix R1.
 * @param ldr The leading dimension of matrix R.
 * @param u On input, the N-element update vector U.  On output,
 *  the rotation sines used to transform R to R1.
 *
 * @return An error code.  The following codes are possible.
 *  - LA_NO_ERROR: No error occurred.  Successful operation.
 *  - LA_INVALID_INPUT_ERROR: Occurs if @p ldr is not correct.
 *  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
 *      there is insufficient memory available.
 */
int la_cholesky_rank1_update_cmplx(int n, double complex *r, int ldr, 
    double complex *u);

/**
 * Computes the rank 1 downdate to a Cholesky factored matrix (upper
 * triangular).
 *
 * @param n The dimension of the matrix.
 * @param r On input, the N-by-N upper triangular matrix R.  On
 *  output, the updated matrix R1.
 * @param ldr The leading dimension of matrix R.
 * @param u On input, the N-element update vector U.  On output,
 *  the rotation sines used to transform R to R1.
 *
 * @return An error code.  The following codes are possible.
 *  - LA_NO_ERROR: No error occurred.  Successful operation.
 *  - LA_INVALID_INPUT_ERROR: Occurs if @p ldr is not correct.
 *  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
 *      there is insufficient memory available.
 *  - LA_MATRIX_FORMAT_ERROR: Occurs if the downdated matrix is not
 *      positive definite.
 */
int la_cholesky_rank1_downdate(int n, double *r, int ldr, double *u);

/**
 * Computes the rank 1 downdate to a Cholesky factored matrix (upper
 * triangular).
 *
 * @param n The dimension of the matrix.
 * @param r On input, the N-by-N upper triangular matrix R.  On
 *  output, the updated matrix R1.
 * @param ldr The leading dimension of matrix R.
 * @param u On input, the N-element update vector U.  On output,
 *  the rotation sines used to transform R to R1.
 *
 * @return An error code.  The following codes are possible.
 *  - LA_NO_ERROR: No error occurred.  Successful operation.
 *  - LA_INVALID_INPUT_ERROR: Occurs if @p ldr is not correct.
 *  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
 *      there is insufficient memory available.
 *  - LA_MATRIX_FORMAT_ERROR: Occurs if the downdated matrix is not
 *      positive definite.
 */
int la_cholesky_rank1_downdate_cmplx(int n, double complex *r, int ldr,
    double complex *u);

/**
 * Computes the singular value decomposition of a matrix A.  The
 *  SVD is defined as: A = U * S * V**T, where U is an M-by-M orthogonal
 *  matrix, S is an M-by-N diagonal matrix, and V is an N-by-N orthogonal
 *  matrix.
 *
 * @param m The number of rows in the matrix.
 * @param n The number of columns in the matrix.
 * @param a On input, the M-by-N matrix to factor.  The matrix is
 *  overwritten on output.
 * @param lda The leading dimension of matrix A.
 * @param s A MIN(M, N)-element array containing the singular values
 *  of @p a sorted in descending order.
 * @param u An M-by-M matrix where the orthogonal U matrix will be
 *  written.
 * @param ldu The leading dimension of matrix U.
 * @param vt An N-by-N matrix where the transpose of the right 
 *  singular vector matrix V.
 * @param ldv The leading dimension of matrix V.
 *
 * @return An error code.  The following codes are possible.
 *  - LA_NO_ERROR: No error occurred.  Successful operation.
 *  - LA_INVALID_INPUT_ERROR: Occurs if @p lda, @p ldu, or @p ldv is not 
 *      correct.
 *  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
 *      there is insufficient memory available.
 *  - LA_CONVERGENCE_ERROR: Occurs as a warning if the QR iteration process
 *      could not converge to a zero value.
 */
int la_svd(int m, int n, double *a, int lda, double *s, double *u, int ldu,
    double *vt, int ldv);

/**
 * Computes the singular value decomposition of a matrix A.  The
 *  SVD is defined as: A = U * S * V**T, where U is an M-by-M orthogonal
 *  matrix, S is an M-by-N diagonal matrix, and V is an N-by-N orthogonal
 *  matrix.
 *
 * @param m The number of rows in the matrix.
 * @param n The number of columns in the matrix.
 * @param a On input, the M-by-N matrix to factor.  The matrix is
 *  overwritten on output.
 * @param lda The leading dimension of matrix A.
 * @param s A MIN(M, N)-element array containing the singular values
 *  of @p a sorted in descending order.
 * @param u An M-by-M matrix where the orthogonal U matrix will be
 *  written.
 * @param ldu The leading dimension of matrix U.
 * @param vt An N-by-N matrix where the transpose of the right 
 *  singular vector matrix V.
 * @param ldv The leading dimension of matrix V.
 *
 * @return An error code.  The following codes are possible.
 *  - LA_NO_ERROR: No error occurred.  Successful operation.
 *  - LA_INVALID_INPUT_ERROR: Occurs if @p lda, @p ldu, or @p ldv is not 
 *      correct.
 *  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
 *      there is insufficient memory available.
 *  - LA_CONVERGENCE_ERROR: Occurs as a warning if the QR iteration process
 *      could not converge to a zero value.
 */
int la_svd_cmplx(int m, int n, double complex *a, int lda, double *s,
    double complex *u, int ldu, double complex *vt, int ldv);

/**
 * Solves one of the matrix equations: op(A) * X = alpha * B, or
 * X * op(A) = alpha * B, where A is a triangular matrix.
 *
 * @param lside Set to true to solve op(A) * X = alpha * B; else, set to
 *  false to solve X * op(A) = alpha * B.
 * @param upper Set to true if A is an upper triangular matrix; else,
 *  set to false if A is a lower triangular matrix.
 * @param trans Set to true if op(A) = A**T; else, set to false if
 *  op(A) = A.
 * @param nounit Set to true if A is not a unit-diagonal matrix (ones on
 *  every diagonal element); else, set to false if A is a unit-diagonal
 *  matrix.
 * @param m The number of rows in matrix B.
 * @param n The number of columns in matrix B.
 * @param alpha The scalar multiplier to B.
 * @param a If @p lside is true, the M-by-M triangular matrix on which
 *  to operate; else, if @p lside is false, the N-by-N triangular matrix on
 *  which to operate.
 * @param lda The leading dimension of matrix A.
 * @param b On input, the M-by-N right-hand-side.  On output, the
 *  M-by-N solution.
 * @param ldb The leading dimension of matrix B.
 *
 * @return An error code.  The following codes are possible.
 *  - LA_NO_ERROR: No error occurred.  Successful operation.
 *  - LA_INVALID_INPUT_ERROR: Occurs if @p lda, or @p ldb is not correct.
 */
int la_solve_tri_mtx(bool lside, bool upper, bool trans, bool nounit, int m,
    int n, double alpha, const double *a, int lda, double *b, int ldb);

/**
 * Solves one of the matrix equations: op(A) * X = alpha * B, or
 * X * op(A) = alpha * B, where A is a triangular matrix.
 *
 * @param lside Set to true to solve op(A) * X = alpha * B; else, set to
 *  false to solve X * op(A) = alpha * B.
 * @param upper Set to true if A is an upper triangular matrix; else,
 *  set to false if A is a lower triangular matrix.
 * @param trans Set to true if op(A) = A**H; else, set to false if
 *  op(A) = A.
 * @param nounit Set to true if A is not a unit-diagonal matrix (ones on
 *  every diagonal element); else, set to false if A is a unit-diagonal
 *  matrix.
 * @param m The number of rows in matrix B.
 * @param n The number of columns in matrix B.
 * @param alpha The scalar multiplier to B.
 * @param a If @p lside is true, the M-by-M triangular matrix on which
 *  to operate; else, if @p lside is false, the N-by-N triangular matrix on
 *  which to operate.
 * @param lda The leading dimension of matrix A.
 * @param b On input, the M-by-N right-hand-side.  On output, the
 *  M-by-N solution.
 * @param ldb The leading dimension of matrix B.
 *
 * @return An error code.  The following codes are possible.
 *  - LA_NO_ERROR: No error occurred.  Successful operation.
 *  - LA_INVALID_INPUT_ERROR: Occurs if @p lda, or @p ldb is not correct.
 */
int la_solve_tri_mtx_cmplx(bool lside, bool upper, bool trans, bool nounit,
    int m, int n, double complex alpha, const double complex *a, int lda,
    double complex *b, int ldb);

/**
 * Solves a system of LU-factored equations.
 *
 * @param m The number of rows in matrix B.
 * @param n The number of columns in matrix B.
 * @param a The M-by-M LU factored matrix.
 * @param lda The leading dimension of matrix A.
 * @param ipvt The M-element pivot array from the LU factorization.
 * @param b On input, the M-by-N right-hand-side.  On output, the
 *  M-by-N solution.
 * @param ldb The leading dimension of matrix B.
 *
 * @return An error code.  The following codes are possible.
 *  - LA_NO_ERROR: No error occurred.  Successful operation.
 *  - LA_INVALID_INPUT_ERROR: Occurs if @p lda, or @p ldb is not correct.
 */
int la_solve_lu(int m, int n, const double *a, int lda, const int *ipvt,
    double *b, int ldb);

/**
 * Solves a system of LU-factored equations.
 *
 * @param m The number of rows in matrix B.
 * @param n The number of columns in matrix B.
 * @param a The M-by-M LU factored matrix.
 * @param lda The leading dimension of matrix A.
 * @param ipvt The M-element pivot array from the LU factorization.
 * @param b On input, the M-by-N right-hand-side.  On output, the
 *  M-by-N solution.
 * @param ldb The leading dimension of matrix B.
 *
 * @return An error code.  The following codes are possible.
 *  - LA_NO_ERROR: No error occurred.  Successful operation.
 *  - LA_INVALID_INPUT_ERROR: Occurs if @p lda, or @p ldb is not correct.
 */
int la_solve_lu_cmplx(int m, int n, const double complex *a, int lda,
    const int *ipvt, double complex *b, int ldb);

/**
 * Solves a system of M QR-factored equations of N unknowns where
 * M >= N.
 *
 * @param m The number of equations (rows in matrix A).
 * @param n The number of unknowns (columns in matrix A).
 * @param k The number of columns in the right-hand-side matrix.
 * @param a On input, the M-by-N QR factored matrix as returned by
 *  qr_factor.  On output, the contents of this matrix are restored.
 * @param lda The leading dimension of matrix A.
 * @param tau A MIN(M, N)-element array containing the scalar factors of
 *  the elementary reflectors as returned by qr_factor.
 * @param b On input, the M-by-K right-hand-side matrix.  On output,
 *  the first N rows are overwritten by the solution matrix X.
 * @param ldb The leading dimension of matrix B.
 *
 * @return An error code.  The following codes are possible.
 *  - LA_NO_ERROR: No error occurred.  Successful operation.
 *  - LA_INVALID_INPUT_ERROR: Occurs if @p lda, or @p ldb is not correct, or
 *      if @p m is less than @p n.
 *  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
 *      there is insufficient memory available.
 */
int la_solve_qr(int m, int n, int k, double *a, int lda, const double *tau,
    double *b, int ldb);

/**
 * Solves a system of M QR-factored equations of N unknowns where
 * M >= N.
 *
 * @param m The number of equations (rows in matrix A).
 * @param n The number of unknowns (columns in matrix A).
 * @param k The number of columns in the right-hand-side matrix.
 * @param a On input, the M-by-N QR factored matrix as returned by
 *  qr_factor.  On output, the contents of this matrix are restored.
 * @param lda The leading dimension of matrix A.
 * @param tau A MIN(M, N)-element array containing the scalar factors of
 *  the elementary reflectors as returned by qr_factor.
 * @param b On input, the M-by-K right-hand-side matrix.  On output,
 *  the first N rows are overwritten by the solution matrix X.
 * @param ldb The leading dimension of matrix B.
 *
 * @return An error code.  The following codes are possible.
 *  - LA_NO_ERROR: No error occurred.  Successful operation.
 *  - LA_INVALID_INPUT_ERROR: Occurs if @p lda, or @p ldb is not correct, or
 *      if @p m is less than @p n.
 *  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
 *      there is insufficient memory available.
 */
int la_solve_qr_cmplx(int m, int n, int k, double complex *a, int lda, 
    const double complex *tau, double complex *b, int ldb);

/**
 * Solves a system of M QR-factored equations of N unknowns.
 *
 * @param m The number of equations (rows in matrix A).
 * @param n The number of unknowns (columns in matrix A).
 * @param k The number of columns in the right-hand-side matrix.
 * @param a On input, the M-by-N QR factored matrix as returned by
 *  qr_factor.  On output, the contents of this matrix are restored.
 * @param lda The leading dimension of matrix A.
 * @param tau A MIN(M, N)-element array containing the scalar factors of
 *  the elementary reflectors as returned by qr_factor.
 * @param jpvt The N-element array that was used to track the column
 *  pivoting operations in the QR factorization.
 * @param b On input, the M-by-K right-hand-side matrix.  On output,
 *  the first N rows are overwritten by the solution matrix X.
 * @param ldb The leading dimension of matrix B.
 *
 * @return An error code.  The following codes are possible.
 *  - LA_NO_ERROR: No error occurred.  Successful operation.
 *  - LA_INVALID_INPUT_ERROR: Occurs if @p lda, or @p ldb is not correct.
 *  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
 *      there is insufficient memory available.
 */
int la_solve_qr_pvt(int m, int n, int k, double *a, int lda, const double *tau,
    const int *jpvt, double *b, int ldb);

/**
 * Solves a system of M QR-factored equations of N unknowns.
 *
 * @param m The number of equations (rows in matrix A).
 * @param n The number of unknowns (columns in matrix A).
 * @param k The number of columns in the right-hand-side matrix.
 * @param a On input, the M-by-N QR factored matrix as returned by
 *  qr_factor.  On output, the contents of this matrix are restored.
 * @param lda The leading dimension of matrix A.
 * @param tau A MIN(M, N)-element array containing the scalar factors of
 *  the elementary reflectors as returned by qr_factor.
 * @param jpvt The N-element array that was used to track the column
 *  pivoting operations in the QR factorization.
 * @param b On input, the M-by-K right-hand-side matrix.  On output,
 *  the first N rows are overwritten by the solution matrix X.
 * @param ldb The leading dimension of matrix B.
 *
 * @return An error code.  The following codes are possible.
 *  - LA_NO_ERROR: No error occurred.  Successful operation.
 *  - LA_INVALID_INPUT_ERROR: Occurs if @p lda, or @p ldb is not correct.
 *  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
 *      there is insufficient memory available.
 */
int la_solve_qr_cmplx_pvt(int m, int n, int k, double complex *a, int lda,
    const double complex *tau, const int *jpvt, double complex *b, int ldb);

/**
 * Solves a system of Cholesky factored equations.
 *
 * @param upper Set to true if the original matrix A was factored such
 *  that A = U**T * U; else, set to false if the factorization of A was
 *  A = L**T * L.
 * @param m The number of rows in matrix B.
 * @param n The number of columns in matrix B.
 * @param a The M-by-M Cholesky factored matrix.
 * @param lda The leading dimension of matrix A.
 * @param b On input, the M-by-N right-hand-side matrix B.  On
 *  output, the M-by-N solution matrix X.
 * @param ldb The leading dimension of matrix B.
 *
 * @return An error code.  The following codes are possible.
 *  - LA_NO_ERROR: No error occurred.  Successful operation.
 *  - LA_INVALID_INPUT_ERROR: Occurs if @p lda, or @p ldb is not correct.
 */
int la_solve_cholesky(bool upper, int m, int n, const double *a, int lda,
    double *b, int ldb);

/**
 * Solves a system of Cholesky factored equations.
 *
 * @param upper Set to true if the original matrix A was factored such
 *  that A = U**T * U; else, set to false if the factorization of A was
 *  A = L**T * L.
 * @param m The number of rows in matrix B.
 * @param n The number of columns in matrix B.
 * @param a The M-by-M Cholesky factored matrix.
 * @param lda The leading dimension of matrix A.
 * @param b On input, the M-by-N right-hand-side matrix B.  On
 *  output, the M-by-N solution matrix X.
 * @param ldb The leading dimension of matrix B.
 *
 * @return An error code.  The following codes are possible.
 *  - LA_NO_ERROR: No error occurred.  Successful operation.
 *  - LA_INVALID_INPUT_ERROR: Occurs if @p lda, or @p ldb is not correct.
 */
int la_solve_cholesky_cmplx(bool upper, int m, int n, const double complex *a,
    int lda, double complex *b, int ldb);

/**
 * Solves the overdetermined or underdetermined system (A*X = B) of
 * M equations of N unknowns using a QR or LQ factorization of the matrix A.
 * Notice, it is assumed that matrix A has full rank.
 *
 * @param m The number of equations (rows in matrix A).
 * @param n The number of unknowns (columns in matrix A).
 * @param k The number of columns in the right-hand-side matrix.
 * @param a On input, the M-by-N matrix A.  On output, if M >= N,
 *  the QR factorization of A in the form as output by qr_factor; else,
 *  if M < N, the LQ factorization of A.
 * @param lda The leading dimension of matrix A.
 * @param b If M >= N, the M-by-NRHS matrix B.  On output, the first
 *  N rows contain the N-by-NRHS solution matrix X.  If M < N, an
 *  N-by-NRHS matrix with the first M rows containing the matrix B.  On
 *  output, the N-by-NRHS solution matrix X.
 * @param ldb The leading dimension of matrix B.
 *
 * @return An error code.  The following codes are possible.
 *  - LA_NO_ERROR: No error occurred.  Successful operation.
 *  - LA_INVALID_INPUT_ERROR: Occurs if @p lda, or @p ldb is not correct.
 *  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
 *      there is insufficient memory available.
 *  - LA_INVALID_OPERATION_ERROR: Occurs if @p a is not of full rank.
 */
int la_solve_least_squares(int m, int n, int k, double *a, int lda, double *b, 
    int ldb);

/**
 * Solves the overdetermined or underdetermined system (A*X = B) of
 * M equations of N unknowns using a QR or LQ factorization of the matrix A.
 * Notice, it is assumed that matrix A has full rank.
 *
 * @param m The number of equations (rows in matrix A).
 * @param n The number of unknowns (columns in matrix A).
 * @param k The number of columns in the right-hand-side matrix.
 * @param a On input, the M-by-N matrix A.  On output, if M >= N,
 *  the QR factorization of A in the form as output by qr_factor; else,
 *  if M < N, the LQ factorization of A.
 * @param lda The leading dimension of matrix A.
 * @param b If M >= N, the M-by-NRHS matrix B.  On output, the first
 *  N rows contain the N-by-NRHS solution matrix X.  If M < N, an
 *  N-by-NRHS matrix with the first M rows containing the matrix B.  On
 *  output, the N-by-NRHS solution matrix X.
 * @param ldb The leading dimension of matrix B.
 *
 * @return An error code.  The following codes are possible.
 *  - LA_NO_ERROR: No error occurred.  Successful operation.
 *  - LA_INVALID_INPUT_ERROR: Occurs if @p lda, or @p ldb is not correct.
 *  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
 *      there is insufficient memory available.
 *  - LA_INVALID_OPERATION_ERROR: Occurs if @p a is not of full rank.
 */
int la_solve_least_squares_cmplx(int m, int n, int k, double complex *a,
    int lda, double complex *b, int ldb);

/**
 * Computes the inverse of a square matrix.
 *
 * @param n The dimension of matrix A.
 * @param a On input, the N-by-N matrix to invert.  On output, the
 *  inverted matrix.
 * @param lda The leading dimension of matrix A.
 *
 * @return An error code.  The following codes are possible.
 *  - LA_NO_ERROR: No error occurred.  Successful operation.
 *  - LA_INVALID_INPUT_ERROR: Occurs if @p lda is not correct.
 *  - LA_SINGULAR_MATRIX_ERROR: Occurs if the input matrix is singular.
 */
int la_inverse(int n, double *a, int lda);

/**
 * Computes the inverse of a square matrix.
 *
 * @param n The dimension of matrix A.
 * @param a On input, the N-by-N matrix to invert.  On output, the
 *  inverted matrix.
 * @param lda The leading dimension of matrix A.
 *
 * @return An error code.  The following codes are possible.
 *  - LA_NO_ERROR: No error occurred.  Successful operation.
 *  - LA_INVALID_INPUT_ERROR: Occurs if @p lda is not correct.
 *  - LA_SINGULAR_MATRIX_ERROR: Occurs if the input matrix is singular.
 */
int la_inverse_cmplx(int n, double complex *a, int lda);

/**
 * Computes the Moore-Penrose pseudo-inverse of an M-by-N matrix by
 * means of singular value decomposition.
 *
 * @param m The number of rows in the matrix.
 * @parma n The number of columns in the matrix.
 * @param a On input, the M-by-N matrix to invert.  The matrix is
 *  overwritten on output.
 * @param lda The leading dimension of matrix A.
 * @param ainv The N-by-M matrix where the pseudo-inverse of @p a
 *  will be written.
 * @param ldai The leading dimension of matrix AINV.
 *
 * @return An error code.  The following codes are possible.
 *  - LA_NO_ERROR: No error occurred.  Successful operation.
 *  - LA_INVALID_INPUT_ERROR: Occurs if @p lda, or @p ldai is not correct.
 */
int la_pinverse(int m, int n, double *a, int lda, double *ainv, int ldai);

/**
 * Computes the Moore-Penrose pseudo-inverse of an M-by-N matrix by
 * means of singular value decomposition.
 *
 * @param m The number of rows in the matrix.
 * @parma n The number of columns in the matrix.
 * @param a On input, the M-by-N matrix to invert.  The matrix is
 *  overwritten on output.
 * @param lda The leading dimension of matrix A.
 * @param ainv The N-by-M matrix where the pseudo-inverse of @p a
 *  will be written.
 * @param ldai The leading dimension of matrix AINV.
 *
 * @return An error code.  The following codes are possible.
 *  - LA_NO_ERROR: No error occurred.  Successful operation.
 *  - LA_INVALID_INPUT_ERROR: Occurs if @p lda, or @p ldai is not correct.
 */
int la_pinverse_cmplx(int m, int n, double complex *a, int lda,
    double complex *ainv, int ldai);

/**
 * Computes the eigenvalues, and optionally the eigenvectors of a
 * real, symmetric matrix.
 *
 * @param vecs Set to true to compute the eigenvectors as well as the
 *  eigenvalues; else, set to false to just compute the eigenvalues.
 * @param n The dimension of the matrix.
 * @param a On input, the N-by-N symmetric matrix on which to
 *  operate.  On output, and if @p vecs is set to true, the matrix will
 *  contain the eigenvectors (one per column) corresponding to each
 *  eigenvalue in @p vals.  If @p vecs is set to false, the lower triangular
 *  portion of the matrix is overwritten.
 * @param lda The leading dimension of matrix A.
 * @param vals An N-element array that will contain the eigenvalues
 *  sorted into ascending order.
 *
 * @return An error code.  The following codes are possible.
 *  - LA_NO_ERROR: No error occurred.  Successful operation.
 *  - LA_INVALID_INPUT_ERROR: Occurs if @p lda is not correct.
 *  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
 *      there is insufficient memory available.
 *  - LA_CONVERGENCE_ERROR: Occurs if the algorithm failed to converge.
 */
int la_eigen_symm(bool vecs, int n, double *a, int lda, double *vals);

/**
 * Computes the eigenvalues, and optionally the right eigenvectors of
 * a square matrix.
 *
 * @param vecs Set to true to compute the eigenvectors as well as the
 *  eigenvalues; else, set to false to just compute the eigenvalues.
 * @param n The dimension of the matrix.
 * @param a On input, the N-by-N matrix on which to operate.  On
 *  output, the contents of this matrix are overwritten.
 * @param lda The leading dimension of matrix A.
 * @param vals An N-element array containing the eigenvalues of the
 *  matrix.  The eigenvalues are not sorted.
 * @param v An N-by-N matrix where the right eigenvectors will be
 *  written (one per column).
 *
 * @return An error code.  The following codes are possible.
 *  - LA_NO_ERROR: No error occurred.  Successful operation.
 *  - LA_INVALID_INPUT_ERROR: Occurs if @p lda or @p ldv is not correct.
 *  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
 *      there is insufficient memory available.
 *  - LA_CONVERGENCE_ERROR: Occurs if the algorithm failed to converge.
 */
int la_eigen_asymm(bool vecs, int n, double *a, int lda,
    double complex *vals, double complex *v, int ldv);

/**
 * Computes the eigenvalues, and optionally the right eigenvectors of
 * a square matrix assuming the structure of the eigenvalue problem is
 * A*X = lambda*B*X.
 *
 * @param vecs Set to true to compute the eigenvectors as well as the
 *  eigenvalues; else, set to false to just compute the eigenvalues.
 * @param n The dimension of the matrix.
 * @param a On input, the N-by-N matrix A.  On output, the contents
 *  of this matrix are overwritten.
 * @param lda The leading dimension of matrix A.
 * @param b On input, the N-by-N matrix B.  On output, the contents
 *  of this matrix are overwritten.
 * @param ldb The leading dimension of matrix B.
 * @param alpha An N-element array a factor of the eigenvalues.  The 
 *  eigenvalues must be computed as ALPHA / BETA.  This however, is not as
 *  trivial as it seems as it is entirely possible, and likely, that
 *  ALPHA / BETA can overflow or underflow.  With that said, the values in
 *  ALPHA will always be less than and usually comparable with the NORM(A).
 * @param beta An N-element array that contains the denominator used to 
 *  determine the eigenvalues as ALPHA / BETA.  If used, the values in this 
 *  array will always be less than and usually comparable with the NORM(B).
 * @param v An N-by-N matrix where the right eigenvectors will be
 *  written (one per column).
 *
 * @return An error code.  The following codes are possible.
 *  - LA_NO_ERROR: No error occurred.  Successful operation.
 *  - LA_INVALID_INPUT_ERROR: Occurs if @p lda or @p ldv is not correct.
 *  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
 *      there is insufficient memory available.
 *  - LA_CONVERGENCE_ERROR: Occurs if the algorithm failed to converge.
 */
int la_eigen_gen(bool vecs, int n, double *a, int lda, double *b, int ldb,
    double complex *alpha, double *beta, double complex *v, int ldv);

/**
 * Computes the eigenvalues, and optionally the right eigenvectors of
 * a square matrix.
 *
 * @param vecs Set to true to compute the eigenvectors as well as the
 *  eigenvalues; else, set to false to just compute the eigenvalues.
 * @param n The dimension of the matrix.
 * @param a On input, the N-by-N matrix on which to operate.  On
 *  output, the contents of this matrix are overwritten.
 * @param lda The leading dimension of matrix A.
 * @param vals An N-element array containing the eigenvalues of the
 *  matrix.  The eigenvalues are not sorted.
 * @param v An N-by-N matrix where the right eigenvectors will be
 *  written (one per column).
 *
 * @return An error code.  The following codes are possible.
 *  - LA_NO_ERROR: No error occurred.  Successful operation.
 *  - LA_INVALID_INPUT_ERROR: Occurs if @p lda or @p ldv is not correct.
 *  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
 *      there is insufficient memory available.
 *  - LA_CONVERGENCE_ERROR: Occurs if the algorithm failed to converge.
 */
int la_eigen_cmplx(bool vecs, int n, double complex *a, int lda,
    double complex *vals, double complex *v, int ldv);

/**
 * A sorting routine specifically tailored for sorting of eigenvalues
 * and their associated eigenvectors using a quick-sort approach.
 *
 * @param ascend
 * @param n The number of eigenvalues.
 * @param vals On input, an N-element array containing the
 *  eigenvalues.  On output, the sorted eigenvalues.
 * @param vecs On input, an N-by-N matrix containing the
 *  eigenvectors associated with @p vals (one vector per column).  On
 *  output, the sorted eigenvector matrix.
 * @param ldv The leading dimension of @p vecs.
 *
 * @return An error code.  The following codes are possible.
 *  - LA_NO_ERROR: No error occurred.  Successful operation.
 *  - LA_INVALID_INPUT_ERROR: Occurs if @p ldv is not correct.
 *  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
 *      there is insufficient memory available.
 */
int la_sort_eigen(bool ascend, int n, double *vals, double *vecs, int ldv);

/**
 * A sorting routine specifically tailored for sorting of eigenvalues
 * and their associated eigenvectors using a quick-sort approach.
 *
 * @param ascend
 * @param n The number of eigenvalues.
 * @param vals On input, an N-element array containing the
 *  eigenvalues.  On output, the sorted eigenvalues.
 * @param vecs On input, an N-by-N matrix containing the
 *  eigenvectors associated with @p vals (one vector per column).  On
 *  output, the sorted eigenvector matrix.
 * @param ldv The leading dimension of @p vecs.
 *
 * @return An error code.  The following codes are possible.
 *  - LA_NO_ERROR: No error occurred.  Successful operation.
 *  - LA_INVALID_INPUT_ERROR: Occurs if @p ldv is not correct.
 *  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
 *      there is insufficient memory available.
 */
int la_sort_eigen_cmplx(bool ascend, int n, double complex *vals,
    double complex *vecs, int ldv);

#ifdef __cplusplus
}
#endif  // __cplusplus
#endif  // LINALG_H_
