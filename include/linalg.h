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
 * @param opb Set to TLA_RANSPOSE to compute op(B) as a direct transpose of B,
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
 * Computes the rank of a matrix.
 *
 * @param m The number of rows in the matrix.
 * @param n The number of columns in the matrix.
 * @param a The M-by-N matrix.  The matrix is overwritten as part of this
 *  operation.
 * @param lda The leading dimension of matrix A.
 * @param[out] rnk The rank of @p a.
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
 * @param[out] rnk The rank of @p a.
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
 * @param[out] d The determinant of @p a.
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
 * @param[out] d The determinant of @p a.
 *
 * @return An error code.  The following codes are possible.
 *  - LA_NO_ERROR: No error occurred.  Successful operation.
 *  - LA_ARRAY_SIZE_ERROR: Occurs if any of the input arrays are not sized
 *      appropriately.
 *  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
 *      there is insufficient memory available.
 */
int la_det_cmplx(int n, double complex *a, int lda, double complex *d);

#ifdef __cplusplus
}
#endif  // __cplusplus
#endif  // LINALG_H_
