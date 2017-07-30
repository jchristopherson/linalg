! linalg_c_binding.f90

!> @brief \b linalg_c_binding
!!
!! @par Purpose
!! Provides a C friendly interface to the LINALG library.
module linalg_c_binding
    use, intrinsic :: iso_c_binding
    use linalg_constants
    use linalg_core
    use linalg_factor
    use linalg_solve
    use linalg_eigen
    use ferror, only : errors
    use ferror_c_binding, only : errorhandler, get_errorhandler, &
        update_errorhandler
contains
! ******************************************************************************
! LINALG_CORE ROUTINES
! ------------------------------------------------------------------------------
    !> @brief Performs the matrix operation: 
    !!  C = alpha * op(A) * op(B) + beta * C.
    !!
    !! @param[in] transa Set to true if op(A) == A**T; else, set to false if
    !!  op(A) == A.
    !! @param[in] transb Set to true if op(B) == B**T; else, set to false if
    !!  op(B) == B.
    !! @param[in] m The number of rows in matrix C, and the number of rows in
    !!  matrix op(A).
    !! @param[in] n The number of columns in matrix C, and the number of columns
    !!  in matrix op(B).
    !! @param[in] k The number of columns in matrix op(A), and the number of
    !!  rows in the matrix op(B).
    !! @param[in] alpha The scalar multiplier to matrix A.
    !! @param[in] a The M-by-K matrix A.
    !! @param[in] lda The leading dimension of matrix A.  If @p transa is true, 
    !!  this value must be at least MAX(1, K); else, if @p transa is false, this
    !!  value must be at least MAX(1, M).
    !! @param[in] b The K-by-N matrix B.
    !! @param[in] ldb The leading dimension of matrix B.  If @p transb is true,
    !!  this value must be at least MAX(1, N); else, if @p transb is false, this
    !!  value must be at least MAX(1, K).
    !! @param[in] beta The scalar multiplier to matrix C.
    !! @param[in,out] c The M-by-N matrix C.
    subroutine mtx_mult_c(transa, transb, m, n, k, alpha, a, lda, b, ldb, &
            beta, c) bind(C, name = "mtx_mult_")
        ! Arguments
        logical(c_bool), intent(in), value :: transa, transb
        integer(i32), intent(in), value :: m, n, k, lda, ldb
        real(dp), intent(in), value :: alpha, beta
        real(dp), intent(in) :: a(lda,*), b(ldb,*)
        real(dp), intent(inout) :: c(m,n)

        ! Process
        character :: ta, tb
        if (transa) then
            ta = 'T'
        else
            ta = 'N'
        end if
        if (transb) then
            tb = 'T'
        else
            tb = 'N'
        end if
        call DGEMM(ta, tb, m, n, k, alpha, a, lda, b, ldb, beta, c, m)
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Computes the matrix operation: C = alpha * A * op(B) + beta * C.
    !!
    !! @param[in] trans Set to true if op(B) == B**T; else, set to false if
    !!  op(B) == B.
    !! @param[in] m The number of rows in matrix C.
    !! @param[in] n The number of columns in matrix C.
    !! @param[in] alpha The scalar multiplier to matrix A.
    !! @param[in] na The length of @p a.
    !! @param[in] a A MIN(M,P)-element array containing the diagonal elements 
    !!  of matrix A.
    !! @param[in] mb The number of rows in matrix B.
    !! @param[in] nb The number of columns in matrix B.
    !! @param[in] b The LDB-by-TDB matrix B where (LDB = leading dimension of B,
    !!  and TDB = trailing dimension of B):
    !!  - @p trans == true: LDB = N, TDB = P
    !!  - @p trans == false: LDB = P, TDB = N
    !! @param[in] beta The scalar multiplier to matrix C.
    !! @param[in,out] c The M-by-N matrix C.
    !! @param[in,out] err The errorhandler object.  If no error handling is
    !!  desired, simply pass NULL, and errors will be dealt with by the default
    !!  internal error handler.  Possible errors that may be encountered are as
    !!  follows.
    !!  - LA_ARRAY_SIZE_ERROR: Occurs if any of the input arrays are not sized
    !!      appropriately.
    subroutine diag_mtx_mult_c(trans, m, n, alpha, na, a, mb, nb, b, beta, &
            c, err) bind(C, name = "diag_mtx_mult_")
        ! Arguments
        logical(c_bool), intent(in), value :: trans
        integer(i32), intent(in), value :: m, n, na, mb, nb
        real(dp), intent(in), value :: alpha, beta
        real(dp), intent(in) :: a(na), b(mb, nb)
        real(dp), intent(inout) :: c(m, n)
        type(errorhandler), intent(inout) :: err

        ! Local Variables
        class(errors), allocatable :: eptr

        ! Process
        call get_errorhandler(err, eptr)
        if (allocated(eptr)) then
            call diag_mtx_mult(.true., logical(trans), alpha, a, b, beta, c, &
                eptr)
            call update_errorhandler(eptr, err)
        else
            call diag_mtx_mult(.true., logical(trans), alpha, a, b, beta, c)
        end if
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Computes the matrix operation: C = alpha * A * op(B) + beta * C,
    !!  where A and C are complex-valued.
    !!
    !! @param[in] trans Set to true if op(B) == B**T; else, set to false if
    !!  op(B) == B.
    !! @param[in] m The number of rows in matrix C.
    !! @param[in] n The number of columns in matrix C.
    !! @param[in] alpha The scalar multiplier to matrix A.
    !! @param[in] na The length of @p a.
    !! @param[in] a A MIN(M,P)-element array containing the diagonal elements 
    !!  of matrix A.
    !! @param[in] mb The number of rows in matrix B.
    !! @param[in] nb The number of columns in matrix B.
    !! @param[in] b The LDB-by-TDB matrix B where (LDB = leading dimension of B,
    !!  and TDB = trailing dimension of B):
    !!  - @p trans == true: LDB = N, TDB = P
    !!  - @p trans == false: LDB = P, TDB = N
    !! @param[in] beta The scalar multiplier to matrix C.
    !! @param[in,out] c THe M-by-N matrix C.
    !! @param[in,out] err The errorhandler object.  If no error handling is
    !!  desired, simply pass NULL, and errors will be dealt with by the default
    !!  internal error handler.  Possible errors that may be encountered are as
    !!  follows.
    !!  - LA_ARRAY_SIZE_ERROR: Occurs if any of the input arrays are not sized
    !!      appropriately.
    subroutine diag_mtx_mult_cmplx_c(trans, m, n, alpha, na, a, mb, nb, b, &
            beta, c, err) bind(C, name = "diag_mtx_mult_cmplx_")
        ! Arguments
        logical(c_bool), intent(in), value :: trans
        integer(i32), intent(in), value :: m, n, na, mb, nb
        real(dp), intent(in), value :: alpha, beta
        complex(dp), intent(in) :: a(na)
        real(dp), intent(in) :: b(mb, nb)
        complex(dp), intent(inout) :: c(m, n)
        type(errorhandler), intent(inout) :: err

        ! Local Variables
        class(errors), allocatable :: eptr

        ! Process
        call get_errorhandler(err, eptr)
        if (allocated(eptr)) then
            call diag_mtx_mult(.true., logical(trans), alpha, a, b, beta, c, &
                eptr)
            call update_errorhandler(eptr, err)
        else
            call diag_mtx_mult(.true., logical(trans), alpha, a, b, beta, c)
        end if
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Performs the rank-1 update to matrix A such that:
    !! A = alpha * X * Y**T + A, where A is an M-by-N matrix, alpha is a scalar,
    !! X is an M-element array, and N is an N-element array.
    !!
    !! @param[in] m The number of elements in @p x, and the number of rows in
    !!  matrix @p a.
    !! @param[in] n The number of elements in @p y, and the number of columns in
    !!  matrix @p a.
    !! @param[in] alpha The scalar multiplier.
    !! @param[in] x An M-element array.
    !! @param[in] y An N-element array.
    !! @param[in,out] a On input, the M-by-N matrix to update.  On output, the
    !!  updated M-by-N matrix.
    !!
    !! @par Notes
    !! This routine is based upon the BLAS routine DGER.
    subroutine rank1_update_c(m, n, alpha, x, y, a) &
            bind(C, name = "rank1_update_")
        ! Arguments
        real(dp), intent(in), value :: alpha
        integer(i32), intent(in), value :: m, n
        real(dp), intent(in) :: x(m), y(n)
        real(dp), intent(inout) :: a(m,n)

        ! Process
        call rank1_update(alpha, x, y, a)
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Computes the trace of a matrix (the sum of the main diagonal
    !! elements).
    !!
    !! @param[in] m The number of rows in the matrix.
    !! @param[in] n The number of columns in the matrix.
    !! @param[in] x The matrix on which to operate.
    !!
    !! @return The trace of @p x.
    pure function trace_c(m, n, x) result(y) bind(C, name = "trace_")
        ! Arguments
        integer(i32), intent(in), value :: m, n
        real(dp), intent(in) :: x(m,n)
        real(dp) :: y

        ! Process
        y = trace(x)
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the rank of a matrix.
    !!
    !! @param[in] m The number of rows in the matrix.
    !! @param[in] n The number of columns in the matrix.
    !! @param[in,out] a On input, the M-by-N matrix of interest.  On output, the
    !!  contents of the matrix are overwritten.
    !! @param[in,out] err The errorhandler object.  If no error handling is
    !!  desired, simply pass NULL, and errors will be dealt with by the default
    !!  internal error handler.  Possible errors that may be encountered are as
    !!  follows.
    !!  - LA_ARRAY_SIZE_ERROR: Occurs if any of the input arrays are not sized
    !!      appropriately.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
    !!      there is insufficient memory available.
    !!  - LA_CONVERGENCE_ERROR: Occurs as a warning if the QR iteration process
    !!      could not converge to a zero value.
    function mtx_rank_c(m, n, a, err) result(rnk) bind(C, name = "mtx_rank_")
        ! Arguments
        integer(i32), intent(in), value :: m, n
        real(dp), intent(inout) :: a(m,n)
        type(errorhandler), intent(inout) :: err
        integer(i32) :: rnk

        ! Local Variables
        class(errors), allocatable :: eptr

        ! Process
        call get_errorhandler(err, eptr)
        if (allocated(eptr)) then
            rnk = mtx_rank(a, err = eptr)
            call update_errorhandler(eptr, err)
        else
            rnk = mtx_rank(a)
        end if
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the determinant of a square matrix.
    !!
    !! @param[in] n The dimension of the matrix.
    !! @param[in,out] a On input, the N-by-N matrix on which to operate.  On
    !! output the contents are overwritten by the LU factorization of the
    !! original matrix.
    !! @param[in,out] err The errorhandler object.  If no error handling is
    !!  desired, simply pass NULL, and errors will be dealt with by the default
    !!  internal error handler.  Possible errors that may be encountered are as
    !!  follows.
    !!  - LA_ARRAY_SIZE_ERROR: Occurs if the input matrix is not square.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
    !!      there is insufficient memory available.
    function det_c(n, a, err) result(x) bind(C, name = "det_")
        ! Arguments
        integer(i32), intent(in), value :: n
        real(dp), intent(inout) :: a(n,n)
        type(errorhandler), intent(inout) :: err
        real(dp) :: x

        ! Local Variables
        class(errors), allocatable :: eptr

        ! Process
        call get_errorhandler(err, eptr)
        if (allocated(eptr)) then
            x = det(a, err = eptr)
            call update_errorhandler(eptr, err)
        else
            x = det(a)
        end if
    end function

! ------------------------------------------------------------------------------
    !> @brief Swaps the contents of two arrays.
    !!
    !! @param[in] n The number of elements either array.
    !! @param[in,out] x One of the N-element arrays.
    !! @param[in,out] y The other N-element array.
    subroutine swap_c(n, x, y) bind(C, name = "swap_")
        ! Arguments
        integer(i32), intent(in), value :: n
        real(dp), intent(inout) :: x(n), y(n)

        ! Process
        call swap(x, y)
    end subroutine

! ------------------------------------------------------------------------------
    !@ brief Computes the triangular matrix operation: 
    !! B = alpha * A**T * A + beta * B, or B = alpha * A * A**T + beta * B, 
    !! where A is a triangular matrix.
    !!
    !! @param[in] upper Set to true if matrix A is upper triangular, and 
    !!  B = alpha * A**T * A + beta * B is to be calculated; else, set to false
    !!  if A is lower triangular, and B = alpha * A * A**T + beta * B is to
    !!  be computed.
    !! @param[in] n The size of the matrix.
    !! @param[in] alpha A scalar multiplier.
    !! @param[in] a The N-by-N triangular matrix.  Notice, if @p upper is true
    !!  only the upper triangular portion of this matrix is referenced; else,
    !!  if @p upper is false, only the lower triangular portion of this matrix
    !!  is referenced.
    !! @param[in] beta A scalar multiplier.
    !! @param[in,out] b On input, the N-by-N matrix B.  On output, the N-by-N
    !!  solution matrix.
    !! @param[in,out] err The errorhandler object.  If no error handling is
    !!  desired, simply pass NULL, and errors will be dealt with by the default
    !!  internal error handler.  Possible errors that may be encountered are as
    !!  follows.
    !!  - LA_ARRAY_SIZE_ERROR: Occurs if any of the input arrays are not sized
    !!      appropriately.
    subroutine tri_mtx_mult_c(upper, n, alpha, a, beta, b, err) &
            bind(C, name = "tri_mtx_mult_")
        ! Arguments
        logical(c_bool), intent(in), value :: upper
        integer(i32), intent(in), value :: n
        real(dp), intent(in), value :: alpha, beta
        real(dp), intent(in) :: a(n,n)
        real(dp), intent(inout) :: b(n,n)
        type(errorhandler), intent(inout) :: err

        ! Local Variables
        class(errors), allocatable :: eptr

        ! Process
        call get_errorhandler(err, eptr)
        if (allocated(eptr)) then
            call tri_mtx_mult(logical(upper), alpha, a, beta, b, eptr)
            call update_errorhandler(eptr, err)
        else
            call tri_mtx_mult(logical(upper), alpha, a, beta, b)
        end if
    end subroutine

! ******************************************************************************
! LINALG_FACTOR ROUTINES
! ------------------------------------------------------------------------------
    !> @brief Computes the LU factorization of an M-by-N matrix.
    !!
    !! @param[in] m The number of rows in the matrix.
    !! @param[in] n The number of columns in the matrix.
    !! @param[in,out] a On input, the M-by-N matrix on which to operate.  On
    !! output, the LU factored matrix in the form [L\\U] where the unit diagonal
    !! elements of L are not stored.
    !! @param[out] ipvt An MIN(M, N)-element array used to track row-pivot
    !!  operations.  The array stored pivot information such that row I is
    !!  interchanged with row IPVT(I).
    !! @param[in,out] err The errorhandler object.  If no error handling is
    !!  desired, simply pass NULL, and errors will be dealt with by the default
    !!  internal error handler.  Possible errors that may be encountered are as
    !!  follows.
    !!  - LA_ARRAY_SIZE_ERROR: Occurs if the pivot array is not sized 
    !!      appropriately.
    !!  - LA_SINGULAR_MATRIX_ERROR: Occurs as a warning if @p a is found to be
    !!      singular.
    subroutine lu_factor_c(m, n, a, ipvt, err) bind(C, name = "lu_factor_")
        ! Arguments
        integer(i32), intent(in), value :: m, n
        real(dp), intent(inout) :: a(m,n)
        integer(i32), intent(out) :: ipvt(min(m,n))
        type(errorhandler), intent(inout) :: err

        ! Local Variables
        class(errors), allocatable :: eptr

        ! Process
        call get_errorhandler(err, eptr)
        if (allocated(eptr)) then
            call lu_factor(a, ipvt, eptr)
            call update_errorhandler(eptr, err)
        else
            call lu_factor(a, ipvt)
        end if
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Extracts the L, U, and P matrices from the output of the
    !! lu_factor routine.
    !!
    !! @param[in] n The dimension of the original matrix.
    !! @param[in,out] lu On input, the N-by-N matrix as output by
    !!  lu_factor.  On output, the N-by-N lower triangular matrix L.
    !! @param[in] ipvt The N-element pivot array as output by
    !!  lu_factor.
    !! @param[out] u An N-by-N matrix where the U matrix will be written.
    !! @param[out] p An N-by-N matrix where the row permutation matrix will be
    !!  written.
    !!
    !! @par Remarks
    !! This routine allows extraction of the actual "L", "U", and "P" matrices
    !! of the decomposition.  To use these matrices to solve the system A*X = B,
    !! the following approach is used.
    !!
    !! 1. First, solve the linear system: L*Y = P*B for Y.
    !! 2. Second, solve the linear system: U*X = Y for X.
    !!
    !! Notice, as both L and U are triangular in structure, the above equations
    !! can be solved by forward and backward substitution.
    !!
    !! @par See Also
    !! - [Wikipedia](https://en.wikipedia.org/wiki/LU_decomposition)
    !! - [Wolfram MathWorld](http://mathworld.wolfram.com/LUDecomposition.html)
    subroutine form_lu_c(n, lu, ipvt, u, p) bind(C, name = "form_lu_")
        ! Arguments
        integer(i32), intent(in), value :: n
        real(dp), intent(inout) :: lu(n,n)
        integer(i32), intent(in) :: ipvt(n)
        real(dp), intent(out) :: u(n,n), p(n,n)

        ! Process
        call form_lu(lu, ipvt, u, p)
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Computes the QR factorization of an M-by-N matrix without
    !! pivoting.
    !!
    !! @param[in] m The number of rows in the matrix.
    !! @param[in] n The number of columns in the matrix.
    !! @param[in,out] a On input, the M-by-N matrix to factor.  On output, the
    !!  elements on and above the diagonal contain the MIN(M, N)-by-N upper
    !!  trapezoidal matrix R (R is upper triangular if M >= N).  The elements
    !!  below the diagonal, along with the array @p tau, represent the
    !!  orthogonal matrix Q as a product of elementary reflectors.
    !! @param[out] tau A MIN(M, N)-element array used to store the scalar
    !!  factors of the elementary reflectors.
    !! @param[in,out] err The errorhandler object.  If no error handling is
    !!  desired, simply pass NULL, and errors will be dealt with by the default
    !!  internal error handler.  Possible errors that may be encountered are as
    !!  follows.
    !!  - LA_ARRAY_SIZE_ERROR: Occurs if the scalar factor array is not sized 
    !!      appropriately.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
    !!      there is insufficient memory available.
    subroutine qr_factor_c(m, n, a, tau, err) bind(C, name = "qr_factor_")
        ! Arguments
        integer(i32), intent(in), value :: m, n
        real(dp), intent(inout) :: a(m,n)
        real(dp), intent(out) :: tau(min(m,n))
        type(errorhandler), intent(inout) :: err

        ! Local Variables
        class(errors), allocatable :: eptr

        ! Process
        call get_errorhandler(err, eptr)
        if (allocated(eptr)) then
            call qr_factor(a, tau, err = eptr)
            call update_errorhandler(eptr, err)
        else
            call qr_factor(a, tau)
        end if
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Computes the QR factorization of an M-by-N matrix with column
    !! pivoting such that A * P = Q * R.
    !!
    !! @param[in] m The number of rows in the matrix.
    !! @param[in] n The number of columns in the matrix.
    !! @param[in,out] a On input, the M-by-N matrix to factor.  On output, the
    !!  elements on and above the diagonal contain the MIN(M, N)-by-N upper
    !!  trapezoidal matrix R (R is upper triangular if M >= N).  The elements
    !!  below the diagonal, along with the array @p tau, represent the
    !!  orthogonal matrix Q as a product of elementary reflectors.
    !! @param[out] tau A MIN(M, N)-element array used to store the scalar
    !!  factors of the elementary reflectors.
    !! @param[in,out] jpvt On input, an N-element array that if JPVT(I) .ne. 0,
    !!  the I-th column of A is permuted to the front of A * P; if JPVT(I) = 0,
    !!  the I-th column of A is a free column.  On output, if JPVT(I) = K, then
    !!  the I-th column of A * P was the K-th column of A.
    !! @param[in,out] err The errorhandler object.  If no error handling is
    !!  desired, simply pass NULL, and errors will be dealt with by the default
    !!  internal error handler.  Possible errors that may be encountered are as
    !!  follows.
    !!  - LA_ARRAY_SIZE_ERROR: Occurs if the scalar factor array is not sized 
    !!      appropriately.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
    !!      there is insufficient memory available.
    subroutine qr_factor_pivot_c(m, n, a, tau, jpvt, err) &
            bind(C, name = "qr_factor_pivot_")
        ! Arguments
        integer(i32), intent(in), value :: m, n
        real(dp), intent(inout) :: a(m,n)
        real(dp), intent(out) :: tau(min(m,n))
        integer(i32), intent(inout) :: jpvt(n)
        type(errorhandler), intent(inout) :: err

        ! Local Variables
        class(errors), allocatable :: eptr

        ! Process
        call get_errorhandler(err, eptr)
        if (allocated(eptr)) then
            call qr_factor(a, tau, jpvt, err = eptr)
            call update_errorhandler(eptr, err)
        else
            call qr_factor(a, tau, jpvt)
        end if
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Forms the full M-by-M orthogonal matrix Q from the elementary
    !! reflectors returned by the base QR factorization algorithm.
    !!
    !! @param[in] m The number of rows in the original matrix.
    !! @param[in] n The number of columns in the original matrix.
    !! @param[in,out] r On input, an M-by-N matrix where the elements below the
    !!  diagonal contain the elementary reflectors generated from the QR
    !!  factorization.  On and above the diagonal, the matrix contains the
    !!  matrix R.  On output, the elements below the diagonal are zeroed such
    !!  that the remaining matrix is simply the M-by-N matrix R.
    !! @param[in] tau A MIN(M, N)-element array containing the scalar factors of
    !!  each elementary reflector defined in @p r.
    !! @param[out] q An M-by-M matrix where the full orthogonal matrix Q will be
    !!  written.  In the event that M > N, Q may be supplied as M-by-N, and
    !!  therefore only return the useful submatrix Q1 (Q = [Q1, Q2]) as the
    !!  factorization can be written as Q * R = [Q1, Q2] * [R1; 0].
    !! @param[in,out] err The errorhandler object.  If no error handling is
    !!  desired, simply pass NULL, and errors will be dealt with by the default
    !!  internal error handler.  Possible errors that may be encountered are as
    !!  follows.
    !!  - LA_ARRAY_SIZE_ERROR: Occurs if the scalar factor array is not sized 
    !!      appropriately.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
    !!      there is insufficient memory available.
    subroutine form_qr_c(m, n, r, tau, q, err) bind(C, name = "form_qr_")
        ! Arguments
        integer(i32), intent(in), value :: m, n
        real(dp), intent(inout) :: r(m,n)
        real(dp), intent(in) :: tau(min(m,n))
        real(dp), intent(out) :: q(m,m)
        type(errorhandler), intent(inout) :: err

        ! Local Variables
        class(errors), allocatable :: eptr

        ! Process
        call get_errorhandler(err, eptr)
        if (allocated(eptr)) then
            call form_qr(r, tau, q, err = eptr)
            call update_errorhandler(eptr, err)
        else
            call form_qr(r, tau, q)
        end if
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Forms the full M-by-M orthogonal matrix Q from the elementary
    !! reflectors returned by the base QR factorization algorithm.
    !!
    !! @param[in] m The number of rows in the original matrix.
    !! @param[in] n The number of columns in the original matrix.
    !! @param[in,out] r On input, an M-by-N matrix where the elements below the
    !!  diagonal contain the elementary reflectors generated from the QR
    !!  factorization.  On and above the diagonal, the matrix contains the
    !!  matrix R.  On output, the elements below the diagonal are zeroed such
    !!  that the remaining matrix is simply the M-by-N matrix R.
    !! @param[in] tau A MIN(M, N)-element array containing the scalar factors of
    !!  each elementary reflector defined in @p r.
    !! @param[in] pvt An N-element column pivot array as returned by the QR
    !!  factorization.
    !! @param[out] q An M-by-M matrix where the full orthogonal matrix Q will be
    !!  written.  In the event that M > N, Q may be supplied as M-by-N, and
    !!  therefore only return the useful submatrix Q1 (Q = [Q1, Q2]) as the
    !!  factorization can be written as Q * R = [Q1, Q2] * [R1; 0].
    !! @param[out] p An N-by-N matrix where the pivot matrix will be written.
    !! @param[in,out] err The errorhandler object.  If no error handling is
    !!  desired, simply pass NULL, and errors will be dealt with by the default
    !!  internal error handler.  Possible errors that may be encountered are as
    !!  follows.
    !!  - LA_ARRAY_SIZE_ERROR: Occurs if the scalar factor array is not sized 
    !!      appropriately.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
    !!      there is insufficient memory available.
    subroutine form_qr_pivot_c(m, n, r, tau, pvt, q, p, err) &
            bind(C, name = "form_qr_pivot_")
        ! Arguments
        integer(i32), intent(in), value :: m, n
        real(dp), intent(inout) :: r(m,n)
        real(dp), intent(in) :: tau(min(m,n))
        integer(i32), intent(in) :: pvt(n)
        real(dp), intent(out) :: q(m,m), p(n,n)
        type(errorhandler), intent(inout) :: err

        ! Local Variables
        class(errors), allocatable :: eptr

        ! Process
        call get_errorhandler(err, eptr)
        if (allocated(eptr)) then
            call form_qr(r, tau, pvt, q, p, err = eptr)
            call update_errorhandler(eptr, err)
        else
            call form_qr(r, tau, pvt, q, p)
        end if
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Multiplies a general matrix by the orthogonal matrix Q from a QR
    !! factorization such that: C = op(Q) * C.
    !!
    !! @param[in] trans Set to true to apply Q**T; else, set to false.
    !! @param[in] m The number of rows in the matrix @p c.
    !! @param[in] n The number of columns in the matrix @p c.
    !! @param[in] q On input, an M-by-M matrix containing the elementary
    !!  reflectors output from the QR factorization.    Notice, the contents of 
    !!  this matrix are restored on exit.
    !!  that the remaining matrix is simply the M-by-N matrix R.
    !! @param[in] tau A MIN(M,N)-element array containing the scalar factors of 
    !!  each elementary reflector defined in @p a.
    !! @param[in,out] c On input, the M-by-N matrix C.  On output, the product
    !!  of the orthogonal matrix Q and the original matrix C.
    !! @param[in,out] err The errorhandler object.  If no error handling is
    !!  desired, simply pass NULL, and errors will be dealt with by the default
    !!  internal error handler.  Possible errors that may be encountered are as
    !!  follows.
    !!  - LA_ARRAY_SIZE_ERROR: Occurs if the scalar factor array is not sized 
    !!      appropriately.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
    !!      there is insufficient memory available.
    subroutine mult_qr_c(trans, m, n, q, tau, c, err) &
            bind(C, name = "mult_qr_")
        ! Arguments
        logical(c_bool), intent(in), value :: trans
        integer(i32), intent(in), value :: m, n
        real(dp), intent(inout) :: q(m,m), c(m,n)
        real(dp), intent(in) :: tau(min(m,n))
        type(errorhandler), intent(inout) :: err

        ! Local Variables
        class(errors), allocatable :: eptr

        ! Process
        call get_errorhandler(err, eptr)
        if (allocated(eptr)) then
            call mult_qr(.true., logical(trans), q, tau, c, err = eptr)
            call update_errorhandler(eptr, err)
        else
            call mult_qr(.true., logical(trans), q, tau, c)
        end if
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Computes the rank 1 update to an M-by-N QR factored matrix A
    !! (M >= N) where A = Q * R, and A1 = A + U * V**T such that A1 = Q1 * R1.
    !!
    !! @param[in] m The number of rows in the original matrix.
    !! @param[in] n The number of columns in the original matrix.
    !! @param[in,out] q On input, the original M-by-M orthogonal matrix Q.  On
    !!  output, the updated matrix Q1.
    !! @param[in,out] r On input, the M-by-N matrix R.  On output, the updated
    !!  matrix R1.
    !! @param[in,out] u On input, the M-element U update vector.  On output,
    !!  the original content of the array is overwritten.
    !! @param[in,out] v On input, the N-element V update vector.  On output,
    !!  the original content of the array is overwritten.
    !! @param[in,out] err The errorhandler object.  If no error handling is
    !!  desired, simply pass NULL, and errors will be dealt with by the default
    !!  internal error handler.  Possible errors that may be encountered are as
    !!  follows.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
    !!      there is insufficient memory available.
    subroutine qr_rank1_update_c(m, n, q, r, u, v, err) &
            bind(C, name = "qr_rank1_update_")
        ! Arguments
        integer(i32), intent(in), value :: m, n
        real(dp), intent(inout) :: q(m,m), r(m,n), u(m), v(n)
        type(errorhandler), intent(inout) :: err

        ! Local Variables
        class(errors), allocatable :: eptr

        ! Process
        call get_errorhandler(err, eptr)
        if (allocated(eptr)) then
            call qr_rank1_update(q, r, u, v, err = eptr)
            call update_errorhandler(eptr, err)
        else
            call qr_rank1_update(q, r, u, v)
        endif
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Computes the Cholesky factorization of a symmetric, positive
    !! definite matrix.
    !!
    !! @param[in] n The dimension of the matrix.
    !! @param[in,out] a On input, the N-by-N matrix to factor.  On output, the
    !!  factored matrix is returned in either the upper or lower triangular
    !!  portion of the matrix, dependent upon the value of @p upper.
    !! @param[in] upper An optional input that, if specified, provides control
    !!  over whether the factorization is computed as A = U**T * U (set to
    !!  true), or as A = L * L**T (set to false).  The default value is true
    !!  such that A = U**T * U.
    !! @param[in,out] err The errorhandler object.  If no error handling is
    !!  desired, simply pass NULL, and errors will be dealt with by the default
    !!  internal error handler.  Possible errors that may be encountered are as
    !!  follows.
    !!  - LA_MATRIX_FORMAT_ERROR: Occurs if @p a is not positive definite.
    subroutine cholesky_factor_c(n, a, upper, err) &
            bind(C, name = "cholesky_factor_")
        ! Arguments
        integer(i32), intent(in), value :: n
        real(dp), intent(inout) :: a(n,n)
        logical(c_bool), intent(in), value :: upper
        type(errorhandler), intent(inout) :: err

        ! Local Variables
        class(errors), allocatable :: eptr

        ! Process
        call get_errorhandler(err, eptr)
        if (allocated(eptr)) then
            call cholesky_factor(a, logical(upper), eptr)
            call update_errorhandler(eptr, err)
        else
            call cholesky_factor(a, logical(upper))
        end if
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Computes the rank 1 update to a Cholesky factored matrix (upper
    !! triangular).
    !!
    !! @param[in] n The dimension of the matrix.
    !! @param[in,out] r On input, the N-by-N upper triangular matrix R.  On
    !!  output, the updated matrix R1.
    !! @param[in,out] u On input, the N-element update vector U.  On output,
    !!  the rotation sines used to transform R to R1.
    !! @param[in,out] err The errorhandler object.  If no error handling is
    !!  desired, simply pass NULL, and errors will be dealt with by the default
    !!  internal error handler.  Possible errors that may be encountered are as
    !!  follows.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
    !!      there is insufficient memory available.
    subroutine cholesky_rank1_update_c(n, r, u, err) &
            bind(C, name = "cholesky_rank1_update_")
        ! Arguments
        integer(i32), intent(in), value :: n
        real(dp), intent(inout) :: r(n,n), u(n)
        type(errorhandler), intent(inout) :: err

        ! Local Variables
        class(errors), allocatable :: eptr

        ! Process
        call get_errorhandler(err, eptr)
        if (allocated(eptr)) then
            call cholesky_rank1_update(r, u, err = eptr)
            call update_errorhandler(eptr, err)
        else
            call cholesky_rank1_update(r, u)
        end if
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Computes the rank 1 downdate to a Cholesky factored matrix (upper
    !! triangular).
    !!
    !! @param[in] n The dimension of the matrix.
    !! @param[in,out] r On input, the N-by-N upper triangular matrix R.  On
    !!  output, the updated matrix R1.
    !! @param[in,out] u On input, the N-element update vector U.  On output,
    !!  the rotation sines used to transform R to R1.
    !! @param[in,out] err The errorhandler object.  If no error handling is
    !!  desired, simply pass NULL, and errors will be dealt with by the default
    !!  internal error handler.  Possible errors that may be encountered are as
    !!  follows.
    !!  - LA_ARRAY_SIZE_ERROR: Occurs if any of the input array sizes are 
    !!      incorrect.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
    !!      there is insufficient memory available.
    !!  - LA_MATRIX_FORMAT_ERROR: Occurs if the downdated matrix is not 
    !!      positive definite.
    !!  - LA_SINGULAR_MATRIX_ERROR: Occurs if @p r is singular.
    subroutine cholesky_rank1_downdate_c(n, r, u, err) &
            bind(C, name = "cholesky_rank1_downdate_")
        ! Arguments
        integer(i32), intent(in), value :: n
        real(dp), intent(inout) :: r(n,n), u(n)
        type(errorhandler), intent(inout) :: err

        ! Local Variables
        class(errors), allocatable :: eptr

        ! Process
        call get_errorhandler(err, eptr)
        if (allocated(eptr)) then
            call cholesky_rank1_downdate(r, u, err = eptr)
            call update_errorhandler(eptr, err)
        else
            call cholesky_rank1_downdate(r, u)
        end if
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Factors an upper trapezoidal matrix by means of orthogonal
    !! transformations such that A = R * Z = (R 0) * Z.  Z is an orthogonal
    !! matrix of dimension N-by-N, and R is an M-by-M upper triangular
    !! matrix.
    !!
    !! @param[in] m The number of rows in the original matrix.
    !! @param[in] n The number of columns in the original matrix.
    !! @param[in,out] a On input, the M-by-N upper trapezoidal matrix to factor.
    !!  On output, the leading M-by-M upper triangular part of the matrix
    !!  contains the upper triangular matrix R, and elements N-L+1 to N of the
    !!  first M rows of A, with the array @p tau, represent the orthogonal
    !!  matrix Z as a product of M elementary reflectors.
    !! @param[out] tau An M-element array used to store the scalar
    !!  factors of the elementary reflectors.
    !! @param[in,out] err The errorhandler object.  If no error handling is
    !!  desired, simply pass NULL, and errors will be dealt with by the default
    !!  internal error handler.  Possible errors that may be encountered are as
    !!  follows.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
    !!      there is insufficient memory available.
    subroutine rz_factor_c(m, n, a, tau, err) bind(C, name = "rz_factor_")
        ! Arguments
        integer(i32), intent(in), value :: m, n
        real(dp), intent(inout) :: a(m,n)
        real(dp), intent(out) :: tau(m)
        type(errorhandler), intent(inout) :: err

        ! Local Variables
        class(errors), allocatable :: eptr

        ! Process
        call get_errorhandler(err, eptr)
        if (allocated(eptr)) then
            call rz_factor(a, tau, err = eptr)
            call update_errorhandler(eptr, err)
        else
            call rz_factor(a, tau)
        end if
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Multiplies a general matrix by the orthogonal matrix Z from an
    !! RZ factorization such that: C = op(Z) * C.
    !!
    !! @param[in] trans Set to true to apply Z**T; else, set to false.
    !! @param[in] m The number of rows in the matrix @p c.
    !! @param[in] n The number of columns in the matrix @p c.
    !! @param[in] l The number of columns in matrix @p a containing the
    !!  meaningful part of the Householder vectors (M >= L >= 0).
    !! @param[in,out] a On input, the M-by-M matrix Z as output by @p rz_factor.
    !!  The matrix is used as in-place storage during execution; however, the
    !!  contents of the matrix are restored on exit.
    !! @param[in] tau An M-element array containing the scalar factors of the
    !!  elementary reflectors found in @p a.
    !! @param[in,out] c On input, the M-by-N matrix C.  On output, the product
    !!  of the orthogonal matrix Z and the original matrix C.
    !! @param[in,out] err The errorhandler object.  If no error handling is
    !!  desired, simply pass NULL, and errors will be dealt with by the default
    !!  internal error handler.  Possible errors that may be encountered are as
    !!  follows.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
    !!      there is insufficient memory available.
    subroutine mult_rz_c(trans, m, n, l, a, tau, c, err) &
            bind(C, name = "mult_rz_")
        ! Arguments
        logical(c_bool), intent(in), value :: trans
        integer(i32), intent(in), value :: m, n, l
        real(dp), intent(inout) :: a(m,m), c(m,n)
        real(dp), intent(in) :: tau(m)
        type(errorhandler), intent(inout) :: err

        ! Local Variables
        class(errors), allocatable :: eptr

        ! Process
        call get_errorhandler(err, eptr)
        if (allocated(eptr)) then
            call mult_rz(.true., logical(trans), l, a, tau, c, err = eptr)
            call update_errorhandler(eptr, err)
        else
            call mult_rz(.true., logical(trans), l, a, tau, c)
        end if
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Computes the singular value decomposition of a matrix A.  The
    !!  SVD is defined as: A = U * S * V**T, where U is an M-by-M orthogonal
    !!  matrix, S is an M-by-N diagonal matrix, and V is an N-by-N orthogonal
    !!  matrix.
    !!
    !! @param[in] m The number of rows in the original matrix.
    !! @param[in] n The number of columns in the original matrix.
    !! @param[in,out] a On input, the M-by-N matrix to factor.  The matrix is
    !!  overwritten on output.
    !!  that the remaining matrix is simply the M-by-N matrix R.
    !! @param[out] s A MIN(M, N)-element array containing the singular values
    !!  of @p a sorted in descending order.
    !! @param[out] u An M-by-M matrix that on output contains the left singular
    !!  vectors (matrix U in the decomposition: A = U * S * V**T)
    !! @param[out] vt An N-by-N matrix that on output contains the right
    !!  singular vectors (matrix V**T in the decomposition: A = U * S * V**T).
    !! @param[in,out] err The errorhandler object.  If no error handling is
    !!  desired, simply pass NULL, and errors will be dealt with by the default
    !!  internal error handler.  Possible errors that may be encountered are as
    !!  follows.
    !!  - LA_ARRAY_SIZE_ERROR: Occurs if the singular value array is not sized 
    !!      appropriately.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
    !!      there is insufficient memory available.
    !!  - LA_CONVERGENCE_ERROR: Occurs as a warning if the QR iteration process
    !!      could not converge to a zero value.
    subroutine svd_c(m, n, a, s, u, vt, err) bind(C, name = "svd_")
        ! Arguments
        integer(i32), intent(in), value :: m, n
        real(dp), intent(inout) :: a(m,n)
        real(dp), intent(out) :: s(min(m,n)), u(m,m), vt(n,n)
        type(errorhandler), intent(inout) :: err

        ! Local Variables
        class(errors), allocatable :: eptr

        ! Process
        call get_errorhandler(err, eptr)
        if (allocated(eptr)) then
            call svd(a, s, u, vt, err = eptr)
            call update_errorhandler(eptr, err)
        else
            call svd(a, s, u, vt)
        end if
    end subroutine

! ******************************************************************************
! LINALG_SOLVE ROUTINES
! ------------------------------------------------------------------------------
    !> @brief Solves one of the matrix equations: op(A) * X = alpha * B, where 
    !! A is a triangular matrix.
    !!
    !! @param[in] upper Set to true if A is an upper triangular matrix; else,
    !!  set to false if A is a lower triangular matrix.
    !! @param[in] trans Set to true if op(A) = A**T; else, set to false if
    !!  op(A) = A.
    !! @param[in] nounit Set to true if A is not a unit-diagonal matrix (ones on
    !!  every diagonal element); else, set to false if A is a unit-diagonal
    !!  matrix.
    !! @param[in] n The dimension of the triangular matrix @p a.
    !! @param[in] nrhs The number of right-hand-side vectors (number of columns
    !!  in matrix @p b).
    !! @param[in] alpha The scalar multiplier to B.
    !! @param[in] a N-by-N triangular matrix on which to operate.
    !! @param[in,out] b On input, the N-by-NRHS right-hand-side.  On output, the
    !!  N-by-NRHS solution.
    subroutine solve_tri_mtx_c(upper, trans, nounit, n, nrhs, alpha, a, b) &
            bind(C, name = "solve_triangular_system_")
        ! Arguments
        logical(c_bool), intent(in), value :: upper, trans, nounit
        integer(i32), intent(in), value :: n, nrhs
        real(dp), intent(in), value :: alpha
        real(dp), intent(in) :: a(n,n)
        real(dp), intent(inout) :: b(n,nrhs)

        ! Process
        call solve_triangular_system(.true., logical(upper), logical(trans), &
            logical(nounit), alpha, a, b)
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Solves a system of LU-factored equations.
    !!
    !! @param[in] n The dimension of the original matrix @p a.
    !! @param[in] nrhs The number of right-hand-side vectors (number of columns
    !!  in matrix @p b).
    !! @param[in] a The N-by-N LU factored matrix as output by lu_factor.
    !! @param[in] ipvt The N-element pivot array as output by lu_factor.
    !! @param[in,out] b On input, the N-by-NRHS right-hand-side matrix.  On
    !!  output, the N-by-NRHS solution matrix.
    subroutine solve_lu_c(n, nrhs, a, ipvt, b) bind(C, name = "solve_lu_")
        ! Arguments
        integer(i32), intent(in), value :: n, nrhs
        real(dp), intent(in) :: a(n,n)
        integer(i32), intent(in) :: ipvt(n)
        real(dp), intent(inout) :: b(n,nrhs)

        ! Process
        call solve_lu(a, ipvt, b)
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Solves a system of M QR-factored equations of N unknowns where
    !! M >= N.
    !!
    !! @param[in] m The number of rows in the original coefficient matrix A.
    !! @param[in] n The number of columns in the original coefficient matrix A.
    !! @param[in] nrhs The number of right-hand-side vectors (number of columns
    !!  in matrix @p b).
    !! @param[in] a On input, the M-by-N QR factored matrix as returned by
    !!  qr_factor.  On output, the contents of this matrix are restored.
    !!  Notice, M must be greater than or equal to N.
    !! @param[in] tau A MIN(M, N)-element array containing the scalar factors of
    !!  the elementary reflectors as returned by qr_factor.
    !! @param[in] b On input, the M-by-NRHS right-hand-side matrix.  On output,
    !!  the first N columns are overwritten by the solution matrix X.
    !! @param[in,out] err The errorhandler object.  If no error handling is
    !!  desired, simply pass NULL, and errors will be dealt with by the default
    !!  internal error handler.  Possible errors that may be encountered are as
    !!  follows.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
    !!      there is insufficient memory available.
    subroutine solve_qr_c(m, n, nrhs, a, tau, b, err) &
            bind(C, name = "solve_qr_")
        ! Arguments
        integer(i32), intent(in), value :: m, n, nrhs
        real(dp), intent(inout) :: a(m,n), b(m,nrhs)
        real(dp), intent(in) :: tau(min(m,n))
        type(errorhandler), intent(inout) :: err

        ! Local Variables
        class(errors), allocatable :: eptr

        ! Process
        call get_errorhandler(err, eptr)
        if (allocated(eptr)) then
            call solve_qr(a, tau, b, err = eptr)
            call update_errorhandler(eptr, err)
        else
            call solve_qr(a, tau, b)
        end if
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Solves a system of M QR-factored equations of N unknowns where the
    !! QR factorization made use of column pivoting.
    !!
    !! @param[in] m The number of rows in the original coefficient matrix A.
    !! @param[in] n The number of columns in the original coefficient matrix A.
    !! @param[in] nrhs The number of right-hand-side vectors (number of columns
    !!  in matrix @p b).
    !! @param[in] a On input, the M-by-N QR factored matrix as returned by
    !!  qr_factor.  On output, the contents of this matrix are altered.
    !! @param[in] tau A MIN(M, N)-element array containing the scalar factors of
    !!  the elementary reflectors as returned by qr_factor.
    !! @param[in] jpvt An N-element array, as output by qr_factor, used to
    !!  track the column pivots.
    !! @param[in] b On input, the MAX(M, N)-by-NRHS matrix where the first M
    !!  rows contain the right-hand-side matrix B.  On output, the first N rows
    !!  are overwritten by the solution matrix X.
    !! @param[in,out] err The errorhandler object.  If no error handling is
    !!  desired, simply pass NULL, and errors will be dealt with by the default
    !!  internal error handler.  Possible errors that may be encountered are as
    !!  follows.
    !!  - LA_ARRAY_SIZE_ERROR: Occurs if any of the input arrays are not sized 
    !!      appropriately.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
    !!      there is insufficient memory available.
    subroutine solve_qr_pivot_c(m, n, nrhs, a, tau, jpvt, b, err) &
            bind(C, name = "solve_qr_pivot_")
        ! Arguments
        integer(i32), intent(in), value :: m, n, nrhs
        real(dp), intent(inout) :: a(m,n), b(max(m,n),nrhs)
        real(dp), intent(in) :: tau(min(m,n))
        integer(i32), intent(in) :: jpvt(n)
        type(errorhandler), intent(inout) :: err

        ! Local Variables
        class(errors), allocatable :: eptr

        ! Process
        call get_errorhandler(err, eptr)
        if (allocated(eptr)) then
            call solve_qr(a, tau, jpvt, b, err = eptr)
            call update_errorhandler(eptr, err)
        else
            call solve_qr(a, tau, jpvt, b)
        end if
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Solves a system of Cholesky factored equations.
    !!
    !! @param[in] upper Set to true if the original matrix A was factored such
    !!  that A = U**T * U; else, set to false if the factorization of A was
    !!  A = L**T * L.
    !! @param[in] n The dimension of the original matrix @p a.
    !! @param[in] nrhs The number of right-hand-side vectors (number of columns
    !!  in matrix @p b).
    !! @param[in] a The N-by-N Cholesky factored matrix.
    !! @param[in,out] b On input, the N-by-NRHS right-hand-side matrix B.  On
    !!  output, the solution matrix X.
    subroutine solve_cholesky_c(upper, n, nrhs, a, b) &
            bind(C, name = "solve_cholesky_")
        ! Arguments
        logical(c_bool), intent(in), value :: upper
        integer(i32), intent(in), value :: n, nrhs
        real(dp), intent(in) :: a(n,n)
        real(dp), intent(inout) :: b(n,nrhs)

        ! Process
        call solve_cholesky(logical(upper), a, b)
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Computes the inverse of a square matrix.
    !!
    !! @param[in] n The dimension of the matrix.
    !! @param[in,out] a On input, the N-by-N matrix to invert.  On output, the
    !!  inverted matrix.
    !! @param[in,out] err The errorhandler object.  If no error handling is
    !!  desired, simply pass NULL, and errors will be dealt with by the default
    !!  internal error handler.  Possible errors that may be encountered are as
    !!  follows.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
    !!      there is insufficient memory available.
    !!  - LA_SINGULAR_MATRIX_ERROR: Occurs if the input matrix is singular.
    subroutine mtx_inverse_c(n, a, err) bind(C, name = "mtx_inverse_")
        ! Arguments
        integer(i32), intent(in), value :: n
        real(dp), intent(inout) :: a(n,n)
        type(errorhandler), intent(inout) :: err

        ! Local Variables
        class(errors), allocatable :: eptr

        ! Process
        call get_errorhandler(err, eptr)
        if (allocated(eptr)) then
            call mtx_inverse(a, err = eptr)
            call update_errorhandler(eptr, err)
        else
            call mtx_inverse(a)
        end if
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Computes the Moore-Penrose pseudo-inverse of a M-by-N matrix
    !! using the singular value decomposition of the matrix.
    !!
    !! @param[in] m The number of rows in the matrix to invert.
    !! @param[in] n The number of columns in the matrix to invert.
    !! @param[in,out] a On input, the M-by-N matrix to invert.  The matrix is
    !!  overwritten on output.
    !! @param[out] ainv The N-by-M matrix where the pseudo-inverse of @p a
    !!  will be written.
    !! @param[in,out] err The errorhandler object.  If no error handling is
    !!  desired, simply pass NULL, and errors will be dealt with by the default
    !!  internal error handler.  Possible errors that may be encountered are as
    !!  follows.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
    !!      there is insufficient memory available.
    !!  - LA_CONVERGENCE_ERROR: Occurs as a warning if the QR iteration process
    !!      could not converge to a zero value.
    subroutine mtx_pinverse_c(m, n, a, ainv, err) &
            bind(C, name = "mtx_pinverse_")
        ! Arguments
        integer(i32), intent(in), value :: m, n
        real(dp), intent(inout) :: a(m,n)
        real(dp), intent(out) :: ainv(n,m)
        type(errorhandler), intent(inout) :: err

        ! Local Variables
        class(errors), allocatable :: eptr

        ! Process
        call get_errorhandler(err, eptr)
        if (allocated(eptr)) then
            call mtx_pinverse(a, ainv, err = eptr)
            call update_errorhandler(eptr, err)
        else
            call mtx_pinverse(a, ainv)
        end if
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Solves the overdetermined or underdetermined system (A*X = B) of
    !! M equations of N unknowns using a QR or LQ factorization of the matrix A.
    !! Notice, it is assumed that matrix A has full rank.
    !!
    !! @param[in] m The number of rows in the original coefficient matrix A.
    !! @param[in] n The number of columns in the original coefficient matrix A.
    !! @param[in] nrhs The number of right-hand-side vectors (number of columns
    !!  in matrix @p b).
    !! @param[in,out] a On input, the M-by-N matrix A.  On output, the matrix
    !!  is overwritten by the details of its complete orthogonal factorization.
    !! @param[in,out] b If M >= N, the M-by-NRHS matrix B.  On output, the first
    !!  N rows contain the N-by-NRHS solution matrix X.  If M < N, an
    !!  N-by-NRHS matrix with the first M rows containing the matrix B.  On
    !!  output, the N-by-NRHS solution matrix X.
    !! @param[in,out] err The errorhandler object.  If no error handling is
    !!  desired, simply pass NULL, and errors will be dealt with by the default
    !!  internal error handler.  Possible errors that may be encountered are as
    !!  follows.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
    !!      there is insufficient memory available.
    subroutine solve_least_squares_c(m, n, nrhs, a, b, err) &
            bind(C, name = "solve_least_squares_")
        ! Arguments
        integer(i32), intent(in), value :: m, n, nrhs
        real(dp), intent(inout) :: a(m, n), b(max(m,n), nrhs)
        type(errorhandler), intent(inout) :: err

        ! Local Variables
        class(errors), allocatable :: eptr

        ! Process
        call get_errorhandler(err, eptr)
        if (allocated(eptr)) then
            call solve_least_squares(a, b, err = eptr)
            call update_errorhandler(eptr, err)
        else
            call solve_least_squares(a, b)
        end if
    end subroutine

! ******************************************************************************
! LINALG_EIGEN ROUTINES
! ------------------------------------------------------------------------------
    !> @brief Computes the eigenvalues, and optionally the eigenvectors of a
    !! real, symmetric matrix.
    !!
    !! @param[in] n The dimension of the matrix.
    !! @param[in] vecs Set to true to compute the eigenvectors as well as the
    !!  eigenvalues; else, set to false to just compute the eigenvalues.
    !! @param[in,out] a On input, the N-by-N symmetric matrix on which to
    !!  operate.  On output, and if @p vecs is set to true, the matrix will
    !!  contain the eigenvectors (one per column) corresponding to each
    !!  eigenvalue in @p vals.  If @p vecs is set to false, the lower triangular
    !!  portion of the matrix is overwritten.
    !! @param[out] vals An N-element array that will contain the eigenvalues
    !!  sorted into ascending order.
    !! @param[in,out] err The errorhandler object.  If no error handling is
    !!  desired, simply pass NULL, and errors will be dealt with by the default
    !!  internal error handler.  Possible errors that may be encountered are as
    !!  follows.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
    !!      there is insufficient memory available.
    !!  - LA_CONVERGENCE_ERROR: Occurs if the algorithm failed to converge.
    subroutine eigen_symm_c(n, vecs, a, vals, err) bind(C, name = "eigen_symm_")
        ! Arguments
        integer(i32), intent(in), value :: n
        logical(c_bool), intent(in), value :: vecs
        real(dp), intent(inout) :: a(n,n)
        real(dp), intent(out) :: vals(n)
        type(errorhandler), intent(inout) :: err

        ! Local Variables
        class(errors), allocatable :: eptr

        ! Process
        call get_errorhandler(err, eptr)
        if (allocated(eptr)) then
            call eigen(logical(vecs), a, vals, err = eptr)
            call update_errorhandler(eptr, err)
        else
            call eigen(logical(vecs), a, vals)
        end if
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Computes the eigenvalues, and the right eigenvectors of a square 
    !!  matrix.
    !!
    !! @param[in] n The dimension of the matrix.
    !! @param[in,out] a On input, the N-by-N matrix on which to operate.  On
    !!  output, the contents of this matrix are overwritten.
    !! @param[out] vals An N-element array containing the eigenvalues of the
    !!  matrix on output.  The eigenvalues are not sorted.
    !! @param[out] vecs An N-by-N matrix containing the right eigenvectors 
    !!  (one per column) on output.
    !! @param[in,out] err The errorhandler object.  If no error handling is
    !!  desired, simply pass NULL, and errors will be dealt with by the default
    !!  internal error handler.  Possible errors that may be encountered are as
    !!  follows.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
    !!      there is insufficient memory available.
    !!  - LA_CONVERGENCE_ERROR: Occurs if the algorithm failed to converge.
    subroutine eigen_asymm_c(n, a, vals, vecs, err) &
            bind(C, name = "eigen_asymm_")
        ! Arguments
        integer(i32), intent(in), value :: n
        real(dp), intent(inout) :: a(n,n)
        complex(dp), intent(out) :: vals(n), vecs(n,n)
        type(errorhandler), intent(inout) :: err

        ! Local Variables
        class(errors), allocatable :: eptr

        ! Process
        call get_errorhandler(err, eptr)
        if (allocated(eptr)) then
            call eigen(a, vals, vecs, err = eptr)
            call update_errorhandler(eptr, err)
        else
            call eigen(a, vals, vecs)
        end if
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Computes the eigenvalues, and optionally the right eigenvectors of
    !! a square matrix assuming the structure of the eigenvalue problem is
    !! A*X = lambda*B*X.
    !!
    !! @param[in] n The dimension of the matrix.
    !! @param[in,out] a On input, the N-by-N matrix A.  On output, the contents
    !!  of this matrix are overwritten.
    !! @param[in,out] b On input, the N-by-N matrix B.  On output, the contents
    !!  of this matrix are overwritten.
    !! @param[out] alpha An N-element array that, on output, contains the 
    !!  numerator of the eigenvalue ration ALPHA / BETA.  Computation of this
    !!  ratio isn't necessarily as trivial as it seems as it is entirely 
    !!  possible, and likely, that ALPHA / BETA can overflow or underflow.  With
    !!  that said, the values in ALPHA will always be less than and usually 
    !!  comparable with the NORM(A).
    !! @param[out] beta An N-element array that, on output, contains the
    !!  denominator used to determine the eigenvalues as ALPHA / BETA.  The 
    !!  values in this array will always be less than and usually comparable
    !!  with the NORM(B).
    !! @param[out] vecs An N-by-N matrix containing the right eigenvectors 
    !!  (one per column) on output.
    !! @param[in,out] err The errorhandler object.  If no error handling is
    !!  desired, simply pass NULL, and errors will be dealt with by the default
    !!  internal error handler.  Possible errors that may be encountered are as
    !!  follows.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
    !!      there is insufficient memory available.
    !!  - LA_CONVERGENCE_ERROR: Occurs if the algorithm failed to converge.
    subroutine eigen_gen_c(n, a, b, alpha, beta, vecs, err) &
            bind(C, name = "eigen_gen_")
        ! Arguments
        integer(i32), intent(in), value :: n
        real(dp), intent(inout) :: a(n, n), b(n, n)
        complex(dp), intent(out) :: alpha(n), vecs(n,n)
        real(dp), intent(out) :: beta(n)
        type(errorhandler), intent(inout) :: err

        ! Local Variables
        class(errors), allocatable :: eptr

        ! Process
        call get_errorhandler(err, eptr)
        if (allocated(eptr)) then
            call eigen(a, b, alpha, beta = beta, vecs = vecs, err = eptr)
            call update_errorhandler(eptr, err)
        else
            call eigen(a, b, alpha, beta = beta, vecs = vecs)
        end if
    end subroutine

! ------------------------------------------------------------------------------
end module
