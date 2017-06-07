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
    use ferror, only : errors
contains
! ******************************************************************************
! LINALG_CORE ROUTINES
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
            bind(C, name = "rank1_update")
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
    pure function trace_c(m, n, x) result(y) bind(C, name = "trace")
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
    !! @param[in] err A pointer to the C error handler object.  If no error
    !!  handling is desired, simply pass NULL, and errors will be dealt with
    !!  by the default internal error handler.  Possible errors that may be
    !!  encountered are as follows.
    !!  - LA_ARRAY_SIZE_ERROR: Occurs if any of the input arrays are not sized
    !!      appropriately.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
    !!      there is insufficient memory available.
    !!  - LA_CONVERGENCE_ERROR: Occurs as a warning if the QR iteration process
    !!      could not converge to a zero value.
    !!
    !! @par See Also
    !! - [Wolfram MathWorld](http://mathworld.wolfram.com/MatrixRank.html)
    function mtx_rank_c(m, n, a, err) result(rnk) bind(C, name = "mtx_rank")
        ! Arguments
        integer(i32), intent(in), value :: m, n
        real(dp), intent(inout) :: a(m,n)
        type(c_ptr), intent(in), value :: err
        integer(i32) :: rnk

        ! Local Variables
        type(errors), pointer :: eptr
        logical :: useError

        ! Set up error handling
        if (c_associated(err)) then
            call c_f_pointer(err, eptr)
            useError = .true.
        else
            useError = .false.
        end if

        ! Process
        if (useError) then
            rnk = mtx_rank(a, err = eptr)
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
    !! @param[in] err A pointer to the C error handler object.  If no error
    !!  handling is desired, simply pass NULL, and errors will be dealt with
    !!  by the default internal error handler.  Possible errors that may be
    !!  encountered are as follows.
    !!  - LA_ARRAY_SIZE_ERROR: Occurs if the input matrix is not square.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
    !!      there is insufficient memory available.
    function det_c(n, a, err) result(x) bind(C, name = "det")
        ! Arguments
        integer(i32), intent(in), value :: n
        real(dp), intent(inout) :: a(n,n)
        type(c_ptr), intent(in), value :: err
        real(dp) :: x

        ! Local Variables
        type(errors), pointer :: eptr
        logical :: useError

        ! Set up error handling
        if (c_associated(err)) then
            call c_f_pointer(err, eptr)
            useError = .true.
        else
            useError = .false.
        end if

        ! Process
        if (useError) then
            x = det(a, err = eptr)
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
    subroutine swap_c(n, x, y) bind(C, name = "swap")
        ! Arguments
        integer(i32), intent(in), value :: n
        real(dp), intent(inout) :: x(n), y(n)

        ! Process
        call swap(x, y)
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
    !! @param[in] ni The number of elements in the pivot array @p ipvt.  This
    !!  value must be equal to MIN(M, N).
    !! @param[out] ipvt An MIN(M, N)-element array used to track row-pivot
    !!  operations.  The array stored pivot information such that row I is
    !!  interchanged with row IPVT(I).
    !! @param[in] err A pointer to the C error handler object.  If no error
    !!  handling is desired, simply pass NULL, and errors will be dealt with
    !!  by the default internal error handler.  Possible errors that may be
    !!  encountered are as follows.
    !!  - LA_ARRAY_SIZE_ERROR: Occurs if the pivot array is not sized 
    !!      appropriately.
    !!  - LA_SINGULAR_MATRIX_ERROR: Occurs as a warning if @p a is found to be
    !!      singular.
    subroutine lu_factor_c(m, n, a, ni, ipvt, err) bind(C, name = "lu_factor")
        ! Arguments
        integer(i32), intent(in), value :: m, n, ni
        real(dp), intent(inout) :: a(m,n)
        integer(i32), intent(out) :: ipvt(ipvt)
        type(c_ptr), intent(in), value :: err

        ! Local Variables
        type(errors), pointer :: eptr
        logical :: useError

        ! Set up error handling
        if (c_associated(err)) then
            call c_f_pointer(err, eptr)
            useError = .true.
        else
            useError = .false.
        end if

        ! Process
        if (useError) then
            call lu_factor(a, ipvt, eptr)
        else
            call lu_factor(a, ipvt)
        end if
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Extracts the L, U, and P matrices from the output of the
    !! @ref lu_factor routine.
    !!
    !! @param[in] n The dimension of the original matrix.
    !! @param[in,out] lu On input, the N-by-N matrix as output by
    !!  @ref lu_factor.  On output, the N-by-N lower triangular matrix L.
    !! @param[in] ipvt The N-element pivot array as output by
    !!  @ref lu_factor.
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
    subroutine form_lu_c(n, lu, ipvt, u, p) bind(C, name = "form_lu")
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
    !! @param[in] nt The number of elements in the scalar factor array @p tau.
    !!  This value must be equal to MIN(M, N).
    !! @param[out] tau A MIN(M, N)-element array used to store the scalar
    !!  factors of the elementary reflectors.
    !! @param[in] err A pointer to the C error handler object.  If no error
    !!  handling is desired, simply pass NULL, and errors will be dealt with
    !!  by the default internal error handler.  Possible errors that may be
    !!  encountered are as follows.
    !!  - LA_ARRAY_SIZE_ERROR: Occurs if the scalar factor array is not sized 
    !!      appropriately.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
    !!      there is insufficient memory available.
    subroutine qr_factor_c(m, n, a, nt, tau, err) bind(C, name = "qr_factor")
        ! Arguments
        integer(i32), intent(in), value :: m, n, nt
        real(dp), intent(inout) :: a(m,n)
        real(dp), intent(out) :: tau(nt)
        type(c_ptr), intent(in), value :: err

        ! Local Variables
        type(errors), pointer :: eptr
        logical :: useError

        ! Set up error handling
        if (c_associated(err)) then
            call c_f_pointer(err, eptr)
            useError = .true.
        else
            useError = .false.
        end if

        ! Process
        if (useError) then
            call qr_factor(a, tau, err = eptr)
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
    !! @param[in] nt The number of elements in the scalar factor array @p tau.
    !!  This value must be equal to MIN(M, N).
    !! @param[out] tau A MIN(M, N)-element array used to store the scalar
    !!  factors of the elementary reflectors.
    !! @param[in,out] jpvt On input, an N-element array that if JPVT(I) .ne. 0,
    !!  the I-th column of A is permuted to the front of A * P; if JPVT(I) = 0,
    !!  the I-th column of A is a free column.  On output, if JPVT(I) = K, then
    !!  the I-th column of A * P was the K-th column of A.
    !! @param[in] err A pointer to the C error handler object.  If no error
    !!  handling is desired, simply pass NULL, and errors will be dealt with
    !!  by the default internal error handler.  Possible errors that may be
    !!  encountered are as follows.
    !!  - LA_ARRAY_SIZE_ERROR: Occurs if the scalar factor array is not sized 
    !!      appropriately.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
    !!      there is insufficient memory available.
    subroutine qr_factor_pivot_c(m, n, a, nt, tau, jpvt, err) &
            bind(C, name = "qr_factor_pivot")
        ! Arguments
        integer(i32), intent(in), value :: m, n, nt
        real(dp), intent(inout) :: a(m,n)
        real(dp), intent(out) :: tau(nt)
        integer(i32), intent(inout) :: jpvt(n)
        type(c_ptr), intent(in), value :: err

        ! Local Variables
        type(errors), pointer :: eptr
        logical :: useError

        ! Set up error handling
        if (c_associated(err)) then
            call c_f_pointer(err, eptr)
            useError = .true.
        else
            useError = .false.
        end if

        ! Process
        if (useError) then
            call qr_factor(a, tau, jpvt, err = eptr)
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
    !! @param[in] nt The number of elements in the scalar factor array @p tau.
    !!  This value must be equal to MIN(M, N).
    !! @param[in] tau A MIN(M, N)-element array containing the scalar factors of
    !!  each elementary reflector defined in @p h.
    !! @param[out] q An M-by-M matrix where the full orthogonal matrix Q will be
    !!  written.  In the event that M > N, Q may be supplied as M-by-N, and
    !!  therefore only return the useful submatrix Q1 (Q = [Q1, Q2]) as the
    !!  factorization can be written as Q * R = [Q1, Q2] * [R1; 0].
    !! @param[in] err A pointer to the C error handler object.  If no error
    !!  handling is desired, simply pass NULL, and errors will be dealt with
    !!  by the default internal error handler.  Possible errors that may be
    !!  encountered are as follows.
    !!  - LA_ARRAY_SIZE_ERROR: Occurs if the scalar factor array is not sized 
    !!      appropriately.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
    !!      there is insufficient memory available.
    subroutine form_qr_c(m, n, r, nt, tau, q, err) bind(C, name = "form_qr")
        ! Arguments
        integer(i32), intent(in), value :: m, n, nt
        real(dp), intent(inout) :: r(m,n)
        real(dp), intent(in) :: tau(nt)
        real(dp), intent(out) :: q(m,m)
        type(c_ptr), intent(in), value :: err

        ! Local Variables
        type(errors), pointer :: eptr
        logical :: useError

        ! Set up error handling
        if (c_associated(err)) then
            call c_f_pointer(err, eptr)
            useError = .true.
        else
            useError = .false.
        end if

        ! Process
        if (useError) then
            call form_qr(r, tau, q, err = eptr)
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
    !! @param[in] nt The number of elements in the scalar factor array @p tau.
    !!  This value must be equal to MIN(M, N).
    !! @param[in] tau A MIN(M, N)-element array containing the scalar factors of
    !!  each elementary reflector defined in @p h.
    !! @param[in] pvt An N-element column pivot array as returned by the QR
    !!  factorization.
    !! @param[out] q An M-by-M matrix where the full orthogonal matrix Q will be
    !!  written.  In the event that M > N, Q may be supplied as M-by-N, and
    !!  therefore only return the useful submatrix Q1 (Q = [Q1, Q2]) as the
    !!  factorization can be written as Q * R = [Q1, Q2] * [R1; 0].
    !! @param[out] p An N-by-N matrix where the pivot matrix will be written.
    !! @param[in] err A pointer to the C error handler object.  If no error
    !!  handling is desired, simply pass NULL, and errors will be dealt with
    !!  by the default internal error handler.  Possible errors that may be
    !!  encountered are as follows.
    !!  - LA_ARRAY_SIZE_ERROR: Occurs if the scalar factor array is not sized 
    !!      appropriately.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
    !!      there is insufficient memory available.
    subroutine form_qr_pivot_c(m, n, r, nt, tau, pvt, q, p, err) &
            bind(C, name = "form_qr_pivot")
        ! Arguments
        integer(i32), intent(in), value :: m, n, nt
        real(dp), intent(inout) :: r(m,n)
        real(dp), intent(in) :: tau(nt)
        integer(i32), intent(in) :: pvt(n)
        real(dp), intent(out) :: q(m,m), p(n,n)
        type(c_ptr), intent(in), value :: err

        ! Local Variables
        type(errors), pointer :: eptr
        logical :: useError

        ! Set up error handling
        if (c_associated(err)) then
            call c_f_pointer(err, eptr)
            useError = .true.
        else
            useError = .false.
        end if

        ! Process
        if (useError) then
            call form_qr(r, tau, pvt, q, err = eptr)
        else
            call form_qr(r, tau, pvt, q)
        end if
    end subroutine

! ------------------------------------------------------------------------------

end module
