! linalg_c_binding.f90

!> @brief \b linalg_c_binding
!!
!! @par Purpose
!! Provides a C friendly interface to the LINALG library.
module linalg_c_binding
    use, intrinsic :: iso_c_binding
    use linalg_constants
    use linalg_core
    use ferror, only : errors
contains
! ******************************************************************************
! LINALG_CORE ROUTINES
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

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
    !
    ! subroutine diag_mtx_mult_c(lside, trans, alpha, a, b, beta, c)
    ! end subroutine

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
    !!  - LA_ARRAY_SIZE_ERROR: Occurs if any of the input arrays are not sized
    !!      appropriately.
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

end module
