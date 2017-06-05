! linalg_c_binding.f90

!> @brief \b linalg_c_binding
!!
!! @par Purpose
!! Provides a C friendly interface to the LINALG library.
module linalg_c_binding
    use, intrinsic :: iso_c_binding
    use linalg_constants
    use linalg_core
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
            bind(c, name = "rank1_update")
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
    pure function trace_c(m, n, x) result(y) bind(c, name = "trace")
        ! Arguments
        integer(i32), intent(in), value :: m, n
        real(dp), intent(in) :: x(m,n)
        real(dp) :: y

        ! Process
        y = trace(x)
    end function


! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

end module
