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
    !! @param[in] alpha The scalar multiplier.
    !! @param[in] x An M-element array.
    !! @param[in] y An N-element array.
    !! @param[in,out] a On input, the M-by-N matrix to update.  On output, the
    !!  updated M-by-N matrix.
    !!
    !! @par Notes
    !! This routine is based upon the BLAS routine DGER.
    subroutine rank1_update_c(alpha, m, n, x, y, a) &
            bind(c, name = "rank1_update")
        ! Arguments
        real(dp), intent(in), value :: alpha
        integer(i32), intent(in), value :: m, n
        real(dp), intent(in) :: x(m), y(n)
        real(dp), intent(inout) :: a(m,n)

        ! Process
        call rank1_update(alpha, x, y, a)
    end subroutine

end module
