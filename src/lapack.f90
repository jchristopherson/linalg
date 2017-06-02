! lapack.f90

!> @brief \b lapack
!!
!! @par Purpose
!! Provides interfaces to various LAPACK routines.
module lapack
    implicit none

! ******************************************************************************
! LAPACK FUNCTION INTERFACES
! ------------------------------------------------------------------------------
    interface
        function DLAMCH(cmach) result(x)
            use linalg_constants, only : dp
            character, intent(in) :: cmach
            real(dp) :: x
        end function
    end interface

end module
