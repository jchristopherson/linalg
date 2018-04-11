! linalg_constants.f90

!> @brief \b linalg_constants
!!
!! @par Purpose
!! Provides a set of constants and error flags for the library.
module linalg_constants
    use, intrinsic :: iso_fortran_env, only : int32, real64
    implicit none

! ******************************************************************************
! ERROR FLAGS
! ------------------------------------------------------------------------------
    !> An error flag denoting an invalid input.
    integer, parameter :: LA_INVALID_INPUT_ERROR = 101
    !> An error flag denoting an improperly sized array.
    integer, parameter :: LA_ARRAY_SIZE_ERROR = 102
    !> An error flag denoting a singular matrix.
    integer, parameter :: LA_SINGULAR_MATRIX_ERROR = 103
    !> An error flag denoting an issue with the matrix format.
    integer, parameter :: LA_MATRIX_FORMAT_ERROR = 104
    !> An error flag denoting that there is insufficient memory available.
    integer, parameter :: LA_OUT_OF_MEMORY_ERROR = 105
    !> An error flag denoting a convergence failure.
    integer, parameter :: LA_CONVERGENCE_ERROR = 106
    !> An error resulting from an invalid operation.
    integer, parameter :: LA_INVALID_OPERATION_ERROR = 107
end module
