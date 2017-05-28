! error_flags.f90

!> @brief \b error_flags
!!
!! @par Purpose
!! A set of error flags.
module error_flags
    implicit none

    !> An error flag denoting an invalid input.
    integer, parameter :: LA_INVALID_INPUT_ERROR = 101
    !> An error flag denoting an improperly sized array.
    integer, parameter :: LA_ARRAY_SIZE_ERROR = 102
    !> An error flag denoting a singular matrix.
    integer, parameter :: LA_SINGULAR_MATRIX_ERROR = 103
end module
