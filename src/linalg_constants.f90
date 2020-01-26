! linalg_constants.f90

!> @brief Provides a set of constants and error flags for the library.
module linalg_constants
    use, intrinsic :: iso_fortran_env, only : int32
    implicit none

! ******************************************************************************
! CONSTANTS
! ------------------------------------------------------------------------------
    !> @brief Defines no operation should be performed on the matrix.
    integer(int32), parameter :: NO_OPERATION = 0
    !> @brief Defines a transpose operation.
    integer(int32), parameter :: TRANSPOSE = 1
    !> @brief Defines a Hermitian transpose operation for a complex-valued matrix.
    integer(int32), parameter :: HERMITIAN_TRANSPOSE = 2

! ******************************************************************************
! ERROR FLAGS
! ------------------------------------------------------------------------------
    !> A flag denoting no error condition.
    integer(int32), parameter :: LA_NO_ERROR = 0
    !> An error flag denoting an invalid input.
    integer(int32), parameter :: LA_INVALID_INPUT_ERROR = 101
    !> An error flag denoting an improperly sized array.
    integer(int32), parameter :: LA_ARRAY_SIZE_ERROR = 102
    !> An error flag denoting a singular matrix.
    integer(int32), parameter :: LA_SINGULAR_MATRIX_ERROR = 103
    !> An error flag denoting an issue with the matrix format.
    integer(int32), parameter :: LA_MATRIX_FORMAT_ERROR = 104
    !> An error flag denoting that there is insufficient memory available.
    integer(int32), parameter :: LA_OUT_OF_MEMORY_ERROR = 105
    !> An error flag denoting a convergence failure.
    integer(int32), parameter :: LA_CONVERGENCE_ERROR = 106
    !> An error resulting from an invalid operation.
    integer(int32), parameter :: LA_INVALID_OPERATION_ERROR = 107
end module
