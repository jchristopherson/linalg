module linalg_errors
    use iso_fortran_env
    use ferror
    implicit none

    integer(int32), parameter :: LA_NO_ERROR = 0
        !! An flag denoting no error condition.
    integer(int32), parameter :: LA_INVALID_INPUT_ERROR = 101
        !! An error flag denoting an invalid input.
    integer(int32), parameter :: LA_ARRAY_SIZE_ERROR = 102
        !! An error flag denoting an improperly sized array.
    integer(int32), parameter :: LA_SINGULAR_MATRIX_ERROR = 103
        !! An error flag denoting a singular matrix.
    integer(int32), parameter :: LA_MATRIX_FORMAT_ERROR = 104
        !! An error flag denoting an issue with the matrix format.
    integer(int32), parameter :: LA_OUT_OF_MEMORY_ERROR = 105
        !! An error flag denoting that there is insufficient memory available.
    integer(int32), parameter :: LA_CONVERGENCE_ERROR = 106
        !! An error flag denoting a convergence failure.
    integer(int32), parameter :: LA_INVALID_OPERATION_ERROR = 107
        !! An error resulting from an invalid operation.

contains
    subroutine report_memory_error(fcn, err, flag)
        !! Reports a memory allocation error.
        character(len=*), intent(in) :: fcn
            !! The name of the function that failed.
        class(errors), intent(inout) :: err
            !! The error object to be updated.
        integer(int32), intent(in) :: flag
            !! The error flag.

        ! Local Variables
        character(len = 256) :: msg

        ! Construct the error message
        write(msg, 100) "Memory allocation failed in ", fcn, &
            " with code ", flag, "."
        call err%report_error(fcn, trim(msg), LA_OUT_OF_MEMORY_ERROR)

        ! Formatting
    100 format(A, A, A, I0, A)
    end subroutine

    subroutine report_array_size_error(fcn, err, name, expected, actual)
        !! Reports an array size error.
        character(len=*), intent(in) :: fcn
            !! The name of the function that failed.
        class(errors), intent(inout) :: err
            !! The error object to be updated.
        character(len=*), intent(in) :: name
            !! The name of the array.
        integer(int32), intent(in) :: expected
            !! The expected size of the array.
        integer(int32), intent(in) :: actual
            !! The actual size of the array.

        ! Local Variables
        character(len = 256) :: msg

        ! Construct the error message
        write(msg, 100) "Expected array size of ", expected, &
            " but received ", actual, " in ", fcn, "for array ", name, " ."
        call err%report_error(fcn, trim(msg), LA_ARRAY_SIZE_ERROR)

        ! Formatting
    100 format(A, I0, A, I0, A, A, A, A, A)
    end subroutine

    subroutine report_matrix_size_error(fcn, err, name, expectedRows, &
        expectedCols, actualRows, actualCols)
        !! Reports a matrix size error.
        character(len=*), intent(in) :: fcn
            !! The name of the function that failed.
        class(errors), intent(inout) :: err
            !! The error object to be updated.
        character(len=*), intent(in) :: name
            !! The name of the matrix.
        integer(int32), intent(in) :: expectedRows
            !! The expected number of rows in the matrix.
        integer(int32), intent(in) :: expectedCols
            !! The expected number of columns in the matrix.
        integer(int32), intent(in) :: actualRows
            !! The actual number of rows in the matrix.
        integer(int32), intent(in) :: actualCols
            !! The actual number of columns in the matrix.

        ! Local Variables
        character(len = 256) :: msg

        ! Construct the error message
        write(msg, 100) "Expected matrix size of ", expectedRows, " x ", &
            expectedCols, " but received ", actualRows, " x ", actualCols, &
            " in ", fcn, "for matrix ", name, " ."
        call err%report_error(fcn, trim(msg), LA_ARRAY_SIZE_ERROR)

        ! Formatting
    100 format(A, I0, A, I0, A, I0, A, I0, A, A, A, A, A)
    end subroutine

    subroutine report_inner_matrix_dimension_error(fcn, err, name1, name2, &
        expected, actual)
        !! Reports an inner matrix dimension error.
        character(len=*), intent(in) :: fcn
            !! The name of the function that failed.
        class(errors), intent(inout) :: err
            !! The error object to be updated.
        character(len=*), intent(in) :: name1
            !! The name of the first matrix.
        character(len=*), intent(in) :: name2
            !! The name of the second matrix.
        integer(int32), intent(in) :: expected
            !! The expected inner dimension.
        integer(int32), intent(in) :: actual
            !! The actual inner dimension.

        ! Local Variables
        character(len = 256) :: msg

        ! Construct the error message
        write(msg, 100) "Expected inner matrix dimension of ", expected, &
            " but received ", actual, " in ", fcn, "for matrices ", name1, &
            " and ", name2, " ."
        call err%report_error(fcn, trim(msg), LA_ARRAY_SIZE_ERROR)

        ! Formatting
    100 format(A, I0, A, I0, A, I0, A, A, A, A, A)
    end subroutine

    subroutine report_square_matrix_error(fcn, err, name, expectedSize, &
        actualRows, actualCols)
        !! Reports an error where a square matrix was expected but a non-square 
        !! matrix was provided.
        character(len=*), intent(in) :: fcn
            !! The name of the function that failed.
        class(errors), intent(inout) :: err
            !! The error object to be updated.
        character(len=*), intent(in) :: name
            !! The name of the matrix.
        integer(int32), intent(in) :: expectedSize
            !! The expected size of the square matrix.
        integer(int32), intent(in) :: actualRows
            !! The actual number of rows in the matrix.
        integer(int32), intent(in) :: actualCols
            !! The actual number of columns in the matrix.

        ! Local Variables
        character(len = 256) :: msg

        ! Construct the error message
        write(msg, 100) "Expected square matrix of size ", expectedSize, &
            " but received ", actualRows, " x ", actualCols, " in ", fcn, &
            "for matrix ", name, " ."
        call err%report_error(fcn, trim(msg), LA_ARRAY_SIZE_ERROR)

        ! Formatting
    100 format(A, I0, A, I0, A, I0, A, I0, A, A, A, A, A)
    end subroutine

    subroutine report_singular_matrix_warning(fcn, err, row)
        !! Reports a singular matrix error.
        character(len=*), intent(in) :: fcn
            !! The name of the function that failed.
        class(errors), intent(inout) :: err
            !! The error object to be updated.
        integer(int32), intent(in) :: row
            !! The row index where the singularity issue was first encountered.

        ! Local Variables
        character(len = 256) :: msg

        ! Write the error message
        write(msg, 100) &
            "A singular matrix was encountered with the issue found at row ", &
            row, " of the matrix."
        call err%report_warning(fcn, trim(msg), LA_SINGULAR_MATRIX_ERROR)

        ! Formatting
    100 format(A, I0, A)
    end subroutine
end module