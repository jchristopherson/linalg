module linalg_types
    use iso_fortran_env
    implicit none
    private
    public :: matrix_factorization
    public :: size

    type, abstract :: matrix_factorization
        !! Defines a matrix factorization container.
    contains
        procedure(get_logical_value), deferred, public, pass :: is_factored
            !! Checks to see if the matrix has been factored.
        procedure(get_matrix_size), deferred, public, pass :: size
            !! Gets the size of the factored matrix.
        procedure(matrix_factor), deferred, public, pass :: factor
            !! Factors the matrix.
        procedure(matrix_solve_multiple_rhs), deferred, public, pass :: solve_multiple_rhs
            !! Solves the factored system with multiple right-hand-side vectors.
        procedure(matrix_solve_single_rhs), deferred, public, pass :: solve_single_rhs
            !! Solves the factored system with a single right-hand-side vector.
        generic, public :: solve => solve_multiple_rhs, solve_single_rhs
            !! Solves the factored system.
    end type

    abstract interface
        pure function get_logical_value(this) result(rst)
            !! Retrieves a logical parameter value.
            import matrix_factorization
            class(matrix_factorization), intent(in) :: this
                !! The [[matrix_factorization]] object.
            logical :: rst
                !! The logical result.
        end function

        pure function get_matrix_size(this, dim) result(rst)
            !! Retrieves the size of the requested matrix dimension.
            import matrix_factorization
            class(matrix_factorization), intent(in) :: this
                !! The [[matrix_factorization]] object.
            integer, intent(in) :: dim
                !! The dimension to query.
            integer :: rst
                !! The requested dimension size.
        end function

        subroutine matrix_factor(this, flag)
            !! Factors the matrix.
            import matrix_factorization
            class(matrix_factorization), intent(inout) :: this
                !! The [[matrix_factorization]] object.
            integer, intent(out), optional :: flag
                !! An error flag.
        end subroutine

        function matrix_solve_multiple_rhs(this, b, flag) result(rst)
            !! Solves the linear system \(A X = B \) where \(B\)
            !! is a matrix defining multiple right-hand-side vectors.
            use iso_fortran_env, only : real64
            import matrix_factorization
            class(matrix_factorization), intent(inout) :: this
                !! The [[matrix_factorization]] object.
            real(real64), intent(in), dimension(:,:) :: b
                !! The right-hand-side matrix.
            integer, intent(out), optional :: flag
                !! An error flag.
            real(real64), allocatable, dimension(:,:) :: rst
                !! The solution matrix.
        end function

        function matrix_solve_single_rhs(this, b, flag) result(rst)
            !! Solves the linear system \(A X = B \) where \(B\)
            !! is a matrix defining a single right-hand-side vector.
            use iso_fortran_env, only : real64
            import matrix_factorization
            class(matrix_factorization), intent(inout) :: this
                !! The [[matrix_factorization]] object.
            real(real64), intent(in), dimension(:) :: b
                !! The right-hand-side vector.
            integer, intent(out), optional :: flag
                !! An error flag.
            real(real64), allocatable, dimension(:) :: rst
                !! The solution vector.
        end function
    end interface

    interface size
        module procedure :: mf_size
    end interface

contains

pure function mf_size(x, dim) result(rst)
    class(matrix_factorization), intent(in) :: x
    integer, intent(in) :: dim
    integer :: rst
    rst = x%size(dim)
end function

end module