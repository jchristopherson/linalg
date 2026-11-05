module linalg_linear_operators
    !! Defines storage-independent real operators and prepared linear solvers.
    !! ARPACK uses operator applications for matrix products and prepared
    !! solvers for generalized and shift-invert eigenvalue problems.
    use iso_fortran_env, only : real64
    implicit none
    private
    public :: linear_operator
    public :: linear_solver

    type, abstract :: linear_operator
        !! Defines a real-valued, possibly rectangular matrix-free operator.
        !! Concurrent applications require read-only operator data and separate
        !! output and scratch storage for each call. INTENT(IN) alone does not
        !! guarantee thread safety for pointer targets or external libraries.
    contains
        procedure(linear_operator_size_query), deferred, public, pass :: size
            !! Performs a size query on the operator.
        procedure(apply_linear_operator), deferred, public, pass :: apply
            !! Applies the linear operator.
    end type

    type, abstract :: linear_solver
        !! Defines a prepared solver for a square real-valued system.
        !! Preparation and factorization take place outside SOLVE. Each call
        !! must own its scratch storage to permit concurrent applications.
    contains
        procedure(linear_solver_size_query), deferred, public, pass :: size
            !! Queries the dimensions of the system being solved.
        procedure(solve_linear_system), deferred, public, pass :: solve
            !! Solves the prepared system without modifying the solver.
    end type

    abstract interface
        pure function linear_operator_size_query(this, dim) result(rst)
            !! Queries the number of rows or columns of the operator.
            import linear_operator
            class(linear_operator), intent(in) :: this
                !! The [[linear_operator]] object.
            integer, intent(in) :: dim
                !! The dimension to query: 1 for rows, 2 for columns.
            integer :: rst
                !! The nonnegative extent of the requested dimension.
        end function

        subroutine apply_linear_operator(this, x, y)
            !! Computes Y = A * X without allocating the output vector.
            !! X and Y must have SIZE(2) and SIZE(1) elements, respectively,
            !! and must not overlap. Any scratch storage must be call-local.
            use iso_fortran_env, only : real64
            import linear_operator
            class(linear_operator), intent(in) :: this
                !! The [[linear_operator]] object.
            real(real64), intent(in), dimension(:) :: x
                !! The input vector, with THIS%SIZE(2) elements.
            real(real64), intent(out), dimension(:) :: y
                !! Caller-owned output, with THIS%SIZE(1) elements. All
                !! elements must be assigned, without reading previous values.
        end subroutine

        pure function linear_solver_size_query(this, dim) result(rst)
            !! Queries the dimensions of the prepared coefficient matrix.
            import linear_solver
            class(linear_solver), intent(in) :: this
                !! The prepared linear solver, which is not modified.
            integer, intent(in) :: dim
                !! The dimension to query: 1 for rows, 2 for columns.
            integer :: rst
                !! The nonnegative extent. Both dimensions must agree for
                !! the square systems required by the ARPACK interface.
        end function

        subroutine solve_linear_system(this, rhs, x)
            !! Solves C * X = RHS for the prepared coefficient matrix C.
            !! RHS and X must have SIZE(1) and SIZE(2) elements, respectively,
            !! and must not overlap. Inner solves must be accurate enough for
            !! the outer iteration; a preconditioner alone is not sufficient.
            import linear_solver, real64
            class(linear_solver), intent(in) :: this
                !! The prepared solver. Coefficients and factors are read-only;
                !! implementations must not hide mutable shared scratch state.
            real(real64), intent(in) :: rhs(:)
                !! Right-hand side, with THIS%SIZE(1) elements.
            real(real64), intent(out) :: x(:)
                !! Caller-owned solution, with THIS%SIZE(2) elements. All
                !! elements must be assigned before returning.
        end subroutine
    end interface
end module