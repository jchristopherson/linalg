! test_eigen.f90
! SPDX-License-Identifier: MIT

! Tests the eigenvalue/eigenvector routines
module test_eigen
    use, intrinsic :: iso_fortran_env, only : real64
    use linalg
    use test_core
    use fortran_test_helper
    implicit none
    private :: block_test_operator, block_test_solver

    type, extends(linear_operator) :: block_test_operator
        !! Matrix-free diagonal operator with an optional trailing 2-by-2 block.
        real(real64), allocatable :: diagonal(:)
            !! Diagonal coefficients; allocation determines both dimensions.
        real(real64) :: coupling = 0.0d0
            !! Negative upper and positive lower off-diagonal block entries.
    contains
        procedure :: size => block_operator_size
            !! Returns the number of rows or columns.
        procedure :: apply => block_operator_apply
            !! Applies the diagonal and trailing block without sparse storage.
    end type

    type, extends(linear_solver) :: block_test_solver
        !! Exact prepared solver for the same diagonal/trailing-block structure.
        real(real64), allocatable :: diagonal(:)
            !! Nonzero diagonal coefficients of the prepared system.
        real(real64) :: coupling = 0.0d0
            !! Negative upper and positive lower off-diagonal block entries.
    contains
        procedure :: size => block_solver_size
            !! Returns the square coefficient-matrix dimensions.
        procedure :: solve => block_solver_solve
            !! Solves the diagonal system and optional trailing block exactly.
    end type
contains
    pure function block_operator_size(this, dim) result(rst)
        !! Queries the dimensions of the matrix-free test operator.
        class(block_test_operator), intent(in) :: this
            !! Initialized operator with allocated diagonal coefficients.
        integer, intent(in) :: dim
            !! Dimension to query: 1 for rows, 2 for columns.
        integer :: rst
            !! Matrix order, or zero for an unsupported dimension.
        rst = 0
        if (dim == 1 .or. dim == 2) rst = size(this%diagonal)
    end function

    subroutine block_operator_apply(this, x, y)
        !! Computes the diagonal product and optional trailing-block product.
        class(block_test_operator), intent(in) :: this
            !! Initialized read-only operator.
        real(real64), intent(in) :: x(:)
            !! Input vector with THIS%SIZE(2) elements; must not overlap Y.
        real(real64), intent(out) :: y(:)
            !! Caller-owned product vector with THIS%SIZE(1) elements.
        integer :: order

        order = size(this%diagonal)
        y = this%diagonal * x
        y(order-1) = y(order-1) - this%coupling * x(order)
        y(order) = y(order) + this%coupling * x(order-1)
    end subroutine

    pure function block_solver_size(this, dim) result(rst)
        !! Queries the dimensions of the prepared test system.
        class(block_test_solver), intent(in) :: this
            !! Initialized solver with allocated diagonal coefficients.
        integer, intent(in) :: dim
            !! Dimension to query: 1 for rows, 2 for columns.
        integer :: rst
            !! Matrix order, or zero for an unsupported dimension.
        rst = 0
        if (dim == 1 .or. dim == 2) rst = size(this%diagonal)
    end function

    subroutine block_solver_solve(this, rhs, x)
        !! Solves a diagonal system with an optional nonsingular trailing block.
        class(block_test_solver), intent(in) :: this
            !! Prepared read-only coefficients of the system to solve.
        real(real64), intent(in) :: rhs(:)
            !! Right-hand side with THIS%SIZE(1) elements; must not overlap X.
        real(real64), intent(out) :: x(:)
            !! Caller-owned solution vector with THIS%SIZE(2) elements.
        integer :: order
        real(real64) :: determinant

        order = size(this%diagonal)
        x = rhs / this%diagonal
        if (this%coupling == 0.0d0) return
        determinant = this%diagonal(order-1) * this%diagonal(order) + &
            this%coupling**2
        x(order-1) = (this%diagonal(order) * rhs(order-1) + &
            this%coupling * rhs(order)) / determinant
        x(order) = (this%diagonal(order-1) * rhs(order) - &
            this%coupling * rhs(order-1)) / determinant
    end subroutine

    function test_operator_eigen() result(rst)
        !! Checks standard, generalized, and shift-invert matrix-free eigenpairs
        !! against independent dense products and known selected eigenvalues.
        logical :: rst
            !! True when all eigenvalues, residuals, and values-only calls agree.
        integer, parameter :: order = 20, count = 3
        real(real64), parameter :: tolerance = 1.0d-7
        type(block_test_operator) :: op, mass
        type(block_test_solver) :: solver
        real(real64) :: matrix(order,order), mass_matrix(order,order)
        real(real64) :: expected(count), shift
        real(real64), allocatable :: vals(:), vecs(:,:), values_only(:)
        complex(real64), allocatable :: cvals(:), cvecs(:,:), cvalues_only(:)
        integer :: index, mode

        rst = .true.
        op%diagonal = [(real(index, real64), index = 1, order)]
        mass%diagonal = 1.0d0 + 0.01d0 * op%diagonal
        matrix = 0.0d0
        mass_matrix = 0.0d0
        do index = 1, order
            matrix(index,index) = op%diagonal(index)
            mass_matrix(index,index) = mass%diagonal(index)
        end do

        call eigen(op, count, vals, vecs, ncv = 12, tol = 1.0d-10)
        rst = assert(size(vals), count) .and. rst
        rst = assert(matmul(matrix, vecs), &
            vecs * spread(vals, 1, order), tol = tolerance) .and. rst
        expected = op%diagonal(order-count+1:order)
        do index = 1, count
            rst = minval(abs(vals - expected(index))) < tolerance .and. rst
        end do
        call eigen(op, count, values_only, ncv = 12, tol = 1.0d-10)
        rst = assert(vals, values_only, tol = tolerance) .and. rst

        do mode = 1, 2
            if (mode == 1) then
                solver%diagonal = mass%diagonal
                expected = op%diagonal(order-count+1:order) / &
                    mass%diagonal(order-count+1:order)
                call eigen(op, mass, solver, count, vals, vecs, &
                    ncv = 12, tol = 1.0d-10)
            else
                shift = 7.4d0
                solver%diagonal = op%diagonal - shift * mass%diagonal
                expected = op%diagonal(7:9) / mass%diagonal(7:9)
                call eigen(op, mass, solver, count, vals, vecs, &
                    sigma = shift, ncv = 12, tol = 1.0d-10)
            end if
            rst = assert(size(vals), count) .and. rst
            rst = assert(matmul(matrix, vecs), &
                matmul(mass_matrix, vecs) * spread(vals, 1, order), &
                tol = tolerance) .and. rst
            do index = 1, count
                rst = minval(abs(vals - expected(index))) < tolerance .and. rst
            end do
        end do

        op%diagonal(order) = op%diagonal(order-1)
        op%coupling = 30.0d0
        matrix(order,order) = op%diagonal(order)
        matrix(order-1,order) = -op%coupling
        matrix(order,order-1) = op%coupling
        call eigen(op, 1, cvals, cvecs, ncv = 12, tol = 1.0d-10)
        rst = assert(size(cvals), 2) .and. rst
        rst = assert(matmul(matrix, cvecs), &
            cvecs * spread(cvals, 1, order), tol = tolerance) .and. rst
        rst = minval(abs(cvals - cmplx(19.0d0, 30.0d0, real64))) &
            < tolerance .and. rst
        rst = minval(abs(cvals - cmplx(19.0d0, -30.0d0, real64))) &
            < tolerance .and. rst
        call eigen(op, 1, cvalues_only, ncv = 12, tol = 1.0d-10)
        rst = assert(cvals, cvalues_only, tol = tolerance) .and. rst

        mass%diagonal = 2.0d0
        mass_matrix = 0.0d0
        do index = 1, order
            mass_matrix(index,index) = mass%diagonal(index)
        end do
        do mode = 1, 2
            if (mode == 1) then
                solver%diagonal = mass%diagonal
                solver%coupling = 0.0d0
                call eigen(op, mass, solver, 1, cvals, cvecs, &
                    ncv = 12, tol = 1.0d-10)
                rst = assert(size(cvals), 2) .and. rst
                rst = minval(abs(cvals - cmplx(9.5d0, 15.0d0, real64))) &
                    < tolerance .and. rst
                rst = minval(abs(cvals - cmplx(9.5d0, -15.0d0, real64))) &
                    < tolerance .and. rst
            else
                shift = 4.2d0
                solver%diagonal = op%diagonal - shift * mass%diagonal
                solver%coupling = op%coupling
                expected = [3.5d0, 4.0d0, 4.5d0]
                call eigen(op, mass, solver, count, cvals, cvecs, &
                    sigma = shift, ncv = 12, tol = 1.0d-10)
                rst = assert(size(cvals), count) .and. rst
                do index = 1, count
                    rst = minval(abs(cvals - expected(index))) &
                        < tolerance .and. rst
                end do
            end if
            rst = assert(matmul(matrix, cvecs), &
                matmul(mass_matrix, cvecs) * spread(cvals, 1, order), &
                tol = tolerance) .and. rst
        end do
        if (.not.rst) print '(A)', 'Test Failed: Matrix-Free ARPACK Eigenpairs'
    end function

! ******************************************************************************
! EIGEN ANALYSIS TESTS
! ------------------------------------------------------------------------------
    function test_eigen_symm() result(rst)
        ! Parameters
        integer, parameter :: n = 100

        ! Local Variables
        real(real64), dimension(n, n) :: a, x, y
        real(real64), allocatable, dimension(:) :: vals
        real(real64), allocatable, dimension(:,:) :: vecs
        logical :: rst

        ! Initialization
        rst = .true.
        call create_random_array(a, mtype = SYMMETRIC_MATRIX)

        ! Compute the eigenvalues and eigenvectors of A
        call eigen(a, vals, vecs)

        ! Compute vecs * vals, where vals is a diagonal matrix
        call diag_mtx_mult(.false., .false., 1.0d0, vals, vecs, 0.0d0, x)

        ! Test
        y = matmul(a, vecs)
        if (.not.is_mtx_equal(x, y, REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Symmetric Eigen Values"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_eigen_asymm() result(rst)
        ! Parameters
        integer, parameter :: n = 100

        ! Local Variables
        real(real64), dimension(n, n) :: a
        complex(real64), dimension(n, n) :: vmtx, x, y
        complex(real64), allocatable, dimension(:) :: vals, vals1
        complex(real64), allocatable, dimension(:,:) :: vecs
        integer :: i
        logical :: rst

        ! Initialization
        rst = .true.
        call create_random_array(a)
        vmtx = cmplx(0.0d0, 0.0d0, real64)

        ! Compute the eigenvalues and eigenvectors of A
        call eigen(a, vals, rvecs = vecs)

        ! Compute vecs * vals, where vals is a diagonal matrix
        do i = 1, n
            vmtx(i,i) = vals(i)
        end do
        x = matmul(vecs, vmtx)

        ! Test 1
        y = matmul(a, vecs)
        if (.not.assert(x, y, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Asymmetric Eigen Values Test 1"
        end if

        ! Compute just the eigenvalues
        call eigen(a, vals1)
        if (.not.assert(vals, vals1, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Asymmetric Eigen Values Test 2"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_eigen_asymm_cmplx() result(rst)
        ! Parameters
        integer, parameter :: n = 100

        ! Local Variables
        complex(real64), dimension(n, n) :: a, vmtx, x, y
        complex(real64), allocatable, dimension(:) :: vals, vals1
        complex(real64), allocatable, dimension(:,:) :: vecs
        integer :: i
        logical :: rst

        ! Initialization
        rst = .true.
        call create_random_array(a, mtype = SYMMETRIC_MATRIX)
        vmtx = cmplx(0.0d0, 0.0d0, real64)

        ! Compute the eigenvalues and eigenvectors of A
        call eigen(a, vals, vecs)

        ! Compute vecs * vals, where vals is a diagonal matrix
        do i = 1, n
            vmtx(i,i) = vals(i)
        end do
        x = matmul(vecs, vmtx)

        ! Test 1
        y = matmul(a, vecs)
        if (.not.assert(x, y, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Complex-Valued Asymmetric Eigen Values Test 1"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_eigen_gen() result(rst)
        ! Parameters
        integer, parameter :: n = 100

        ! Local Variables
        real(real64), dimension(n, n) :: a, b
        complex(real64), dimension(n, n) :: x, y
        complex(real64), allocatable, dimension(:) :: vals, vals2
        complex(real64), allocatable, dimension(:,:) :: vecs
        logical :: rst

        ! Initialization
        rst = .true.
        call create_random_array(a)
        call create_random_array(b)

        ! Test 1
        call eigen(a, b, vals, rvecs = vecs)

        ! Compute X = A * VECS
        x = matmul(a, vecs)

        ! And Y = B * (VECS * LAMBDA)
        call diag_mtx_mult(.false., LA_NO_OPERATION, 1.0d0, vals, vecs, 0.0d0, y)
        y = matmul(b, y)

        ! Check
        if (.not.assert(x, y, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Generalized Eigen Values Test 1"
        end if

        ! Test 2 - Eigenvalues Only
        call eigen(a, b, vals2)
        if (.not.assert(vals, vals2, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Generalized Eigen Values Test 2"
        end if
    end function

! ------------------------------------------------------------------------------
end module
