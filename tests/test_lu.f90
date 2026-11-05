! test_lu.f90
! SPDX-License-Identifier: MIT

! Tests the LU factorization routines
module test_lu
    use, intrinsic :: iso_fortran_env, only : real64
    use test_core
    use linalg
    use linalg_errors, only : LA_NO_ERROR, LA_MATRIX_NOT_SQUARE_ERROR, &
        LA_NOT_INITIALIZED_ERROR, LA_INVALID_OPERATION_ERROR, &
        LA_MATRIX_DIMENSION_ERROR, LA_INVALID_INPUT_ERROR, &
        LA_MATRIX_FORMAT_ERROR, LA_SINGULAR_MATRIX_ERROR, LA_CONVERGENCE_ERROR
    use fortran_test_helper
contains
! ******************************************************************************
! LU FACTORIZATION TEST
! ------------------------------------------------------------------------------
    function test_lu_factor() result(rst)
        ! Parameters
        integer, parameter :: n = 75

        ! Local Variables
        real(real64), dimension(n, n) :: a
        real(real64), allocatable, dimension(:,:) :: l, u, p
        logical :: rst

        ! Initialization
        rst = .true.
        call create_random_array(a)

        ! Compute the factorization
        call lu_factor(a, l = l, u = u, p = p)

        ! Determine if P * A = L * U
        if (.not.assert(matmul(p, a), matmul(l, u), tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: LU Factorization Test"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_lu_solve() result(rst)
        ! Parameters
        integer, parameter :: n = 75
        integer, parameter :: nrhs = 20

        ! Local Variables
        real(real64), dimension(n, n) :: a
        real(real64), dimension(n, nrhs) :: b, x
        integer, allocatable :: ipvt(:)
        real(real64), allocatable :: lu(:,:)
        logical :: rst

        ! Initialization
        rst = .true.
        call create_random_array(a)
        call create_random_array(b)

        ! Factor A
        call lu_factor(a, ipvt = ipvt, lu = lu)

        ! Solve for X
        x = solve_lu(lu, ipvt, b)

        ! Test by determining if A * X = B
        if (.not.assert(matmul(a, x), b, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: LU Factorization & Solution Test"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_lu_factor_cmplx() result(rst)
        ! Parameters
        integer, parameter :: n = 75

        ! Local Variables
        complex(real64), dimension(n, n) :: a
        complex(real64), allocatable, dimension(:,:) :: l, u
        real(real64), allocatable, dimension(:,:) :: p
        logical :: rst

        ! Initialization
        rst = .true.
        call create_random_array(a)

        ! Compute the factorization
        call lu_factor(a, l = l, u = u, p = p)

        ! Determine if P * A = L * U
        if (.not.assert(matmul(p, a), matmul(l, u), tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Complex-Valued LU Factorization Test"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_lu_solve_cmplx() result(rst)
        ! Parameters
        integer, parameter :: n = 75
        integer, parameter :: nrhs = 20

        ! Local Variables
        complex(real64), dimension(n, n) :: a
        complex(real64), dimension(n, nrhs) :: b, x
        integer :: i, j
        integer, allocatable :: ipvt(:)
        complex(real64), allocatable :: lu(:,:)
        logical :: rst

        ! Initialization
        rst = .true.
        call create_random_array(a)
        call create_random_array(b)

        ! Factor A
        call lu_factor(a, ipvt = ipvt, lu = lu)

        ! Solve for X
        x = solve_lu(lu, ipvt, b)

        ! Test by determining if A * X = B
        if (.not.assert(matmul(a, x), b, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Complex-Valued LU Factorization & Solution Test"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_lu_factorization_type() result(rst)
        integer, parameter :: n = 3
        real(real64), dimension(n, n), parameter :: a = reshape([ &
            0.0d0, 1.0d0, 0.0d0, &
            2.0d0, 0.0d0, 1.0d0, &
            0.0d0, 1.0d0, 1.0d0], [n, n])
        real(real64), dimension(n, 2), parameter :: x_expected = reshape([ &
            1.0d0, 2.0d0, 3.0d0, &
            4.0d0, 5.0d0, 6.0d0], [n, 2])

        type(lu_factorization) :: fac, fac_form, fac_det, fac_inverse
        type(lu_factorization) :: fac_multiple_rhs, fac_single_rhs, fac_rect
        real(real64), allocatable, dimension(:,:) :: l, u, p, x_multiple, ainv
        real(real64), dimension(n, n) :: identity_mtx
        real(real64), allocatable, dimension(:) :: x_single
        real(real64), dimension(n, 2) :: b_multiple
        real(real64), dimension(n) :: b_single
        real(real64) :: det
        integer :: flag
        logical :: rst

        rst = .true.
        b_multiple = matmul(a, x_expected)
        b_single = matmul(a, x_expected(:,1))
        identity_mtx = identity(n)

        fac = lu_factorization(a)
        if (fac%size(1) /= n .or. fac%size(2) /= n .or. fac%is_factored()) then
            rst = .false.
            print '(A)', "Test Failed: LU Factorization Type Initialization"
        end if

        call fac%factor(flag)
        if (flag /= LA_NO_ERROR .or. .not.fac%is_factored()) then
            rst = .false.
            print '(A)', "Test Failed: LU Factorization Type Factor"
        end if
        call fac%factor(flag)
        if (flag /= LA_NO_ERROR .or. .not.fac%is_factored()) then
            rst = .false.
            print '(A)', "Test Failed: LU Factorization Type Repeated Factor"
        end if

        fac_form = lu_factorization(a)
        call fac_form%get_factored_form(l, u, p, flag)
        if (flag /= LA_NO_ERROR .or. .not.fac_form%is_factored()) then
            rst = .false.
            print '(A)', "Test Failed: LU Factorization Type Factored Form Status"
        else if (.not.is_mtx_equal(matmul(p, a), matmul(l, u), REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: LU Factorization Type Factored Form"
        end if
        fac_det = lu_factorization(a)
        det = fac_det%determinant(flag)
        if (flag /= LA_NO_ERROR .or. abs(det + 2.0d0) > REAL64_TOL .or. &
            .not.fac_det%is_factored()) then
            rst = .false.
            print '(A)', "Test Failed: LU Factorization Type Determinant"
        end if

        fac_inverse = lu_factorization(a)
        call fac_inverse%inverse(ainv, flag)
        if (flag /= LA_NO_ERROR) then
            rst = .false.
            print '(A)', "Test Failed: LU Factorization Type Inverse Status"
        else if (.not.is_mtx_equal(matmul(a, ainv), identity_mtx, REAL64_TOL) .or. &
            .not.is_mtx_equal(matmul(ainv, a), identity_mtx, REAL64_TOL) .or. &
            .not.fac_inverse%is_factored()) then
            rst = .false.
            print '(A)', "Test Failed: LU Factorization Type Inverse"
        end if

        fac_multiple_rhs = lu_factorization(a)
        x_multiple = fac_multiple_rhs%solve(b_multiple, flag)
        if (flag /= LA_NO_ERROR .or. &
            .not.is_mtx_equal(x_multiple, x_expected, REAL64_TOL) .or. &
            .not.fac_multiple_rhs%is_factored()) then
            rst = .false.
            print '(A)', "Test Failed: LU Factorization Type Multiple RHS Solve"
        end if

        fac_single_rhs = lu_factorization(a)
        x_single = fac_single_rhs%solve(b_single, flag)
        if (flag /= LA_NO_ERROR .or. &
            .not.is_mtx_equal(x_single, x_expected(:,1), REAL64_TOL) .or. &
            .not.fac_single_rhs%is_factored()) then
            rst = .false.
            print '(A)', "Test Failed: LU Factorization Type Single RHS Solve"
        end if

        fac_rect = lu_factorization(a(:,1:2))
        if (fac_rect%size(1) /= n .or. fac_rect%size(2) /= n - 1) then
            rst = .false.
            print '(A)', "Test Failed: LU Factorization Type Rectangular Size"
        end if
        det = fac_rect%determinant(flag)
        if (flag /= LA_MATRIX_NOT_SQUARE_ERROR) then
            rst = .false.
            print '(A)', "Test Failed: LU Factorization Type Non-Square Determinant Validation"
        end if
        call fac_rect%inverse(ainv, flag)
        if (flag /= LA_MATRIX_NOT_SQUARE_ERROR .or. allocated(ainv)) then
            rst = .false.
            print '(A)', "Test Failed: LU Factorization Type Non-Square Inverse Validation"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_sparse_lu_factorization_type() result(rst)
        !! Verifies complete solves, not merely incomplete-factor application,
        !! together with lifecycle, storage validation, and error guarantees.
        logical :: rst
            !! True when every solution, residual, and status check succeeds.
        integer, parameter :: order = 8
        real(real64), parameter :: tolerance = 1.0d-10
        real(real64) :: matrix(order,order), expected(order,3), rhs(order,3)
        real(real64) :: in_place(order,3), vector(order), prepared(order)
        real(real64) :: limited_rhs(order,2), limited_copy(order,2)
        real(real64), allocatable :: solution(:), solutions(:,:)
        type(csr_matrix) :: csr, malformed
        type(msr_matrix) :: factors
        integer, allocatable :: upper_rows(:)
        type(sparse_lu_factorization) :: fac, lazy, invalid, limited, scalar
        integer :: index, flag

        rst = .true.
        matrix = 0.0d0
        do index = 1, order
            matrix(index,index) = 4.0d0
            if (index < order) then
                matrix(index,index+1) = -1.0d0
                matrix(index+1,index) = -1.0d0
            end if
        end do
        expected(:,1) = [(real(index, real64), index = 1, order)]
        expected(:,2) = 0.0d0
        expected(:,3) = 1.0d0 - expected(:,1) / 3.0d0
        rhs = matmul(matrix, expected)
        csr = dense_to_csr(matrix)

        rst = fac%size(1) == 0 .and. fac%size(2) == 0 .and. rst
        rst = .not.fac%is_factored() .and. rst
        call fac%factor(flag)
        rst = flag == LA_NOT_INITIALIZED_ERROR .and. rst
        vector = rhs(:,1)
        call fac%solve_in_place(vector, flag)
        rst = flag == LA_NOT_INITIALIZED_ERROR .and. rst
        rst = assert(vector, rhs(:,1)) .and. rst

        fac = sparse_lu_factorization(csr, droptol = 1.0d0, &
            tol = tolerance, im = order)
        rst = fac%size(1) == order .and. fac%size(2) == order .and. rst
        rst = size(fac, 1) == order .and. .not.fac%is_factored() .and. rst
        call fac%solve_prepared(rhs(:,1), prepared, flag)
        rst = flag == LA_INVALID_OPERATION_ERROR .and. rst
        rst = all(prepared == 0.0d0) .and. rst
        csr%values = 0.0d0
        call fac%factor(flag)
        rst = flag == LA_NO_ERROR .and. fac%is_factored() .and. rst
        call fac%factor(flag)
        rst = flag == LA_NO_ERROR .and. fac%is_factored() .and. rst
        call fac%get_factored_form(factors, upper_rows, flag)
        rst = flag == LA_NO_ERROR .and. rst
        solution = solve_lu(factors, upper_rows, rhs(:,1))
        rst = norm2(matmul(matrix, solution) - rhs(:,1)) > 1.0d-3 .and. rst
        call fac%solve_prepared(rhs(:,1), prepared, flag)
        rst = flag == LA_NO_ERROR .and. rst
        rst = assert(prepared, expected(:,1), tol = 1.0d-8) .and. rst

        solutions = fac%solve(rhs, flag)
        rst = flag == LA_NO_ERROR .and. rst
        rst = assert(solutions, expected, tol = 1.0d-8) .and. rst
        rst = assert(matmul(matrix, solutions), rhs, tol = 1.0d-8) .and. rst
        rst = assert(rhs, matmul(matrix, expected)) .and. rst
        vector = rhs(:,1)
        call fac%solve_in_place(vector, flag)
        rst = flag == LA_NO_ERROR .and. rst
        rst = assert(vector, expected(:,1), tol = 1.0d-8) .and. rst
        in_place = rhs
        call fac%solve_in_place(in_place, flag)
        rst = flag == LA_NO_ERROR .and. rst
        rst = assert(in_place, expected, tol = 1.0d-8) .and. rst
        solution = fac%solve(rhs(:,1))
        rst = assert(solution, expected(:,1), tol = 1.0d-8) .and. rst
        call fac%solve_prepared(rhs(:,1), prepared(:order-1), flag)
        rst = flag == LA_MATRIX_DIMENSION_ERROR .and. rst
        vector = rhs(:,1)
        call fac%solve_in_place(vector(:order-1), flag)
        rst = flag == LA_MATRIX_DIMENSION_ERROR .and. rst
        rst = assert(vector, rhs(:,1)) .and. rst

        csr = dense_to_csr(matrix)
        lazy = sparse_lu_factorization(csr, droptol = 1.0d0, tol = tolerance)
        solution = lazy%solve(rhs(:,1), flag)
        rst = flag == LA_NO_ERROR .and. lazy%is_factored() .and. rst
        rst = assert(solution, expected(:,1), tol = 1.0d-8) .and. rst
        call lazy%update(dense_to_csr(2.0d0 * matrix), tol = tolerance)
        rst = .not.lazy%is_factored() .and. rst
        call lazy%solve_prepared(rhs(:,1), prepared, flag)
        rst = flag == LA_INVALID_OPERATION_ERROR .and. rst
        solution = lazy%solve(rhs(:,1), flag)
        rst = flag == LA_NO_ERROR .and. rst
        rst = assert(solution, 0.5d0 * expected(:,1), tol = 1.0d-8) .and. rst

        limited = sparse_lu_factorization(csr, droptol = 1.0d0, &
            tol = tolerance, im = 1, maxits = 1)
        vector = rhs(:,1)
        call limited%solve_in_place(vector, flag)
        rst = flag == LA_CONVERGENCE_ERROR .and. rst
        rst = assert(vector, rhs(:,1)) .and. rst
        limited_rhs(:,1) = sin([(real(index, real64), index = 1, order)] * &
            acos(-1.0d0) / real(order+1, real64))
        limited_rhs(:,1) = matmul(matrix, limited_rhs(:,1))
        limited_rhs(:,2) = rhs(:,1)
        limited_copy = limited_rhs
        call limited%solve_in_place(limited_rhs, flag)
        rst = flag == LA_CONVERGENCE_ERROR .and. rst
        rst = assert(limited_rhs, limited_copy) .and. rst

        invalid = sparse_lu_factorization(dense_to_csr(matrix(:,:order-1)))
        call invalid%factor(flag)
        rst = flag == LA_MATRIX_NOT_SQUARE_ERROR .and. rst
        invalid = sparse_lu_factorization(csr, tol = 0.0d0)
        call invalid%factor(flag)
        rst = flag == LA_INVALID_INPUT_ERROR .and. rst
        invalid = sparse_lu_factorization(csr, droptol = -1.0d0)
        call invalid%factor(flag)
        rst = flag == LA_INVALID_INPUT_ERROR .and. rst
        invalid = sparse_lu_factorization(csr, im = 51)
        call invalid%factor(flag)
        rst = flag == LA_INVALID_INPUT_ERROR .and. rst
        invalid = sparse_lu_factorization(csr, maxits = 0)
        call invalid%factor(flag)
        rst = flag == LA_INVALID_INPUT_ERROR .and. rst
        malformed = csr
        malformed%column_indices(1) = order + 1
        invalid = sparse_lu_factorization(malformed)
        call invalid%factor(flag)
        rst = flag == LA_MATRIX_FORMAT_ERROR .and. rst
        malformed = csr
        malformed%row_indices(1) = 0
        invalid = sparse_lu_factorization(malformed)
        call invalid%factor(flag)
        rst = flag == LA_MATRIX_FORMAT_ERROR .and. rst
        malformed = csr
        malformed%column_indices(2) = malformed%column_indices(1)
        invalid = sparse_lu_factorization(malformed)
        call invalid%factor(flag)
        rst = flag == LA_MATRIX_FORMAT_ERROR .and. rst
        malformed = csr
        malformed%values = 0.0d0
        invalid = sparse_lu_factorization(malformed)
        call invalid%factor(flag)
        rst = flag == LA_SINGULAR_MATRIX_ERROR .and. .not.invalid%is_factored() .and. rst

        scalar = sparse_lu_factorization(dense_to_csr(reshape([2.0d0], [1,1])))
        solution = scalar%solve([6.0d0], flag)
        rst = flag == LA_NO_ERROR .and. rst
        rst = assert(solution, [3.0d0], tol = 1.0d-12) .and. rst
        call scalar%get_factored_form(factors, upper_rows, flag)
        rst = flag == LA_NO_ERROR .and. size(factors%values) == 2 .and. rst
        if (.not.rst) print '(A)', 'Test Failed: Sparse LU Factorization Type'
    end function

! ------------------------------------------------------------------------------
end module
