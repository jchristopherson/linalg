! test_lu.f90
! SPDX-License-Identifier: MIT

! Tests the LU factorization routines
module test_lu
    use, intrinsic :: iso_fortran_env, only : int32, real64
    use test_core
    use linalg
    use linalg_errors, only : LA_NO_ERROR, LA_MATRIX_NOT_SQUARE_ERROR
    use fortran_test_helper
contains
! ******************************************************************************
! LU FACTORIZATION TEST
! ------------------------------------------------------------------------------
    function test_lu_factor() result(rst)
        ! Parameters
        integer(int32), parameter :: n = 75

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
        integer(int32), parameter :: n = 75
        integer(int32), parameter :: nrhs = 20

        ! Local Variables
        real(real64), dimension(n, n) :: a
        real(real64), dimension(n, nrhs) :: b, x
        integer(int32), allocatable :: ipvt(:)
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
        integer(int32), parameter :: n = 75

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
        integer(int32), parameter :: n = 75
        integer(int32), parameter :: nrhs = 20

        ! Local Variables
        complex(real64), dimension(n, n) :: a
        complex(real64), dimension(n, nrhs) :: b, x
        integer(int32) :: i, j
        integer(int32), allocatable :: ipvt(:)
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
        integer(int32), parameter :: n = 3
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
        integer(int32) :: flag
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
end module
