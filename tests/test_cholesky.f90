! test_cholesky.f90
! SPDX-License-Identifier: MIT

! Tests the Cholesky factorization/solution operations
module test_cholesky
    use, intrinsic :: iso_fortran_env, only : int32, real64
    use test_core
    use linalg
    use linalg_errors, only : LA_NO_ERROR, LA_MATRIX_NOT_SQUARE_ERROR, &
        LA_POSITIVE_DEFINITE_ERROR
    use fortran_test_helper
    implicit none
contains
! ******************************************************************************
! CHOLESKY FACTORIZATION TESTS
! ------------------------------------------------------------------------------
    function test_cholesky_factor() result(rst)
        ! Parameters
        integer(int32), parameter :: n = 100
        integer(int32), parameter :: nrhs = 20

        ! Local Variables
        real(real64), dimension(n, n) :: a, u, l
        real(real64), dimension(n, nrhs) :: b, x
        logical :: rst

        ! Initialization
        rst = .true.
        call create_random_array(a, mtype = POSITIVE_DEFINITE_MATRIX)
        call create_random_array(b)

        ! Test 1: A = L * L**T
        l = cholesky_factor(a, .false.)
        if (.not.assert(a, matmul(l, transpose(l)), tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Cholesky Factorization Test 1"
        end if

        ! Test 2: A = U**T * U
        u = cholesky_factor(a, .true.)
        if (.not.assert(a, matmul(transpose(u), u), tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Cholesky Factorization Test 2"
        end if

        ! Test 3: Solve L*L**T * X = B
        x = solve_cholesky(.false., l, b)
        if (.not.assert(matmul(a, x), b, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Cholesky Factorization Test 3"
        end if

        ! Test 4: Solve U**T * U * X = B
        x = solve_cholesky(.true., u, b)
        if (.not.assert(matmul(a, x), b, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Cholesky Factorization Test 4"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_cholesky_factorization_type() result(rst)
        integer(int32), parameter :: n = 3
        real(real64), dimension(n, n), parameter :: a = reshape([ &
            4.0d0, 1.0d0, 1.0d0, &
            1.0d0, 3.0d0, 0.0d0, &
            1.0d0, 0.0d0, 2.0d0], [n, n])
        real(real64), dimension(n, 2), parameter :: x_expected = reshape([ &
            1.0d0, 2.0d0, 3.0d0, &
            -1.0d0, 0.0d0, 2.0d0], [n, 2])
        real(real64), dimension(2, 3), parameter :: a_nonsquare = reshape([ &
            1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0, 6.0d0], [2, 3])
        real(real64), dimension(2, 2), parameter :: a_not_positive_definite = reshape([ &
            1.0d0, 2.0d0, 2.0d0, 1.0d0], [2, 2])

        type(cholesky_factorization) :: fac, fac_lower, fac_single, fac_multiple_in_place
        type(cholesky_factorization) :: fac_single_in_place, fac_invalid, fac_not_positive_definite
        real(real64), allocatable, dimension(:,:) :: r, x_multiple
        real(real64), allocatable, dimension(:) :: x_single
        real(real64), dimension(n, 2) :: b_multiple, b_multiple_in_place
        real(real64), dimension(n) :: b_single, b_single_in_place
        integer(int32) :: flag
        logical :: rst

        rst = .true.
        b_multiple = matmul(a, x_expected)
        b_single = matmul(a, x_expected(:,1))

        fac = cholesky_factorization(a, flag = flag)
        if (flag /= LA_NO_ERROR .or. fac%size(1) /= n .or. &
            fac%size(2) /= n .or. fac%is_factored() .or. .not.fac%is_upper()) then
            rst = .false.
            print '(A)', "Test Failed: Cholesky Factorization Type Initialization"
        end if

        call fac%get_factored_form(r, flag)
        if (flag /= LA_NO_ERROR .or. .not.fac%is_factored()) then
            rst = .false.
            print '(A)', "Test Failed: Cholesky Factorization Type Upper Form Status"
        else if (.not.is_mtx_equal(a, matmul(transpose(r), r), REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Cholesky Factorization Type Upper Form"
        end if

        fac_lower = cholesky_factorization(a, upper = .false., flag = flag)
        if (flag /= LA_NO_ERROR .or. fac_lower%is_upper()) then
            rst = .false.
            print '(A)', "Test Failed: Cholesky Factorization Type Lower Initialization"
        end if
        call fac_lower%factor(flag)
        if (flag /= LA_NO_ERROR .or. .not.fac_lower%is_factored()) then
            rst = .false.
            print '(A)', "Test Failed: Cholesky Factorization Type Lower Factor"
        else
            call fac_lower%get_factored_form(r, flag)
            if (flag /= LA_NO_ERROR) then
                rst = .false.
                print '(A)', "Test Failed: Cholesky Factorization Type Lower Form Status"
            else if (.not.is_mtx_equal(a, matmul(r, transpose(r)), REAL64_TOL)) then
                rst = .false.
                print '(A)', "Test Failed: Cholesky Factorization Type Lower Form"
            end if
        end if

        x_multiple = fac%solve(b_multiple, flag)
        if (flag /= LA_NO_ERROR) then
            rst = .false.
            print '(A)', "Test Failed: Cholesky Factorization Type Multiple RHS Solve Status"
        else if (.not.is_mtx_equal(x_multiple, x_expected, REAL64_TOL) .or. &
            .not.fac%is_factored()) then
            rst = .false.
            print '(A)', "Test Failed: Cholesky Factorization Type Multiple RHS Solve"
        end if

        fac_single = cholesky_factorization(a)
        x_single = fac_single%solve(b_single, flag)
        if (flag /= LA_NO_ERROR) then
            rst = .false.
            print '(A)', "Test Failed: Cholesky Factorization Type Single RHS Solve Status"
        else if (.not.is_mtx_equal(x_single, x_expected(:,1), REAL64_TOL) .or. &
            .not.fac_single%is_factored()) then
            rst = .false.
            print '(A)', "Test Failed: Cholesky Factorization Type Single RHS Solve"
        end if

        fac_multiple_in_place = cholesky_factorization(a)
        b_multiple_in_place = b_multiple
        call fac_multiple_in_place%solve_in_place(b_multiple_in_place, flag)
        if (flag /= LA_NO_ERROR .or. &
            .not.is_mtx_equal(b_multiple_in_place, x_expected, REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Cholesky Factorization Type Multiple RHS In-Place Solve"
        end if

        fac_single_in_place = cholesky_factorization(a)
        b_single_in_place = b_single
        call fac_single_in_place%solve_in_place(b_single_in_place, flag)
        if (flag /= LA_NO_ERROR .or. &
            .not.is_mtx_equal(b_single_in_place, x_expected(:,1), REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Cholesky Factorization Type Single RHS In-Place Solve"
        end if

        fac_invalid = cholesky_factorization(a_nonsquare, flag = flag)
        if (flag /= LA_MATRIX_NOT_SQUARE_ERROR .or. fac_invalid%size(1) /= 0 .or. &
            fac_invalid%size(2) /= 0) then
            rst = .false.
            print '(A)', "Test Failed: Cholesky Factorization Type Non-Square Validation"
        end if

        fac_not_positive_definite = cholesky_factorization(a_not_positive_definite)
        call fac_not_positive_definite%factor(flag)
        if (flag /= LA_POSITIVE_DEFINITE_ERROR .or. &
            fac_not_positive_definite%is_factored()) then
            rst = .false.
            print '(A)', "Test Failed: Cholesky Factorization Type Positive-Definite Validation"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_cholesky_rank1_update() result(rst)
        ! Parameters
        integer(int32), parameter :: n = 100

        ! Local Variables
        real(real64), dimension(n, n) :: a, r
        real(real64), dimension(n) :: u
        logical :: rst

        ! Initialization
        rst = .true.
        call create_random_array(a, mtype = POSITIVE_DEFINITE_MATRIX)
        call create_random_array(u)

        ! Compute the Cholesky factorization of the original matrix
        r = cholesky_factor(a)

        ! Update the original matrix
        call rank1_update(1.0d0, u, u, a)

        ! Update the factored matrix
        call cholesky_rank1_update(r, u)

        ! Test
        if (.not.assert(a, matmul(transpose(r), r), tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Cholesky Rank 1 Update Test 1"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_cholesky_rank1_downdate() result(rst)
        ! Parameters
        integer(int32), parameter :: n = 100

        ! Local Variables
        real(real64), dimension(n, n) :: a, r
        real(real64), dimension(n) :: u
        logical :: rst

        ! Initialization
        rst = .true.
        call create_random_array(a, mtype = POSITIVE_DEFINITE_MATRIX)
        call create_random_array(u)

        ! Start with a positive definite matrix, and then update it
        call rank1_update(1.0d0, u, u, a)

        ! Compute the Cholesky factorization of the original matrix
        r = cholesky_factor(a)

        ! Update the original matrix: A = A - u * u**T
        call rank1_update(-1.0d0, u, u, a)

        ! Update the factored matrix
        call cholesky_rank1_downdate(r, u)

        ! Test
        if (.not.is_mtx_equal(a, matmul(transpose(r), r), tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Cholesky Rank 1 Downdate Test 1"
        end if
    end function

! ******************************************************************************
! COMPLEX VALUED ROUTINES
! ------------------------------------------------------------------------------
    function test_cholesky_factor_cmplx() result(rst)
        ! Parameters
        integer(int32), parameter :: n = 10
        integer(int32), parameter :: nrhs = 20

        ! Local Variables
        complex(real64), dimension(n, n) :: a, u, l
        complex(real64), dimension(n, nrhs) :: b, x
        logical :: rst

        ! Initialization
        rst = .true.
        call create_random_array(a, mtype = POSITIVE_DEFINITE_MATRIX)
        call create_random_array(b)

        ! Test 1: A = L * L**T
        l = cholesky_factor(a, .false.)
        if (.not.assert(a, matmul(l, conjg(transpose(l))), tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Complex Cholesky Factorization Test 1"
        end if

        ! Test 2: A = U**T * U
        u = cholesky_factor(a, .true.)
        if (.not.assert(a, matmul(conjg(transpose(u)), u), tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Complex Cholesky Factorization Test 2"
        end if

        ! Test 3: Solve L*L**T * X = B
        x = solve_cholesky(.false., l, b)
        if (.not.assert(matmul(a, x), b, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Complex Cholesky Factorization Test 3"
        end if

        ! Test 4: Solve U**T * U * X = B
        x = solve_cholesky(.true., u, b)
        if (.not.assert(matmul(a, x), b, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Complex Cholesky Factorization Test 4"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_cholesky_rank1_update_cmplx() result(rst)
        ! Parameters
        integer(int32), parameter :: n = 100
        complex(real64), parameter :: one = (1.0d0, 0.0d0)

        ! Local Variables
        complex(real64), dimension(n, n) :: a, r
        complex(real64), dimension(n) :: u
        logical :: rst

        ! Initialization
        rst = .true.
        call create_random_array(a, mtype = POSITIVE_DEFINITE_MATRIX)
        call create_random_array(u)

        ! Compute the Cholesky factorization of the original matrix
        r = cholesky_factor(a)

        ! Update the original matrix
        call rank1_update(one, u, u, a)

        ! Update the factored matrix
        call cholesky_rank1_update(r, u)

        ! Test
        if (.not.assert(a, matmul(conjg(transpose(r)), r), tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Complex-Valued Cholesky Rank 1 Update Test 1"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_cholesky_rank1_downdate_cmplx() result(rst)
        ! Parameters
        integer(int32), parameter :: n = 100
        complex(real64), parameter :: one = (1.0d0, 0.0d0)

        ! Local Variables
        complex(real64), dimension(n, n) :: a, r
        complex(real64), dimension(n) :: u
        logical :: rst

        ! Initialization
        rst = .true.
        call create_random_array(a, mtype = POSITIVE_DEFINITE_MATRIX)
        call create_random_array(u)

        ! Start with a positive definite matrix, and then update it
        call rank1_update(one, u, u, a)

        ! Compute the Cholesky factorization of the original matrix
        r = cholesky_factor(a)

        ! Update the original matrix: A = A - u * u**T
        call rank1_update(-one, u, u, a)

        ! Update the factored matrix
        call cholesky_rank1_downdate(r, u)

        ! Test
        if (.not.is_mtx_equal(a, matmul(conjg(transpose(r)), r), tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Complex-Valued Cholesky Rank 1 Downdate Test 1"
        end if
    end function

! ------------------------------------------------------------------------------
end module
