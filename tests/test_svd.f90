! test_svd.f90
! SPDX-License-Identifier: MIT

! Tests for SVD operations.
module test_svd_ops
    use, intrinsic :: iso_fortran_env, only : real64
    use test_core
    use linalg
    use linalg_errors, only : LA_NO_ERROR, LA_NOT_INITIALIZED_ERROR, &
        LA_MATRIX_DIMENSION_ERROR
    use fortran_test_helper
    implicit none
contains
! ------------------------------------------------------------------------------
    function test_svd() result(rst)
        ! Parameters
        integer, parameter :: m = 60
        integer, parameter :: n = 60

        ! Local Variables
        real(real64), dimension(m, n) :: a, a1
        real(real64), allocatable, dimension(:,:) :: u1, vt1
        real(real64), allocatable, dimension(:) :: s1
        logical :: rst

        ! Initialization
        rst = .true.
        call create_random_array(a)

        ! Compute the full SVD of A
        call svd(a, s = s1, u = u1, vt = vt1)

        ! Ensure A = U * S * V**T
        call diag_mtx_mult(.false., 1.0d0, s1, u1)  ! U1 = U1 * S1
        a1 = matmul(u1, vt1)

        ! Test
        if (.not.assert(a, a1, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Singular Value Decomposition"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_svd_cmplx() result(rst)
        ! Parameters
        integer, parameter :: m = 60
        integer, parameter :: n = 60
        complex(real64), parameter :: one = (1.0d0, 0.0d0)

        ! Local Variables
        complex(real64), dimension(m, n) :: a, a1
        complex(real64), allocatable, dimension(:,:) :: u1, vt1
        real(real64), allocatable, dimension(:) :: s1
        logical :: rst

        ! Initialization
        rst = .true.
        call create_random_array(a)

        ! Compute the full SVD of A
        call svd(a, s = s1, u = u1, vt = vt1)

        ! Ensure A = U * S * V**T
        call diag_mtx_mult(.false., one, s1, u1)  ! U1 = U1 * S1
        a1 = matmul(u1, vt1)

        ! Test
        if (.not.assert(a, a1, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Complex-Valued Singular Value Decomposition"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_svd_od() result(rst)
        ! Parameters
        integer, parameter :: m = 60
        integer, parameter :: n = 50

        ! Local Variables
        real(real64), dimension(m, n) :: a, a1, t1
        real(real64), allocatable, dimension(:,:) :: u1, vt1
        real(real64), allocatable, dimension(:) :: s1
        logical :: rst

        ! Initialization
        rst = .true.
        call create_random_array(a)

        ! Compute the full SVD of A
        call svd(a, s = s1, u = u1, vt = vt1)

        ! Ensure A = U * S * V**T
        t1 = u1(:,1:n)
        call diag_mtx_mult(.false., 1.0d0, s1, t1)  ! T1 = U1 * S1
        a1 = matmul(t1, vt1)

        ! Test
        if (.not.assert(a, a1, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Overdetermined Singular Value Decomposition, Test 1"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_svd_od_cmplx() result(rst)
        ! Parameters
        integer, parameter :: m = 60
        integer, parameter :: n = 50
        complex(real64), parameter :: one = (1.0d0, 0.0d0)

        ! Local Variables
        complex(real64), dimension(m, n) :: a, a1, t1
        complex(real64), allocatable, dimension(:,:) :: u1, vt1
        real(real64), allocatable, dimension(:) :: s1
        logical :: rst

        ! Initialization
        rst = .true.
        call create_random_array(a)

        ! Compute the full SVD of A
        call svd(a, s = s1, u = u1, vt = vt1)

        ! Ensure A = U * S * V**T
        t1 = u1(:,1:n)
        call diag_mtx_mult(.false., one, s1, t1)  ! T1 = U1 * S1
        a1 = matmul(t1, vt1)

        ! Test
        if (.not.assert(a, a1, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Complex-Valued Overdetermined Singular Value Decomposition, Test 1"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_svd_ud() result(rst)
        ! Parameters
        integer, parameter :: m = 50
        integer, parameter :: n = 60

        ! Local Variables
        real(real64), dimension(m, n) :: a, a1, t1
        real(real64), allocatable :: u1(:,:), vt1(:,:), s1(:)
        logical :: rst

        ! Initialization
        rst = .true.
        call create_random_array(a)
        t1 = 0.0d0

        ! Compute the full SVD of A
        call svd(a, s = s1, u = u1, vt = vt1)

        ! Ensure A = U * S * V**T
        t1(:,1:m) = u1
        call diag_mtx_mult(.false., 1.0d0, s1, t1)  ! T1 = U1 * S1
        a1 = matmul(t1, vt1)

        ! Test
        if (.not.assert(a, a1, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Underdetermined Singular Value Decomposition, Test 1"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_svd_ud_cmplx() result(rst)
        ! Parameters
        integer, parameter :: m = 50
        integer, parameter :: n = 60
        complex(real64), parameter :: zero = (0.0d0, 0.0d0)
        complex(real64), parameter :: one = (1.0d0, 0.0d0)

        ! Local Variables
        complex(real64), dimension(m, n) :: a, t1, a1
        complex(real64), allocatable, dimension(:,:) :: u1, vt1
        real(real64), allocatable, dimension(:) :: s1
        logical :: rst

        ! Initialization
        rst = .true.
        call create_random_array(a)
        t1 = zero

        ! Compute the full SVD of A
        call svd(a, s = s1, u = u1, vt = vt1)

        ! Ensure A = U * S * V**T
        t1(:,1:m) = u1
        call diag_mtx_mult(.false., one, s1, t1)  ! T1 = U1 * S1
        a1 = matmul(t1, vt1)

        ! Test
        if (.not.assert(a, a1, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Complex-Valued Underdetermined Singular Value Decomposition, Test 1"
        end if
    end function

! ******************************************************************************
! SVD FACTORIZATION TYPE TESTS
! ------------------------------------------------------------------------------
    function test_svd_factorization_type() result(rst)
        ! Parameters
        integer, parameter :: n = 6
        integer, parameter :: nrhs = 3

        ! Local Variables
        real(real64), dimension(n, n) :: a, t
        real(real64), dimension(n, nrhs) :: x_expected, b_multiple, b_copy
        real(real64), dimension(n) :: b_single, b_single_copy
        real(real64), allocatable, dimension(:,:) :: u, vt, x_multiple
        real(real64), allocatable, dimension(:) :: s, x_single
        type(svd_factorization) :: fac, fac_form, fac_single, fac_multiple
        type(svd_factorization) :: fac_inplace
        integer :: i, flag
        logical :: rst

        ! Initialization
        rst = .true.
        call create_random_array(a)
        do i = 1, n
            a(i,i) = a(i,i) + real(n, real64)
        end do
        call create_random_array(x_expected)
        b_multiple = matmul(a, x_expected)
        b_single = matmul(a, x_expected(:,1))

        ! Initialization status
        fac = svd_factorization(a)
        if (fac%size(1) /= n .or. fac%size(2) /= n .or. fac%is_factored()) then
            rst = .false.
            print '(A)', "Test Failed: SVD Factorization Type Initialization"
        end if

        ! Factored form: A = U * S * V**T
        fac_form = svd_factorization(a)
        call fac_form%get_factored_form(u, s, vt, flag)
        if (flag /= LA_NO_ERROR .or. .not.fac_form%is_factored()) then
            rst = .false.
            print '(A)', "Test Failed: SVD Factorization Type Factored Form Status"
        else
            t = u
            do i = 1, n
                t(:,i) = t(:,i) * s(i)
            end do
            if (.not.assert(a, matmul(t, vt), tol = REAL64_TOL)) then
                rst = .false.
                print '(A)', "Test Failed: SVD Factorization Type Factored Form"
            end if
            if (.not.all(s(1:n-1) >= s(2:n))) then
                rst = .false.
                print '(A)', "Test Failed: SVD Factorization Type Singular Value Order"
            end if
        end if

        ! Single right-hand-side solve (allocating)
        fac_single = svd_factorization(a)
        x_single = fac_single%solve(b_single, flag)
        if (flag /= LA_NO_ERROR .or. .not.fac_single%is_factored()) then
            rst = .false.
            print '(A)', "Test Failed: SVD Factorization Type Single RHS Status"
        else if (.not.assert(x_single, x_expected(:,1), tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: SVD Factorization Type Single RHS"
        end if

        ! Multiple right-hand-side solve (allocating)
        fac_multiple = svd_factorization(a)
        x_multiple = fac_multiple%solve(b_multiple, flag)
        if (flag /= LA_NO_ERROR .or. .not.fac_multiple%is_factored()) then
            rst = .false.
            print '(A)', "Test Failed: SVD Factorization Type Multiple RHS Status"
        else if (.not.assert(x_multiple, x_expected, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: SVD Factorization Type Multiple RHS"
        end if

        ! Single right-hand-side solve (in-place, already factored)
        b_single_copy = b_single
        call fac_single%solve_in_place(b_single_copy, flag)
        if (flag /= LA_NO_ERROR) then
            rst = .false.
            print '(A)', "Test Failed: SVD Factorization Type In-Place Single RHS Status"
        else if (.not.assert(b_single_copy, x_expected(:,1), tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: SVD Factorization Type In-Place Single RHS"
        end if

        ! Multiple right-hand-side solve (in-place, automatically factored)
        fac_inplace = svd_factorization(a)
        b_copy = b_multiple
        call fac_inplace%solve_in_place(b_copy, flag)
        if (flag /= LA_NO_ERROR .or. .not.fac_inplace%is_factored()) then
            rst = .false.
            print '(A)', "Test Failed: SVD Factorization Type In-Place Multiple RHS Status"
        else if (.not.assert(b_copy, x_expected, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: SVD Factorization Type In-Place Multiple RHS"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_svd_factorization_skinny() result(rst)
        ! Parameters
        integer, parameter :: m = 8
        integer, parameter :: n = 5

        ! Local Variables
        real(real64), dimension(m, n) :: a
        real(real64), dimension(m) :: b, b_inplace
        real(real64), dimension(n) :: zero_n
        real(real64), dimension(n, m) :: wide
        real(real64), dimension(n) :: b_wide
        real(real64), allocatable, dimension(:,:) :: u, vt, t
        real(real64), allocatable, dimension(:) :: s, x, x_wide, x_ref
        type(svd_factorization) :: fac_skinny, fac_full, fac_wide
        integer :: i, flag
        logical :: rst

        ! Initialization
        rst = .true.
        call create_random_array(a)
        zero_n = 0.0d0

        ! Skinny U and V: A = U * S * V**T with U(M,N) and V**T(N,N)
        fac_skinny = svd_factorization(a, skinny_u = .true., skinny_v = .true.)
        call fac_skinny%get_factored_form(u, s, vt, flag)
        if (flag /= LA_NO_ERROR .or. any(shape(u) /= [m, n]) .or. &
            size(s) /= n .or. any(shape(vt) /= [n, n])) then
            rst = .false.
            print '(A)', "Test Failed: SVD Factorization Type Skinny Shapes"
        else
            t = u
            do i = 1, n
                t(:,i) = t(:,i) * s(i)
            end do
            if (.not.assert(a, matmul(t, vt), tol = REAL64_TOL)) then
                rst = .false.
                print '(A)', "Test Failed: SVD Factorization Type Skinny Reconstruction"
            end if
        end if

        ! Full U and V: A = U(:,1:N) * S * V**T
        fac_full = svd_factorization(a)
        call fac_full%get_factored_form(u, s, vt, flag)
        if (flag /= LA_NO_ERROR .or. any(shape(u) /= [m, m]) .or. &
            any(shape(vt) /= [n, n])) then
            rst = .false.
            print '(A)', "Test Failed: SVD Factorization Type Full Shapes"
        else
            t = u(:,1:n)
            do i = 1, n
                t(:,i) = t(:,i) * s(i)
            end do
            if (.not.assert(a, matmul(t, vt), tol = REAL64_TOL)) then
                rst = .false.
                print '(A)', "Test Failed: SVD Factorization Type Full Reconstruction"
            end if
        end if

        ! Overdetermined least-squares solution: A**T * (A * X - B) = 0
        call create_random_array(b)
        x = fac_full%solve(b, flag)
        if (flag /= LA_NO_ERROR .or. size(x) /= n) then
            rst = .false.
            print '(A)', "Test Failed: SVD Factorization Type Skinny Solve Status"
        else if (.not.assert(matmul(transpose(a), matmul(a, x) - b), &
            zero_n, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: SVD Factorization Type Skinny Solve"
        end if

        ! In-place solve: B has MAX(M, N) rows; the solution is in the first N
        b_inplace = 0.0d0
        b_inplace(1:m) = b
        call fac_full%solve_in_place(b_inplace, flag)
        if (flag /= LA_NO_ERROR) then
            rst = .false.
            print '(A)', "Test Failed: SVD Factorization Type Skinny In-Place Status"
        else if (.not.assert(b_inplace(1:n), x, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: SVD Factorization Type Skinny In-Place"
        end if

        ! Underdetermined (wide) system: minimum-norm solution X = V S^+ U**T B
        wide = transpose(a)
        b_wide = 1.0d0
        fac_wide = svd_factorization(wide, skinny_u = .true., skinny_v = .true.)
        x_wide = fac_wide%solve(b_wide, flag)
        if (flag /= LA_NO_ERROR .or. size(x_wide) /= m) then
            rst = .false.
            print '(A)', "Test Failed: SVD Factorization Type Wide Solve Status"
        else if (.not.assert(matmul(wide, x_wide), b_wide, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: SVD Factorization Type Wide Solve"
        else
            call fac_wide%get_factored_form(u, s, vt, flag)
            x_ref = matmul(transpose(vt), matmul(transpose(u), b_wide) / s)
            if (.not.assert(x_wide, x_ref, tol = REAL64_TOL)) then
                rst = .false.
                print '(A)', "Test Failed: SVD Factorization Type Wide Minimum Norm"
            end if
        end if
    end function

! ------------------------------------------------------------------------------
    function test_svd_factorization_errors() result(rst)
        ! Parameters
        integer, parameter :: n = 4

        ! Local Variables
        real(real64), dimension(n, n) :: a, a_singular
        real(real64), dimension(n) :: b, zero_n
        real(real64), dimension(n - 1) :: b_bad
        real(real64), allocatable, dimension(:) :: x
        type(svd_factorization) :: empty, fac, fac_singular
        integer :: flag, i
        logical :: rst

        ! Initialization
        rst = .true.
        call create_random_array(a)
        do i = 1, n
            a(i,i) = a(i,i) + real(n, real64)
        end do
        call create_random_array(b)
        b_bad = 1.0d0
        zero_n = 0.0d0

        ! Uninitialized object
        call empty%factor(flag)
        if (flag /= LA_NOT_INITIALIZED_ERROR) then
            rst = .false.
            print '(A)', "Test Failed: SVD Factorization Type Uninitialized Factor"
        end if
        x = empty%solve(b, flag)
        if (flag /= LA_NOT_INITIALIZED_ERROR) then
            rst = .false.
            print '(A)', "Test Failed: SVD Factorization Type Uninitialized Solve"
        end if

        ! Dimension mismatch
        fac = svd_factorization(a)
        x = fac%solve(b_bad, flag)
        if (flag /= LA_MATRIX_DIMENSION_ERROR) then
            rst = .false.
            print '(A)', "Test Failed: SVD Factorization Type Dimension Mismatch"
        end if

        ! Rank-deficient matrix (zero column): minimum-norm solution with a
        ! zero component in the null direction
        a_singular = a
        a_singular(:,2) = 0.0d0
        fac_singular = svd_factorization(a_singular)
        x = fac_singular%solve(b, flag)
        if (flag /= LA_NO_ERROR .or. abs(x(2)) > REAL64_TOL) then
            rst = .false.
            print '(A)', "Test Failed: SVD Factorization Type Rank-Deficient Solve"
        else if (.not.assert(matmul(transpose(a_singular), &
            matmul(a_singular, x) - b), zero_n, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: SVD Factorization Type Rank-Deficient Least Squares"
        end if

        ! Updating the object destroys the existing factorization
        call fac%factor(flag)
        if (flag /= LA_NO_ERROR .or. .not.fac%is_factored()) then
            rst = .false.
            print '(A)', "Test Failed: SVD Factorization Type Update Pre-Factor"
        end if
        call fac%update(a)
        if (fac%is_factored()) then
            rst = .false.
            print '(A)', "Test Failed: SVD Factorization Type Update"
        end if
    end function

! ------------------------------------------------------------------------------
end module
