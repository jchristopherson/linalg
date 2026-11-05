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

! ------------------------------------------------------------------------------
    function test_svd_optional_outputs() result(rst)
        integer, parameter :: m = 3, n = 2
        real(real64), parameter :: zero = 0.0_real64
        real(real64) :: a(m,n), expected(n)
        complex(real64) :: ac(m,n)
        real(real64), allocatable :: s(:)
        real(real64), allocatable :: u(:,:), vt(:,:)
        complex(real64), allocatable :: uc(:,:), vtc(:,:)
        logical :: rst

        rst = .true.
        a = reshape([3.0_real64, zero, zero, zero, 2.0_real64, zero], [m,n])
        ac = cmplx(a, 0.0_real64, kind = real64)
        ac(1,1) = cmplx(0.0_real64, 3.0_real64, kind = real64)
        expected = [3.0_real64, 2.0_real64]

        ! Exercise each optional-output combination of the real interface.
        call svd(a)
        call svd(a, s = s)
        if (.not.assert(s, expected, tol = REAL64_TOL)) rst = .false.
        call svd(a, u = u)
        if (any(shape(u) /= [m,m])) rst = .false.
        call svd(a, vt = vt)
        if (any(shape(vt) /= [n,n])) rst = .false.
        call svd(a, s = s, u = u)
        if (.not.assert(s, expected, tol = REAL64_TOL) .or. &
            any(shape(u) /= [m,m])) rst = .false.
        call svd(a, s = s, vt = vt)
        if (.not.assert(s, expected, tol = REAL64_TOL) .or. &
            any(shape(vt) /= [n,n])) rst = .false.
        call svd(a, u = u, vt = vt)
        if (any(shape(u) /= [m,m]) .or. any(shape(vt) /= [n,n])) rst = .false.

        ! Exercise the same optional-output branches for the complex interface.
        call svd(ac)
        call svd(ac, s = s)
        if (.not.assert(s, expected, tol = REAL64_TOL)) rst = .false.
        call svd(ac, u = uc)
        if (any(shape(uc) /= [m,m])) rst = .false.
        call svd(ac, vt = vtc)
        if (any(shape(vtc) /= [n,n])) rst = .false.
        call svd(ac, s = s, u = uc)
        if (.not.assert(s, expected, tol = REAL64_TOL) .or. &
            any(shape(uc) /= [m,m])) rst = .false.
        call svd(ac, s = s, vt = vtc)
        if (.not.assert(s, expected, tol = REAL64_TOL) .or. &
            any(shape(vtc) /= [n,n])) rst = .false.
        call svd(ac, u = uc, vt = vtc)
        if (any(shape(uc) /= [m,m]) .or. any(shape(vtc) /= [n,n])) rst = .false.

        if (.not.rst) print '(A)', "Test Failed: SVD Optional Outputs"
    end function

! ------------------------------------------------------------------------------
    function test_svd_complex_invariants() result(rst)
        integer, parameter :: m = 3, n = 2
        complex(real64), parameter :: zero = (0.0_real64, 0.0_real64)
        complex(real64) :: a(m,n), reconstructed(m,n)
        complex(real64), allocatable :: u(:,:), vt(:,:), scaled(:,:)
        real(real64), allocatable :: s(:)
        real(real64) :: expected(n)
        logical :: rst

        rst = .true.
        a = zero
        a(1,1) = (0.0_real64, 3.0_real64)
        a(2,2) = (2.0_real64, 0.0_real64)
        expected = [3.0_real64, 2.0_real64]
        call svd(a, s = s, u = u, vt = vt)

        if (.not.assert(s, expected, tol = REAL64_TOL)) rst = .false.
        if (.not.all(s(1:n-1) >= s(2:n))) rst = .false.
        scaled = u(:,1:n)
        scaled(:,1) = scaled(:,1) * s(1)
        scaled(:,2) = scaled(:,2) * s(2)
        reconstructed = matmul(scaled, vt)
        if (.not.assert(a, reconstructed, tol = REAL64_TOL)) rst = .false.
        if (.not.check_complex_orthogonality(u, vt)) rst = .false.

        if (.not.rst) print '(A)', "Test Failed: Complex SVD Invariants"
    end function

! ------------------------------------------------------------------------------
    function test_svd_edge_cases() result(rst)
        integer, parameter :: n = 3
        real(real64) :: a(n,n), b(n), expected(n)
        real(real64), allocatable :: s(:), u(:,:), vt(:,:), x(:)
        type(svd_factorization) :: fac
        integer :: flag
        logical :: rst

        rst = .true.
        a = 0.0_real64
        a(1,1) = 4.0_real64
        a(2,2) = 2.0_real64
        call svd(a, s = s, u = u, vt = vt)
        expected = [4.0_real64, 2.0_real64, 0.0_real64]
        if (.not.assert(s, expected, tol = REAL64_TOL)) rst = .false.
        if (.not.check_real_orthogonality(u, vt)) rst = .false.

        b = [1.0_real64, 2.0_real64, 7.0_real64]
        fac = svd_factorization(a)
        x = fac%solve(b, flag)
        expected = [0.25_real64, 1.0_real64, 0.0_real64]
        if (flag /= LA_NO_ERROR .or. &
            .not.assert(x, expected, tol = REAL64_TOL)) rst = .false.

        ! Values below the documented relative rank threshold are discarded.
        a = 0.0_real64
        a(1,1) = 1.0_real64
        a(2,2) = epsilon(1.0_real64)
        fac = svd_factorization(a)
        b = 1.0_real64
        x = fac%solve(b, flag)
        expected = [1.0_real64, 0.0_real64, 0.0_real64]
        if (flag /= LA_NO_ERROR .or. &
            .not.assert(x, expected, tol = REAL64_TOL)) rst = .false.

        ! A nonempty all-zero matrix has a zero minimum-norm solution.
        a = 0.0_real64
        fac = svd_factorization(a)
        x = fac%solve(b, flag)
        if (flag /= LA_NO_ERROR .or. any(x /= 0.0_real64)) rst = .false.

        if (.not.rst) print '(A)', "Test Failed: SVD Edge Cases"
    end function

! ------------------------------------------------------------------------------
    function test_svd_skinny_flag_combinations() result(rst)
        real(real64) :: tall(4,2), wide(2,4)
        logical :: rst, skinny_u, skinny_v
        integer :: i

        rst = .true.
        tall = reshape([1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64, &
            2.0_real64, -1.0_real64, 0.5_real64, 3.0_real64], [4,2])
        wide = transpose(tall)
        do i = 0, 3
            skinny_u = mod(i, 2) == 1
            skinny_v = i >= 2
            call check_svd_factorization_shape(tall, skinny_u, skinny_v, rst)
            call check_svd_factorization_shape(wide, skinny_u, skinny_v, rst)
        end do
        if (.not.rst) print '(A)', "Test Failed: SVD Skinny Flag Combinations"
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
    subroutine check_svd_factorization_shape(a, skinny_u, skinny_v, passed)
        real(real64), intent(in) :: a(:,:)
        logical, intent(in) :: skinny_u, skinny_v
        logical, intent(inout) :: passed
        type(svd_factorization) :: fac
        real(real64), allocatable :: u(:,:), vt(:,:), s(:), scaled(:,:)
        integer :: m, n, mn, flag, i

        m = size(a, 1)
        n = size(a, 2)
        mn = min(m, n)
        fac = svd_factorization(a, skinny_u = skinny_u, skinny_v = skinny_v)
        call fac%get_factored_form(u, s, vt, flag)
        if (flag /= LA_NO_ERROR) then
            passed = .false.
            return
        end if
        if (.not.allocated(u) .or. .not.allocated(s) .or. .not.allocated(vt)) then
            passed = .false.
            return
        end if
        if (any(shape(u) /= [m, merge(mn, m, skinny_u)]) .or. &
            any(shape(vt) /= [merge(mn, n, skinny_v), n]) .or. size(s) /= mn) then
            passed = .false.
            return
        end if

        scaled = u(:,1:mn)
        do i = 1, mn
            scaled(:,i) = scaled(:,i) * s(i)
        end do
        if (.not.assert(a, matmul(scaled, vt(1:mn,:)), tol = REAL64_TOL)) passed = .false.
        if (.not.all(s(1:mn-1) >= s(2:mn))) passed = .false.
        if (.not.check_real_orthogonality(u(:,1:mn), vt(1:mn,:))) passed = .false.
    end subroutine

    ! ------------------------------------------------------------------------------
    function check_real_orthogonality(u, vt) result(passed)
        real(real64), intent(in) :: u(:,:), vt(:,:)
        real(real64) :: gram(size(u,2),size(u,2)), identity(size(u,2),size(u,2))
        logical :: passed
        integer :: i

        identity = 0.0_real64
        do i = 1, size(identity,1)
            identity(i,i) = 1.0_real64
        end do
        gram = matmul(transpose(u), u)
        passed = assert(gram, identity, tol = REAL64_TOL)
        gram = matmul(vt, transpose(vt))
        passed = passed .and. assert(gram, identity, tol = REAL64_TOL)
    end function

    ! ------------------------------------------------------------------------------
    function check_complex_orthogonality(u, vt) result(passed)
        complex(real64), intent(in) :: u(:,:), vt(:,:)
        complex(real64) :: gram_u(size(u,2),size(u,2))
        complex(real64) :: gram_v(size(vt,1),size(vt,1))
        complex(real64) :: identity_u(size(u,2),size(u,2))
        complex(real64) :: identity_v(size(vt,1),size(vt,1))
        logical :: passed
        integer :: i

        identity_u = (0.0_real64, 0.0_real64)
        do i = 1, size(identity_u,1)
            identity_u(i,i) = (1.0_real64, 0.0_real64)
        end do
        identity_v = (0.0_real64, 0.0_real64)
        do i = 1, size(identity_v,1)
            identity_v(i,i) = (1.0_real64, 0.0_real64)
        end do
        gram_u = matmul(transpose(conjg(u)), u)
        passed = assert(gram_u, identity_u, tol = REAL64_TOL)
        gram_v = matmul(vt, transpose(conjg(vt)))
        passed = passed .and. assert(gram_v, identity_v, tol = REAL64_TOL)
    end function

    ! ------------------------------------------------------------------------------
    end module
