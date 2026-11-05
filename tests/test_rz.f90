! test_rz.f90
! SPDX-License-Identifier: MIT
module test_rz
    use, intrinsic :: iso_fortran_env, only : real64
    use test_core
    use linalg
    use fortran_test_helper
    implicit none
contains
! ------------------------------------------------------------------------------
    function test_rz_factor() result(rst)
        ! Parameters
        integer, parameter :: m = 10
        integer, parameter :: n = 15

        ! Local Variables
        real(real64) :: a(m, n)
        real(real64), allocatable, dimension(:) :: tau
        real(real64), allocatable, dimension(:,:) :: r, rz, a1
        logical :: rst

        ! Initialization
        rst = .true.
        call random_number(a)
        a = extract_upper_triangular(a) ! RZ is for upper trapezoidal matrices

        ! Compute the factorization
        call rz_factor(a, tau, rz)

        ! Extract R & zero out any columns beyond M if N > M
        r = extract_upper_triangular(rz)
        if (n > m) then
            r(:,m+1:) = 0.0d0
        end if

        ! Compute R * Z
        a1 = mult_rz(.false., .false., n - m, rz, tau, r)

        ! Test that R * Z = A
        if (.not.assert(a1, a, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: RZ Test 1"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_rz_factor_cmplx() result(rst)
        ! Parameters
        integer, parameter :: m = 10
        integer, parameter :: n = 15

        ! Local Variables
        complex(real64) :: a(m, n)
        complex(real64), allocatable, dimension(:) :: tau
        complex(real64), allocatable, dimension(:,:) :: r, rz, a1
        logical :: rst

        ! Initialization
        rst = .true.
        call create_random_array(a)
        a = extract_upper_triangular(a) ! RZ is for upper trapezoidal matrices

        ! Compute the factorization
        call rz_factor(a, tau, rz)

        ! Extract R & zero out any columns beyond M if N > M
        r = extract_upper_triangular(rz)
        if (n > m) then
            r(:,m+1:) = (0.0d0, 0.0d0)
        end if

        ! Compute R * Z
        a1 = mult_rz(.false., .false., n - m, rz, tau, r)

        ! Test that R * Z = A
        if (.not.assert(a1, a, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Complex-Valued RZ Test 1"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_rz_square() result(rst)
        integer, parameter :: n = 3
        real(real64) :: a(n,n)
        complex(real64) :: ac(n,n)
        real(real64), allocatable :: tau(:), rz(:,:), r(:,:), reconstructed(:,:)
        complex(real64), allocatable :: ctau(:), crz(:,:), cr(:,:), &
            creconstructed(:,:)
        logical :: rst

        rst = .true.
        a = reshape([2.0_real64, 0.0_real64, 0.0_real64, &
            -1.0_real64, 3.0_real64, 0.0_real64, &
            4.0_real64, 5.0_real64, 6.0_real64], [n,n])
        call rz_factor(a, tau, rz)
        r = extract_upper_triangular(rz)
        reconstructed = mult_rz(.false., .false., 0, rz, tau, r)
        if (.not.assert(reconstructed, a, tol = REAL64_TOL)) rst = .false.

        ac = cmplx(a, -0.5_real64 * a, kind = real64)
        call rz_factor(ac, ctau, crz)
        cr = extract_upper_triangular(crz)
        creconstructed = mult_rz(.false., .false., 0, crz, ctau, cr)
        if (.not.assert(creconstructed, ac, tol = REAL64_TOL)) rst = .false.

        if (.not.rst) print '(A)', "Test Failed: Square RZ Factorization"
    end function

! ------------------------------------------------------------------------------
end module
