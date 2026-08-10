module test_rz
    use, intrinsic :: iso_fortran_env, only : int32, real64
    use test_core
    use linalg
    use fortran_test_helper
    implicit none
contains
! ------------------------------------------------------------------------------
    function test_rz_factor() result(rst)
        ! Parameters
        integer(int32), parameter :: m = 10
        integer(int32), parameter :: n = 15

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
        integer(int32), parameter :: m = 10
        integer(int32), parameter :: n = 15

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
end module