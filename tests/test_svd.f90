! test_svd.f90

! Tests for SVD operations.
module test_svd_ops
    use, intrinsic :: iso_fortran_env, only : int32, real64
    use test_core
    use linalg
    implicit none
contains
! ------------------------------------------------------------------------------
    function test_svd() result(rst)
        ! Parameters
        integer(int32), parameter :: m = 60
        integer(int32), parameter :: n = 60

        ! Local Variables
        real(real64), dimension(m, n) :: a, a1
        real(real64), dimension(m, m) :: u1
        real(real64), dimension(n, n) :: vt1
        real(real64), dimension(n) :: s1
        logical :: rst

        ! Initialization
        rst = .true.
        call random_number(a)
        a1 = a

        ! Compute the full SVD of A
        call svd(a1, s1, u1, vt1)

        ! Ensure A = U * S * V**T
        call diag_mtx_mult(.false., 1.0d0, s1, u1)  ! U1 = U1 * S1
        a1 = matmul(u1, vt1)

        ! Test
        if (.not.is_mtx_equal(a, a1, REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Singular Value Decomposition"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_svd_od() result(rst)
        ! Parameters
        integer(int32), parameter :: m = 60
        integer(int32), parameter :: n = 50

        ! Local Variables
        real(real64), dimension(m, n) :: a, a1, t1, u2
        real(real64), dimension(m, m) :: u1
        real(real64), dimension(n, n) :: vt1
        real(real64), dimension(n) :: s1
        logical :: rst

        ! Initialization
        rst = .true.
        call random_number(a)
        a1 = a

        ! Compute the full SVD of A
        call svd(a1, s1, u1, vt1)

        ! Ensure A = U * S * V**T
        t1 = u1(:,1:n)
        call diag_mtx_mult(.false., 1.0d0, s1, t1)  ! T1 = U1 * S1
        a1 = matmul(t1, vt1)

        ! Test
        if (.not.is_mtx_equal(a, a1, REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Overdetermined Singular Value Decomposition, Test 1"
        end if

        ! Compute the partial SVD - Incomplete U
        a1 = a
        call svd(a1, s1, u2, vt1)

        ! Ensure A = U * S * V**T
        call diag_mtx_mult(.false., 1.0d0, s1, u2) ! U2 = U2 * S1
        a1 = matmul(u2, vt1)

        ! Test
        if (.not.is_mtx_equal(a, a1, REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Overdetermined Singular Value Decomposition, Test 2"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_svd_ud() result(rst)
        ! Parameters
        integer(int32), parameter :: m = 50
        integer(int32), parameter :: n = 60

        ! Local Variables
        real(real64), dimension(m, n) :: a, a1, t1
        real(real64), dimension(m, m) :: u1
        real(real64), dimension(n, n) :: vt1
        real(real64), dimension(m) :: s1
        logical :: rst

        ! Initialization
        rst = .true.
        call random_number(a)
        a1 = a
        t1 = 0.0d0

        ! Compute the full SVD of A
        call svd(a1, s1, u1, vt1)

        ! Ensure A = U * S * V**T
        t1(:,1:m) = u1
        call diag_mtx_mult(.false., 1.0d0, s1, t1)  ! T1 = U1 * S1
        a1 = matmul(t1, vt1)

        ! Test
        if (.not.is_mtx_equal(a, a1, REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Underdetermined Singular Value Decomposition, Test 1"
        end if
    end function


end module
