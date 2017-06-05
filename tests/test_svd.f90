! test_svd.f90

! Tests for SVD operations.
module test_svd_ops
    use linalg_constants
    use test_core
    use linalg_factor, only : svd
    use linalg_core, only : diag_mtx_mult
    implicit none
contains
! ------------------------------------------------------------------------------
    subroutine test_svd()
        ! Parameters
        integer(i32), parameter :: m = 60
        integer(i32), parameter :: n = 60
        real(dp), parameter :: tol = 1.0d-8

        ! Local Variables
        real(dp), dimension(m, n) :: a, a1
        real(dp), dimension(m, m) :: u1
        real(dp), dimension(n, n) :: vt1
        real(dp), dimension(n) :: s1
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
        if (.not.is_mtx_equal(a, a1, tol)) then
            rst = .false.
            print '(A)', "Test Failed: Singular Value Decomposition"
        end if
        if (rst) print '(A)', "Test Passed: Singular Value Decomposition"
    end subroutine

! ------------------------------------------------------------------------------
    subroutine test_svd_od()
        ! Parameters
        integer(i32), parameter :: m = 60
        integer(i32), parameter :: n = 50
        real(dp), parameter :: tol = 1.0d-8

        ! Local Variables
        real(dp), dimension(m, n) :: a, a1, t1, u2
        real(dp), dimension(m, m) :: u1
        real(dp), dimension(n, n) :: vt1
        real(dp), dimension(n) :: s1
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
        if (.not.is_mtx_equal(a, a1, tol)) then
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
        if (.not.is_mtx_equal(a, a1, tol)) then
            rst = .false.
            print '(A)', "Test Failed: Overdetermined Singular Value Decomposition, Test 2"
        end if
        if (rst) print '(A)', "Test Passed: Overdetermined Singular Value Decomposition"
    end subroutine

! ------------------------------------------------------------------------------
    subroutine test_svd_ud()
        ! Parameters
        integer(i32), parameter :: m = 50
        integer(i32), parameter :: n = 60
        real(dp), parameter :: tol = 1.0d-8

        ! Local Variables
        real(dp), dimension(m, n) :: a, a1, t1
        real(dp), dimension(m, m) :: u1
        real(dp), dimension(n, n) :: vt1
        real(dp), dimension(m) :: s1
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
        if (.not.is_mtx_equal(a, a1, tol)) then
            rst = .false.
            print '(A)', "Test Failed: Underdetermined Singular Value Decomposition, Test 1"
        end if
        if (rst) print '(A)', "Test Passed: Underdetermined Singular Value Decomposition"
    end subroutine


end module
