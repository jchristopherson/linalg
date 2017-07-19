! test_misc.f90

! Tests miscellaneous routines.
module test_misc
    use linalg_constants
    use linalg_core
    use test_core
    implicit none
contains
! ******************************************************************************
! DIAGONAL MATRIX MULTIPLICATION TEST
! ------------------------------------------------------------------------------
    subroutine test_diagonal_mtx_mult()
        ! Parameters
        integer(i32), parameter :: m = 30
        integer(i32), parameter :: n = 30
        integer(i32), parameter :: k = 30
        real(dp), parameter :: tol = 1.0d-8
        real(dp), parameter :: alpha = 0.5d0
        real(dp), parameter :: beta = 0.25d0

        ! Local Variables
        integer(i32) :: i
        real(dp), dimension(m, n) :: c1, ans1
        real(dp), dimension(k, n) :: b1, b2, ans2
        real(dp), dimension(m, k) :: d1
        real(dp), dimension(k) :: d1v
        logical :: rst

        ! Initialization
        rst = .true.
        call random_number(c1)
        call random_number(b1)
        call random_number(d1v)
        d1 = 0.0d0
        do i = 1, k
            d1(i,i) = d1v(i)
        end do

        ! Compute C1 = D1 * B1 + C1
        ans1 = alpha * matmul(d1, b1) + beta * c1
        call diag_mtx_mult(.true., .false., alpha, d1v, b1, beta, c1)

        ! Test
        if (.not.is_mtx_equal(ans1, c1, tol)) then
            rst = .false.
            print '(A)', "Test Failed: Diagonal Matrix Multiply Test 1"
        end if

        ! Compute C1 = D1 * B1**T + C1
        ans1 = alpha * matmul(d1, transpose(b1)) + beta * c1
        call diag_mtx_mult(.true., .true., alpha, d1v, b1, beta, c1)

        ! Test
        if (.not.is_mtx_equal(ans1, c1, tol)) then
            rst = .false.
            print '(A)', "Test Failed: Diagonal Matrix Multiply Test 2"
        end if

        ! Compute B2 = D1 * B2
        ans2 = alpha * matmul(d1, b2)
        call diag_mtx_mult(.true., alpha, d1v, b2)

        ! Test
        if (.not.is_mtx_equal(ans2, b2, tol)) then
            rst = .false.
            print '(A)', "Test Failed: Diagonal Matrix Multiply Test 3"
        end if
        if (rst) print '(A)', "Test Passed: Diagonal Matrix Multiply"
    end subroutine

! ******************************************************************************
! RANK 1 UPDATE TEST
! ------------------------------------------------------------------------------
    subroutine test_rank1_update()
        ! Parameters
        integer(i32), parameter :: m = 50
        integer(i32), parameter :: n = 20
        real(dp), parameter :: alpha = 0.5d0
        real(dp), parameter :: tol = 1.0d-12

        ! Local Variables
        real(dp), dimension(m, n) :: a, a1, b
        real(dp), dimension(m, 1) :: x
        real(dp), dimension(n, 1) :: y

        ! Initialization
        call random_number(a)
        call random_number(x)
        call random_number(y)
        a1 = a

        ! Define the solution
        b = alpha * matmul(x, transpose(y)) + a

        ! Perform the operation
        call rank1_update(alpha, x(:,1), y(:,1), a1)

        ! Compare the results
        if (.not.is_mtx_equal(a1, b, tol)) then
            print '(A)', "Test Failed: Rank 1 Update"
        else
            print '(A)', "Test Passed: Rank 1 Update"
        end if
    end subroutine

! ******************************************************************************
! MATRIX RANK TESTS
! ------------------------------------------------------------------------------
    ! REF: http://www.mathworks.com/help/matlab/ref/pinv.html?s_tid=srchtitle
    subroutine test_rank()
        ! Local Variables
        real(dp), dimension(8, 6) :: a

        ! Initialization
        a = reshape([64, 9, 17, 40, 32, 41, 49, 8, 2, 55, 47, 26, 34, 23, 15, &
            58, 3, 54, 46, 27, 35, 22, 14, 59, 61, 12, 20, 37, 29, 44, 52, 5, &
            60, 13, 21, 36, 28, 45, 53, 4, 6, 51, 43, 30, 38, 19, 11, 62], &
            [8, 6])

        ! The rank of A should be 3
        if (mtx_rank(a) /= 3) then
            print '(A)', "Test Failed: Matrix Rank"
        else
            print '(A)', "Test Passed: Matrix Rank"
        end if
    end subroutine

! ******************************************************************************
! TRIANGULAR MATRIX MULTIPLICATION TESTS
! ------------------------------------------------------------------------------
    subroutine test_tri_mtx_mult_1()
        ! Parameters
        integer(i32), parameter :: n = 100
        real(dp), parameter :: alpha = 1.0d0
        real(dp), parameter :: beta = -3.0d0
        real(dp), parameter :: tol = 1.0d-8

        ! Local Variables
        logical :: check
        integer(i32) :: j
        real(dp) :: a(n,n), b(n,n), bans(n,n)

        ! Initialization
        check = .true.
        call random_number(a)
        do j = 1, n
            a(j+1:n,j) = 0.0d0
        end do

        ! Test1 (beta = 0)
        call tri_mtx_mult(.true., alpha, a, 0.0d0, b)
        bans = alpha * matmul(transpose(a), a)
        if (.not.is_mtx_equal(b, bans, tol)) then
            check = .false.
            print '(A)', "Test Failed: Triangular Matrix Update - Test 1A"
        end if

        if (check) &
            print '(A)', "Test Passed: Triangular Matrix Update - Test 1A"
    end subroutine
end module
