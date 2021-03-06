! test_misc.f90

! Tests miscellaneous routines.
module test_misc
    use, intrinsic :: iso_fortran_env, only : int32, real64
    use linalg_core
    use test_core
    implicit none
contains
! ******************************************************************************
! DIAGONAL MATRIX MULTIPLICATION TEST
! ------------------------------------------------------------------------------
    function test_diagonal_mtx_mult() result(rst)
        ! Parameters
        integer(int32), parameter :: m = 30
        integer(int32), parameter :: n = 30
        integer(int32), parameter :: k = 30
        real(real64), parameter :: tol = 1.0d-8
        real(real64), parameter :: alpha = 0.5d0
        real(real64), parameter :: beta = 0.25d0

        ! Local Variables
        integer(int32) :: i
        real(real64), dimension(m, n) :: c1, ans1
        real(real64), dimension(k, n) :: b1, b2, ans2
        real(real64), dimension(m, k) :: d1
        real(real64), dimension(k) :: d1v
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
    end function

! ******************************************************************************
! RANK 1 UPDATE TEST
! ------------------------------------------------------------------------------
    function test_rank1_update() result(rst)
        ! Parameters
        integer(int32), parameter :: m = 50
        integer(int32), parameter :: n = 20
        real(real64), parameter :: alpha = 0.5d0
        real(real64), parameter :: tol = 1.0d-12

        ! Local Variables
        real(real64), dimension(m, n) :: a, a1, b
        real(real64), dimension(m, 1) :: x
        real(real64), dimension(n, 1) :: y
        logical :: rst

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
            rst = .false.
        else
            rst = .true.
        end if
    end function

! ******************************************************************************
! MATRIX RANK TESTS
! ------------------------------------------------------------------------------
    ! REF: http://www.mathworks.com/help/matlab/ref/pinv.html?s_tid=srchtitle
    function test_rank() result(rst)
        ! Local Variables
        real(real64), dimension(8, 6) :: a
        logical :: rst

        ! Initialization
        a = reshape([64, 9, 17, 40, 32, 41, 49, 8, 2, 55, 47, 26, 34, 23, 15, &
            58, 3, 54, 46, 27, 35, 22, 14, 59, 61, 12, 20, 37, 29, 44, 52, 5, &
            60, 13, 21, 36, 28, 45, 53, 4, 6, 51, 43, 30, 38, 19, 11, 62], &
            [8, 6])

        ! The rank of A should be 3
        if (mtx_rank(a) /= 3) then
            print '(A)', "Test Failed: Matrix Rank"
            rst = .false.
        else
            rst = .true.
        end if
    end function

! ******************************************************************************
! TRIANGULAR MATRIX MULTIPLICATION TESTS
! ------------------------------------------------------------------------------
    function test_tri_mtx_mult_1() result(rst)
        ! Parameters
        integer(int32), parameter :: n = 100
        real(real64), parameter :: alpha = 1.5d0
        real(real64), parameter :: beta = -3.0d0
        real(real64), parameter :: tol = 1.0d-8

        ! Local Variables
        logical :: check, rst
        integer(int32) :: j
        real(real64) :: a(n,n), b(n,n), bans(n,n)

        ! Initialization
        check = .true.
        call random_number(a)
        do j = 1, n
            a(j+1:n,j) = 0.0d0
        end do

        ! Test 1 (beta = 0)
        call tri_mtx_mult(.true., alpha, a, 0.0d0, b)
        bans = alpha * matmul(transpose(a), a)
        if (.not.is_mtx_equal(b, bans, tol)) then
            check = .false.
            print '(A)', "Test Failed: Triangular Matrix Update - Test 1A"
        end if
        rst = check
        
        ! Test 2 (beta /= 0)
        check = .true.
        call tri_mtx_mult(.true., alpha, a, beta, b)
        bans = alpha * matmul(transpose(a), a) + beta * bans
        if (.not.is_mtx_equal(b, bans, tol)) then
            check = .false.
            rst = .false.
            print '(A)', "Test Failed: Triangular Matrix Update - Test 1B"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_tri_mtx_mult_2() result(rst)
        ! Parameters
        integer(int32), parameter :: n = 100
        real(real64), parameter :: alpha = 1.5d0
        real(real64), parameter :: beta = -3.0d0
        real(real64), parameter :: tol = 1.0d-8

        ! Local Variables
        logical :: check, rst
        integer(int32) :: j
        real(real64) :: a(n,n), b(n,n), bans(n,n)

        ! Initialization
        check = .true.
        call random_number(a)
        do j = 2, n
            a(1:j-1,j) = 0.0d0
        end do

        ! Test 1 (beta = 0)
        call tri_mtx_mult(.false., alpha, a, 0.0d0, b)
        bans = alpha * matmul(a, transpose(a))
        if (.not.is_mtx_equal(b, bans, tol)) then
            check = .false.
            print '(A)', "Test Failed: Triangular Matrix Update - Test 2A"
        end if
        rst = check
        
        ! Test 2 (beta /= 0)
        check = .true.
        call tri_mtx_mult(.false., alpha, a, beta, b)
        bans = alpha * matmul(a, transpose(a)) + beta * bans
        if (.not.is_mtx_equal(b, bans, tol)) then
            check = .false.
            rst = .false.
            print '(A)', "Test Failed: Triangular Matrix Update - Test 2B"
        end if
    end function

! ******************************************************************************
! TEST MATRIX MULTIPLICATION
! ------------------------------------------------------------------------------
    !
    function test_mtx_mult_1() result(rst)
        ! Arguments
        logical :: rst

        ! Local Variables
        integer(int32), parameter :: m = 50
        integer(int32), parameter :: n = 20
        integer(int32), parameter :: k = 30
        real(real64), parameter :: alpha = 1.0d0
        real(real64), parameter :: beta = 0.0d0
        real(real64), parameter :: tol = 1.0d-8
        real(real64) :: a(m, k), b(n, k), c(m, n), ans(m, n)

        ! Initialization
        rst = .true.
        call random_number(a)
        call random_number(b)
        call random_number(c)

        ! Compute the solution
        ans = alpha * matmul(a, transpose(b)) + beta * c

        ! Test
        call mtx_mult(.false., .true., alpha, a, b, beta, c)
        if (.not.is_mtx_equal(c, ans, tol)) then
            rst = .false.
            print '(A)', "Test Failed: Matrix Multiplication - Test 1"
        end if
    end function

! ******************************************************************************
! TEST TRIANGULAR SYSTEM SOLUTIONS
! ------------------------------------------------------------------------------
    function test_tri_mtx_solve_1() result(rst)
        ! Arguments
        logical :: rst

        ! Local Variables
        integer(int32), parameter :: n = 200
        integer(int32), parameter :: nrhs = 20
        real(real64), parameter :: alpha = 1.0d0
        real(real64), parameter :: tol = 1.0d-8
        integer(int32) :: j
        real(real64) :: a(n,n), b1(n,nrhs), x1(n,nrhs), check1(n,nrhs)

        ! Initialization - upper triangular systems
        rst = .true.
        call random_number(a)
        do j = 1, n
            a(j+1:n,j) = 0.0d0
            a(j,j) = 2.0d0  ! Make sure we don't have too small of diagonal
        end do
        call random_number(b1)
        x1 = b1

        ! Compute the solution to A * X1 = B1
        call solve_triangular_system(.true., .true., .false., .true., &
            alpha, a, x1)

        ! Verify that A * X1 = B1
        check1 = matmul(a, x1)
        if (.not.is_mtx_equal(b1, check1, tol)) then
            rst = .false.
            print '(A)', "Test Failed: Tri Matrix Solve - Test 1A"
        end if
    end function

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end module
