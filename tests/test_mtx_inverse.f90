! test_mtx_inverse.f90

! Tests matrix inversion routines
module test_mtx_inverse
    use linalg_constants
    use test_core
    use linalg_solve, only : mtx_inverse, mtx_pinverse
    implicit none
contains
! ******************************************************************************
! PSEUDO-INVERSE TESTS
! ------------------------------------------------------------------------------
    subroutine test_pinv()
        ! Parameters
        integer(i32), parameter :: m = 60
        integer(i32), parameter :: n = 60
        integer(i32), parameter :: nrhs = 20
        real(dp), parameter :: tol = 1.0d-8

        ! Local Variables
        real(dp), dimension(m, n) :: a, a1
        real(dp), dimension(n, m) :: ainv
        real(dp), dimension(m, nrhs) :: b
        real(dp), dimension(n, nrhs) :: x
        logical :: rst

        ! Initialization
        rst = .true.
        call random_number(a)
        call random_number(b)
        a1 = a

        ! Compute the inverse
        call mtx_pinverse(a1, ainv)

        ! Compute X = inv(A) * B
        x = matmul(ainv, b)

        ! Test: A * X = B
        if (.not.is_mtx_equal(matmul(a, x), b, tol)) then
            rst = .false.
            print '(A)', "Test Failed: Pseudo-Inverse Test 1"
        end if
        if (rst) print '(A)', "Test Passed: Psuedo-Inverse"
    end subroutine

! ------------------------------------------------------------------------------
    ! REF: http://www.mathworks.com/help/matlab/ref/pinv.html?s_tid=srchtitle
    subroutine test_pinv_od()
        ! Parameters
        real(dp), parameter :: tol = 1.0d-8

        ! Local Variables
        real(dp), dimension(8, 6) :: a, a1
        real(dp), dimension(6, 8) :: ainv
        real(dp), dimension(8) :: b
        real(dp), dimension(6) :: x, x1
        logical :: rst

        ! Initialization
        rst = .true.
        a = reshape(&
            [64, 9, 17, 40, 32, 41, 49, 8, 2, 55, 47, 26, 34, 23, 15, 58, &
            3, 54, 46, 27, 35, 22, 14, 59, 61, 12, 20, 37, 29, 44, 52, 5, 60, &
            13, 21, 36, 28, 45, 53, 4, 6, 51, 43, 30, 38, 19, 11, 62], [8, 6])
        b = 260.0d0
        x = [1.15384615384615d0, 1.46153846153846d0, 1.38461538461539d0, &
            1.38461538461538d0, 1.46153846153846d0, 1.15384615384615d0]
        a1 = a

        ! Compute the pseudoinverse of A
        call mtx_pinverse(a1, ainv)

        ! Compute X = AINV*B
        x1 = matmul(ainv, b)

        ! Test
        if (.not.is_mtx_equal(x1, x, tol)) then
            rst = .false.
            print '(A)', "Test Failed: Overdetermined Pseudo-Inverse Test 1"
        end if
        if (rst) print '(A)', "Test Passed: Overdetermined Pseudo-Inverse"
    end subroutine

! ------------------------------------------------------------------------------
    subroutine test_inv()
        ! Parameters
        integer(i32), parameter :: n = 100
        real(dp), parameter :: tol = 1.0d-8

        ! Local Variables
        real(dp), dimension(n, n) :: a, a1, ainv
        logical :: rst

        ! Initialization
        rst = .true.
        call random_number(a)
        a1 = a

        ! Compute the inverse of A
        call mtx_inverse(a1)

        ! Compute the inverse using the already tested pseudo-inverse
        call mtx_pinverse(a, ainv)

        ! Test
        if (.not.is_mtx_equal(ainv, a1, tol)) then
            rst = .false.
            print '(A)', "Test Failed: Matrix Inverse Test 1"
        end if
        if (rst) print '(A)', "Test Passed: Matrix Inverse"
    end subroutine


end module
