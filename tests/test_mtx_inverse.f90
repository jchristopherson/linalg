! test_mtx_inverse.f90

! Tests matrix inversion routines
module test_mtx_inverse
    use, intrinsic :: iso_fortran_env, only : int32, real64
    use test_core
    use linalg_core
    implicit none
contains
! ******************************************************************************
! PSEUDO-INVERSE TESTS
! ------------------------------------------------------------------------------
    function test_pinv() result(rst)
        ! Parameters
        integer(int32), parameter :: m = 60
        integer(int32), parameter :: n = 60
        integer(int32), parameter :: nrhs = 20
        real(real64), parameter :: tol = 1.0d-8

        ! Local Variables
        real(real64), dimension(m, n) :: a, a1
        real(real64), dimension(n, m) :: ainv
        real(real64), dimension(m, nrhs) :: b
        real(real64), dimension(n, nrhs) :: x
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
    end function

! ------------------------------------------------------------------------------
    ! REF: http://www.mathworks.com/help/matlab/ref/pinv.html?s_tid=srchtitle
    function test_pinv_od() result(rst)
        ! Parameters
        real(real64), parameter :: tol = 1.0d-8

        ! Local Variables
        real(real64), dimension(8, 6) :: a, a1
        real(real64), dimension(6, 8) :: ainv
        real(real64), dimension(8) :: b
        real(real64), dimension(6) :: x, x1
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
    end function

! ------------------------------------------------------------------------------
    function test_inv() result(rst)
        ! Parameters
        integer(int32), parameter :: n = 100
        real(real64), parameter :: tol = 1.0d-8

        ! Local Variables
        real(real64), dimension(n, n) :: a, a1, ainv
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
    end function


end module
