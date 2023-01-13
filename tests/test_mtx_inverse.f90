! test_mtx_inverse.f90

! Tests matrix inversion routines
module test_mtx_inverse
    use, intrinsic :: iso_fortran_env, only : int32, real64
    use test_core
    use linalg
    use fortran_test_helper
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

        ! Local Variables
        real(real64), dimension(m, n) :: a, a1
        real(real64), dimension(n, m) :: ainv
        real(real64), dimension(m, nrhs) :: b
        real(real64), dimension(n, nrhs) :: x
        logical :: rst

        ! Initialization
        rst = .true.
        call create_random_array(a)
        call create_random_array(b)
        a1 = a

        ! Compute the inverse
        call mtx_pinverse(a1, ainv)

        ! Compute X = inv(A) * B
        x = matmul(ainv, b)

        ! Test: A * X = B
        if (.not.assert(matmul(a, x), b, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Pseudo-Inverse Test 1"
        end if
    end function

! ------------------------------------------------------------------------------
    ! REF: http://www.mathworks.com/help/matlab/ref/pinv.html?s_tid=srchtitle
    function test_pinv_od() result(rst)
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
        if (.not.assert(x1, x, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Overdetermined Pseudo-Inverse Test 1"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_inv() result(rst)
        ! Parameters
        integer(int32), parameter :: n = 100

        ! Local Variables
        real(real64), dimension(n, n) :: a, a1, ainv
        logical :: rst

        ! Initialization
        rst = .true.
        call create_random_array(a)
        a1 = a

        ! Compute the inverse of A
        call mtx_inverse(a1)

        ! Compute the inverse using the already tested pseudo-inverse
        call mtx_pinverse(a, ainv)

        ! Test
        if (.not.assert(ainv, a1, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Matrix Inverse Test 1"
        end if
    end function
! ------------------------------------------------------------------------------
    function test_pinv_cmplx() result(rst)
        ! Parameters
        integer(int32), parameter :: m = 60
        integer(int32), parameter :: n = 60
        integer(int32), parameter :: nrhs = 20

        ! Local Variables
        real(real64), dimension(m, n) :: ar, ai
        real(real64), dimension(m, nrhs) :: br, bi
        complex(real64), dimension(m, n) :: a, a1
        complex(real64), dimension(n, m) :: ainv
        complex(real64), dimension(m, nrhs) :: b
        complex(real64), dimension(n, nrhs) :: x
        logical :: rst

        ! Initialization
        rst = .true.
        call create_random_array(ar)
        call create_random_array(ai)
        a = cmplx(ar, ai, real64)
        a1 = a

        call create_random_array(br)
        call create_random_array(bi)
        b = cmplx(br, bi, real64)

        ! Compute the inverse
        call mtx_pinverse(a1, ainv)

        ! Compute X = inv(A) * B
        x = matmul(ainv, b)

        ! Test: A * X = B
        if (.not.assert(matmul(a, x), b, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Complex Pseudo-Inverse Test 1"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_pinv_od_cmplx() result(rst)
        ! Parameters
        integer(int32), parameter :: m = 80
        integer(int32), parameter :: n = 60

        ! Parameters
        complex(real64), parameter :: zero = (0.0d0, 0.0d0)
        complex(real64), parameter :: one = (1.0d0, 0.0d0)

        ! Local Variables
        real(real64), dimension(m, n) :: ar, ai
        complex(real64), dimension(m, n) :: a, a1
        complex(real64), dimension(n, m) :: ainv
        complex(real64), dimension(n, n) :: identity, check
        logical :: rst
        integer(int32) :: i

        ! Initialization
        rst = .true.
        call create_random_array(ar)
        call create_random_array(ai)
        a = cmplx(ar, ai, real64)
        a1 = a

        identity = zero
        do i = 1, size(identity, 1)
            identity(i,i) = one
        end do
        
        ! Compute the inverse
        call mtx_pinverse(a1, ainv)

        ! Compute A+ * A - should = I
        check = matmul(ainv, a)
        if (.not.assert(check, identity, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Complex Pseudo-Inverse Test 2"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_pinv_ud_cmplx() result(rst)
        ! Parameters
        integer(int32), parameter :: m = 60
        integer(int32), parameter :: n = 80

        ! Parameters
        complex(real64), parameter :: zero = (0.0d0, 0.0d0)
        complex(real64), parameter :: one = (1.0d0, 0.0d0)

        ! Local Variables
        real(real64), dimension(m, n) :: ar, ai
        complex(real64), dimension(m, n) :: a, a1
        complex(real64), dimension(n, m) :: ainv
        complex(real64), dimension(m, m) :: identity, check
        logical :: rst
        integer(int32) :: i

        ! Initialization
        rst = .true.
        call create_random_array(ar)
        call create_random_array(ai)
        a = cmplx(ar, ai, real64)
        a1 = a

        identity = zero
        do i = 1, size(identity, 1)
            identity(i,i) = one
        end do
        
        ! Compute the inverse
        call mtx_pinverse(a1, ainv)

        ! Compute A * A+ - should = I
        check = matmul(a, ainv)
        if (.not.assert(check, identity, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Complex Pseudo-Inverse Test 3"
        end if
    end function

! ------------------------------------------------------------------------------
end module
