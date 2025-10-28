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
    function test_pinv_od() result(rst)
        ! Parameters
        integer(int32), parameter :: m = 80
        integer(int32), parameter :: n = 60

        ! Local Variables
        real(real64), dimension(m, n) :: a, a1
        real(real64), dimension(n, m) :: ainv
        real(real64), dimension(n, n) :: identity, check
        logical :: rst
        integer(int32) :: i

        ! Initialization
        rst = .true.
        call create_random_array(a)
        a1 = a

        identity = 0.0d0
        do i = 1, size(identity, 1)
            identity(i,i) = 1.0d0
        end do

        ! Compute the inverse
        call mtx_pinverse(a1, ainv)

        ! Compute A+ * A - should = I
        check = matmul(ainv, a)
        if (.not.assert(check, identity, tol = REAL64_TOL)) then
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
        complex(real64), dimension(m, n) :: a, a1
        complex(real64), dimension(n, m) :: ainv
        complex(real64), dimension(m, nrhs) :: b
        complex(real64), dimension(n, nrhs) :: x
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
        complex(real64), dimension(m, n) :: a, a1
        complex(real64), dimension(n, m) :: ainv
        complex(real64), dimension(n, n) :: identity, check
        logical :: rst
        integer(int32) :: i

        ! Initialization
        rst = .true.
        call create_random_array(a)
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
        complex(real64), dimension(m, n) :: a, a1
        complex(real64), dimension(n, m) :: ainv
        complex(real64), dimension(m, m) :: identity, check
        logical :: rst
        integer(int32) :: i

        ! Initialization
        rst = .true.
        call create_random_array(a)
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
    function test_inv_cmplx() result(rst)
        ! Parameters
        integer(int32), parameter :: n = 100

        ! Local Variables
        complex(real64), dimension(n, n) :: a, a1, ainv
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
            print '(A)', "Test Failed: Complex-Valued Matrix Inverse Test 1"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_inverse_pure() result(rst)
        use linear_algebra

        ! Arguments
        logical :: rst

        ! Parameters & Variables
        integer(int32), parameter :: n = 100
        integer(int32), parameter :: nrhs = 20
        real(real64) :: a(n, n), b(n, nrhs), x(n, nrhs), ainv(n, n)

        ! Initialization
        rst = .true.
        call create_random_array(a)
        call create_random_array(b)

        ! Test
        ainv = inverse(a)
        x = matmul(ainv, b)
        if (.not.assert(matmul(a, x), b, REAL64_TOL)) then
            rst = .false.
            print '(A)', "TEST FAILED: test_inverse_pure"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_pinverse_pure_1() result(rst)
        use linear_algebra

        ! Arguments
        logical :: rst

        ! Parameters & Variables
        integer(int32), parameter :: n = 100
        integer(int32), parameter :: nrhs = 20
        real(real64) :: a(n, n), b(n, nrhs), x(n, nrhs), ainv(n, n)

        ! Initialization
        rst = .true.
        call create_random_array(a)
        call create_random_array(b)

        ! Test
        ainv = pinverse(a)
        x = matmul(ainv, b)
        if (.not.assert(matmul(a, x), b, REAL64_TOL)) then
            rst = .false.
            print '(A)', "TEST FAILED: test_pinverse_pure_1"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_pinverse_pure_2() result(rst)
        use linear_algebra

        ! Arguments
        logical :: rst

        ! Parameters & Variables
        integer(int32), parameter :: m = 100
        integer(int32), parameter :: n = 70
        real(real64) :: a(m, n), ainv(n, m), ident(n, n)

        ! Initialization
        rst = .true.
        call create_random_array(a)
        ident = identity(n)

        ! Compute the inverse
        ainv = pinverse(a)

        ! Compute pinv(A) * A, and that should equal I
        if (.not.assert(matmul(ainv, a), ident, REAL64_TOL)) then
            rst = .false.
            print '(A)', "TEST FAILED: test_pinverse_pure_2"
        end if
    end function

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end module
