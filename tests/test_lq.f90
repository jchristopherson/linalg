! test_lq.f90

! Tests for LQ factorization/solution operations.
module test_lq
    use iso_fortran_env
    use test_core
    use linalg
    use fortran_test_helper
    implicit none
contains
! ******************************************************************************
! LQ FACTORIZATION TESTS
! ------------------------------------------------------------------------------
    function test_lq_factor() result(rst)
        ! Parameters
        integer(int32), parameter :: m = 50
        integer(int32), parameter :: n = 50

        ! Local Variables
        real(real64) :: a(m, n)
        real(real64), allocatable :: q(:,:), l(:,:)
        logical :: rst

        ! Initialization
        rst = .true.
        call create_random_array(a)

        ! Compute the LQ factorization of A
        call lq_factor(a, l = l, q = q)

        ! Perform the check L*Q = A
        if (.not.assert(matmul(l, q), a, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: LQ Factorization Test 1"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_lq_factor_ud() result(rst)
        ! Parameters
        integer(int32), parameter :: m = 50
        integer(int32), parameter :: n = 60

        ! Local Variables
        real(real64) :: a(m, n)
        real(real64), allocatable :: q(:,:), l(:,:)
        logical :: rst

        ! Initialization
        rst = .true.
        call create_random_array(a)

        ! Compute the LQ factorization of A
        call lq_factor(a, l = l, q = q)

        ! Perform the check
        if (.not.assert(matmul(l, q), a, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Underdetermined LQ Factorization Test 1"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_lq_factor_cmplx() result(rst)
        ! Parameters
        integer(int32), parameter :: m = 50
        integer(int32), parameter :: n = 50

        ! Local Variables
        complex(real64) :: a(m, n)
        complex(real64), allocatable :: l(:,:), q(:,:)
        logical :: rst

        ! Initialization
        rst = .true.
        call create_random_array(a)

        ! Compute the LQ factorization of A
        call lq_factor(a, l = l, q = q)

        ! Perform the check
        if (.not.assert(matmul(l, q), a, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Complex LQ Factorization Test 1"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_lq_factor_ud_cmplx() result(rst)
        ! Parameters
        integer(int32), parameter :: m = 50
        integer(int32), parameter :: n = 60

        ! Local Variables
        complex(real64) :: a(m, n)
        complex(real64), allocatable :: l(:,:), q(:,:)
        logical :: rst

        ! Initialization
        rst = .true.
        call create_random_array(a)

        ! Compute the LQ factorization of A
        call lq_factor(a, l = l, q = q)

        ! Perform the check
        if (.not.assert(matmul(l, q), a, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Underdetermined Complex LQ Factorization Test 1"
        end if
    end function

! ******************************************************************************
! LQ MULTIPLICATION TEST
! ------------------------------------------------------------------------------
    function test_lq_mult() result(rst)
        ! Parameters
        integer(int32), parameter :: m = 50
        integer(int32), parameter :: n = 50

        ! Local Variables
        real(real64) :: a(m, n), c1(n, n), c2(n, n), ans(n, n), c3(n), &
            c4(n), ans2(n)
        real(real64), allocatable :: tau(:), lq(:,:), q(:,:)
        logical :: rst

        ! Initialization
        rst = .true.
        call create_random_array(a)
        call create_random_array(c1)
        call create_random_array(c3)
        c2 = c1
        c4 = c3

        ! Compute the LQ factorization of A
        call lq_factor(a, tau = tau, lq = lq, q = q)

        ! Compute C = Q * C
        c1 = mult_lq(.true., .false., lq, tau, c1)

        ! Compute the answer
        ans = matmul(q, c2)

        ! Test
        if (.not.assert(c1, ans, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: LQ Multiplication Test 1"
        end if

        ! Vector RHS
        c3 = mult_lq(.false., lq, tau, c3)

        ! Compute the answer
        ans2 = matmul(q, c4)

        ! Test
        if (.not.assert(c3, ans2, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: LQ Multiplication Test 2"
        end if

        ! ----------
        ! Q**T

        ! Compute C = Q**T * C
        c1 = c2
        c1 = mult_lq(.true., .true., lq, tau, c1)

        ! Compute the answer
        call mtx_mult(.true., .false., 1.0d0, q, c2, 0.0d0, ans)

        ! Test
        if (.not.assert(c1, ans, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: LQ Multiplication Test 3"
        end if

        ! Vector RHS
        c3 = c4
        c3 = mult_lq(.true., lq, tau, c3)

        ! Compute the answer
        call mtx_mult(.true., 1.0d0, q, c4, 0.0d0, ans2)

        ! Test
        if (.not.assert(c3, ans2, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: LQ Multiplication Test 4"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_lq_mult_ud() result(rst)
        ! Parameters
        integer(int32), parameter :: m = 50
        integer(int32), parameter :: n = 60

        ! Local Variables
        real(real64) :: a(m, n), c1(n, n), c2(n, n), ans(n, n), c3(n), &
            c4(n), ans2(n)
        real(real64), allocatable :: tau(:), lq(:,:), q(:,:)
        logical :: rst

        ! Initialization
        rst = .true.
        call create_random_array(a)
        call create_random_array(c1)
        call create_random_array(c3)
        c2 = c1
        c4 = c3

        ! Compute the LQ factorization of A
        call lq_factor(a, tau = tau, lq = lq, q = q)

        ! Compute C = Q * C
        c1 = mult_lq(.true., .false., lq, tau, c1)

        ! Compute the answer
        ans = matmul(q, c2)

        ! Test
        if (.not.assert(c1, ans, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Underdetermined LQ Multiplication Test 1"
        end if

        ! Vector RHS
        c3 = mult_lq(.false., lq, tau, c3)

        ! Compute the answer
        ans2 = matmul(q, c4)

        ! Test
        if (.not.assert(c3, ans2, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Underdetermined LQ Multiplication Test 2"
        end if

        ! ----------
        ! Q**T

        ! Compute C = Q**T * C
        c1 = c2
        c1 = mult_lq(.true., .true., lq, tau, c1)

        ! Compute the answer
        ans = matmul(transpose(q), c2)

        ! Test
        if (.not.assert(c1, ans, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: LQ Multiplication Test 3"
        end if

        ! Vector RHS
        c3 = c4
        c3 = mult_lq(.true., lq, tau, c3)

        ! Compute the answer
        ans2 = matmul(transpose(q), c4)

        ! Test
        if (.not.assert(c3, ans2, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: LQ Multiplication Test 4"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_lq_mult_cmplx() result(rst)
        ! Parameters
        integer(int32), parameter :: m = 50
        integer(int32), parameter :: n = 50
        complex(real64), parameter :: zero = (0.0d0, 0.0d0)
        complex(real64), parameter :: one = (1.0d0, 0.0d0)

        ! Local Variables
        complex(real64) :: a(m, n), c1(n, n), c2(n, n), ans(n, n), c3(n), &
            c4(n), ans2(n)
        complex(real64), allocatable :: tau(:), lq(:,:), q(:,:)
        logical :: rst

        ! Initialization
        rst = .true.
        call create_random_array(a)
        call create_random_array(c1)
        c3 = c1(:,1)
        c2 = c1
        c4 = c3

        ! Compute the LQ factorization of A
        call lq_factor(a, tau = tau, lq = lq, q = q)

        ! Compute C = Q * C
        c1 = mult_lq(.true., .false., lq, tau, c1)

        ! Compute the answer
        ans = matmul(q, c2)

        ! Test
        if (.not.assert(c1, ans, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Complex LQ Multiplication Test 1"
        end if

        ! Vector RHS
        c3 = mult_lq(.false., lq, tau, c3)

        ! Compute the answer
        ans2 = matmul(q, c4)

        ! Test
        if (.not.assert(c3, ans2, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Complex LQ Multiplication Test 2"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_lq_mult_cmplx_ud() result(rst)
        ! Parameters
        integer(int32), parameter :: m = 50
        integer(int32), parameter :: n = 60
        complex(real64), parameter :: zero = (0.0d0, 0.0d0)
        complex(real64), parameter :: one = (1.0d0, 0.0d0)

        ! Local Variables
        complex(real64) :: a(m, n), c1(n, n), c2(n, n), ans(n, n), c3(n), &
            c4(n), ans2(n)
        complex(real64), allocatable :: lq(:,:), tau(:), q(:,:)
        logical :: rst

        ! Initialization
        rst = .true.
        call create_random_array(a)
        call create_random_array(c1)
        c3 = c1(:,1)
        c2 = c1
        c4 = c3

        ! Compute the LQ factorization of A
        call lq_factor(a, tau = tau, lq = lq, q = q)

        ! Compute C = Q * C
        c1 = mult_lq(.true., .false., lq, tau, c1)

        ! Compute the answer
        ans = matmul(q, c2)

        ! Test
        if (.not.assert(c1, ans, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Underdetermined Complex LQ Multiplication Test 1"
        end if

        ! Vector RHS
        c3 = mult_lq(.false., lq, tau, c3)

        ! Compute the answer
        ans2 = matmul(q, c4)

        ! Test
        if (.not.assert(c3, ans2, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Underdetermined Complex LQ Multiplication Test 2"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_lq_mult_right() result(rst)
        ! Parameters
        integer(int32), parameter :: m = 50
        integer(int32), parameter :: n = 50

        ! Local Variables
        logical :: rst
        real(real64) :: a(m, n), c1(m, n), c2(m, n), ans(m, n)
        real(real64), allocatable :: tau(:), lq(:,:), q(:,:)

        ! Initialization
        rst = .true.
        call create_random_array(a)
        call create_random_array(c1)
        c2 = c1

        ! Compute the LQ factorization
        call lq_factor(a, tau = tau, lq = lq, q = q)

        ! Compute C = C * Q
        c1 = mult_lq(.false., .false., lq, tau, c1)

        ! Compute the answer
        ans = matmul(c2, q)

        ! Test
        if (.not.assert(c1, ans, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: LQ Right Multiplication Test 1"
        end if

        ! Transpose
        c1 = c2
        c1 = mult_lq(.false., .true., lq, tau, c1)

        ! Compute the answer: C = C * Q**T
        call mtx_mult(.false., .true., 1.0d0, c2, q, 0.0d0, ans)

        ! Test
        if (.not.assert(c1, ans, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: LQ Right Multiplication Test 2"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_lq_mult_right_cmplx() result(rst)
        ! Parameters
        integer(int32), parameter :: m = 50
        integer(int32), parameter :: n = 50
        complex(real64), parameter :: zero = (0.0d0, 0.0d0)
        complex(real64), parameter :: one = (1.0d0, 0.0d0)

        ! Local Variables
        logical :: rst
        complex(real64) :: a(m, n), c1(m, n), c2(m, n), ans(m, n)
        complex(real64), allocatable :: tau(:), lq(:,:), q(:,:)

        ! Initialization
        rst = .true.
        call create_random_array(a)
        call create_random_array(c1)
        c2 = c1

        ! Compute the LQ factorization
        call lq_factor(a, tau = tau, lq = lq, q = q)

        ! Compute C = C * Q
        c1 = mult_lq(.false., .false., lq, tau, c1)

        ! Compute the answer
        ans = matmul(c2, q)

        ! Test
        if (.not.assert(c1, ans, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Complex LQ Right Multiplication Test 1"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_lq_mult_right_ud() result(rst)
        ! Parameters
        integer(int32), parameter :: m = 50
        integer(int32), parameter :: n = 60

        ! Local Variables
        logical :: rst
        real(real64) :: a(m, n), c1(m, n), c2(m, n), ans(m, n)
        real(real64), allocatable :: tau(:), lq(:,:), q(:,:)

        ! Initialization
        rst = .true.
        call create_random_array(a)
        call create_random_array(c1)
        c2 = c1

        ! Compute the LQ factorization
        call lq_factor(a, tau = tau, lq = lq, q = q)

        ! Compute C = C * Q
        c1 = mult_lq(.false., .false., lq, tau, c1)

        ! Compute the answer
        ans = matmul(c2, q)

        ! Test
        if (.not.assert(c1, ans, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Underdetermined LQ Right Multiplication Test 1"
        end if

        ! Transpose
        c1 = c2
        c1 = mult_lq(.false., .true., lq, tau, c1)

        ! Compute the answer: C = C * Q**T
        call mtx_mult(.false., .true., 1.0d0, c2, q, 0.0d0, ans)

        ! Test
        if (.not.assert(c1, ans, REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Underdetermined LQ Right Multiplication Test 2"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_lq_mult_right_cmplx_ud() result(rst)
        ! Parameters
        integer(int32), parameter :: m = 50
        integer(int32), parameter :: n = 60
        complex(real64), parameter :: zero = (0.0d0, 0.0d0)
        complex(real64), parameter :: one = (1.0d0, 0.0d0)

        ! Local Variables
        logical :: rst
        complex(real64) :: a(m, n), c1(m, n), c2(m, n), ans(m, n)
        complex(real64), allocatable :: tau(:), lq(:,:), q(:,:)

        ! Initialization
        rst = .true.
        call create_random_array(a)
        call create_random_array(c1)
        c2 = c1

        ! Compute the LQ factorization
        call lq_factor(a, tau = tau, lq = lq, q = q)

        ! Compute C = C * Q
        c1 = mult_lq(.false., .false., lq, tau, c1)

        ! Compute the answer
        ans = matmul(c2, q)

        ! Test
        if (.not.assert(c1, ans, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Underdetermined Complex LQ Right Multiplication Test 1"
        end if
    end function

! ------------------------------------------------------------------------------
end module