! test_lq.f90

! Tests for LQ factorization/solution operations.
module test_lq
    use iso_fortran_env
    use test_core
    use linalg
    implicit none
contains
! ******************************************************************************
! LQ FACTORIZATION TESTS
! ------------------------------------------------------------------------------
    function test_lq_factor() result(rst)
        ! Parameters
        integer(int32), parameter :: m = 50
        integer(int32), parameter :: n = 50
        real(real64), parameter :: tol = 1.0d-8

        ! Local Variables
        real(real64) :: a(m, n), aref(m, n), tau(m), q(m, n)
        logical :: rst

        ! Initialization
        rst = .true.
        call random_number(a)
        aref = a

        ! Compute the LQ factorization of A
        call lq_factor(a, tau)

        ! Extract L and Q and check that L * Q = A
        call form_lq(a, tau, q)

        ! Perform the check
        if (.not.is_mtx_equal(matmul(a, q), aref, tol)) then
            rst = .false.
            print '(A)', "Test Failed: LQ Factorization Test 1"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_lq_factor_ud() result(rst)
        ! Parameters
        integer(int32), parameter :: m = 50
        integer(int32), parameter :: n = 60
        real(real64), parameter :: tol = 1.0d-8

        ! Local Variables
        real(real64) :: a(m, n), aref(m, n), tau(m), q(m, n), temp(m, n)
        logical :: rst

        ! Initialization
        rst = .true.
        call random_number(a)
        aref = a

        ! Compute the LQ factorization of A
        call lq_factor(a, tau)

        ! Extract L and Q and check that L * Q = A
        call form_lq(a, tau, q)

        ! Perform the check
        if (.not.is_mtx_equal(matmul(a(:,1:m), q), aref, tol)) then
            rst = .false.
            print '(A)', "Test Failed: Underdetermined LQ Factorization Test 1"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_lq_factor_cmplx() result(rst)
        ! Parameters
        integer(int32), parameter :: m = 50
        integer(int32), parameter :: n = 50
        real(real64), parameter :: tol = 1.0d-8

        ! Local Variables
        real(real64) :: ar(m, n), ai(m, n)
        complex(real64) :: a(m, n), aref(m, n), tau(m), q(m, n)
        logical :: rst

        ! Initialization
        rst = .true.
        call random_number(ar)
        call random_number(ai)
        a = cmplx(ar, ai, real64)
        aref = a

        ! Compute the LQ factorization of A
        call lq_factor(a, tau)

        ! Extract L and Q and check that L * Q = A
        call form_lq(a, tau, q)

        ! Perform the check
        if (.not.is_mtx_equal(matmul(a, q), aref, tol)) then
            rst = .false.
            print '(A)', "Test Failed: Complex LQ Factorization Test 1"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_lq_factor_cmplx_ud() result(rst)
        ! Parameters
        integer(int32), parameter :: m = 50
        integer(int32), parameter :: n = 60
        real(real64), parameter :: tol = 1.0d-8

        ! Local Variables
        real(real64) :: ar(m, n), ai(m, n)
        complex(real64) :: a(m, n), aref(m, n), tau(m), q(m, n)
        logical :: rst

        ! Initialization
        rst = .true.
        call random_number(ar)
        call random_number(ai)
        a = cmplx(ar, ai, real64)
        aref = a

        ! Compute the LQ factorization of A
        call lq_factor(a, tau)

        ! Extract L and Q and check that L * Q = A
        call form_lq(a, tau, q)

        ! Perform the check
        if (.not.is_mtx_equal(matmul(a(:,1:m), q), aref, tol)) then
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
        real(real64), parameter :: tol = 1.0d-8

        ! Local Variables
        real(real64) :: a(m, n), l(m, n), tau(m), q(m, n), c1(m, n), c2(m, n), &
            ans(m, n), c3(m), c4(m), ans2(m)
        logical :: rst

        ! Initialization
        rst = .true.
        call random_number(a)
        call random_number(c1)
        call random_number(c3)
        c2 = c1
        c4 = c3

        ! Compute the LQ factorization of A
        call lq_factor(a, tau)
        l = a

        ! Extract L and Q and check that L * Q = A
        call form_lq(l, tau, q)

        ! Compute C = Q * C
        call mult_lq(.true., .false., a, tau, c1)

        ! Compute the answer
        ans = matmul(q, c2)

        ! Test
        if (.not.is_mtx_equal(c1, ans, tol)) then
            rst = .false.
            print '(A)', "Test Failed: LQ Multiplication Test 1"
        end if

        ! Vector RHS
        call mult_lq(.false., a, tau, c3)

        ! Compute the answer
        ans2 = matmul(q, c4)

        ! Test
        if (.not.is_mtx_equal(c3, ans2, tol)) then
            rst = .false.
            print '(A)', "Test Failed: LQ Multiplication Test 2"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_lq_mult_ud() result(rst)
        ! Parameters
        integer(int32), parameter :: m = 50
        integer(int32), parameter :: n = 60
        real(real64), parameter :: tol = 1.0d-8

        ! Local Variables
        real(real64) :: a(m, n), l(m, n), tau(m), q(m, n), c1(m, n), c2(m, n), &
            ans(m, n), c3(m), c4(m), ans2(m)
        logical :: rst

        ! Initialization
        rst = .true.
        call random_number(a)
        call random_number(c1)
        call random_number(c3)
        c2 = c1
        c4 = c3

        ! Compute the LQ factorization of A
        call lq_factor(a, tau)
        l = a

        ! Extract L and Q and check that L * Q = A
        call form_lq(l, tau, q)

        ! Compute C = Q * C
        call mult_lq(.true., .false., a(:,1:m), tau, c1)

        ! Compute the answer
        ans = matmul(q(:,1:m), c2)

        ! Test
        if (.not.is_mtx_equal(c1, ans, tol)) then
            rst = .false.
            print '(A)', "Test Failed: Underdetermined LQ Multiplication Test 1"
        end if

        ! Vector RHS
        call mult_lq(.false., a(:,1:m), tau, c3)

        ! Compute the answer
        ans2 = matmul(q(:,1:m), c4)

        ! Test
        if (.not.is_mtx_equal(c3, ans2, tol)) then
            rst = .false.
            print '(A)', "Test Failed: Underdetermined LQ Multiplication Test 2"
        end if
    end function

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end module