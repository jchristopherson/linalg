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
        real(real64) :: a(m, n), aref(m, n), tau(m), q(n, n)
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
        real(real64) :: a(m, n), aref(m, n), tau(m), q(n, n), temp(m, n)
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
        complex(real64) :: a(m, n), aref(m, n), tau(m), q(n, n)
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
        complex(real64) :: a(m, n), aref(m, n), tau(m), q(n, n)
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
        real(real64) :: a(m, n), l(m, n), tau(m), q(n, n), c1(n, n), c2(n, n), &
            ans(n, n), c3(n), c4(n), ans2(n)
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

        ! Extract L and Q
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

        ! ----------
        ! Q**T

        ! Compute C = Q**T * C
        c1 = c2
        call mult_lq(.true., .true., a, tau, c1)

        ! Compute the answer
        call mtx_mult(.true., .false., 1.0d0, q, c2, 0.0d0, ans)

        ! Test
        if (.not.is_mtx_equal(c1, ans, tol)) then
            rst = .false.
            print '(A)', "Test Failed: LQ Multiplication Test 3"
        end if

        ! Vector RHS
        c3 = c4
        call mult_lq(.true., a, tau, c3)

        ! Compute the answer
        call mtx_mult(.true., 1.0d0, q, c4, 0.0d0, ans2)

        ! Test
        if (.not.is_mtx_equal(c3, ans2, tol)) then
            rst = .false.
            print '(A)', "Test Failed: LQ Multiplication Test 4"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_lq_mult_ud() result(rst)
        ! Parameters
        integer(int32), parameter :: m = 50
        integer(int32), parameter :: n = 60
        real(real64), parameter :: tol = 1.0d-8

        ! Local Variables
        real(real64) :: a(m, n), l(m, n), tau(m), q(n, n), c1(n, n), c2(n, n), &
            ans(n, n), c3(n), c4(n), ans2(n)
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

        ! Extract L and Q
        call form_lq(l, tau, q)

        ! Compute C = Q * C
        call mult_lq(.true., .false., a, tau, c1)

        ! Compute the answer
        ans = matmul(q, c2)

        ! Test
        if (.not.is_mtx_equal(c1, ans, tol)) then
            rst = .false.
            print '(A)', "Test Failed: Underdetermined LQ Multiplication Test 1"
        end if

        ! Vector RHS
        call mult_lq(.false., a, tau, c3)

        ! Compute the answer
        ans2 = matmul(q, c4)

        ! Test
        if (.not.is_mtx_equal(c3, ans2, tol)) then
            rst = .false.
            print '(A)', "Test Failed: Underdetermined LQ Multiplication Test 2"
        end if

        ! ----------
        ! Q**T

        ! Compute C = Q**T * C
        c1 = c2
        call mult_lq(.true., .true., a, tau, c1)

        ! Compute the answer
        ans = matmul(transpose(q), c2)

        ! Test
        if (.not.is_mtx_equal(c1, ans, tol)) then
            rst = .false.
            print '(A)', "Test Failed: LQ Multiplication Test 3"
        end if

        ! Vector RHS
        c3 = c4
        call mult_lq(.true., a, tau, c3)

        ! Compute the answer
        ans2 = matmul(transpose(q), c4)

        ! Test
        if (.not.is_mtx_equal(c3, ans2, tol)) then
            rst = .false.
            print '(A)', "Test Failed: LQ Multiplication Test 4"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_lq_mult_cmplx() result(rst)
        ! Parameters
        integer(int32), parameter :: m = 50
        integer(int32), parameter :: n = 50
        real(real64), parameter :: tol = 1.0d-8
        complex(real64), parameter :: zero = (0.0d0, 0.0d0)
        complex(real64), parameter :: one = (1.0d0, 0.0d0)

        ! Local Variables
        complex(real64) :: a(m, n), l(m, n), tau(m), q(n, n), c1(n, n), &
            c2(n, n), ans(n, n), c3(n), c4(n), ans2(n)
        real(real64) :: ar(m, n), ai(m, n), cr(n, n), ci(n, n)
        logical :: rst

        ! Initialization
        rst = .true.
        call random_number(ar)
        call random_number(ai)
        a = cmplx(ar, ai, real64)
        call random_number(cr)
        call random_number(ci)
        c1 = cmplx(cr, ci, real64)
        c3 = c1(:,1)
        c2 = c1
        c4 = c3

        ! Compute the LQ factorization of A
        call lq_factor(a, tau)
        l = a

        ! Extract L and Q
        call form_lq(l, tau, q)

        ! Compute C = Q * C
        call mult_lq(.true., .false., a, tau, c1)

        ! Compute the answer
        ans = matmul(q, c2)

        ! Test
        if (.not.is_mtx_equal(c1, ans, tol)) then
            rst = .false.
            print '(A)', "Test Failed: Complex LQ Multiplication Test 1"
        end if

        ! Vector RHS
        call mult_lq(.false., a, tau, c3)

        ! Compute the answer
        ans2 = matmul(q, c4)

        ! Test
        if (.not.is_mtx_equal(c3, ans2, tol)) then
            rst = .false.
            print '(A)', "Test Failed: Complex LQ Multiplication Test 2"
        end if

        ! ----------
        ! Q**H

        ! Compute C = Q**H * C
        ! c1 = c2
        ! call mult_lq(.true., .true., a, tau, c1)

        ! ! Compute the answer
        ! call mtx_mult(LA_HERMITIAN_TRANSPOSE, LA_NO_OPERATION, one, q, c2, zero, ans)

        ! ! Test
        ! if (.not.is_mtx_equal(c1, ans, tol)) then
        !     rst = .false.
        !     print '(A)', "Test Failed: Complex LQ Multiplication Test 3"
        ! end if

        ! ! Vector RHS
        ! c3 = c4
        ! call mult_lq(.true., a, tau, c3)

        ! ! Compute the answer
        ! call mtx_mult(LA_HERMITIAN_TRANSPOSE, one, q, c4, zero, ans2)

        ! ! Test
        ! if (.not.is_mtx_equal(c3, ans2, tol)) then
        !     rst = .false.
        !     print '(A)', "Test Failed: Complex LQ Multiplication Test 4"
        ! end if
    end function

! ------------------------------------------------------------------------------
    function test_lq_mult_cmplx_ud() result(rst)
        ! Parameters
        integer(int32), parameter :: m = 50
        integer(int32), parameter :: n = 60
        real(real64), parameter :: tol = 1.0d-8
        complex(real64), parameter :: zero = (0.0d0, 0.0d0)
        complex(real64), parameter :: one = (1.0d0, 0.0d0)

        ! Local Variables
        complex(real64) :: a(m, n), l(m, n), tau(m), q(n, n), c1(n, n), &
            c2(n, n), ans(n, n), c3(n), c4(n), ans2(n)
        real(real64) :: ar(m, n), ai(m, n), cr(n, n), ci(n, n)
        logical :: rst

        ! Initialization
        rst = .true.
        call random_number(ar)
        call random_number(ai)
        a = cmplx(ar, ai, real64)
        call random_number(cr)
        call random_number(ci)
        c1 = cmplx(cr, ci, real64)
        c3 = c1(:,1)
        c2 = c1
        c4 = c3

        ! Compute the LQ factorization of A
        call lq_factor(a, tau)
        l = a

        ! Extract L and Q
        call form_lq(l, tau, q)

        ! Compute C = Q * C
        call mult_lq(.true., .false., a, tau, c1)

        ! Compute the answer
        ans = matmul(q, c2)

        ! Test
        if (.not.is_mtx_equal(c1, ans, tol)) then
            rst = .false.
            print '(A)', "Test Failed: Underdetermined Complex LQ Multiplication Test 1"
        end if

        ! Vector RHS
        call mult_lq(.false., a, tau, c3)

        ! Compute the answer
        ans2 = matmul(q, c4)

        ! Test
        if (.not.is_mtx_equal(c3, ans2, tol)) then
            rst = .false.
            print '(A)', "Test Failed: Underdetermined Complex LQ Multiplication Test 2"
        end if

        ! ----------
        ! Q**H

        ! Compute C = Q**H * C
        ! c1 = c2
        ! call mult_lq(.true., .true., a, tau, c1)

        ! ! Compute the answer
        ! call mtx_mult(LA_HERMITIAN_TRANSPOSE, LA_NO_OPERATION, one, q, c2, zero, ans)

        ! ! Test
        ! if (.not.is_mtx_equal(c1, ans, tol)) then
        !     rst = .false.
        !     print '(A)', "Test Failed: Underdetermined Complex LQ Multiplication Test 3"
        ! end if

        ! ! Vector RHS
        ! c3 = c4
        ! call mult_lq(.true., a, tau, c3)

        ! ! Compute the answer
        ! call mtx_mult(LA_HERMITIAN_TRANSPOSE, one, q, c4, zero, ans2)

        ! ! Test
        ! if (.not.is_mtx_equal(c3, ans2, tol)) then
        !     rst = .false.
        !     print '(A)', "Test Failed: Underdetermined Complex LQ Multiplication Test 4"
        ! end if
    end function

! ------------------------------------------------------------------------------
    function test_lq_mult_right() result(rst)
        ! Parameters
        integer(int32), parameter :: m = 50
        integer(int32), parameter :: n = 50
        real(real64), parameter :: tol = 1.0d-8

        ! Local Variables
        logical :: rst
        real(real64) :: a(m, n), l(m, n), tau(m), q(n, n), c1(m, n), c2(m, n), &
            ans(m, n)

        ! Initialization
        rst = .true.
        call random_number(a)
        call random_number(c1)
        c2 = c1

        ! Compute the LQ factorization
        call lq_factor(a, tau)
        l = a

        ! Extract L & Q
        call form_lq(l, tau, q)

        ! Compute C = C * Q
        call mult_lq(.false., .false., a, tau, c1)

        ! Compute the answer
        ans = matmul(c2, q)

        ! Test
        if (.not.is_mtx_equal(c1, ans, tol)) then
            rst = .false.
            print '(A)', "Test Failed: LQ Right Multiplication Test 1"
        end if

        ! Transpose
        c1 = c2
        call mult_lq(.false., .true., a, tau, c1)

        ! Compute the answer: C = C * Q**T
        call mtx_mult(.false., .true., 1.0d0, c2, q, 0.0d0, ans)

        ! Test
        if (.not.is_mtx_equal(c1, ans, tol)) then
            rst = .false.
            print '(A)', "Test Failed: LQ Right Multiplication Test 2"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_lq_mult_right_cmplx() result(rst)
        ! Parameters
        integer(int32), parameter :: m = 50
        integer(int32), parameter :: n = 50
        real(real64), parameter :: tol = 1.0d-8
        complex(real64), parameter :: zero = (0.0d0, 0.0d0)
        complex(real64), parameter :: one = (1.0d0, 0.0d0)

        ! Local Variables
        logical :: rst
        complex(real64) :: a(m, n), l(m, n), tau(m), q(n, n), c1(m, n), &
            c2(m, n), ans(m, n)
        real(real64) :: ar(m, n), ai(m, n), cr(m, n), ci(m, n)

        ! Initialization
        rst = .true.
        call random_number(ar)
        call random_number(ai)
        call random_number(cr)
        call random_number(ci)
        a = cmplx(ar, ai, real64)
        c1 = cmplx(cr, ci, real64)
        c2 = c1

        ! Compute the LQ factorization
        call lq_factor(a, tau)
        l = a

        ! Extract L & Q
        call form_lq(l, tau, q)

        ! Compute C = C * Q
        call mult_lq(.false., .false., a, tau, c1)

        ! Compute the answer
        ans = matmul(c2, q)

        ! Test
        if (.not.is_mtx_equal(c1, ans, tol)) then
            rst = .false.
            print '(A)', "Test Failed: Complex LQ Right Multiplication Test 1"
        end if

        ! Transpose
        ! c1 = c2
        ! call mult_lq(.false., .true., a, tau, c1)

        ! ! Compute the answer: C = C * Q**H
        ! call mtx_mult(LA_NO_OPERATION, LA_HERMITIAN_TRANSPOSE, one, c2, q, &
        !     zero, ans)

        ! ! Test
        ! if (.not.is_mtx_equal(c1, ans, tol)) then
        !     rst = .false.
        !     print '(A)', "Test Failed: Complex LQ Right Multiplication Test 2"
        ! end if
    end function

! ------------------------------------------------------------------------------
    function test_lq_mult_right_ud() result(rst)
        ! Parameters
        integer(int32), parameter :: m = 50
        integer(int32), parameter :: n = 60
        real(real64), parameter :: tol = 1.0d-8

        ! Local Variables
        logical :: rst
        real(real64) :: a(m, n), l(m, n), tau(m), q(n, n), c1(m, n), c2(m, n), &
            ans(m, n)

        ! Initialization
        rst = .true.
        call random_number(a)
        call random_number(c1)
        c2 = c1

        ! Compute the LQ factorization
        call lq_factor(a, tau)
        l = a

        ! Extract L & Q
        call form_lq(l, tau, q)

        ! Compute C = C * Q
        call mult_lq(.false., .false., a, tau, c1)

        ! Compute the answer
        ans = matmul(c2, q)

        ! Test
        if (.not.is_mtx_equal(c1, ans, tol)) then
            rst = .false.
            print '(A)', "Test Failed: Underdetermined LQ Right Multiplication Test 1"
        end if

        ! Transpose
        c1 = c2
        call mult_lq(.false., .true., a, tau, c1)

        ! Compute the answer: C = C * Q**T
        call mtx_mult(.false., .true., 1.0d0, c2, q, 0.0d0, ans)

        ! Test
        if (.not.is_mtx_equal(c1, ans, tol)) then
            rst = .false.
            print '(A)', "Test Failed: Underdetermined LQ Right Multiplication Test 2"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_lq_mult_right_cmplx_ud() result(rst)
        ! Parameters
        integer(int32), parameter :: m = 50
        integer(int32), parameter :: n = 60
        real(real64), parameter :: tol = 1.0d-8
        complex(real64), parameter :: zero = (0.0d0, 0.0d0)
        complex(real64), parameter :: one = (1.0d0, 0.0d0)

        ! Local Variables
        logical :: rst
        complex(real64) :: a(m, n), l(m, n), tau(m), q(n, n), c1(m, n), &
            c2(m, n), ans(m, n)
        real(real64) :: ar(m, n), ai(m, n), cr(m, n), ci(m, n)

        ! Initialization
        rst = .true.
        call random_number(ar)
        call random_number(ai)
        call random_number(cr)
        call random_number(ci)
        a = cmplx(ar, ai, real64)
        c1 = cmplx(cr, ci, real64)
        c2 = c1

        ! Compute the LQ factorization
        call lq_factor(a, tau)
        l = a

        ! Extract L & Q
        call form_lq(l, tau, q)

        ! Compute C = C * Q
        call mult_lq(.false., .false., a, tau, c1)

        ! Compute the answer
        ans = matmul(c2, q)

        ! Test
        if (.not.is_mtx_equal(c1, ans, tol)) then
            rst = .false.
            print '(A)', "Test Failed: Underdetermined Complex LQ Right Multiplication Test 1"
        end if

        ! Transpose
        ! c1 = c2
        ! call mult_lq(.false., .true., a, tau, c1)

        ! ! Compute the answer: C = C * Q**H
        ! call mtx_mult(LA_NO_OPERATION, LA_HERMITIAN_TRANSPOSE, one, c2, q, &
        !     zero, ans)

        ! ! Test
        ! if (.not.is_mtx_equal(c1, ans, tol)) then
        !     rst = .false.
        !     print '(A)', "Test Failed: Underdetermined Complex LQ Right Multiplication Test 2"
        ! end if
    end function

! ------------------------------------------------------------------------------
end module