! test_cholesky.f90

! Tests the Cholesky factorization/solution operations
module test_cholesky
    use, intrinsic :: iso_fortran_env, only : int32, real64
    use test_core
    use linalg
    use fortran_test_helper
    implicit none
contains
! ******************************************************************************
! CHOLESKY FACTORIZATION TESTS
! ------------------------------------------------------------------------------
    function test_cholesky_factor() result(rst)
        ! Parameters
        integer(int32), parameter :: n = 100
        integer(int32), parameter :: nrhs = 20

        ! Local Variables
        real(real64), dimension(n, n) :: a, u, l
        real(real64), dimension(n, nrhs) :: b, b1, b2
        logical :: rst

        ! Initialization
        rst = .true.
        call create_random_array(a, mtype = POSITIVE_DEFINITE_MATRIX)
        call create_random_array(b)
        u = a
        l = a
        b1 = b
        b2 = b

        ! Test 1: A = L * L**T
        call cholesky_factor(l, .false.)
        if (.not.assert(a, matmul(l, transpose(l)), tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Cholesky Factorization Test 1"
        end if

        ! Test 2: A = U**T * U
        call cholesky_factor(u, .true.)
        if (.not.assert(a, matmul(transpose(u), u), tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Cholesky Factorization Test 2"
        end if

        ! Test 3: Solve L*L**T * X = B
        call solve_cholesky(.false., l, b1)
        if (.not.assert(matmul(a, b1), b, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Cholesky Factorization Test 3"
        end if

        ! Test 4: Solve U**T * U * X = B
        call solve_cholesky(.true., u, b2)
        if (.not.assert(matmul(a, b2), b, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Cholesky Factorization Test 4"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_cholesky_rank1_update() result(rst)
        ! Parameters
        integer(int32), parameter :: n = 100

        ! Local Variables
        real(real64), dimension(n, n) :: a, r
        real(real64), dimension(n) :: u
        logical :: rst

        ! Initialization
        rst = .true.
        call create_random_array(a, mtype = POSITIVE_DEFINITE_MATRIX)
        call create_random_array(u)
        r = a

        ! Compute the Cholesky factorization of the original matrix
        call cholesky_factor(r)

        ! Update the original matrix
        call rank1_update(1.0d0, u, u, a)

        ! Update the factored matrix
        call cholesky_rank1_update(r, u)

        ! Test
        if (.not.assert(a, matmul(transpose(r), r), tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Cholesky Rank 1 Update Test 1"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_cholesky_rank1_downdate() result(rst)
        ! Parameters
        integer(int32), parameter :: n = 100

        ! Local Variables
        real(real64), dimension(n, n) :: a, r
        real(real64), dimension(n) :: u
        logical :: rst

        ! Initialization
        rst = .true.
        call create_random_array(a, mtype = POSITIVE_DEFINITE_MATRIX)
        call create_random_array(u)

        ! Start with a positive definite matrix, and then update it
        call rank1_update(1.0d0, u, u, a)
        r = a

        ! Compute the Cholesky factorization of the original matrix
        call cholesky_factor(r)

        ! Update the original matrix: A = A - u * u**T
        call rank1_update(-1.0d0, u, u, a)

        ! Update the factored matrix
        call cholesky_rank1_downdate(r, u)

        ! Test
        if (.not.is_mtx_equal(a, matmul(transpose(r), r), tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Cholesky Rank 1 Downdate Test 1"
        end if
    end function

! ******************************************************************************
! COMPLEX VALUED ROUTINES
! ------------------------------------------------------------------------------
    function test_cholesky_factor_cmplx() result(rst)
        ! Parameters
        integer(int32), parameter :: n = 10
        integer(int32), parameter :: nrhs = 20

        ! Local Variables
        complex(real64), dimension(n, n) :: a, u, l
        complex(real64), dimension(n, nrhs) :: b, b1, b2
        logical :: rst

        ! Initialization
        rst = .true.
        call create_random_array( &
            a, &
            mtype = POSITIVE_DEFINITE_MATRIX)
        call create_random_array(b)
        u = a
        l = a
        b1 = b
        b2 = b

        ! Test 1: A = L * L**T
        call cholesky_factor(l, .false.)
        if (.not.assert(a, matmul(l, transpose(l)), tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Complex Cholesky Factorization Test 1"
        end if

        ! Test 2: A = U**T * U
        call cholesky_factor(u, .true.)
        if (.not.assert(a, matmul(transpose(u), u), tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Complex Cholesky Factorization Test 2"
        end if

        ! Test 3: Solve L*L**T * X = B
        call solve_cholesky(.false., l, b1)
        if (.not.assert(matmul(a, b1), b, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Complex Cholesky Factorization Test 3"
        end if

        ! Test 4: Solve U**T * U * X = B
        call solve_cholesky(.true., u, b2)
        if (.not.assert(matmul(a, b2), b, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Complex Cholesky Factorization Test 4"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_cholesky_rank1_update_cmplx() result(rst)
        ! Parameters
        integer(int32), parameter :: n = 100
        complex(real64), parameter :: one = (1.0d0, 0.0d0)

        ! Local Variables
        complex(real64), dimension(n, n) :: a, r
        complex(real64), dimension(n) :: u
        logical :: rst

        ! Initialization
        rst = .true.
        call create_random_array(a, mtype = POSITIVE_DEFINITE_MATRIX)
        call create_random_array(u)
        r = a

        ! Compute the Cholesky factorization of the original matrix
        call cholesky_factor(r)

        ! Update the original matrix
        call rank1_update(one, u, u, a)

        ! Update the factored matrix
        call cholesky_rank1_update(r, u)

        ! Test
        if (.not.assert(a, matmul(transpose(r), r), tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Complex-Valued Cholesky Rank 1 Update Test 1"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_cholesky_rank1_downdate_cmplx() result(rst)
        ! Parameters
        integer(int32), parameter :: n = 100
        complex(real64), parameter :: one = (1.0d0, 0.0d0)

        ! Local Variables
        complex(real64), dimension(n, n) :: a, r
        complex(real64), dimension(n) :: u
        logical :: rst

        ! Initialization
        rst = .true.
        call create_random_array(a, mtype = POSITIVE_DEFINITE_MATRIX)
        call create_random_array(u)

        ! Start with a positive definite matrix, and then update it
        call rank1_update(one, u, u, a)
        r = a

        ! Compute the Cholesky factorization of the original matrix
        call cholesky_factor(r)

        ! Update the original matrix: A = A - u * u**T
        call rank1_update(-one, u, u, a)

        ! Update the factored matrix
        call cholesky_rank1_downdate(r, u)

        ! Test
        if (.not.is_mtx_equal(a, matmul(transpose(r), r), tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Complex-Valued Cholesky Rank 1 Downdate Test 1"
        end if
    end function

! ------------------------------------------------------------------------------
end module
