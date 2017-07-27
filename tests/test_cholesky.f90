! test_cholesky.f90

! Tests the Cholesky factorization/solution operations
module test_cholesky
    use linalg_constants
    use test_core
    use linalg_core, only : rank1_update
    use linalg_factor, only : cholesky_factor, cholesky_rank1_update, &
        cholesky_rank1_downdate
    use linalg_solve, only : solve_cholesky
    implicit none
contains
! ******************************************************************************
! CHOLESKY FACTORIZATION TESTS
! ------------------------------------------------------------------------------
    function test_cholesky_factor() result(rst)
        ! Parameters
        integer(i32), parameter :: n = 100
        integer(i32), parameter :: nrhs = 20
        real(dp), parameter :: tol = 1.0d-8

        ! Local Variables
        real(dp), dimension(n, n) :: a, a1, u, l
        real(dp), dimension(n, nrhs) :: b, b1, b2
        logical :: rst

        ! Initialization
        rst = .true.
        call random_number(a1)
        call random_number(b)
        a = matmul(a1, transpose(a1))
        u = a
        l = a
        b1 = b
        b2 = b

        ! Test 1: A = L * L**T
        call cholesky_factor(l, .false.)
        if (.not.is_mtx_equal(a, matmul(l, transpose(l)), tol)) then
            rst = .false.
            print '(A)', "Test Failed: Cholesky Factorization Test 1"
        end if

        ! Test 2: A = U**T * U
        call cholesky_factor(u, .true.)
        if (.not.is_mtx_equal(a, matmul(transpose(u), u), tol)) then
            rst = .false.
            print '(A)', "Test Failed: Cholesky Factorization Test 2"
        end if

        ! Test 3: Solve L*L**T * X = B
        call solve_cholesky(.false., l, b1)
        if (.not.is_mtx_equal(matmul(a, b1), b, tol)) then
            rst = .false.
            print '(A)', "Test Failed: Cholesky Factorization Test 3"
        end if

        ! Test 4: Solve U**T * U * X = B
        call solve_cholesky(.true., u, b2)
        if (.not.is_mtx_equal(matmul(a, b2), b, tol)) then
            rst = .false.
            print '(A)', "Test Failed: Cholesky Factorization Test 4"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_cholesky_rank1_update() result(rst)
        ! Parameters
        integer(i32), parameter :: n = 100
        real(dp), parameter :: tol = 1.0d-8

        ! Local Variables
        real(dp), dimension(n, n) :: a, a1, r
        real(dp), dimension(n) :: u
        logical :: rst

        ! Initialization
        rst = .true.
        call random_number(a1)
        call random_number(u)
        a = matmul(a1, transpose(a1))
        r = a

        ! Compute the Cholesky factorization of the original matrix
        call cholesky_factor(r)

        ! Update the original matrix
        call rank1_update(1.0d0, u, u, a)

        ! Update the factored matrix
        call cholesky_rank1_update(r, u)

        ! Test
        if (.not.is_mtx_equal(a, matmul(transpose(r), r), tol)) then
            rst = .false.
            print '(A)', "Test Failed: Cholesky Rank 1 Update Test 1"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_cholesky_rank1_downdate() result(rst)
        ! Parameters
        integer(i32), parameter :: n = 100
        real(dp), parameter :: tol = 1.0d-8

        ! Local Variables
        real(dp), dimension(n, n) :: a, a1, r
        real(dp), dimension(n) :: u
        logical :: rst

        ! Initialization
        rst = .true.
        call random_number(a1)
        call random_number(u)
        a = matmul(a1, transpose(a1))

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
        if (.not.is_mtx_equal(a, matmul(transpose(r), r), tol)) then
            rst = .false.
            print '(A)', "Test Failed: Cholesky Rank 1 Downdate Test 1"
        end if
    end function

end module
