! test_immutable.f90

module test_immutable
    use, intrinsic :: iso_fortran_env, only : int32, real64
    use test_core
    use linalg_immutable
    implicit none
contains
! ******************************************************************************
! LU FACTORIZATION
! ------------------------------------------------------------------------------
    function test_im_lu_factor() result(rst)
        ! Parameters
        integer(int32), parameter :: n = 75
        real(real64), parameter :: tol = 1.0d-8

        ! Local Variables
        real(real64), dimension(n, n) :: a
        type(lu_results) :: lu_rst
        logical :: rst

        ! Initialization
        rst = .true.
        call random_number(a) 

        ! Compute the factorization
        lu_rst = mat_lu(a)

        ! Ensuure P * A = L * U
        if (.not.is_mtx_equal(matmul(lu_rst%p, a), &
                matmul(lu_rst%l, lu_rst%u), tol)) then
            rst = .false.
            print '(A)', "Test Failed: Immutable LU Factorization Test"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_im_lu_solve() result(rst)
        ! Parameters
        integer(int32), parameter :: n = 75
        integer(int32), parameter :: nrhs = 20
        real(real64), parameter :: tol = 1.0d-8

        ! Local Variables
        real(real64) :: a(n, n), b(n, nrhs), x(n, nrhs), y(n, nrhs)
        type(lu_results) :: lu_rst
        logical :: rst

        ! Initialization
        rst = .true.
        call random_number(a)
        call random_number(b)

        ! Compute the factorization
        lu_rst = mat_lu(a)

        ! Solve L * Y = P * B, for Y
        y = mat_solve_lower_tri(lu_rst%l, matmul(lu_rst%p, b))

        ! Solve U * X = Y, for X
        x = mat_solve_upper_tri(lu_rst%u, y)

        ! Ensure A * X = B
        if (.not.is_mtx_equal(matmul(a, x), b, tol)) then
            rst = .false.
            print '(A)', "Test Failed: Immutable LU Solution Test"
        end if
    end function

end module
