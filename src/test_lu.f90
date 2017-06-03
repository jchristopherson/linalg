! test_lu.f90

! Tests the LU factorization routines
module test_lu
    use linalg_constants
    use test_core
    use linalg_factor, only : lu_factor, form_lu
    use linalg_solve, only : solve_lu
contains
! ******************************************************************************
! LU FACTORIZATION TEST
! ------------------------------------------------------------------------------
    subroutine test_lu_factor()
        ! Parameters
        integer(i32), parameter :: n = 75
        real(dp), parameter :: tol = 1.0d-8

        ! Local Variables
        real(dp), dimension(n, n) :: a, a1, l, u, p
        integer(i32), dimension(n) :: ipvt
        logical :: rst

        ! Initialization
        rst = .true.
        call random_number(a)
        a1 = a  ! Forces us to keep a copy of the original matrix

        ! Compute the factorization
        call lu_factor(a1, ipvt)

        ! Extract L, U, and P to determine if P * A = L * U
        l = a1
        call form_lu(l, ipvt, u, p)
        if (.not.is_mtx_equal(matmul(p, a), matmul(l, u), tol)) then
            rst = .false.
            print '(A)', "Test Failed: LU Factorization Test 2"
        end if
        if (rst) print '(A)', "Test Passed: LU Factorization Test"
    end subroutine

end module
