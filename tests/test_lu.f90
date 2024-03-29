! test_lu.f90

! Tests the LU factorization routines
module test_lu
    use, intrinsic :: iso_fortran_env, only : int32, real64
    use test_core
    use linalg
    use fortran_test_helper
contains
! ******************************************************************************
! LU FACTORIZATION TEST
! ------------------------------------------------------------------------------
    function test_lu_factor() result(rst)
        ! Parameters
        integer(int32), parameter :: n = 75

        ! Local Variables
        real(real64), dimension(n, n) :: a, a1, l, u, p
        integer(int32), dimension(n) :: ipvt
        logical :: rst

        ! Initialization
        rst = .true.
        call create_random_array(a)
        a1 = a  ! Forces us to keep a copy of the original matrix

        ! Compute the factorization
        call lu_factor(a1, ipvt)

        ! Extract L, U, and P to determine if P * A = L * U
        l = a1
        call form_lu(l, ipvt, u, p)
        if (.not.assert(matmul(p, a), matmul(l, u), tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: LU Factorization Test"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_lu_solve() result(rst)
        ! Parameters
        integer(int32), parameter :: n = 75
        integer(int32), parameter :: nrhs = 20

        ! Local Variables
        real(real64), dimension(n, n) :: a, a1
        real(real64), dimension(n, nrhs) :: b, x
        integer(int32), dimension(n) :: ipvt
        logical :: rst

        ! Initialization
        rst = .true.
        call create_random_array(a)
        call create_random_array(b)
        a1 = a
        x = b

        ! Factor A
        call lu_factor(a1, ipvt)

        ! Solve for X
        call solve_lu(a1, ipvt, x)

        ! Test by determining if A * X = B
        if (.not.assert(matmul(a, x), b, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: LU Factorization & Solution Test"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_lu_factor_cmplx() result(rst)
        ! Parameters
        integer(int32), parameter :: n = 75

        ! Local Variables
        complex(real64), dimension(n, n) :: a, a1, l, u
        real(real64), dimension(n, n) :: p
        integer(int32) :: i, j
        integer(int32), dimension(n) :: ipvt
        logical :: rst

        ! Initialization
        rst = .true.
        call create_random_array(a)
        a1 = a

        ! Compute the factorization
        call lu_factor(a1, ipvt)

        ! Extract L, U, and P to determine if P * A = L * U
        l = a1
        call form_lu(l, ipvt, u, p)
        if (.not.assert(matmul(p, a), matmul(l, u), tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Complex-Valued LU Factorization Test"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_lu_solve_cmplx() result(rst)
        ! Parameters
        integer(int32), parameter :: n = 75
        integer(int32), parameter :: nrhs = 20

        ! Local Variables
        complex(real64), dimension(n, n) :: a, a1
        complex(real64), dimension(n, nrhs) :: b, x
        integer(int32) :: i, j
        integer(int32), dimension(n) :: ipvt
        logical :: rst

        ! Initialization
        rst = .true.
        call create_random_array(a)
        call create_random_array(b)
        a1 = a
        x = b

        ! Factor A
        call lu_factor(a1, ipvt)

        ! Solve for X
        call solve_lu(a1, ipvt, x)

        ! Test by determining if A * X = B
        if (.not.assert(matmul(a, x), b, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Complex-Valued LU Factorization & Solution Test"
        end if
    end function

end module
