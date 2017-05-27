! lu.f90

!> @brief \b lu
!!
!! @par Purpose
!! Provides a set of routines for solving systems of equations using LU 
!! factorization.
module lu
    use, intrinsic :: iso_fortran_env, only : int32, real64
    implicit none
    private
    public :: lu_factor

! ******************************************************************************
! NUMERIC TYPE CONSTANTS
! ------------------------------------------------------------------------------
    !> @brief Defines a double-precision (64-bit) floating-point type.
    integer, parameter :: dp = real64
    !> @brief Defines a 32-bit signed integer type.
    integer, parameter :: i32 = int32

contains
! ------------------------------------------------------------------------------
    !> @brief Computes the LU factorization of an M-by-N matrix.
    !!
    !! @param[in,out] a On input, the M-by-N matrix on which to operate.  On
    !! output, the LU factored matrix in the form [L\\U] where the unit diagonal
    !! elements of L are not stored.
    !! @param[out] ipvt An MIN(M, N)-element array used to track row-pivot
    !!  operations.  The array stored pivot information such that row I is
    !!  interchanged with row IPVT(I).
    !!
    !! @par Usage
    !! To solve a system of N equations of N unknowns using LU factorization,
    !! the following code will suffice.
    !! @code{.f90}
    !! ! Solve the system: A*X = B, where A is an N-by-N matrix, and B and X are
    !! ! N-by-NRHS in size.
    !!
    !! ! Variables
    !! real(dp), dimension(n, n) :: a
    !! real(dp), dimension(n, nrhs) :: b
    !!
    !! ! Define the array used to track row pivots.
    !! integer(i32), dimension(n) :: pvt
    !!
    !! ! Initialize A and B...
    !!
    !! ! Compute the LU factorization of A.  On output, A contains [L\U].
    !! call lu_factor(a, pvt)
    !!
    !! ! Solve A*X = B for X - Note: X overwrites B.
    !! call solve_lu(a, pvt, b)
    !! @endcode
    !!
    !! @par Notes
    !! This routine utilizes the LAPACK routine DGETRF.
    !!
    !! @par See Also
    !! - [Wikipedia](https://en.wikipedia.org/wiki/LU_decomposition)
    !! - [Wolfram MathWorld](http://mathworld.wolfram.com/LUDecomposition.html)
    subroutine lu_factor(a, ipvt)
        ! Arguments
        real(dp), intent(inout), dimension(:,:) :: a
        integer(i32), intent(out), dimension(:) :: ipvt

        ! Parameters
        real(dp), parameter :: zero = 0.0d0
        real(dp), parameter :: one = 1.0d0

        ! Local Variables
        integer(i32) :: m, n, mn, flag

        ! Initialization
        m = size(a, 1)
        n = size(a, 2)
        mn = min(m, n)

        ! Input Check
        flag = 0
        if (size(ipvt) /= mn) then
            ! ERROR: IPVT not sized correctly
        end if

        ! Compute the LU factorization by calling the LAPACK routine DGETRF
        call DGETRF(m, n, a, m, ipvt, flag)

        ! If flag > 0, the matrix is singular.  Notice, flag should not be
        ! able to be < 0 as we've already verrified inputs prior to making the
        ! call to LAPACK
        if (flag > 0) then
            ! ERROR: Singular matrix
        end if
    end subroutine

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end module
