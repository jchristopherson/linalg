! lu.f90

!> @brief \b lu
!!
!! @par Purpose
!! Provides a set of routines for solving systems of equations using LU 
!! factorization.
module lu
    use, intrinsic :: iso_fortran_env, only : int32, real64
    use ferror, only : errors
    use error_flags
    implicit none
    private
    public :: lu_factor
    public :: solve_lu
    public :: form_lu

! ******************************************************************************
! NUMERIC TYPE CONSTANTS
! ------------------------------------------------------------------------------
    !> @brief Defines a double-precision (64-bit) floating-point type.
    integer, parameter :: dp = real64
    !> @brief Defines a 32-bit signed integer type.
    integer, parameter :: i32 = int32

! ******************************************************************************
! INTERFACES
! ------------------------------------------------------------------------------
    !> @brief Solves a system of LU-factored equations.
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
    !! @par See Also
    !! - [Wikipedia](https://en.wikipedia.org/wiki/LU_decomposition)
    !! - [Wolfram MathWorld](http://mathworld.wolfram.com/LUDecomposition.html)
    interface solve_lu
        module procedure :: solve_lu_mtx
        module procedure :: solve_lu_vec
    end interface

! ------------------------------------------------------------------------------
    !> @brief Extracts the L and U matrices from the condensed [L\\U] storage 
    !! format used by the @ref lu_factor.
    interface form_lu
        module procedure :: form_lu_all
        module procedure :: form_lu_only
    end interface

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
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - LA_ARRAY_SIZE_ERROR: Occurs if @p ipvt is not sized appropriately.
    !!  - LA_SINGULAR_MATRIX_ERROR: Occurs as a warning if @p a is found to be
    !!      singular.
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
    subroutine lu_factor(a, ipvt, err)
        ! Arguments
        real(dp), intent(inout), dimension(:,:) :: a
        integer(i32), intent(out), dimension(:) :: ipvt
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        real(dp), parameter :: zero = 0.0d0
        real(dp), parameter :: one = 1.0d0

        ! Local Variables
        integer(i32) :: m, n, mn, flag
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg

        ! Initialization
        m = size(a, 1)
        n = size(a, 2)
        mn = min(m, n)
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        flag = 0
        if (size(ipvt) /= mn) then
            ! ERROR: IPVT not sized correctly
            call errmgr%report_error("lu_factor", &
                "Incorrectly sized input array IPVT, argument 2.", &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Compute the LU factorization by calling the LAPACK routine DGETRF
        call DGETRF(m, n, a, m, ipvt, flag)

        ! If flag > 0, the matrix is singular.  Notice, flag should not be
        ! able to be < 0 as we've already verrified inputs prior to making the
        ! call to LAPACK
        if (flag > 0) then
            ! WARNING: Singular matrix
            write(errmsg, '(AI0A)') &
                "Singular matrix encountered (row ", flag, ")"
            call errmgr%report_warning("lu_factor", trim(errmsg), &
                LA_SINGULAR_MATRIX_ERROR)
        end if
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Solves a system of LU-factored equations.
    !!
    !! @param[in] a The N-by-N LU factored matrix as output by @ref lu_factor.
    !! @param[in] ipvt The N-element pivot array as output by @ref lu_factor.
    !! @param[in,out] b On input, the N-by-NRHS right-hand-side matrix.  On
    !!  output, the N-by-NRHS solution matrix.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - LA_ARRAY_SIZE_ERROR: Occurs if any of the input array sizes are 
    !!      incorrect.
    !!
    !! @par Notes
    !! The routine is based upon the LAPACK routine DGETRS.
    subroutine solve_lu_mtx(a, ipvt, b, err)
        ! Arguments
        real(dp), intent(in), dimension(:,:) :: a
        integer(i32), intent(in), dimension(:) :: ipvt
        real(dp), intent(inout), dimension(:,:) :: b
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        integer(i32) :: n, nrhs, flag
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg

        ! Initialization
        n = size(a, 1)
        nrhs = size(b, 2)
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        flag = 0
        if (size(a, 2) /= n) then
            flag = 1
        else if (size(ipvt) /= n) then
            flag = 2
        else if (size(b, 1) /= n) then
            flag = 3
        end if
        if (flag /= 0) then
            ! One of the input arrays is not sized correctly
            write(errmsg, '(AI0A)') "Input number ", flag, &
                " is not sized correctly"
            call errmgr%report_error("solve_lu_mtx", errmsg, &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Call DGETRS
        call DGETRS("N", n, nrhs, a, n, ipvt, b, n, flag)
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Solves a system of LU-factored equations.
    !!
    !! @param[in] a The N-by-N LU factored matrix as output by @ref lu_factor.
    !! @param[in] ipvt The N-element pivot array as output by @ref lu_factor.
    !! @param[in,out] b On input, the N-element right-hand-side array.  On
    !!  output, the N-element solution array.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - LA_ARRAY_SIZE_ERROR: Occurs if any of the input array sizes are 
    !!      incorrect.
    !!
    !! @par Notes
    !! The routine is based upon the LAPACK routine DGETRS.
    subroutine solve_lu_vec(a, ipvt, b, err)
        ! Arguments
        real(dp), intent(in), dimension(:,:) :: a
        integer(i32), intent(in), dimension(:) :: ipvt
        real(dp), intent(inout), dimension(:) :: b
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        integer(i32) :: n, flag
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg

        ! Initialization
        n = size(a, 1)
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        flag = 0
        if (size(a, 2) /= n) then
            flag = 1
        else if (size(ipvt) /= n) then
            flag = 2
        else if (size(b) /= n) then
            flag = 3
        end if
        if (flag /= 0) then
            ! One of the input arrays is not sized correctly
            write(errmsg, '(AI0A)') "Input number ", flag, &
                " is not sized correctly"
            call errmgr%report_error("solve_lu_vec", errmsg, &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Call DGETRS
        call DGETRS("N", n, 1, a, n, ipvt, b, n, flag)
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Extracts the L, U, and P matrices from the output of the
    !! @ref lu_factor routine.
    !!
    !! @param[in,out] lu On input, the N-by-N matrix as output by
    !!  @ref lu_factor.  On output, the N-by-N lower triangular matrix L.
    !! @param[in] ipvt The N-element pivot array as output by
    !!  @ref lu_factor.
    !! @param[out] u An N-by-N matrix where the U matrix will be written.
    !! @param[out] p An N-by-N matrix where the row permutation matrix will be
    !!  written.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - LA_ARRAY_SIZE_ERROR: Occurs if any of the input array sizes are 
    !!      incorrect.
    !!
    !! @par Remarks
    !! This routine allows extraction of the actual "L", "U", and "P" matrices
    !! of the decomposition.  To use these matrices to solve the system A*X = B,
    !! the following approach is used.
    !!
    !! 1. First, solve the linear system: L*Y = P*B for Y.
    !! 2. Second, solve the linear system: U*X = Y for X.
    !!
    !! Notice, as both L and U are triangular in structure, the above equations
    !! can be solved by forward and backward substitution.
    !!
    !! @par See Also
    !! - [Wikipedia](https://en.wikipedia.org/wiki/LU_decomposition)
    !! - [Wolfram MathWorld](http://mathworld.wolfram.com/LUDecomposition.html)
    subroutine form_lu_all(lu, ipvt, u, p, err)
        ! Arguments
        real(dp), intent(inout), dimension(:,:) :: lu
        integer(i32), intent(in), dimension(:) :: ipvt
        real(dp), intent(out), dimension(:,:) :: u, p
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        integer(i32) :: j, jp, n, flag
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg

        ! Parameters
        real(dp), parameter :: zero = 0.0d0
        real(dp), parameter :: one = 1.0d0

        ! Initialization
        n = size(lu, 1)
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        flag = 0
        if (size(lu, 2) /= n) then
            flag = 1
        else if (size(ipvt) /= n) then
            flag = 2
        else if (size(u, 1) /= n .or. size(u, 2) /= n) then
            flag = 3
        else if (size(p, 1) /= n .or. size(p, 2) /= n) then
            flag = 4
        end if
        if (flag /= 0) then
            ! One of the input arrays is not sized correctly
            write(errmsg, '(AI0A)') "Input number ", flag, &
                " is not sized correctly"
            call errmgr%report_error("form_lu_all", errmsg, &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Ensure P starts off as an identity matrix
        call dlaset('A', n, n, zero, one, p, n)

        ! Process
        do j = 1, n
            ! Define the pivot matrix
            jp = ipvt(j)
            if (j /= jp) call DSWAP(n, p(j,1:n), n, p(jp,1:n), n)

            ! Build L and U
            u(1:j,j) = lu(1:j,j)
            u(j+1:n,j) = zero

            if (j > 1) lu(1:j-1,j) = zero
            lu(j,j) = one
        end do
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Extracts the L, and U matrices from the output of the
    !! @ref lu_factor routine.
    !!
    !! @param[in,out] lu On input, the N-by-N matrix as output by
    !!  @ref lu_factor.  On output, the N-by-N lower triangular matrix L.
    !! @param[out] u An N-by-N matrix where the U matrix will be written.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - LA_ARRAY_SIZE_ERROR: Occurs if any of the input array sizes are 
    !!      incorrect.
    subroutine form_lu_only(lu, u, err)
        ! Arguments
        real(dp), intent(inout), dimension(:,:) :: lu
        real(dp), intent(out), dimension(:,:) :: u
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        integer(i32) :: j, n, flag
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg

        ! Parameters
        real(dp), parameter :: zero = 0.0d0
        real(dp), parameter :: one = 1.0d0

        ! Initialization
        n = size(lu, 1)
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        flag = 0
        if (size(lu, 2) /= n) then
            flag = 2
        else if (size(u, 1) /= n .or. size(u, 2) /= n) then
            flag = 3
        end if
        if (flag /= 0) then
            ! One of the input arrays is not sized correctly
            write(errmsg, '(AI0A)') "Input number ", flag, &
                " is not sized correctly"
            call errmgr%report_error("form_lu_only", errmsg, &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Process
        do j = 1, n
            ! Build L and U
            u(1:j,j) = lu(1:j,j)
            u(j+1:n,j) = zero

            if (j > 1) lu(1:j-1,j) = zero
            lu(j,j) = one
        end do
    end subroutine

! ------------------------------------------------------------------------------

end module
