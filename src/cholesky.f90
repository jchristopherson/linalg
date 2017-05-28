! cholesky.f90

!> @brief \b cholesky
!!
!! @par Purpose
!! Provides a set of routines for solving systems of equations using Cholesky 
!! factorization.
module cholesky
    use ferror, only : errors
    use linalg_constants
    implicit none
    private
    public :: cholesky_factor
    public :: solve_cholesky
    public :: cholesky_rank1_update

! ******************************************************************************
! INTERFACES
! ------------------------------------------------------------------------------
    !> @brief Solves a system of Cholesky factored equations.
    interface solve_cholesky
        module procedure :: solve_cholesky_mtx
        module procedure :: solve_cholesky_vec
    end interface


contains
! ------------------------------------------------------------------------------
    !> @brief Computes the Cholesky factorization of a symmetric, positive
    !! definite matrix.
    !!
    !! @param[in,out] a On input, the N-by-N matrix to factor.  On output, the
    !!  factored matrix is returned in either the upper or lower triangular
    !!  portion of the matrix, dependent upon the value of @p upper.
    !! @param[in] upper An optional input that, if specified, provides control
    !!  over whether the factorization is computed as A = U**T * U (set to
    !!  true), or as A = L * L**T (set to false).  The default value is true
    !!  such that A = U**T * U.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - LA_ARRAY_SIZE_ERROR: Occurs if @p a is not square.
    !!  - LA_MATRIX_FORMAT_ERROR: Occurs if @p a is not positive definite.
    !!
    !! @par Usage
    !! To solve a system of N equations of N unknowns using Cholesky 
    !! factorization, the following code will suffice.  Notice, the system of
    !! equations must be positive definite.
    !! @code{.f90}
    !! ! Solve the system: A*X = B, where A is an N-by-N matrix, and B and X are
    !! ! N-by-NRHS in size.
    !!
    !! ! Variables
    !! real(dp), dimension(n, n) :: a
    !! real(dp), dimension(n, nrhs) :: b
    !! logical :: upper
    !!
    !! ! Initialize A and B...
    !!
    !! ! Specify that we're using the upper portion of A (remember positive 
    !! ! definite matrices are symmetric)
    !! upper = .true.
    !!
    !! ! Compute the factorization of A.
    !! call cholesky_factor(a, upper)
    !!
    !! ! Solve A*X = B for X - Note: X overwrites B.
    !! call solve_cholesky(upper, a, b)
    !! @endcode
    !!
    !! @par Notes
    !! This routine utilizes the LAPACK routine DPOTRF.
    subroutine cholesky_factor(a, upper, err)
        ! Arguments
        real(dp), intent(inout), dimension(:,:) :: a
        logical, intent(in), optional :: upper
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        real(dp), parameter :: zero = 0.0d0

        ! Local Variables
        character :: uplo
        integer(i32) :: i, n, flag
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg

        ! Initialization
        n = size(a, 1)
        if (present(upper)) then
            if (upper) then
                uplo = 'U'
            else
                uplo = 'L'
            end if
        else
            uplo = 'U'
        end if
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        if (size(a, 2) /= n) then
            ! ERROR: A must be square
            call errmgr%report_error("cholesky_factor", &
                "The input matrix must be square.", LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Process
        call DPOTRF(uplo, n, a, n, flag)
        if (flag > 0) then
            ! ERROR: Matrix is not positive definite
            write(errmsg, '(AI0A)') "The leading minor of order ", flag, &
                " is not positive definite."
            call errmgr%report_error("cholesky_factor", trim(errmsg), &
                LA_MATRIX_FORMAT_ERROR)
        end if

        ! Zero out the non-used upper or lower diagonal
        if (uplo == 'U') then
            ! Zero out the lower
            do i = 1, n - 1
                a(i+1:n,i) = zero
            end do
        else
            ! Zero out the upper
            do i = 2, n
                a(1:i-1,i) = zero
            end do
        end if
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Solves a system of Cholesky factored equations.
    !!
    !! @param[in] upper Set to true if the original matrix A was factored such
    !!  that A = U**T * U; else, set to false if the factorization of A was
    !!  A = L**T * L.
    !! @param[in] a The N-by-N Cholesky factored matrix.
    !! @param[in,out] b On input, the N-by-NRHS right-hand-side matrix B.  On
    !!  output, the solution matrix X.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - LA_ARRAY_SIZE_ERROR: Occurs if any of the input array sizes are 
    !!      incorrect.
    !!
    !! @par Notes
    !! This routine utilizes the LAPACK routine DPOTRS.
    subroutine solve_cholesky_mtx(upper, a, b, err)
        ! Arguments
        logical, intent(in) :: upper
        real(dp), intent(in), dimension(:,:) :: a
        real(dp), intent(inout), dimension(:,:) :: b
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        character :: uplo
        integer(i32) :: n, nrhs, flag
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg

        ! Initialization
        n = size(a, 1)
        nrhs = size(b, 2)
        if (upper) then
            uplo = 'U'
        else
            uplo = 'L'
        end if
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        flag = 0
        if (size(a, 2) /= n) then
            flag = 2
        else if (size(b, 1) /= n) then
            flag = 3
        end if
        if (flag /= 0) then
            write(errmsg, '(AI0A)') "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("solve_cholesky_mtx", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Process
        call DPOTRS(uplo, n, nrhs, a, n, b, n, flag)
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Solves a system of Cholesky factored equations.
    !!
    !! @param[in] upper Set to true if the original matrix A was factored such
    !!  that A = U**T * U; else, set to false if the factorization of A was
    !!  A = L**T * L.
    !! @param[in] a The N-by-N Cholesky factored matrix.
    !! @param[in,out] b On input, the N-element right-hand-side vector B.  On
    !!  output, the solution vector X.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - LA_ARRAY_SIZE_ERROR: Occurs if any of the input array sizes are 
    !!      incorrect.
    !!
    !! @par Notes
    !! This routine utilizes the LAPACK routine DPOTRS.
    subroutine solve_cholesky_vec(upper, a, b, err)
        ! Arguments
        logical, intent(in) :: upper
        real(dp), intent(in), dimension(:,:) :: a
        real(dp), intent(inout), dimension(:) :: b
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        character :: uplo
        integer(i32) :: n, flag
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg

        ! Initialization
        n = size(a, 1)
        if (upper) then
            uplo = 'U'
        else
            uplo = 'L'
        end if
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        flag = 0
        if (size(a, 2) /= n) then
            flag = 2
        else if (size(b) /= n) then
            flag = 3
        end if
        if (flag /= 0) then
            write(errmsg, '(AI0A)') "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("solve_cholesky_vec", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Process
        call DPOTRS(uplo, n, 1, a, n, b, n, flag)
    end subroutine

! ******************************************************************************
! FACTORIZATION UPDATES
! ------------------------------------------------------------------------------
    !> @brief Computes the rank 1 update to a Cholesky factored matrix (upper
    !! triangular).
    !!
    !! @param[in,out] r On input, the N-by-N upper triangular matrix R.  On
    !!  output, the updated matrix R1.
    !! @param[in,out] u On input, the N-element update vector U.  On output,
    !!  the rotation sines used to transform R to R1.
    !! @param[out] work An optional argument that if supplied prevents local
    !!  memory allocation.  If provided, the array must have at least N
    !!  elements.  Additionally, this workspace array is used to contain the
    !!  rotation cosines used to transform R to R1.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - LA_ARRAY_SIZE_ERROR: Occurs if any of the input array sizes are 
    !!      incorrect.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
    !!      there is insufficient memory available.
    !!
    !! @par Notes
    !! This routine utilizes the QRUPDATE routine DCH1UP.
    !!
    !! @par See Also
    !! [Source](https://sourceforge.net/projects/qrupdate/)
    subroutine cholesky_rank1_update(r, u, work, err)
        ! Arguments
        real(dp), intent(inout), dimension(:,:) :: r
        real(dp), intent(inout), dimension(:) :: u
        real(dp), intent(out), pointer, optional, dimension(:) :: work
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        integer(i32) :: n, lwork, istat, flag
        real(dp), pointer, dimension(:) :: wptr
        real(dp), allocatable, target, dimension(:) :: wrk
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg

        ! Initialization
        n = size(r, 1)
        lwork = n
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        flag = 0
        if (size(r, 2) /= n) then
            flag = 1
        else if (size(u) /= n) then
            flag = 2
        end if
        if (flag /= 0) then
            write(errmsg, '(AI0A)') "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("cholesky_rank1_update", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (.not.associated(work)) then
                allocate(wrk(lwork), stat = istat)
                if (istat /= 0) then
                    call errmgr%report_error("cholesky_rank1_update", &
                        "Insufficient memory available.", &
                        LA_OUT_OF_MEMORY_ERROR)
                    return
                end if
                wptr => wrk
            else
                if (size(work) < lwork) then
                    call errmgr%report_error("cholesky_rank1_update", &
                        "The workspace array is too short.", &
                        LA_ARRAY_SIZE_ERROR)
                    return
                end if
                wptr => work(1:lwork)
            end if
        else
            allocate(wrk(lwork), stat = istat)
            if (istat /= 0) then
                call errmgr%report_error("cholesky_rank1_update", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if

        ! Process
        call DCH1UP(n, r, n, u, wptr)

        ! End
        if (allocated(wrk)) deallocate(wrk)
    end subroutine

! ------------------------------------------------------------------------------
end module
