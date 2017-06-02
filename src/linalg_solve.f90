! linalg_solve.f90

!> @brief \b linalg_solve
!!
!! @par Purpose
!! Provides a set of routines for solving systems of linear equations.
module linalg_solve
    use ferror, only : errors
    use lapack
    use linalg_constants
    use linalg_factor, only : rz_factor, mult_rz, mult_qr
    use linalg_core, only : solve_triangular_system, mtx_mult
    implicit none
    private
    public :: solve_lu
    public :: solve_qr
    public :: solve_cholesky
    public :: mtx_inverse
    public :: mtx_pinverse
    public :: least_squares_solve
    public :: least_squares_solve_full
    public :: least_squares_solve_svd

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
    !> @brief Solves a system of M QR-factored equations of N unknowns.
    !!
    !! @par See Also
    !! - [Wikipedia](https://en.wikipedia.org/wiki/QR_decomposition)
    !! - [LAPACK Users Manual](http://netlib.org/lapack/lug/node39.html)
    interface solve_qr
        module procedure :: solve_qr_no_pivot_mtx
        module procedure :: solve_qr_no_pivot_vec
        module procedure :: solve_qr_pivot_mtx
        module procedure :: solve_qr_pivot_vec
    end interface

! ------------------------------------------------------------------------------
    !> @brief Solves a system of Cholesky factored equations.
    interface solve_cholesky
        module procedure :: solve_cholesky_mtx
        module procedure :: solve_cholesky_vec
    end interface

! ------------------------------------------------------------------------------
    !> @brief Solves the overdetermined or underdetermined system (A*X = B) of
    !! M equations of N unknowns.
    interface least_squares_solve
        module procedure :: least_squares_solve_mtx
        module procedure :: least_squares_solve_vec
        module procedure :: least_squares_solve_mtx_1
    end interface

! ------------------------------------------------------------------------------
    !> @brief Solves the overdetermined or underdetermined system (A*X = B) of
    !! M equations of N unknowns, but uses a full orthogonal factorization of
    !! the system.
    interface least_squares_solve_full
        module procedure :: least_squares_solve_mtx_pvt
        module procedure :: least_squares_solve_vec_pvt
    end interface

! ------------------------------------------------------------------------------
    !> @brief Solves the overdetermined or underdetermined system (A*X = B) of
    !! M equations of N unknowns using a singular value decomposition of
    !! matrix A.
    interface least_squares_solve_svd
        module procedure :: least_squares_solve_mtx_svd
        module procedure :: least_squares_solve_vec_svd
    end interface


contains
! ******************************************************************************
! LU SOLUTION
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
                " is not sized correctly."
            call errmgr%report_error("solve_lu_mtx", trim(errmsg), &
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
                " is not sized correctly."
            call errmgr%report_error("solve_lu_vec", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Call DGETRS
        call DGETRS("N", n, 1, a, n, ipvt, b, n, flag)
    end subroutine

! ******************************************************************************
! QR SOLUTION
! ------------------------------------------------------------------------------
    !> @brief Solves a system of M QR-factored equations of N unknowns where
    !! M >= N.
    !!
    !! @param[in] a On input, the M-by-N QR factored matrix as returned by
    !!  @ref qr_factor.  On output, the contents of this matrix are restored.
    !!  Notice, M must be greater than or equal to N.
    !! @param[in] tau A MIN(M, N)-element array containing the scalar factors of
    !!  the elementary reflectors as returned by @ref qr_factor.
    !! @param[in] b On input, the M-by-NRHS right-hand-side matrix.  On output,
    !!  the first N columns are overwritten by the solution matrix X.
    !! @param[out] work An optional input, that if provided, prevents any local
    !!  memory allocation.  If not provided, the memory required is allocated
    !!  within.  If provided, the length of the array must be at least
    !!  @p olwork.
    !! @param[out] olwork An optional output used to determine workspace size.
    !!  If supplied, the routine determines the optimal size for @p work, and
    !!  returns without performing any actual calculations.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - LA_ARRAY_SIZE_ERROR: Occurs if any of the input arrays are not sized
    !!      appropriately.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
    !!      there is insufficient memory available.
    !!
    !! @par Notes
    !! This routine is based upon a subset of the LAPACK routine DGELS.
    subroutine solve_qr_no_pivot_mtx(a, tau, b, work, olwork, err)
        ! Arguments
        real(dp), intent(inout), dimension(:,:) :: a, b
        real(dp), intent(in), dimension(:) :: tau
        real(dp), intent(out), target, optional, dimension(:) :: work
        integer(i32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        real(dp), parameter :: one = 1.0d0

        ! Local Variables
        integer(i32) :: m, n, nrhs, k, lwork, flag, istat
        real(dp), pointer, dimension(:) :: wptr
        real(dp), allocatable, target, dimension(:) :: wrk
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg

        ! Initialization
        m = size(a, 1)
        n = size(a, 2)
        nrhs = size(b, 2)
        k = min(m, n)
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        flag = 0
        if (m < n) then
            flag = 1
        else if (size(tau) /= k) then
            flag = 2
        else if (size(b, 1) /= m) then
            flag = 3
        end if
        if (flag /= 0) then
            ! ERROR: One of the input arrays is not sized correctly
            write(errmsg, '(AI0A)') "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("solve_qr_no_pivot_mtx", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Workspace Query
        call mult_qr(.true., .true., a, tau, b, olwork = lwork)
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (size(work) < lwork) then
                ! ERROR: WORK not sized correctly
                call errmgr%report_error("solve_qr_no_pivot_mtx", &
                    "Incorrectly sized input array WORK, argument 4.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lwork)
        else
            allocate(wrk(lwork), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("solve_qr_no_pivot_mtx", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if

        ! Compute Q**T * B, and store in B
        call mult_qr(.true., .true., a, tau, b, wptr)

        ! Solve the triangular system: A(1:N,1:N)*X = B(1:N,:)
        call solve_triangular_system(.true., .true., .false., .true., one, &
            a(1:n,1:n), b(1:n,:))
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Solves a system of M QR-factored equations of N unknowns where
    !! M >= N.
    !!
    !! @param[in] a On input, the M-by-N QR factored matrix as returned by
    !!  @ref qr_factor.  On output, the contents of this matrix are restored.
    !!  Notice, M must be greater than or equal to N.
    !! @param[in] tau A MIN(M, N)-element array containing the scalar factors of
    !!  the elementary reflectors as returned by @ref qr_factor.
    !! @param[in] b On input, the M-element right-hand-side vector.  On output,
    !!  the first N elements are overwritten by the solution vector X.
    !! @param[out] work An optional input, that if provided, prevents any local
    !!  memory allocation.  If not provided, the memory required is allocated
    !!  within.  If provided, the length of the array must be at least
    !!  @p olwork.
    !! @param[out] olwork An optional output used to determine workspace size.
    !!  If supplied, the routine determines the optimal size for @p work, and
    !!  returns without performing any actual calculations.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - LA_ARRAY_SIZE_ERROR: Occurs if any of the input arrays are not sized
    !!      appropriately.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
    !!      there is insufficient memory available.
    !!
    !! @par Notes
    !! This routine is based upon a subset of the LAPACK routine DGELS.
    subroutine solve_qr_no_pivot_vec(a, tau, b, work, olwork, err)
        ! Arguments
        real(dp), intent(inout), dimension(:,:) :: a
        real(dp), intent(in), dimension(:) :: tau
        real(dp), intent(inout), dimension(:) :: b
        real(dp), intent(out), target, optional, dimension(:) :: work
        integer(i32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        integer(i32) :: m, n, k, flag, lwork, istat
        real(dp), pointer, dimension(:) :: wptr
        real(dp), allocatable, target, dimension(:) :: wrk
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg

        ! Initialization
        m = size(a, 1)
        n = size(a, 2)
        k = min(m, n)
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        flag = 0
        if (m < n) then
            flag = 1
        else if (size(tau) /= k) then
            flag = 2
        else if (size(b) /= m) then
            flag = 3
        end if
        if (flag /= 0) then
            ! ERROR: One of the input arrays is not sized correctly
            write(errmsg, '(AI0A)') "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("solve_qr_no_pivot_vec", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Workspace Query
        call mult_qr(.true., a, tau, b, olwork = lwork)
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (size(work) < lwork) then
                ! ERROR: WORK not sized correctly
                call errmgr%report_error("solve_qr_no_pivot_vec", &
                    "Incorrectly sized input array WORK, argument 4.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lwork)
        else
            allocate(wrk(lwork), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("solve_qr_no_pivot_vec", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if

        ! Compute Q**T * B, and store in B
        call mult_qr(.true., a, tau, b, work = wptr)

        ! Solve the triangular system: A(1:N,1:N)*X = B(1:N)
        call solve_triangular_system(.true., .false., .true., a(1:n,1:n), b)
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Solves a system of M QR-factored equations of N unknowns where the
    !! QR factorization made use of column pivoting.
    !!
    !! @param[in] a On input, the M-by-N QR factored matrix as returned by
    !!  @ref qr_factor.  On output, the contents of this matrix are altered.
    !!  Notice, M must be greater than or equal to N.
    !! @param[in] tau A MIN(M, N)-element array containing the scalar factors of
    !!  the elementary reflectors as returned by @ref qr_factor.
    !! @param[in] jpvt An N-element array, as output by @ref qr_factor, used to
    !!  track the column pivots.
    !! @param[in] b On input, the MAX(M, N)-by-NRHS matrix where the first M
    !!  rows contain the right-hand-side matrix B.  On output, the first N rows
    !!  are overwritten by the solution matrix X.
    !! @param[out] work An optional input, that if provided, prevents any local
    !!  memory allocation.  If not provided, the memory required is allocated
    !!  within.  If provided, the length of the array must be at least
    !!  @p olwork.
    !! @param[out] olwork An optional output used to determine workspace size.
    !!  If supplied, the routine determines the optimal size for @p work, and
    !!  returns without performing any actual calculations.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - LA_ARRAY_SIZE_ERROR: Occurs if any of the input arrays are not sized
    !!      appropriately.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
    !!      there is insufficient memory available.
    !!
    !! @par Notes
    !! This routine is based upon a subset of the LAPACK routine DGELSY.
    subroutine solve_qr_pivot_mtx(a, tau, jpvt, b, work, olwork, err)
        ! Arguments
        real(dp), intent(inout), dimension(:,:) :: a
        real(dp), intent(in), dimension(:) :: tau
        integer(i32), intent(in), dimension(:) :: jpvt
        real(dp), intent(inout), dimension(:,:) :: b
        real(dp), intent(out), target, optional, dimension(:) :: work
        integer(i32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        integer(i32), parameter :: imin = 2
        integer(i32), parameter :: imax = 1
        real(dp), parameter :: zero = 0.0d0
        real(dp), parameter :: one = 1.0d0

        ! Local Variables
        integer(i32) :: i, j, m, n, mn, nrhs, lwork, ismin, ismax, &
            rnk, maxmn, flag, istat, lwork1, lwork2, lwork3
        real(dp) :: rcond, smax, smin, smaxpr, sminpr, s1, c1, s2, c2
        real(dp), pointer, dimension(:) :: wptr, w, tau2
        real(dp), allocatable, target, dimension(:) :: wrk
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg

        ! Initialization
        m = size(a, 1)
        n = size(a, 2)
        mn = min(m, n)
        maxmn = max(m, n)
        nrhs = size(b, 2)
        ismin = mn + 1
        ismax = 2 * mn + 1
        rcond = epsilon(rcond)
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        flag = 0
        if (size(tau) /= mn) then
            flag = 2
        else if (size(jpvt) /= n) then
            flag = 3
        else if (size(b, 1) /= maxmn) then
            flag = 4
        end if
        if (flag /= 0) then
            ! ERROR: One of the input arrays is not sized correctly
            write(errmsg, '(AI0A)') "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("solve_qr_pivot_mtx", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Workspace Query
        call rz_factor(a(1:mn,:), a(1:mn,1), olwork = lwork1)
        call mult_qr(.true., .true., a, tau, b(1:m,:), olwork = lwork2)
        call mult_rz(.true., .true., n, a(1:mn,:), a(1:mn,1), b, &
            olwork = lwork3)
        lwork = max(lwork1, lwork2, lwork3, 2 * mn + 1) + mn
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (size(work) < lwork) then
                ! ERROR: WORK not sized correctly
                call errmgr%report_error("solve_qr_no_pivot_mtx", &
                    "Incorrectly sized input array WORK, argument 5.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lwork)
        else
            allocate(wrk(lwork), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("solve_qr_pivot_mtx", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if

        ! Determine the rank of R11 using an incremental condition estimation
        wptr(ismin) = one
        wptr(ismax) = one
        smax = abs(a(1,1))
        smin = smax
        if (abs(a(1,1)) == zero) then
            rnk = 0
            b(maxmn,nrhs) = zero
            return
        else
            rnk = 1
        end if
        do
            if (rnk < mn) then
                i = rnk + 1
                call DLAIC1(imin, rnk, wptr(ismin:ismin+rnk-1), smin, &
                    a(1:rnk-1,i), a(i,i), sminpr, s1, c1)
                call DLAIC1(imax, rnk, wptr(ismax:ismax+rnk-1), smax, &
                    a(1:rnk-1,i), a(i,i), smaxpr, s2, c2)
                if (smaxpr * rcond <= sminpr) then
                    do i = 1, rnk
                        wptr(ismin+i-1) = s1 * wptr(ismin+i-1)
                        wptr(ismax+i-1) = s2 * wptr(ismax+i-1)
                    end do
                    wptr(ismin+rnk) = c1
                    wptr(ismax+rnk) = c2
                    smin = sminpr
                    smax = smaxpr
                    rnk = rnk + 1
                    cycle
                end if
            end if
            exit
        end do

        ! Partition R = [R11 R12]
        !               [ 0  R22]
        tau2 => wptr(1:rnk)
        w => wptr(rnk+1:lwork)
        if (rnk < n) call rz_factor(a(1:rnk,:), tau2, w)

        ! Compute B(1:m,1:NRHS) = Q**T * B(1:M,1:NRHS)
        call mult_qr(.true., .true., a, tau, b(1:m,:), w)

        ! Solve the triangular system T11 * B(1:rnk,1:nrhs) = B(1:rnk,1:nrhs)
        call solve_triangular_system(.true., .true., .false., .true., one, &
            a(1:rnk,1:rnk), b(1:rnk,:))
        if (n > rnk) b(rnk+1:n,:) = zero

        ! Compute B(1:n,1:nrhs) = Y**T * B(1:n,1:nrhs)
        if (rnk < n) then
            call mult_rz(.true., .true., n - rnk, a(1:rnk,:), tau2, b, w)
        end if

        ! Apply the pivoting: B(1:N,1:NRHS) = P * B(1:N,1:NRHS)
        do j = 1, nrhs
            do i = 1, n
                wptr(jpvt(i)) = b(i,j)
            end do
            b(:,j) = wptr(1:n)
        end do
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Solves a system of M QR-factored equations of N unknowns where the
    !! QR factorization made use of column pivoting.
    !!
    !! @param[in] a On input, the M-by-N QR factored matrix as returned by
    !!  @ref qr_factor.  On output, the contents of this matrix are altered.
    !!  Notice, M must be greater than or equal to N.
    !! @param[in] tau A MIN(M, N)-element array containing the scalar factors of
    !!  the elementary reflectors as returned by @ref qr_factor.
    !! @param[in] jpvt An N-element array, as output by @ref qr_factor, used to
    !!  track the column pivots.
    !! @param[in] b On input, the MAX(M, N)-element array where the first M
    !!  elements contain the right-hand-side vector B.  On output, the first N
    !!  elements are overwritten by the solution vector X.
    !! @param[out] work An optional input, that if provided, prevents any local
    !!  memory allocation.  If not provided, the memory required is allocated
    !!  within.  If provided, the length of the array must be at least
    !!  @p olwork.
    !! @param[out] olwork An optional output used to determine workspace size.
    !!  If supplied, the routine determines the optimal size for @p work, and
    !!  returns without performing any actual calculations.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - LA_ARRAY_SIZE_ERROR: Occurs if any of the input arrays are not sized
    !!      appropriately.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
    !!      there is insufficient memory available.
    !!
    !! @par Notes
    !! This routine is based upon a subset of the LAPACK routine DGELSY.
    subroutine solve_qr_pivot_vec(a, tau, jpvt, b, work, olwork, err)
        ! Arguments
        real(dp), intent(inout), dimension(:,:) :: a
        real(dp), intent(in), dimension(:) :: tau
        integer(i32), intent(in), dimension(:) :: jpvt
        real(dp), intent(inout), dimension(:) :: b
        real(dp), intent(out), target, optional, dimension(:) :: work
        integer(i32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        integer(i32), parameter :: imin = 2
        integer(i32), parameter :: imax = 1
        real(dp), parameter :: zero = 0.0d0
        real(dp), parameter :: one = 1.0d0

        ! Local Variables
        integer(i32) :: i, m, n, mn, lwork, ismin, ismax, rnk, maxmn, flag, &
            istat, lwork1, lwork2
        real(dp) :: rcond, smax, smin, smaxpr, sminpr, s1, c1, s2, c2
        real(dp), pointer, dimension(:) :: wptr, w, tau2
        real(dp), allocatable, target, dimension(:) :: wrk
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg

        ! Initialization
        m = size(a, 1)
        n = size(a, 2)
        mn = min(m, n)
        maxmn = max(m, n)
        ismin = mn + 1
        ismax = 2 * mn + 1
        rcond = epsilon(rcond)
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        flag = 0
        if (size(tau) /= mn) then
            flag = 2
        else if (size(jpvt) /= n) then
            flag = 3
        else if (size(b) /= maxmn) then
            flag = 4
        end if
        if (flag /= 0) then
            ! ERROR: One of the input arrays is not sized correctly
            write(errmsg, '(AI0A)') "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("solve_qr_pivot_vec", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Workspace Query
        call rz_factor(a(1:mn,:), a(1:mn,1), olwork = lwork1)
        call mult_rz(.true., n, a(1:mn,:), a(1:mn,1), b(1:n), olwork = lwork2)
        lwork = max(lwork1, lwork2, 2 * mn + 1) + mn
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (size(work) < lwork) then
                ! ERROR: WORK not sized correctly
                call errmgr%report_error("solve_qr_no_pivot_mtx", &
                    "Incorrectly sized input array WORK, argument 5.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lwork)
        else
            allocate(wrk(lwork), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("solve_qr_pivot_vec", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if

        ! Determine the rank of R11 using an incremental condition estimation
        wptr(ismin) = one
        wptr(ismax) = one
        smax = abs(a(1,1))
        smin = smax
        if (abs(a(1,1)) == zero) then
            rnk = 0
            b(maxmn) = zero
            return
        else
            rnk = 1
        end if
        do
            if (rnk < mn) then
                i = rnk + 1
                call DLAIC1(imin, rnk, wptr(ismin:ismin+rnk-1), smin, &
                    a(1:rnk-1,i), a(i,i), sminpr, s1, c1)
                call DLAIC1(imax, rnk, wptr(ismax:ismax+rnk-1), smax, &
                    a(1:rnk-1,i), a(i,i), smaxpr, s2, c2)
                if (smaxpr * rcond <= sminpr) then
                    do i = 1, rnk
                        wptr(ismin+i-1) = s1 * wptr(ismin+i-1)
                        wptr(ismax+i-1) = s2 * wptr(ismax+i-1)
                    end do
                    wptr(ismin+rnk) = c1
                    wptr(ismax+rnk) = c2
                    smin = sminpr
                    smax = smaxpr
                    rnk = rnk + 1
                    cycle
                end if
            end if
            exit
        end do

        ! Partition R = [R11 R12]
        !               [ 0  R22]
        tau2 => wptr(1:rnk)
        w => wptr(rnk+1:lwork)
        if (rnk < n) call rz_factor(a(1:rnk,:), tau2, w)

        ! Compute B(1:m,1:NRHS) = Q**T * B(1:M,1:NRHS)
        call mult_qr(.true., a, tau, b(1:m))

        ! Solve the triangular system T11 * B(1:rnk) = B(1:rnk)
        call solve_triangular_system(.true., .false., .true., a(1:rnk,1:rnk), &
            b(1:rnk))
        if (n > rnk) b(rnk+1:n) = zero

        ! Compute B(1:n) = Y**T * B(1:n)
        if (rnk < n) then
            call mult_rz(.true., n - rnk, a(1:rnk,:), tau2, b(1:n), w)
        end if

        ! Apply the pivoting: B(1:N) = P * B(1:N)
        do i = 1, n
            wptr(jpvt(i)) = b(i)
        end do
        b = wptr(1:n)
    end subroutine

! ******************************************************************************
! CHOLESKY SOLVE
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
! MATRIX INVERSION ROUTINES
! ------------------------------------------------------------------------------
    !> @brief Computes the inverse of a square matrix.
    !!
    !! @param[in,out] a On input, the N-by-N matrix to invert.  On output, the
    !!  inverted matrix.
    !! @param[out] iwork An optional N-element integer workspace array.
    !! @param[out] work An optional input, that if provided, prevents any local
    !!  memory allocation.  If not provided, the memory required is allocated
    !!  within.  If provided, the length of the array must be at least
    !!  @p olwork.
    !! @param[out] olwork An optional output used to determine workspace size.
    !!  If supplied, the routine determines the optimal size for @p work, and
    !!  returns without performing any actual calculations.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - LA_ARRAY_SIZE_ERROR: Occurs if @p a is not square.  Will also occur if
    !!      incorrectly sized workspace arrays are provided.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
    !!      there is insufficient memory available.
    !!  - LA_SINGULAR_MATRIX_ERROR: Occurs if the input matrix is singular.
    !!
    !! @par Usage
    !! @code {.f90}
    !! ! The following example illustrates how to solve a system of linear
    !! ! equations by matrix inversion.  Notice, this is not a preferred
    !! ! solution technique (use LU factorization instead), but is merely a
    !! ! means of illustrating how to compute the inverse of a square matrix.
    !!
    !! ! Variables
    !! real(dp), dimension(n, n) :: a
    !! real(dp), dimension(n, nrhs) :: b, x
    !!
    !! ! Initialize A and B...
    !!
    !! ! Compute the inverse of A.  The inverse will overwrite the original
    !! ! matrix.
    !! call mtx_inverse(a)
    !!
    !! ! Solve A*X = B as X = inv(A) * B.
    !! x = matmul(a, b)
    !! @endcode
    !!
    !! @par Notes
    !! This routine utilizes the LAPACK routines DGETRF to perform an LU
    !! factorization of the matrix, and DGETRI to invert the LU factored
    !! matrix.
    !!
    !! @par See Also
    !! - [Wikipedia](https://en.wikipedia.org/wiki/Invertible_matrix)
    !! - [Wolfram MathWorld](http://mathworld.wolfram.com/MatrixInverse.html)
    subroutine mtx_inverse(a, iwork, work, olwork, err)
        ! Arguments
        real(dp), intent(inout), dimension(:,:) :: a
        integer(i32), intent(out), target, optional, dimension(:) :: iwork
        real(dp), intent(out), target, optional, dimension(:) :: work
        integer(i32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        integer(i32) :: n, liwork, lwork, istat, flag
        integer(i32), pointer, dimension(:) :: iptr
        integer(i32), allocatable, target, dimension(:) :: iwrk
        real(dp), pointer, dimension(:) :: wptr
        real(dp), allocatable, target, dimension(:) :: wrk
        real(dp), dimension(1) :: temp
        class(errors), pointer :: errmgr
        type(errors), target :: deferr

        ! Initialization
        n = size(a, 1)
        liwork = n
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        if (size(a, 2) /= n) then
            call errmgr%report_error("mtx_inverse", &
                "The matrix must be squre to invert.", LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Workspace Query
        call DGETRI(n, a, n, istat, temp, -1, flag)
        lwork = int(temp(1), i32)
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Workspace Allocation
        if (present(work)) then
            if (size(work) < lwork) then
                ! ERROR: WORK not sized correctly
                call errmgr%report_error("svd", &
                    "Incorrectly sized input array WORK, argument 3.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lwork)
        else
            allocate(wrk(lwork), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("mtx_inverse", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if

        ! Integer Workspace Allocation
        if (present(iwork)) then
            if (size(iwork) < liwork) then
                ! ERROR: IWORK not sized correctly
                call errmgr%report_error("svd", &
                    "Incorrectly sized input array IWORK, argument 2.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            iptr => iwork(1:liwork)
        else
            allocate(iwrk(liwork), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("mtx_inverse", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            iptr => iwrk
        end if

        ! Compute the LU factorization of A
        call DGETRF(n, n, a, n, iptr, flag)

        ! Compute the inverse of the LU factored matrix
        call DGETRI(n, a, n, iptr, wptr, lwork, flag)

        ! Check for a singular matrix
        if (flag > 0) then
            call errmgr%report_error("mtx_inverse", &
                "The matrix is singular; therefore, the inverse could " // &
                "not be computed.", LA_SINGULAR_MATRIX_ERROR)
        end if
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Computes the Moore-Penrose pseudo-inverse of a M-by-N matrix
    !! using the singular value decomposition of the matrix.
    !!
    !! @param[in,out] a On input, the M-by-N matrix to invert.  The matrix is
    !!  overwritten on output.
    !! @param[out] ainv The N-by-M matrix where the pseudo-inverse of @p a
    !!  will be written.
    !! @param[in] tol An optional input, that if supplied, overrides the default
    !!  tolerance on singular values such that singular values less than this
    !!  tolerance are forced to have a reciprocal of zero, as opposed to 1/S(I).
    !!  The default tolerance is: MAX(M, N) * EPS * MAX(S).  If the supplied
    !!  value is less than a value that causes an overflow, the tolerance
    !!  reverts back to its default value, and the operation continues;
    !!  however, a warning message is issued.
    !! @param[out] work An optional input, that if provided, prevents any local
    !!  memory allocation.  If not provided, the memory required is allocated
    !!  within.  If provided, the length of the array must be at least
    !!  @p olwork.
    !! @param[out] olwork An optional output used to determine workspace size.
    !!  If supplied, the routine determines the optimal size for @p work, and
    !!  returns without performing any actual calculations.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - LA_ARRAY_SIZE_ERROR: Occurs if any of the input arrays are not sized
    !!      appropriately.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
    !!      there is insufficient memory available.
    !!  - LA_CONVERGENCE_ERROR: Occurs as a warning if the QR iteration process
    !!      could not converge to a zero value.
    !!
    !! @par Usage
    !! @code{.f90}
    !! ! Use the pseudo-inverse to obtain a least-squares solution to the
    !! ! overdetermined problem A*X = B, where A is an M-by-N matrix (M >= N),
    !! ! B is an M-by-NRHS matrix, and X is an N-by-NRHS matrix.
    !!
    !! ! Variables
    !! real(dp), dimension(m, n) :: a
    !! real(dp), dimension(n, m) :: ainv
    !! real(dp), dimension(m, nrhs) :: b
    !! real(dp), dimension(n, nrhs) :: x
    !!
    !! ! Initialize A, and B...
    !!
    !! ! Compute the pseudo-inverse of A.  Let the subroutine allocate its
    !! ! own workspace array.
    !! call mtx_pinverse(a, ainv)
    !!
    !! ! Compute X = AINV * B to obtain the solution.
    !! x = matmul(ainv, b)
    !! @endcode
    !!
    !! @par See Also
    !! - [Wikipedia](https://en.wikipedia.org/wiki/Moore%E2%80%93Penrose_pseudoinverse)
    !! - [Wolfram MathWorld](http://mathworld.wolfram.com/Moore-PenroseMatrixInverse.html)
    !! - [MathWorks](http://www.mathworks.com/help/matlab/ref/pinv.html?s_tid=srchtitle)
    subroutine mtx_pinverse(a, ainv, tol, work, olwork, err)
        ! Arguments
        real(dp), intent(inout), dimension(:,:) :: a
        real(dp), intent(out), dimension(:,:) :: ainv
        real(dp), intent(in), optional :: tol
        real(dp), intent(out), target, dimension(:), optional :: work
        integer(i32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        real(dp), parameter :: zero = 0.0d0
        real(dp), parameter :: one = 1.0d0

        ! Local Variables
        integer(i32) :: i, m, n, mn, lwork, istat, flag, i1, i2a, i2b, i3a, &
            i3b, i4
        real(dp), pointer, dimension(:) :: s, wptr, w
        real(dp), pointer, dimension(:,:) :: u, vt
        real(dp), allocatable, target, dimension(:) :: wrk
        real(dp), dimension(1) :: temp
        real(dp) :: t, tref, tolcheck
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg

        ! Initialization
        m = size(a, 1)
        n = size(a, 2)
        mn = min(m, n)
        i1 = m * mn
        i2a = i1 + 1
        i2b = i2a + n * n - 1
        i3a = i2b + 1
        i3b = i3a + mn - 1
        i4 = i3b + 1
        tolcheck = dlamch('s')
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        if (size(ainv, 1) /= n .or. size(ainv, 2) /= m) then
            write(errmsg, '(AI0AI0A)') &
                "The output matrix AINV is not sized appropriately.  " // &
                "It is expected to be ", n, "-by-", m, "."
            call errmgr%report_error("mtx_pinverse", errmsg, &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Workspace Query
        call DGESVD('S', 'A', m, n, a, m, a(1:mn,:), a, m, a, n, temp, -1, flag)
        lwork = int(temp(1), i32)
        lwork = lwork + m * mn + n * n + mn
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (size(work) < lwork) then
                ! ERROR: WORK not sized correctly
                call errmgr%report_error("mtx_pinverse", &
                    "Incorrectly sized input array WORK, argument 4.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lwork)
        else
            allocate(wrk(lwork), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("mtx_pinverse", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if
        u(1:m,1:mn) => wptr(1:i1)
        vt(1:n,1:n) => wptr(i2a:i2b)
        s => wptr(i3a:i3b)
        w => wptr(i4:lwork)

        ! Compute the SVD of A
        call DGESVD('S', 'A', m, n, a, m, s, u, m, vt, n, w, size(w), flag)

        ! Check for convergence
        if (flag > 0) then
            write(errmsg, '(I0A)') flag, " superdiagonals could not " // &
                "converge to zero as part of the QR iteration process."
            call errmgr%report_warning("mtx_pinverse", errmsg, &
                LA_CONVERGENCE_ERROR)
            return
        end if

        ! Determine the threshold tolerance for the singular values such that
        ! singular values less than the threshold result in zero when inverted.
        tref = max(m, n) * epsilon(t) * s(1)
        if (present(tol)) then
            t = tol
        else
            t = tref
        end if
        !if (t < safe_denom(t)) then
        if (t < tolcheck) then
            ! The supplied tolerance is too small, simply fall back to the
            ! default, but issue a warning to the user
            t = tref
            ! call errmgr%report_warning("pinverse_1", "The supplied tolerance was " // &
            !     "smaller than a value that would result in an overflow " // &
            !     "condition, or is negative; therefore, the tolerance has " // &
            !     "been reset to its default value.")
        end if

        ! Compute the pseudoinverse such that pinv(A) = V * inv(S) * U**T by
        ! first computing V * inv(S) (result is N-by-M), and store in the first
        ! MN rows of VT in a transposed manner.
        do i = 1, mn
            ! Apply 1 / S(I) to VT(I,:)
            if (s(i) < t) then
                vt(i,:) = zero
            else
                call DRSCL(n, s(i), vt(i,1:n), 1)
            end if
        end do

        ! Compute (VT**T * inv(S)) * U**T
        call mtx_mult(.true., .true., one, vt(1:mn,:), u, zero, ainv)
    end subroutine

! ******************************************************************************
! LEAST SQUARES SOLUTION ROUTINES
! ------------------------------------------------------------------------------
    !> @brief Solves the overdetermined or underdetermined system (A*X = B) of
    !! M equations of N unknowns using a QR or LQ factorization of the matrix A.
    !! Notice, it is assumed that matrix A has full rank.
    !!
    !! @param[in,out] a On input, the M-by-N matrix A.  On output, if M >= N,
    !!  the QR factorization of A in the form as output by @ref qr_factor; else,
    !!  if M < N, the LQ factorization of A in the form as output by
    !!  @ref lq_factor.
    !! @param[in,out] b If M >= N, the M-by-NRHS matrix B.  On output, the first
    !!  N rows contain the N-by-NRHS solution matrix X.  If M < N, an
    !!  N-by-NRHS matrix with the first M rows containing the matrix B.  On
    !!  output, the N-by-NRHS solution matrix X.
    !! @param[out] work An optional input, that if provided, prevents any local
    !!  memory allocation.  If not provided, the memory required is allocated
    !!  within.  If provided, the length of the array must be at least
    !!  @p olwork.
    !! @param[out] olwork An optional output used to determine workspace size.
    !!  If supplied, the routine determines the optimal size for @p work, and
    !!  returns without performing any actual calculations.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - LA_ARRAY_SIZE_ERROR: Occurs if any of the input arrays are not sized
    !!      appropriately.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
    !!      there is insufficient memory available.
    !!  - LA_INVALID_OPERATION_ERROR: Occurs if @p a is not of full rank.
    !!
    !! @par Notes
    !! This routine utilizes the LAPACK routine DGELS.
    subroutine least_squares_solve_mtx(a, b, work, olwork, err)
        ! Arguments
        real(dp), intent(inout), dimension(:,:) :: a, b
        real(dp), intent(out), target, optional, dimension(:) :: work
        integer(i32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        integer(i32) :: m, n, maxmn, nrhs, lwork, istat, flag
        real(dp), pointer, dimension(:) :: wptr
        real(dp), allocatable, target, dimension(:) :: wrk
        real(dp), dimension(1) :: temp
        class(errors), pointer :: errmgr
        type(errors), target :: deferr

        ! Initialization
        m = size(a, 1)
        n = size(a, 2)
        maxmn = max(m, n)
        nrhs = size(b, 2)
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        if (size(b, 1) /= maxmn) then
            call errmgr%report_error("least_squares_solve_mtx", &
                "Input 2 is not sized correctly.", LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Workspace Query
        call DGELS('N', m, n, nrhs, a, m, b, maxmn, temp, -1, flag)
        lwork = int(temp(1), i32)
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (size(work) < lwork) then
                ! ERROR: WORK not sized correctly
                call errmgr%report_error("least_squares_solve_mtx", &
                    "Incorrectly sized input array WORK, argument 3.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lwork)
        else
            allocate(wrk(lwork), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("least_squares_solve_mtx", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if

        ! Process
        call DGELS('N', m, n, nrhs, a, m, b, maxmn, wptr, lwork, flag)
        if (flag > 0) then
            call errmgr%report_error("least_squares_solve_mtx", &
                "The supplied matrix is not of full rank; therefore, " // &
                "the solution could not be computed via this routine.  " // &
                "Try a routine that utilizes column pivoting.", &
                LA_INVALID_OPERATION_ERROR)
        end if
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Solves the overdetermined or underdetermined system (A*X = B) of
    !! M equations of N unknowns using a QR or LQ factorization of the matrix A.
    !! Notice, it is assumed that matrix A has full rank.
    !!
    !! @param[in,out] a On input, the M-by-N matrix A.  On output, if M >= N,
    !!  the QR factorization of A in the form as output by @ref qr_factor; else,
    !!  if M < N, the LQ factorization of A in the form as output by
    !!  @ref lq_factor.
    !! @param[in,out] b If M >= N, the M-element array B.  On output, the first
    !!  N elements contain the N-element solution array X.  If M < N, an
    !!  N-element array with the first M elements containing the array B.  On
    !!  output, the N-element solution array X.
    !! @param[out] work An optional input, that if provided, prevents any local
    !!  memory allocation.  If not provided, the memory required is allocated
    !!  within.  If provided, the length of the array must be at least
    !!  @p olwork.
    !! @param[out] olwork An optional output used to determine workspace size.
    !!  If supplied, the routine determines the optimal size for @p work, and
    !!  returns without performing any actual calculations.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - LA_ARRAY_SIZE_ERROR: Occurs if any of the input arrays are not sized
    !!      appropriately.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
    !!      there is insufficient memory available.
    !!  - LA_INVALID_OPERATION_ERROR: Occurs if @p a is not of full rank.
    !!
    !! @par Notes
    !! This routine utilizes the LAPACK routine DGELS.
    subroutine least_squares_solve_vec(a, b, work, olwork, err)
        ! Arguments
        real(dp), intent(inout), dimension(:,:) :: a
        real(dp), intent(inout), dimension(:) :: b
        real(dp), intent(out), target, optional, dimension(:) :: work
        integer(i32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        integer(i32) :: m, n, maxmn, lwork, istat, flag
        real(dp), pointer, dimension(:) :: wptr
        real(dp), allocatable, target, dimension(:) :: wrk
        real(dp), dimension(1) :: temp
        class(errors), pointer :: errmgr
        type(errors), target :: deferr

        ! Initialization
        m = size(a, 1)
        n = size(a, 2)
        maxmn = max(m, n)
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        if (size(b) /= maxmn) then
            call errmgr%report_error("least_squares_solve_vec", &
                "Input 2 is not sized correctly.", LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Workspace Query
        call DGELS('N', m, n, 1, a, m, b, maxmn, temp, -1, flag)
        lwork = int(temp(1), i32)
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (size(work) < lwork) then
                ! ERROR: WORK not sized correctly
                call errmgr%report_error("least_squares_solve_vec", &
                    "Incorrectly sized input array WORK, argument 3.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lwork)
        else
            allocate(wrk(lwork), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("least_squares_solve_vec", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if

        ! Process
        call DGELS('N', m, n, 1, a, m, b, maxmn, wptr, lwork, flag)
        if (flag > 0) then
            call errmgr%report_error("least_squares_solve_mtx", &
                "The supplied matrix is not of full rank; therefore, " // &
                "the solution could not be computed via this routine.  " // &
                "Try a routine that utilizes column pivoting.", &
                LA_INVALID_OPERATION_ERROR)
        end if
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Solves the overdetermined or underdetermined system (A*X = B) of
    !! M equations of N unknowns using a QR or LQ factorization of the matrix A.
    !! Notice, it is assumed that matrix A has full rank.
    !!
    !! @param[in,out] a On input, the M-by-N matrix A.  On output, if M >= N,
    !!  the QR factorization of A in the form as output by @ref qr_factor; else,
    !!  if M < N, the LQ factorization of A in the form as output by
    !!  @ref lq_factor.
    !! @param[in,out] b On input, the M-by-NRHS matrix B.  On output the
    !!  contents are overwritten.
    !! @param[out] x The N-by-NRHS solution matrix X.
    !! @param[out] work An optional input, that if provided, prevents any local
    !!  memory allocation.  If not provided, the memory required is allocated
    !!  within.  If provided, the length of the array must be at least
    !!  @p olwork.
    !! @param[out] olwork An optional output used to determine workspace size.
    !!  If supplied, the routine determines the optimal size for @p work, and
    !!  returns without performing any actual calculations.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - LA_ARRAY_SIZE_ERROR: Occurs if any of the input arrays are not sized
    !!      appropriately.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
    !!      there is insufficient memory available.
    !!  - LA_INVALID_OPERATION_ERROR: Occurs if @p a is not of full rank.
    !!
    !! @par Notes
    !! This routine utilizes the LAPACK routine DGELS.
    subroutine least_squares_solve_mtx_1(a, b, x, work, olwork, err)
        ! Arguments
        real(dp), intent(inout), dimension(:,:) :: a, b
        real(dp), intent(out), dimension(:,:) :: x
        real(dp), intent(out), target, optional, dimension(:) :: work
        integer(i32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        integer(i32) :: m, n
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg

        ! Initialization
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Workspace Query
        if (present(olwork)) then
            call least_squares_solve(a, b, work, olwork)
            return
        end if

        ! Process
        m = size(a, 1)
        n = size(a, 2)
        if (size(b, 2) /= size(x, 2) .or. size(x, 1) /= n) then
            if (size(b, 2) /= size(x, 2)) then
                write(errmsg, '(AI0AI0A)') &
                    "The number of columns in matrices B and X must " // &
                    "match.  B has ", size(b, 2), " columns, and X has ", &
                    size(x, 2), " columns."
            else
                write(errmsg, '(AI0AI0A)') &
                    "The matrix X was expected to have ", n, &
                    " rows, but was found to have ", size(x, 1), " rows."
            end if

            call errmgr%report_error("least_squares_solve_mtx_1", &
                trim(errmsg), LA_ARRAY_SIZE_ERROR)
            return
        end if
        if (m >= n) then
            call least_squares_solve(a, b, work, err = errmgr)
            x = b(1:n,:)
        else
            x(1:m,:) = b
            call least_squares_solve(a, x, work, err = errmgr)
        end if
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Solves the overdetermined or underdetermined system (A*X = B) of
    !! M equations of N unknowns using a complete orthogonal factorization of
    !! matrix A.
    !!
    !! @param[in,out] a On input, the M-by-N matrix A.  On output, the matrix
    !!  is overwritten by the details of its complete orthogonal factorization.
    !! @param[in,out] b If M >= N, the M-by-NRHS matrix B.  On output, the first
    !!  N rows contain the N-by-NRHS solution matrix X.  If M < N, an
    !!  N-by-NRHS matrix with the first M rows containing the matrix B.  On
    !!  output, the N-by-NRHS solution matrix X.

    !! @param[out] ipvt On input, an N-element array that if IPVT(I) .ne. 0,
    !!  the I-th column of A is permuted to the front of A * P; if IPVT(I) = 0,
    !!  the I-th column of A is a free column.  On output, if IPVT(I) = K, then
    !!  the I-th column of A * P was the K-th column of A.

    !! @param[out] arnk An optional output, that if provided, will return the
    !!  rank of @p a.
    !! @param[out] work An optional input, that if provided, prevents any local
    !!  memory allocation.  If not provided, the memory required is allocated
    !!  within.  If provided, the length of the array must be at least
    !!  @p olwork.
    !! @param[out] olwork An optional output used to determine workspace size.
    !!  If supplied, the routine determines the optimal size for @p work, and
    !!  returns without performing any actual calculations.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - LA_ARRAY_SIZE_ERROR: Occurs if any of the input arrays are not sized
    !!      appropriately.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
    !!      there is insufficient memory available.
    !!
    !! @par Notes
    !! This routine utilizes the LAPACK routine DGELSY.
    subroutine least_squares_solve_mtx_pvt(a, b, ipvt, arnk, work, olwork, err)
        ! Arguments
        real(dp), intent(inout), dimension(:,:) :: a, b
        integer(i32), intent(inout), target, optional, dimension(:) :: ipvt
        integer(i32), intent(out), optional :: arnk
        real(dp), intent(out), target, optional, dimension(:) :: work
        integer(i32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        integer(i32) :: m, n, maxmn, nrhs, lwork, istat, flag, rnk
        real(dp), pointer, dimension(:) :: wptr
        real(dp), allocatable, target, dimension(:) :: wrk
        integer(i32), allocatable, target, dimension(:) :: iwrk
        integer(i32), pointer, dimension(:) :: iptr
        real(dp), dimension(1) :: temp
        integer(i32), dimension(1) :: itemp
        real(dp) :: rc
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg

        ! Initialization
        m = size(a, 1)
        n = size(a, 2)
        maxmn = max(m, n)
        nrhs = size(b, 2)
        rc = epsilon(rc)
        if (present(arnk)) arnk = 0
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        flag = 0
        if (size(b, 1) /= maxmn) then
            flag = 2
        else if (size(ipvt) /= n) then
            flag = 3
        end if
        if (flag /= 0) then
            write(errmsg, '(AI0A)') "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("least_squares_solve_mtx_pvt", &
                trim(errmsg), LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Workspace Query
        call DGELSY(m, n, nrhs, a, m, b, maxmn, itemp, rc, rnk, temp, -1, flag)
        lwork = int(temp(1), i32)
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Local Memory Allocation
        if (present(ipvt)) then
            if (size(ipvt) < n) then
                ! ERROR: IPVT is not big enough
                call errmgr%report_error("least_squares_solve_mtx_pvt", &
                    "Incorrectly sized pivot array, argument 3.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            iptr => ipvt(1:n)
        else
            allocate(iwrk(n), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("least_squares_solve_mtx_pvt", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            iptr => iwrk
            iptr = 0
        end if

        if (present(work)) then
            if (size(work) < lwork) then
                ! ERROR: WORK not sized correctly
                call errmgr%report_error("least_squares_solve_mtx_pvt", &
                    "Incorrectly sized input array WORK, argument 5.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lwork)
        else
            allocate(wrk(lwork), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("least_squares_solve_mtx_pvt", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if

        ! Process
        call DGELSY(m, n, nrhs, a, m, b, maxmn, iptr, rc, rnk, wptr, lwork, &
            flag)
        if (present(arnk)) arnk = rnk
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Solves the overdetermined or underdetermined system (A*X = B) of
    !! M equations of N unknowns using a complete orthogonal factorization of
    !! matrix A.
    !!
    !! @param[in,out] a On input, the M-by-N matrix A.  On output, the matrix
    !!  is overwritten by the details of its complete orthogonal factorization.
    !! @param[in,out] b If M >= N, the M-element array B.  On output, the first
    !!  N elements contain the N-element solution array X.  If M < N, an
    !!  N-element array with the first M elements containing the array B.  On
    !!  output, the N-element solution array X.
    !! @param[out] ipvt On input, an N-element array that if IPVT(I) .ne. 0,
    !!  the I-th column of A is permuted to the front of A * P; if IPVT(I) = 0,
    !!  the I-th column of A is a free column.  On output, if IPVT(I) = K, then
    !!  the I-th column of A * P was the K-th column of A.
    !! @param[out] arnk An optional output, that if provided, will return the
    !!  rank of @p a.
    !! @param[out] work An optional input, that if provided, prevents any local
    !!  memory allocation.  If not provided, the memory required is allocated
    !!  within.  If provided, the length of the array must be at least
    !!  @p olwork.
    !! @param[out] olwork An optional output used to determine workspace size.
    !!  If supplied, the routine determines the optimal size for @p work, and
    !!  returns without performing any actual calculations.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - LA_ARRAY_SIZE_ERROR: Occurs if any of the input arrays are not sized
    !!      appropriately.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
    !!      there is insufficient memory available.
    !!
    !! @par Notes
    !! This routine utilizes the LAPACK routine DGELSY.
    subroutine least_squares_solve_vec_pvt(a, b, ipvt, arnk, work, olwork, err)
        ! Arguments
        real(dp), intent(inout), dimension(:,:) :: a
        real(dp), intent(inout), dimension(:) :: b
        integer(i32), intent(inout), target, optional, dimension(:) :: ipvt
        integer(i32), intent(out), optional :: arnk
        real(dp), intent(out), target, optional, dimension(:) :: work
        integer(i32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        integer(i32) :: m, n, maxmn, lwork, istat, flag, rnk
        real(dp), pointer, dimension(:) :: wptr
        real(dp), allocatable, target, dimension(:) :: wrk
        integer(i32), allocatable, target, dimension(:) :: iwrk
        integer(i32), pointer, dimension(:) :: iptr
        real(dp), dimension(1) :: temp
        integer(i32), dimension(1) :: itemp
        real(dp) :: rc
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg

        ! Initialization
        m = size(a, 1)
        n = size(a, 2)
        maxmn = max(m, n)
        rc = epsilon(rc)
        if (present(arnk)) arnk = 0
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        flag = 0
        if (size(b, 1) /= maxmn) then
            flag = 2
        else if (size(ipvt) /= n) then
            flag = 3
        end if
        if (flag /= 0) then
            write(errmsg, '(AI0A)') "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("least_squares_solve_vec_pvt", &
                trim(errmsg), LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Workspace Query
        call DGELSY(m, n, 1, a, m, b, maxmn, itemp, rc, rnk, temp, -1, flag)
        lwork = int(temp(1), i32)
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Local Memory Allocation
        if (present(ipvt)) then
            if (size(ipvt) < n) then
                ! ERROR: IPVT is not big enough
                call errmgr%report_error("least_squares_solve_mtx_pvt", &
                    "Incorrectly sized pivot array, argument 3.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            iptr => ipvt(1:n)
        else
            allocate(iwrk(n), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("least_squares_solve_mtx_pvt", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            iptr => iwrk
            iptr = 0
        end if

        if (present(work)) then
            if (size(work) < lwork) then
                ! ERROR: WORK not sized correctly
                call errmgr%report_error("least_squares_solve_vec_pvt", &
                    "Incorrectly sized input array WORK, argument 5.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lwork)
        else
            allocate(wrk(lwork), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("least_squares_solve_vec_pvt", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if

        ! Process
        call DGELSY(m, n, 1, a, m, b, maxmn, ipvt, rc, rnk, wptr, lwork, flag)
        if (present(arnk)) arnk = rnk
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Solves the overdetermined or underdetermined system (A*X = B) of
    !! M equations of N unknowns using a singular value decomposition of
    !! matrix A.
    !!
    !! @param[in,out] a On input, the M-by-N matrix A.  On output, the matrix
    !!  is overwritten by the details of its complete orthogonal factorization.
    !! @param[in,out] b If M >= N, the M-by-NRHS matrix B.  On output, the first
    !!  N rows contain the N-by-NRHS solution matrix X.  If M < N, an
    !!  N-by-NRHS matrix with the first M rows containing the matrix B.  On
    !!  output, the N-by-NRHS solution matrix X.
    !! @param[out] arnk An optional output, that if provided, will return the
    !!  rank of @p a.
    !! @param[out] s A MIN(M, N)-element array that on output contains the
    !!  singular values of @p a in descending order.  Notice, the condition
    !!  number of @p a can be determined by S(1) / S(MIN(M, N)).
    !! @param[out] work An optional input, that if provided, prevents any local
    !!  memory allocation.  If not provided, the memory required is allocated
    !!  within.  If provided, the length of the array must be at least
    !!  @p olwork.
    !! @param[out] olwork An optional output used to determine workspace size.
    !!  If supplied, the routine determines the optimal size for @p work, and
    !!  returns without performing any actual calculations.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - LA_ARRAY_SIZE_ERROR: Occurs if any of the input arrays are not sized
    !!      appropriately.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
    !!      there is insufficient memory available.
    !!  - LA_CONVERGENCE_ERROR: Occurs as a warning if the QR iteration process
    !!      could not converge to a zero value.
    !!
    !! @par Notes
    !! This routine utilizes the LAPACK routine DGELSS.
    subroutine least_squares_solve_mtx_svd(a, b, arnk, s, work, olwork, err)
        ! Arguments
        real(dp), intent(inout), dimension(:,:) :: a, b
        real(dp), intent(out), dimension(:) :: s
        integer(i32), intent(out), optional :: arnk
        real(dp), intent(out), target, optional, dimension(:) :: work
        integer(i32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        integer(i32) :: m, n, nrhs, mn, maxmn, istat, flag, lwork, rnk
        real(dp), pointer, dimension(:) :: wptr
        real(dp), allocatable, target, dimension(:) :: wrk
        real(dp), dimension(1) :: temp
        real(dp) :: rcond
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg

        ! Initialization
        m = size(a, 1)
        n = size(a, 2)
        nrhs = size(b, 2)
        mn = min(m, n)
        maxmn = max(m, n)
        rcond = epsilon(rcond)
        if (present(arnk)) arnk = 0
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        flag = 0
        if (size(b, 1) /= maxmn) then
            flag = 2
        else if (size(s) /= mn) then
            flag = 4
        end if
        if (flag /= 0) then
            ! ERROR: One of the input arrays is not sized correctly
            write(errmsg, '(AI0A)') "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("least_squares_solve_mtx_svd", &
                trim(errmsg), LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Workspace Query
        call DGELSS(m, n, nrhs, a, m, b, maxmn, s, rcond, rnk, temp, -1, flag)
        lwork = int(temp(1), i32)
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (size(work) < lwork) then
                ! ERROR: WORK not sized correctly
                call errmgr%report_error("least_squares_solve_mtx_svd", &
                    "Incorrectly sized input array WORK, argument 5.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lwork)
        else
            allocate(wrk(lwork), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("least_squares_solve_mtx_svd", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if

        ! Process
        call DGELSS(m, n, nrhs, a, m, b, maxmn, s, rcond, rnk, wptr, lwork, &
            flag)
        if (present(arnk)) arnk = rnk
        if (flag > 0) then
            write(errmsg, '(I0A)') flag, " superdiagonals could not " // &
                "converge to zero as part of the QR iteration process."
            call errmgr%report_warning("least_squares_solve_mtx_svd", errmsg, &
                LA_CONVERGENCE_ERROR)
        end if
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Solves the overdetermined or underdetermined system (A*X = B) of
    !! M equations of N unknowns using a singular value decomposition of
    !! matrix A.
    !!
    !! @param[in,out] a On input, the M-by-N matrix A.  On output, the matrix
    !!  is overwritten by the details of its complete orthogonal factorization.
    !! @param[in,out] b If M >= N, the M-by-NRHS matrix B.  On output, the first
    !!  N rows contain the N-by-NRHS solution matrix X.  If M < N, an
    !!  N-by-NRHS matrix with the first M rows containing the matrix B.  On
    !!  output, the N-by-NRHS solution matrix X.
    !! @param[out] arnk An optional output, that if provided, will return the
    !!  rank of @p a.
    !! @param[out] s A MIN(M, N)-element array that on output contains the
    !!  singular values of @p a in descending order.  Notice, the condition
    !!  number of @p a can be determined by S(1) / S(MIN(M, N)).
    !! @param[out] work An optional input, that if provided, prevents any local
    !!  memory allocation.  If not provided, the memory required is allocated
    !!  within.  If provided, the length of the array must be at least
    !!  @p olwork.
    !! @param[out] olwork An optional output used to determine workspace size.
    !!  If supplied, the routine determines the optimal size for @p work, and
    !!  returns without performing any actual calculations.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - LA_ARRAY_SIZE_ERROR: Occurs if any of the input arrays are not sized
    !!      appropriately.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
    !!      there is insufficient memory available.
    !!  - LA_CONVERGENCE_ERROR: Occurs as a warning if the QR iteration process
    !!      could not converge to a zero value.
    !!
    !! @par Notes
    !! This routine utilizes the LAPACK routine DGELSS.
    subroutine least_squares_solve_vec_svd(a, b, arnk, s, work, olwork, err)
        ! Arguments
        real(dp), intent(inout), dimension(:,:) :: a
        real(dp), intent(inout), dimension(:) :: b
        real(dp), intent(out), dimension(:) :: s
        integer(i32), intent(out), optional :: arnk
        real(dp), intent(out), target, optional, dimension(:) :: work
        integer(i32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        integer(i32) :: m, n, mn, maxmn, istat, flag, lwork, rnk
        real(dp), pointer, dimension(:) :: wptr
        real(dp), allocatable, target, dimension(:) :: wrk
        real(dp), dimension(1) :: temp
        real(dp) :: rcond
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg

        ! Initialization
        m = size(a, 1)
        n = size(a, 2)
        mn = min(m, n)
        maxmn = max(m, n)
        rcond = epsilon(rcond)
        if (present(arnk)) arnk = 0
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        flag = 0
        if (size(b) /= maxmn) then
            flag = 2
        else if (size(s) /= mn) then
            flag = 4
        end if
        if (flag /= 0) then
            ! ERROR: One of the input arrays is not sized correctly
            write(errmsg, '(AI0A)') "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("least_squares_solve_vec_svd", &
                trim(errmsg), LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Workspace Query
        call DGELSS(m, n, 1, a, m, b, maxmn, s, rcond, rnk, temp, -1, flag)
        lwork = int(temp(1), i32)
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (size(work) < lwork) then
                ! ERROR: WORK not sized correctly
                call errmgr%report_error("least_squares_solve_vec_svd", &
                    "Incorrectly sized input array WORK, argument 5.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lwork)
        else
            allocate(wrk(lwork), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("least_squares_solve_vec_svd", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if

        ! Process
        call DGELSS(m, n, 1, a, m, b, maxmn, s, rcond, rnk, wptr, lwork, flag)
        if (present(arnk)) arnk = rnk
        if (flag > 0) then
            write(errmsg, '(I0A)') flag, " superdiagonals could not " // &
                "converge to zero as part of the QR iteration process."
            call errmgr%report_warning("least_squares_solve_vec_svd", errmsg, &
                LA_CONVERGENCE_ERROR)
        end if
    end subroutine


end module
