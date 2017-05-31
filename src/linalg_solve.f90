! linalg_solve.f90

!> @brief \b linalg_solve
!!
!! @par Purpose
!! Provides a set of routines for solving systems of linear equations.
module linalg_solve
    use ferror, only : errors
    use linalg_constants
    use linalg_factor, only : rz_factor, mult_rz, mult_qr
    use linalg_core, only : solve_triangular_system
    implicit none
    private
    public :: solve_lu
    public :: solve_qr
    public :: solve_cholesky

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
        real(dp), intent(out), pointer, optional, dimension(:) :: work
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
            if (.not.associated(work)) then
                allocate(wrk(lwork), stat = istat)
                if (istat /= 0) then
                    ! ERROR: Out of memory
                    call errmgr%report_error("solve_qr_no_pivot_mtx", &
                        "Insufficient memory available.", &
                        LA_OUT_OF_MEMORY_ERROR)
                    return
                end if
                wptr => wrk
            else
                if (size(work) < lwork) then
                    ! ERROR: WORK not sized correctly
                    call errmgr%report_error("solve_qr_no_pivot_mtx", &
                        "Incorrectly sized input array WORK, argument 4.", &
                        LA_ARRAY_SIZE_ERROR)
                    return
                end if
                wptr => work(1:lwork)
            end if
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
        real(dp), intent(out), pointer, optional, dimension(:) :: work
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
            if (.not.associated(work)) then
                allocate(wrk(lwork), stat = istat)
                if (istat /= 0) then
                    ! ERROR: Out of memory
                    call errmgr%report_error("solve_qr_no_pivot_vec", &
                        "Insufficient memory available.", &
                        LA_OUT_OF_MEMORY_ERROR)
                    return
                end if
                wptr => wrk
            else
                if (size(work) < lwork) then
                    ! ERROR: WORK not sized correctly
                    call errmgr%report_error("solve_qr_no_pivot_vec", &
                        "Incorrectly sized input array WORK, argument 4.", &
                        LA_ARRAY_SIZE_ERROR)
                    return
                end if
                wptr => work(1:lwork)
            end if
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
        real(dp), intent(out), pointer, optional, dimension(:) :: work
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
            if (.not.associated(work)) then
                allocate(wrk(lwork), stat = istat)
                if (istat /= 0) then
                    ! ERROR: Out of memory
                    call errmgr%report_error("solve_qr_pivot_mtx", &
                        "Insufficient memory available.", &
                        LA_OUT_OF_MEMORY_ERROR)
                    return
                end if
                wptr => wrk
            else
                if (size(work) < lwork) then
                    ! ERROR: WORK not sized correctly
                    call errmgr%report_error("solve_qr_no_pivot_mtx", &
                        "Incorrectly sized input array WORK, argument 5.", &
                        LA_ARRAY_SIZE_ERROR)
                    return
                end if
                wptr => work(1:lwork)
            end if
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
        real(dp), intent(out), pointer, optional, dimension(:) :: work
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
            if (.not.associated(work)) then
                allocate(wrk(lwork), stat = istat)
                if (istat /= 0) then
                    ! ERROR: Out of memory
                    call errmgr%report_error("solve_qr_pivot_vec", &
                        "Insufficient memory available.", &
                        LA_OUT_OF_MEMORY_ERROR)
                    return
                end if
                wptr => wrk
            else
                if (size(work) < lwork) then
                    ! ERROR: WORK not sized correctly
                    call errmgr%report_error("solve_qr_no_pivot_mtx", &
                        "Incorrectly sized input array WORK, argument 5.", &
                        LA_ARRAY_SIZE_ERROR)
                    return
                end if
                wptr => work(1:lwork)
            end if
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

end module
