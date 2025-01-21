module linalg_inverse
    use iso_fortran_env, only : int32, real64
    use lapack
    use blas
    use linalg_errors
    implicit none
    private
    public :: mtx_inverse
    public :: mtx_pinverse

    interface mtx_inverse
        module procedure :: mtx_inverse_dbl
        module procedure :: mtx_inverse_cmplx
    end interface

    interface mtx_pinverse
        module procedure :: mtx_pinverse_dbl
        module procedure :: mtx_pinverse_cmplx
    end interface
contains
! ------------------------------------------------------------------------------
subroutine mtx_inverse_dbl(a, iwork, work, olwork, err)
    !! Computes the inverse of a square matrix.
    real(real64), intent(inout), dimension(:,:) :: a
        !! On input, the N-by-N matrix to invert.  On output, the inverted
        !! matrix.
    integer(int32), intent(out), target, optional, dimension(:) :: iwork
        !! An optional N-element integer workspace array.
    real(real64), intent(out), target, optional, dimension(:) :: work
        !! An optional input, that if provided, prevents any local memory 
        !! allocation.  If not provided, the memory required is allocated
        !! within.  If provided, the length of the array must be at least 
        !! olwork.
    integer(int32), intent(out), optional :: olwork
        !! An optional output used to determine workspace size.  If supplied, 
        !! the routine determines the optimal size for work, and returns 
        !! without performing any actual calculations.
    class(errors), intent(inout), optional, target :: err
        !! The error object to be updated.

    ! Local Variables
    integer(int32) :: n, liwork, lwork, istat, flag, itemp(1)
    integer(int32), pointer, dimension(:) :: iptr
    integer(int32), allocatable, target, dimension(:) :: iwrk
    real(real64), pointer, dimension(:) :: wptr
    real(real64), allocatable, target, dimension(:) :: wrk
    real(real64), dimension(1) :: temp
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
        call report_square_matrix_error("mtx_inverse_dbl", errmgr, "a", &
            n, size(a, 1), size(a, 2))
        return
    end if

    ! Workspace Query
    call DGETRI(n, a, n, itemp, temp, -1, flag)
    lwork = int(temp(1), int32)
    if (present(olwork)) then
        olwork = lwork
        return
    end if

    ! Workspace Allocation
    if (present(work)) then
        if (size(work) < lwork) then
            ! ERROR: WORK not sized correctly
            call report_array_size_error("mtx_inverse_dbl", errmgr, "work", &
                lwork, size(work))
            return
        end if
        wptr => work(1:lwork)
    else
        allocate(wrk(lwork), stat = istat)
        if (istat /= 0) then
            call report_memory_error("mtx_inverse_dbl", errmgr, istat)
            return
        end if
        wptr => wrk
    end if

    ! Integer Workspace Allocation
    if (present(iwork)) then
        if (size(iwork) < liwork) then
            call report_array_size_error("mtx_inverse_dbl", errmgr, "iwork", &
                liwork, size(iwork))
            return
        end if
        iptr => iwork(1:liwork)
    else
        allocate(iwrk(liwork), stat = istat)
        if (istat /= 0) then
            call report_memory_error("mtx_inverse_dbl", errmgr, istat)
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
        call errmgr%report_error("mtx_inverse_dbl", &
            "The matrix is singular; therefore, the inverse could " // &
            "not be computed.", LA_SINGULAR_MATRIX_ERROR)
    end if
end subroutine

! ------------------------------------------------------------------------------
subroutine mtx_inverse_cmplx(a, iwork, work, olwork, err)
    !! Computes the inverse of a square matrix.
    complex(real64), intent(inout), dimension(:,:) :: a
        !! On input, the N-by-N matrix to invert.  On output, the inverted
        !! matrix.
    integer(int32), intent(out), target, optional, dimension(:) :: iwork
        !! An optional N-element integer workspace array.
    complex(real64), intent(out), target, optional, dimension(:) :: work
        !! An optional input, that if provided, prevents any local memory 
        !! allocation.  If not provided, the memory required is allocated
        !! within.  If provided, the length of the array must be at least 
        !! olwork.
    integer(int32), intent(out), optional :: olwork
        !! An optional output used to determine workspace size.  If supplied, 
        !! the routine determines the optimal size for work, and returns 
        !! without performing any actual calculations.
    class(errors), intent(inout), optional, target :: err
        !! The error object to be updated.

    ! Local Variables
    integer(int32) :: n, liwork, lwork, istat, flag, itemp(1)
    integer(int32), pointer, dimension(:) :: iptr
    integer(int32), allocatable, target, dimension(:) :: iwrk
    complex(real64), pointer, dimension(:) :: wptr
    complex(real64), allocatable, target, dimension(:) :: wrk
    complex(real64), dimension(1) :: temp
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
        call report_square_matrix_error("mtx_inverse_cmplx", errmgr, "a", &
            n, size(a, 1), size(a, 2))
        return
    end if

    ! Workspace Query
    call ZGETRI(n, a, n, itemp, temp, -1, flag)
    lwork = int(temp(1), int32)
    if (present(olwork)) then
        olwork = lwork
        return
    end if

    ! Workspace Allocation
    if (present(work)) then
        if (size(work) < lwork) then
            call report_array_size_error("mtx_inverse_cmplx", errmgr, "work", &
                lwork, size(work))
            return
        end if
        wptr => work(1:lwork)
    else
        allocate(wrk(lwork), stat = istat)
        if (istat /= 0) then
            call report_memory_error("mtx_inverse_cmplx", errmgr, istat)
            return
        end if
        wptr => wrk
    end if

    ! Integer Workspace Allocation
    if (present(iwork)) then
        if (size(iwork) < liwork) then
            call report_array_size_error("mtx_inverse_cmplx", errmgr, "iwork", &
                liwork, size(iwork))
            return
        end if
        iptr => iwork(1:liwork)
    else
        allocate(iwrk(liwork), stat = istat)
        if (istat /= 0) then
            call report_memory_error("mtx_inverse_cmplx", errmgr, istat)
            return
        end if
        iptr => iwrk
    end if

    ! Compute the LU factorization of A
    call ZGETRF(n, n, a, n, iptr, flag)

    ! Compute the inverse of the LU factored matrix
    call ZGETRI(n, a, n, iptr, wptr, lwork, flag)

    ! Check for a singular matrix
    if (flag > 0) then
        call errmgr%report_error("mtx_inverse_cmplx", &
            "The matrix is singular; therefore, the inverse could " // &
            "not be computed.", LA_SINGULAR_MATRIX_ERROR)
    end if
end subroutine

! ------------------------------------------------------------------------------
subroutine mtx_pinverse_dbl(a, ainv, tol, work, olwork, err)
    !! Computes the Moore-Penrose pseudo-inverse of a M-by-N matrix using the
    !! singular value decomposition of the matrix.
    real(real64), intent(inout), dimension(:,:) :: a
        !! On input, the M-by-N matrix to invert.  The matrix is overwritten 
        !! on output.
    real(real64), intent(out), dimension(:,:) :: ainv
        !! The N-by-M matrix where the pseudo-inverse of \(A\) will be written.
    real(real64), intent(in), optional :: tol
        !! An optional input, that if supplied, overrides the default tolerance 
        !! on singular values such that singular values less than this
        !! tolerance are forced to have a reciprocal of zero, as opposed to 
        !! 1/S(I).  The default tolerance is: MAX(M, N) * EPS * MAX(S).
    real(real64), intent(out), target, dimension(:), optional :: work
        !! An optional input, that if provided, prevents any local memory 
        !! allocation.  If not provided, the memory required is allocated
        !! within.  If provided, the length of the array must be at least 
        !! olwork.
    integer(int32), intent(out), optional :: olwork
        !! An optional output used to determine workspace size.  If supplied, 
        !! the routine determines the optimal size for work, and returns 
        !! without performing any actual calculations.
    class(errors), intent(inout), optional, target :: err
        !! The error object to be updated.

    ! Parameters
    real(real64), parameter :: zero = 0.0d0
    real(real64), parameter :: one = 1.0d0

    ! Local Variables
    integer(int32) :: i, m, n, mn, lwork, istat, flag, i1, i2a, i2b, i3a, &
        i3b, i4, lrwork
    real(real64), pointer, dimension(:) :: s, wptr, w
    real(real64), pointer, dimension(:,:) :: u, vt
    real(real64), allocatable, target, dimension(:) :: wrk
    real(real64), dimension(1) :: temp
    real(real64) :: t, tref, tolcheck, ss
    class(errors), pointer :: errmgr
    type(errors), target :: deferr

    ! Initialization
    m = size(a, 1)
    n = size(a, 2)
    mn = min(m, n)
    i1 = m * mn
    i2a = i1 + 1
    i2b = i2a + n * mn - 1
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
        call report_matrix_size_error("mtx_pinverse_dbl", errmgr, "ainv", &
            n, m, size(ainv, 1), size(ainv, 2))
        return
    end if

    ! Workspace Query
    call DGESVD('S', 'S', m, n, a, m, temp, a, m, a, n, temp, -1, flag)
    lrwork = int(temp(1), int32)
    lwork = lrwork + m * m + n * n + mn
    if (present(olwork)) then
        olwork = lwork
        return
    end if

    ! Local Memory Allocation
    if (present(work)) then
        if (size(work) < lwork) then
            call report_array_size_error("mtx_pinverse_dbl", errmgr, "work", &
                lwork, size(work))
            return
        end if
        wptr => work(1:lwork)
    else
        allocate(wrk(lwork), stat = istat)
        if (istat /= 0) then
            call report_memory_error("mtx_pinverse_dbl", errmgr, istat)
            return
        end if
        wptr => wrk
    end if
    u(1:m,1:mn) => wptr(1:i1)
    vt(1:mn,1:n) => wptr(i2a:i2b)
    s => wptr(i3a:i3b)
    w => wptr(i4:lwork)

    ! Compute the SVD of A
    call DGESVD('S', 'S', m, n, a, m, s, u, m, vt, n, w, lrwork, flag)

    ! Check for convergence
    if (flag > 0) then
        call errmgr%report_error("mtx_pinverse_dbl", &
            "The QR iteration process could not converge.", &
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
    ! first computing inv(S) * U**T 
    do i = 1, mn
        if (s(i) < t) then
            ss = s(i)
        else
            ss = 1.0d0 / s(i)
        end if
        call DSCAL(m, ss, u(:,i), 1)
    end do
    call DGEMM("T", "T", n, m, mn, one, vt, n, u, m, zero, ainv, n)
end subroutine

! ------------------------------------------------------------------------------
subroutine mtx_pinverse_cmplx(a, ainv, tol, work, olwork, rwork, err)
    !! Computes the Moore-Penrose pseudo-inverse of a M-by-N matrix using the
    !! singular value decomposition of the matrix.
    complex(real64), intent(inout), dimension(:,:) :: a
        !! On input, the M-by-N matrix to invert.  The matrix is overwritten 
        !! on output.
    complex(real64), intent(out), dimension(:,:) :: ainv
        !! The N-by-M matrix where the pseudo-inverse of \(A\) will be written.
    real(real64), intent(in), optional :: tol
        !! An optional input, that if supplied, overrides the default tolerance 
        !! on singular values such that singular values less than this
        !! tolerance are forced to have a reciprocal of zero, as opposed to 
        !! 1/S(I).  The default tolerance is: MAX(M, N) * EPS * MAX(S).
    complex(real64), intent(out), target, dimension(:), optional :: work
        !! An optional input, that if provided, prevents any local memory 
        !! allocation.  If not provided, the memory required is allocated
        !! within.  If provided, the length of the array must be at least 
        !! olwork.
    integer(int32), intent(out), optional :: olwork
        !! An optional output used to determine workspace size.  If supplied, 
        !! the routine determines the optimal size for work, and returns 
        !! without performing any actual calculations.
    real(real64), intent(out), target, dimension(:), optional :: rwork
        !! An optional input, that if provided, prevents any local memory 
        !! allocation for real-valued workspaces.  If not provided, the 
        !! memory required is allocated within.  If provided, the length of the 
        !! array must be at least 6 * MIN(M, N).
    class(errors), intent(inout), optional, target :: err
        !! The error object to be updated.

    ! External Function Interfaces
    interface
        function DLAMCH(cmach) result(x)
            use, intrinsic :: iso_fortran_env, only : real64
            character, intent(in) :: cmach
            real(real64) :: x
        end function
    end interface

    ! Parameters
    complex(real64), parameter :: zero = (0.0d0, 0.0d0)
    complex(real64), parameter :: one = (1.0d0, 0.0d0)

    ! Local Variables
    integer(int32) :: i, m, n, mn, lwork, istat, flag, i1, i2a, i2b, i3, &
        lrwork, j, k
    real(real64), pointer, dimension(:) :: s, rwptr, rw
    real(real64), allocatable, target, dimension(:) :: rwrk
    complex(real64), pointer, dimension(:) :: wptr, w
    complex(real64), pointer, dimension(:,:) :: u, vt
    complex(real64), allocatable, target, dimension(:) :: wrk
    complex(real64) :: temp(1), val
    real(real64) :: t, tref, tolcheck, rtemp(1)
    class(errors), pointer :: errmgr
    type(errors), target :: deferr

    ! Initialization
    m = size(a, 1)
    n = size(a, 2)
    mn = min(m, n)
    lrwork = 6 * mn
    i1 = m * mn
    i2a = i1 + 1
    i2b = i2a + n * n - 1
    i3 = i2b + 1
    tolcheck = dlamch('s')
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if

    ! Input Check
    if (size(ainv, 1) /= n .or. size(ainv, 2) /= m) then
        call report_matrix_size_error("mtx_pinverse_cmplx", errmgr, "ainv", &
            n, m, size(ainv, 1), size(ainv, 2))
        return
    end if

    ! Workspace Query
    call ZGESVD('S', 'A', m, n, a, m, rtemp, a, m, a, n, temp, -1, &
        rtemp, flag)
    lwork = int(temp(1), int32)
    lwork = lwork + m * mn + n * n
    if (present(olwork)) then
        olwork = lwork
        return
    end if

    ! Local Memory Allocation
    if (present(work)) then
        if (size(work) < lwork) then
            call report_array_size_error("mtx_pinverse_cmplx", errmgr, "work", &
                lwork, size(work))
            return
        end if
        wptr => work(1:lwork)
    else
        allocate(wrk(lwork), stat = istat)
        if (istat /= 0) then
            call report_memory_error("mtx_pinverse_cmplx", errmgr, istat)
            return
        end if
        wptr => wrk
    end if

    if (present(rwork)) then
        if (size(rwork) < lrwork) then
            call report_array_size_error("mtx_pinverse_cmplx", errmgr, &
                "rwork", lrwork, size(rwork))
            return
        end if
        rwptr => rwork(1:lrwork)
    else
        allocate(rwrk(lrwork), stat = istat)
        if (istat /= 0) then
            call report_memory_error("mtx_pinverse_cmplx", errmgr, istat)
            return
        end if
        rwptr => rwrk
    end if
    u(1:m,1:mn) => wptr(1:i1)
    vt(1:n,1:n) => wptr(i2a:i2b)
    w => wptr(i3:lwork)
    s => rwptr(1:mn)
    rw => rwptr(mn+1:lrwork)

    ! Compute the SVD of A
    call ZGESVD('S', 'A', m, n, a, m, s, u, m, vt, n, w, size(w), rw, flag)

    ! Check for convergence
    if (flag > 0) then
        call errmgr%report_error("mtx_pinverse_cmplx", &
            "The QR iteration process could not converge.", &
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
            ! call recip_mult_array(s(i), vt(i,1:n))
            vt(i,1:n) = conjg(vt(i,1:n)) / s(i)
        end if
    end do

    ! Compute (VT**T * inv(S)) * U**H
    ! ainv = n-by-m
    ! vt is n-by-n
    ! u is m-by-mn such that u**H = mn-by-m
    ! Compute ainv = vt**T * u**H
    do j = 1, m
        do i = 1, n
            val = zero
            do k = 1, mn
                val = val + vt(k,i) * conjg(u(j,k))
            end do
            ainv(i,j) = val
        end do
    end do
end subroutine

! ------------------------------------------------------------------------------
end module