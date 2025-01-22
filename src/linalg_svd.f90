module linalg_svd
    use iso_fortran_env, only : int32, real64
    use lapack
    use linalg_errors
    implicit none
    private
    public :: svd

    interface svd
        module procedure :: svd_dbl
        module procedure :: svd_cmplx
    end interface
contains
! ------------------------------------------------------------------------------
subroutine svd_dbl(a, s, u, vt, work, olwork, err)
    !! Computes the singular value decomposition of an M-by-N matrix \(A\) such 
    !! that \(A = U S V^T\) where \(U\) is an M-by-M orthogonal matrix, \(S\)
    !! is an M-by-N diagonal matrix containing the singular values, and \(V\)
    !! is an N-by-N orthogonal matrix.
    real(real64), intent(inout), dimension(:,:) :: a
        !! On input, the M-by-N matrix to factor.  The matrix is overwritten on 
        !! output.
    real(real64), intent(out), dimension(:) :: s
        !! A MIN(M, N)-element array containing the singular values of a sorted 
        !! in descending order.
    real(real64), intent(out), optional, dimension(:,:) :: u
        !! An optional argument, that if supplied, is used to contain the 
        !! orthogonal matrix \(U\) from the decomposition.  The matrix \(U\) 
        !! contains the left singular vectors, and can be either M-by-M 
        !! (all left singular vectors are computed), or M-by-MIN(M,N) (only the 
        !! first MIN(M, N) left singular vectors are computed).
    real(real64), intent(out), optional, dimension(:,:) :: vt
        !! An optional argument, that if supplied, is used to contain the 
        !! transpose of the N-by-N orthogonal matrix \(V\).  The matrix \(V\) 
        !! contains the right singular vectors.
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
    character :: jobu, jobvt
    integer(int32) :: m, n, mn, istat, lwork, flag
    real(real64), pointer, dimension(:) :: wptr
    real(real64), allocatable, target, dimension(:) :: wrk
    real(real64), dimension(1) :: temp
    class(errors), pointer :: errmgr
    type(errors), target :: deferr

    ! Initialization
    m = size(a, 1)
    n = size(a, 2)
    mn = min(m, n)
    if (present(u)) then
        if (size(u, 2) == m) then
            jobu = 'A'
        else if (size(u, 2) == mn) then
            jobu = 'S'
        end if
    else
        jobu = 'N'
    end if
    if (present(vt)) then
        jobvt = 'A'
    else
        jobvt = 'N'
    end if
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if

    ! Input Check
    if (size(s) /= mn) then
        call report_array_size_error("svd_dbl", errmgr, "s", mn, size(s))
        return
    else if (present(u)) then
        if (size(u, 1) /= m) then
            call report_matrix_size_error("svd_dbl", errmgr, "u", m, m, &
                size(u, 1), size(u, 2))
            return
        end if
        if (size(u, 2) /= m .and. size(u, 2) /= mn) then
            call report_matrix_size_error("svd_dbl", errmgr, "u", m, m, &
                size(u, 1), size(u, 2))
            return
        end if
    else if (present(vt)) then
        if (size(vt, 1) /= n .or. size(vt, 2) /= n) then
            call report_matrix_size_error("svd_dbl", errmgr, "vt", n, n, &
                size(vt, 1), size(vt, 2))
            return
        end if
    end if

    ! Workspace Query
    call DGESVD(jobu, jobvt, m, n, a, m, s, temp, m, temp, n, temp, -1, &
        flag)
    lwork = int(temp(1), int32)
    if (present(olwork)) then
        olwork = lwork
        return
    end if

    ! Local Memory Allocation
    if (present(work)) then
        if (size(work) < lwork) then
            call report_array_size_error("svd_dbl", errmgr, "work", lwork, &
                size(work))
            return
        end if
        wptr => work(1:lwork)
    else
        allocate(wrk(lwork), stat = istat)
        if (istat /= 0) then
            call report_memory_error("svd_dbl", errmgr, istat)
            return
        end if
        wptr => wrk
    end if

    ! Call DGESVD
    if (present(u) .and. present(vt)) then
        call DGESVD(jobu, jobvt, m, n, a, m, s, u, m, vt, n, wptr, lwork, &
            flag)
    else if (present(u) .and. .not.present(vt)) then
        call DGESVD(jobu, jobvt, m, n, a, m, s, u, m, temp, n, wptr, &
            lwork, flag)
    else if (.not.present(u) .and. present(vt)) then
        call DGESVD(jobu, jobvt, m, n, a, m, s, temp, m, vt, n, wptr, &
            lwork, flag)
    else
        call DGESVD(jobu, jobvt, m, n, a, m, s, temp, m, temp, n, wptr, &
            lwork, flag)
    end if

    ! Check for convergence
    if (flag > 0) then
        call errmgr%report_error("svd_dbl", &
            "The QR iteration process could not converge.", &
            LA_CONVERGENCE_ERROR)
        return
    end if
end subroutine

! ------------------------------------------------------------------------------
subroutine svd_cmplx(a, s, u, vt, work, olwork, rwork, err)
    !! Computes the singular value decomposition of an M-by-N matrix \(A\) such 
    !! that \(A = U S V^H\) where \(U\) is an M-by-M orthogonal matrix, \(S\)
    !! is an M-by-N diagonal matrix containing the singular values, and \(V\)
    !! is an N-by-N orthogonal matrix.
    complex(real64), intent(inout), dimension(:,:) :: a
        !! On input, the M-by-N matrix to factor.  The matrix is overwritten on 
        !! output.
    real(real64), intent(out), dimension(:) :: s
        !! A MIN(M, N)-element array containing the singular values of a sorted 
        !! in descending order.
    complex(real64), intent(out), optional, dimension(:,:) :: u
        !! An optional argument, that if supplied, is used to contain the 
        !! orthogonal matrix \(U\) from the decomposition.  The matrix \(U\) 
        !! contains the left singular vectors, and can be either M-by-M 
        !! (all left singular vectors are computed), or M-by-MIN(M,N) (only the 
        !! first MIN(M, N) left singular vectors are computed).
    complex(real64), intent(out), optional, dimension(:,:) :: vt
        !! An optional argument, that if supplied, is used to contain the 
        !! conjugate transpose of the N-by-N orthogonal matrix \(V\).  The 
        !! matrix \(V\) contains the right singular vectors.
    complex(real64), intent(out), target, optional, dimension(:) :: work
        !! An optional input, that if provided, prevents any local memory 
        !! allocation.  If not provided, the memory required is allocated
        !! within.  If provided, the length of the array must be at least 
        !! olwork.
    integer(int32), intent(out), optional :: olwork
        !! An optional output used to determine workspace size.  If supplied, 
        !! the routine determines the optimal size for work, and returns 
        !! without performing any actual calculations.
    real(real64), intent(out), target, optional, dimension(:) :: rwork
        !! An optional input, that if provided, prevents any local memory 
        !! allocation for real-valued workspaces.  If not provided, the memory 
        !! required is allocated within.  If provided, the length of the array 
        !! must be at least 5 * MIN(M, N).
    class(errors), intent(inout), optional, target :: err
        !! The error object to be updated.

    ! Local Variables
    character :: jobu, jobvt
    integer(int32) :: m, n, mn, istat, lwork, flag, lrwork
    complex(real64), pointer, dimension(:) :: wptr
    complex(real64), allocatable, target, dimension(:) :: wrk
    complex(real64), dimension(1) :: temp
    real(real64), dimension(1) :: rtemp
    real(real64), pointer, dimension(:) :: rwptr
    real(real64), allocatable, target, dimension(:) :: rwrk
    class(errors), pointer :: errmgr
    type(errors), target :: deferr

    ! Initialization
    m = size(a, 1)
    n = size(a, 2)
    mn = min(m, n)
    lrwork = 5 * mn
    if (present(u)) then
        if (size(u, 2) == m) then
            jobu = 'A'
        else if (size(u, 2) == mn) then
            jobu = 'S'
        end if
    else
        jobu = 'N'
    end if
    if (present(vt)) then
        jobvt = 'A'
    else
        jobvt = 'N'
    end if
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if

    ! Input Check
    if (size(s) /= mn) then
        call report_array_size_error("svd_cmplx", errmgr, "s", mn, size(s))
        return
    else if (present(u)) then
        if (size(u, 1) /= m) then
            call report_matrix_size_error("svd_cmplx", errmgr, "u", m, m, &
                size(u, 1), size(u, 2))
            return
        end if
        if (size(u, 2) /= m .and. size(u, 2) /= mn) then
            call report_matrix_size_error("svd_cmplx", errmgr, "u", m, m, &
                size(u, 1), size(u, 2))
            return
        end if
    else if (present(vt)) then
        if (size(vt, 1) /= n .or. size(vt, 2) /= n) then
            call report_matrix_size_error("svd_cmplx", errmgr, "vt", n, n, &
                size(vt, 1), size(vt, 2))
            return
        end if
    end if

    ! Workspace Query
    call ZGESVD(jobu, jobvt, m, n, a, m, s, temp, m, temp, n, temp, -1, &
        rtemp, flag)
    lwork = int(temp(1), int32)
    if (present(olwork)) then
        olwork = lwork
        return
    end if

    ! Local Memory Allocation
    if (present(work)) then
        if (size(work) < lwork) then
            call report_array_size_error("svd_cmplx", errmgr, "work", lwork, &
                size(work))
            return
        end if
        wptr => work(1:lwork)
    else
        allocate(wrk(lwork), stat = istat)
        if (istat /= 0) then
            call report_memory_error("svd_cmplx", errmgr, istat)
            return
        end if
        wptr => wrk
    end if

    if (present(rwork)) then
        if (size(rwork) < lrwork) then
            call report_array_size_error("svd_cmplx", errmgr, "rwork", lrwork, &
                size(rwork))
            return
        end if
        rwptr => rwork(1:lrwork)
    else
        allocate(rwrk(lrwork), stat = istat)
        if (istat /= 0) then
            call report_memory_error("svd_cmplx", errmgr, istat)
            return
        end if
        rwptr => rwrk
    end if

    ! Call ZGESVD
    if (present(u) .and. present(vt)) then
        call ZGESVD(jobu, jobvt, m, n, a, m, s, u, m, vt, n, wptr, lwork, &
            rwptr, flag)
    else if (present(u) .and. .not.present(vt)) then
        call ZGESVD(jobu, jobvt, m, n, a, m, s, u, m, temp, n, wptr, &
            lwork, rwptr, flag)
    else if (.not.present(u) .and. present(vt)) then
        call ZGESVD(jobu, jobvt, m, n, a, m, s, temp, m, vt, n, wptr, &
            lwork, rwptr, flag)
    else
        call ZGESVD(jobu, jobvt, m, n, a, m, s, temp, m, temp, n, wptr, &
            lwork, rwptr, flag)
    end if

    ! Check for convergence
    if (flag > 0) then
        call errmgr%report_error("svd_cmplx", &
            "The QR iteration process could not converge.", &
            LA_CONVERGENCE_ERROR)
        return
    end if
end subroutine

! ------------------------------------------------------------------------------
end module