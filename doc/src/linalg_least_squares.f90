module linalg_least_squares
    use iso_fortran_env, only : int32, real64
    use lapack
    use linalg_errors
    implicit none
    private
    public :: solve_least_squares
    public :: solve_least_squares_full
    public :: solve_least_squares_svd

    interface solve_least_squares
        module procedure :: solve_least_squares_mtx
        module procedure :: solve_least_squares_mtx_cmplx
        module procedure :: solve_least_squares_vec
        module procedure :: solve_least_squares_vec_cmplx
    end interface

    interface solve_least_squares_full
        module procedure :: solve_least_squares_mtx_pvt
        module procedure :: solve_least_squares_mtx_pvt_cmplx
        module procedure :: solve_least_squares_vec_pvt
        module procedure :: solve_least_squares_vec_pvt_cmplx
    end interface

    interface solve_least_squares_svd
        module procedure :: solve_least_squares_mtx_svd
        module procedure :: solve_least_squares_mtx_svd_cmplx
        module procedure :: solve_least_squares_vec_svd
        module procedure :: solve_least_squares_vec_svd_cmplx
    end interface
contains
! ------------------------------------------------------------------------------
subroutine solve_least_squares_mtx(a, b, work, olwork, err)
    !! Solves the system of equations \(A X = B\) assuming matrix \(A\) is of 
    !! full rank.
    real(real64), intent(inout), dimension(:,:) :: a 
        !! On input, the M-by-N matrix \(A\).  On output, if M is greater than
        !! or equal to N, the QR factorization of \(A\) in the form provided
        !! by qr_factor; else, if M is less than N, the LQ factorization of
        !! \(A\) as returned by lq_factor.
    real(real64), intent(inout), dimension(:,:) :: b
        !! If the system is overdetermined, the M-by-NRHS matrix \(B\); else,
        !! the matrix should be sized as N-by-NRHS with the first M rows 
        !! containing \(B\).  On output, the first N rows will contain the
        !! solution matrix \(X\).
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
    integer(int32) :: m, n, maxmn, nrhs, lwork, istat, flag
    real(real64), pointer, dimension(:) :: wptr
    real(real64), allocatable, target, dimension(:) :: wrk
    real(real64), dimension(1) :: temp
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
        call report_matrix_size_error("solve_least_squares_mtx", errmgr, &
            "b", maxmn, nrhs, size(b, 1), size(b, 2))
        return
    end if

    ! Workspace Query
    call DGELS('N', m, n, nrhs, a, m, b, maxmn, temp, -1, flag)
    lwork = int(temp(1), int32)
    if (present(olwork)) then
        olwork = lwork
        return
    end if

    ! Local Memory Allocation
    if (present(work)) then
        if (size(work) < lwork) then
            call report_array_size_error("solve_least_squares_mtx", errmgr, &
                "work", lwork, size(work))
            return
        end if
        wptr => work(1:lwork)
    else
        allocate(wrk(lwork), stat = istat)
        if (istat /= 0) then
            call report_memory_error("solve_least_squares_mtx", errmgr, istat)
            return
        end if
        wptr => wrk
    end if

    ! Process
    call DGELS('N', m, n, nrhs, a, m, b, maxmn, wptr, lwork, flag)
    if (flag > 0) then
        call errmgr%report_error("solve_least_squares_mtx", &
            "The supplied matrix is not of full rank; therefore, " // &
            "the solution could not be computed via this routine.  " // &
            "Try a routine that utilizes column pivoting.", &
            LA_INVALID_OPERATION_ERROR)
    end if
end subroutine

! ------------------------------------------------------------------------------
subroutine solve_least_squares_mtx_cmplx(a, b, work, olwork, err)
    !! Solves the system of equations \(A X = B\) assuming matrix \(A\) is of 
    !! full rank.
    complex(real64), intent(inout), dimension(:,:) :: a 
        !! On input, the M-by-N matrix \(A\).  On output, if M is greater than
        !! or equal to N, the QR factorization of \(A\) in the form provided
        !! by qr_factor; else, if M is less than N, the LQ factorization of
        !! \(A\) as returned by lq_factor.
    complex(real64), intent(inout), dimension(:,:) :: b
        !! If the system is overdetermined, the M-by-NRHS matrix \(B\); else,
        !! the matrix should be sized as N-by-NRHS with the first M rows 
        !! containing \(B\).  On output, the first N rows will contain the
        !! solution matrix \(X\).
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
    integer(int32) :: m, n, maxmn, nrhs, lwork, istat, flag
    complex(real64), pointer, dimension(:) :: wptr
    complex(real64), allocatable, target, dimension(:) :: wrk
    complex(real64), dimension(1) :: temp
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
        call report_matrix_size_error("solve_least_squares_mtx_cmplx", &
            errmgr, "b", maxmn, nrhs, size(b, 1), size(b, 2))
        return
    end if

    ! Workspace Query
    call ZGELS('N', m, n, nrhs, a, m, b, maxmn, temp, -1, flag)
    lwork = int(temp(1), int32)
    if (present(olwork)) then
        olwork = lwork
        return
    end if

    ! Local Memory Allocation
    if (present(work)) then
        if (size(work) < lwork) then
            call report_array_size_error("solve_least_squares_mtx_cmplx", &
                errmgr, "work", lwork, size(work))
            return
        end if
        wptr => work(1:lwork)
    else
        allocate(wrk(lwork), stat = istat)
        if (istat /= 0) then
            call report_memory_error("solve_least_squares_mtx_cmplx", errmgr, &
                istat)
            return
        end if
        wptr => wrk
    end if

    ! Process
    call ZGELS('N', m, n, nrhs, a, m, b, maxmn, wptr, lwork, flag)
    if (flag > 0) then
        call errmgr%report_error("solve_least_squares_mtx_cmplx", &
            "The supplied matrix is not of full rank; therefore, " // &
            "the solution could not be computed via this routine.  " // &
            "Try a routine that utilizes column pivoting.", &
            LA_INVALID_OPERATION_ERROR)
    end if
end subroutine

! ------------------------------------------------------------------------------
subroutine solve_least_squares_vec(a, b, work, olwork, err)
    !! Solves the system of equations \(A \vec{x} = \vec{b}\) assuming matrix 
    !! \(A\) is of full rank.
    real(real64), intent(inout), dimension(:,:) :: a
        !! On input, the M-by-N matrix \(A\).  On output, if M is greater than
        !! or equal to N, the QR factorization of \(A\) in the form provided
        !! by qr_factor; else, if M is less than N, the LQ factorization of
        !! \(A\) as returned by lq_factor.
    real(real64), intent(inout), dimension(:) :: b
        !! If the system is overdetermined, the M-element vector \(\vec{b}\);
        !! else, the array should be sized as N-element with the first M
        !! elements containing \(\vec{b}\).  On output, the first N rows will
        !! contain the solution vector \(\vec{x}\).
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
    integer(int32) :: m, n, maxmn, lwork, istat, flag
    real(real64), pointer, dimension(:) :: wptr
    real(real64), allocatable, target, dimension(:) :: wrk
    real(real64), dimension(1) :: temp
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
        call report_array_size_error("solve_least_squares_vec", errmgr, "b", &
            maxmn, size(b))
        return
    end if

    ! Workspace Query
    call DGELS('N', m, n, 1, a, m, b, maxmn, temp, -1, flag)
    lwork = int(temp(1), int32)
    if (present(olwork)) then
        olwork = lwork
        return
    end if

    ! Local Memory Allocation
    if (present(work)) then
        if (size(work) < lwork) then
            call report_array_size_error("solve_least_squares_vec", errmgr, &
                "work", lwork, size(work))
            return
        end if
        wptr => work(1:lwork)
    else
        allocate(wrk(lwork), stat = istat)
        if (istat /= 0) then
            call report_memory_error("solve_least_squares_vec", errmgr, istat)
            return
        end if
        wptr => wrk
    end if

    ! Process
    call DGELS('N', m, n, 1, a, m, b, maxmn, wptr, lwork, flag)
    if (flag > 0) then
        call errmgr%report_error("solve_least_squares_vec", &
            "The supplied matrix is not of full rank; therefore, " // &
            "the solution could not be computed via this routine.  " // &
            "Try a routine that utilizes column pivoting.", &
            LA_INVALID_OPERATION_ERROR)
    end if
end subroutine

! ------------------------------------------------------------------------------
subroutine solve_least_squares_vec_cmplx(a, b, work, olwork, err)
    !! Solves the system of equations \(A \vec{x} = \vec{b}\) assuming matrix 
    !! \(A\) is of full rank.
    complex(real64), intent(inout), dimension(:,:) :: a
        !! On input, the M-by-N matrix \(A\).  On output, if M is greater than
        !! or equal to N, the QR factorization of \(A\) in the form provided
        !! by qr_factor; else, if M is less than N, the LQ factorization of
        !! \(A\) as returned by lq_factor.
    complex(real64), intent(inout), dimension(:) :: b
        !! If the system is overdetermined, the M-element vector \(\vec{b}\);
        !! else, the array should be sized as N-element with the first M
        !! elements containing \(\vec{b}\).  On output, the first N rows will
        !! contain the solution vector \(\vec{x}\).
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
    integer(int32) :: m, n, maxmn, lwork, istat, flag
    complex(real64), pointer, dimension(:) :: wptr
    complex(real64), allocatable, target, dimension(:) :: wrk
    complex(real64), dimension(1) :: temp
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
        call report_array_size_error("solve_least_squares_vec_cmplx", errmgr, &
            "b", maxmn, size(b))
        return
    end if

    ! Workspace Query
    call ZGELS('N', m, n, 1, a, m, b, maxmn, temp, -1, flag)
    lwork = int(temp(1), int32)
    if (present(olwork)) then
        olwork = lwork
        return
    end if

    ! Local Memory Allocation
    if (present(work)) then
        if (size(work) < lwork) then
            call report_array_size_error("solve_least_squares_vec_cmplx", &
                errmgr, "work", lwork, size(work))
            return
        end if
        wptr => work(1:lwork)
    else
        allocate(wrk(lwork), stat = istat)
        if (istat /= 0) then
            call report_memory_error("solve_least_squares_vec_cmplx", errmgr, &
                istat)
            return
        end if
        wptr => wrk
    end if

    ! Process
    call ZGELS('N', m, n, 1, a, m, b, maxmn, wptr, lwork, flag)
    if (flag > 0) then
        call errmgr%report_error("solve_least_squares_mtx_cmplx", &
            "The supplied matrix is not of full rank; therefore, " // &
            "the solution could not be computed via this routine.  " // &
            "Try a routine that utilizes column pivoting.", &
            LA_INVALID_OPERATION_ERROR)
    end if
end subroutine

! ------------------------------------------------------------------------------
subroutine solve_least_squares_mtx_pvt(a, b, ipvt, arnk, work, olwork, err)
    !! Solves the system of equations \(A X = B\) using a full orthogonal
    !! factorization of \(A\).
    real(real64), intent(inout), dimension(:,:) :: a
        !! On input, the M-by-N matrix \(A\).  On output, the matrix is 
        !! overwritten by its orthogonal factorization.
    real(real64), intent(inout), dimension(:,:) :: b
        !! If the system is overdetermined, the M-by-NRHS matrix \(B\); else,
        !! the matrix should be sized as N-by-NRHS with the first M rows 
        !! containing \(B\).  On output, the first N rows will contain the
        !! solution matrix \(X\).
    integer(int32), intent(inout), target, optional, dimension(:) :: ipvt
        !! An optional input that on input, an N-element array that if 
        !! IPVT(I) .ne. 0, the I-th column of A is permuted to the front
        !! of A * P; if IPVT(I) = 0, the I-th column of A is a free column.  On
        !! output, if IPVT(I) = K, then the I-th column of A * P was the K-th
        !! column of A.  If not supplied, memory is allocated internally, and 
        !! IPVT is set to all zeros such that all columns are treated as free.
    integer(int32), intent(out), optional :: arnk
        !! An optional output, that if provided, will return the rank of \(A\).
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
    integer(int32) :: m, n, maxmn, nrhs, lwork, istat, flag, rnk
    real(real64), pointer, dimension(:) :: wptr
    real(real64), allocatable, target, dimension(:) :: wrk
    integer(int32), allocatable, target, dimension(:) :: iwrk
    integer(int32), pointer, dimension(:) :: iptr
    real(real64), dimension(1) :: temp
    integer(int32), dimension(1) :: itemp
    real(real64) :: rc
    class(errors), pointer :: errmgr
    type(errors), target :: deferr

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
    if (size(b, 1) /= maxmn) then
        call report_matrix_size_error("solve_least_squares_mtx_pvt", errmgr, &
            "b", maxmn, nrhs, size(b, 1), size(b, 2))
        return
    end if

    ! Workspace Query
    call DGELSY(m, n, nrhs, a, m, b, maxmn, itemp, rc, rnk, temp, -1, flag)
    lwork = int(temp(1), int32)
    if (present(olwork)) then
        olwork = lwork
        return
    end if

    ! Local Memory Allocation
    if (present(ipvt)) then
        if (size(ipvt) < n) then
            call report_array_size_error("solve_least_squares_mtx_pvt", &
                errmgr, "ipvt", n, size(ipvt))
            return
        end if
        iptr => ipvt(1:n)
    else
        allocate(iwrk(n), stat = istat)
        if (istat /= 0) then
            call report_memory_error("solve_least_squares_mtx_pvt", errmgr, &
                istat)
            return
        end if
        iptr => iwrk
        iptr = 0
    end if

    if (present(work)) then
        if (size(work) < lwork) then
            call report_array_size_error("solve_least_squares_mtx_pvt", &
                errmgr, "work", lwork, size(work))
            return
        end if
        wptr => work(1:lwork)
    else
        allocate(wrk(lwork), stat = istat)
        if (istat /= 0) then
            call report_memory_error("solve_least_squares_mtx_pvt", errmgr, &
                istat)
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
subroutine solve_least_squares_mtx_pvt_cmplx(a, b, ipvt, arnk, &
    work, olwork, rwork, err)
    !! Solves the system of equations \(A X = B\) using a full orthogonal
    !! factorization of \(A\).
    complex(real64), intent(inout), dimension(:,:) :: a
        !! On input, the M-by-N matrix \(A\).  On output, the matrix is 
        !! overwritten by its orthogonal factorization.
    complex(real64), intent(inout), dimension(:,:) :: b
        !! If the system is overdetermined, the M-by-NRHS matrix \(B\); else,
        !! the matrix should be sized as N-by-NRHS with the first M rows 
        !! containing \(B\).  On output, the first N rows will contain the
        !! solution matrix \(X\).
    integer(int32), intent(inout), target, optional, dimension(:) :: ipvt
        !! An optional input that on input, an N-element array that if 
        !! IPVT(I) .ne. 0, the I-th column of A is permuted to the front
        !! of A * P; if IPVT(I) = 0, the I-th column of A is a free column.  On
        !! output, if IPVT(I) = K, then the I-th column of A * P was the K-th
        !! column of A.  If not supplied, memory is allocated internally, and 
        !! IPVT is set to all zeros such that all columns are treated as free.
    integer(int32), intent(out), optional :: arnk
        !! An optional output, that if provided, will return the rank of \(A\).
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
        !! must be at least 2 * N.
    class(errors), intent(inout), optional, target :: err
        !! The error object to be updated.

    ! Local Variables
    integer(int32) :: m, n, maxmn, nrhs, lwork, istat, flag, rnk, lrwork
    complex(real64), pointer, dimension(:) :: wptr
    complex(real64), allocatable, target, dimension(:) :: wrk
    real(real64), pointer, dimension(:) :: rwptr
    real(real64), allocatable, target, dimension(:) :: rwrk
    integer(int32), allocatable, target, dimension(:) :: iwrk
    integer(int32), pointer, dimension(:) :: iptr
    complex(real64), dimension(1) :: temp
    real(real64), dimension(1) :: rtemp
    integer(int32), dimension(1) :: itemp
    real(real64) :: rc
    class(errors), pointer :: errmgr
    type(errors), target :: deferr

    ! Initialization
    m = size(a, 1)
    n = size(a, 2)
    maxmn = max(m, n)
    nrhs = size(b, 2)
    lrwork = 2 * n
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
        call report_matrix_size_error("solve_least_squares_mtx_pvt_cmplx", &
            errmgr, "b", maxmn, nrhs, size(b, 1), size(b, 2))
        return
    end if

    ! Workspace Query
    call ZGELSY(m, n, nrhs, a, m, b, maxmn, itemp, rc, rnk, temp, -1, &
        rtemp, flag)
    lwork = int(temp(1), int32)
    if (present(olwork)) then
        olwork = lwork
        return
    end if

    ! Local Memory Allocation
    if (present(ipvt)) then
        if (size(ipvt) < n) then
            call report_array_size_error("solve_least_squares_mtx_pvt_cmplx", &
                errmgr, "ipvt", n, size(ipvt))
            return
        end if
        iptr => ipvt(1:n)
    else
        allocate(iwrk(n), stat = istat)
        if (istat /= 0) then
            call report_memory_error("solve_least_squares_mtx_pvt_cmplx", &
                errmgr, istat)
            return
        end if
        iptr => iwrk
        iptr = 0
    end if

    if (present(work)) then
        if (size(work) < lwork) then
            call report_array_size_error("solve_least_squares_mtx_pvt_cmplx", &
                errmgr, "work", lwork, size(work))
            return
        end if
        wptr => work(1:lwork)
    else
        allocate(wrk(lwork), stat = istat)
        if (istat /= 0) then
            call report_memory_error("solve_least_squares_mtx_pvt_cmplx", &
                errmgr, istat)
            return
        end if
        wptr => wrk
    end if

    if (present(rwork)) then
        if (size(rwork) < lrwork) then
            call report_array_size_error("solve_least_squares_mtx_pvt_cmplx", &
                errmgr, "rwork", lrwork, size(rwork))
            return
        end if
        rwptr => rwork(1:lrwork)
    else
        allocate(rwrk(lrwork), stat = istat)
        if (istat /= 0) then
            call report_memory_error("solve_least_squares_mtx_pvt_cmplx", &
                errmgr, istat)
            return
        end if
        rwptr => rwrk
    end if

    ! Process
    call ZGELSY(m, n, nrhs, a, m, b, maxmn, iptr, rc, rnk, wptr, lwork, &
        rwptr, flag)
    if (present(arnk)) arnk = rnk
end subroutine

! ------------------------------------------------------------------------------
subroutine solve_least_squares_vec_pvt(a, b, ipvt, arnk, work, olwork, err)
    !! Solves the system of equations \(A \vec{x} = \vec{b}\) using a full 
    !! orthogonal factorization of \(A\).
    real(real64), intent(inout), dimension(:,:) :: a
        !! On input, the M-by-N matrix \(A\).  On output, the matrix is 
        !! overwritten by its orthogonal factorization.
    real(real64), intent(inout), dimension(:) :: b
        !! If the system is overdetermined, the M-element vector \(\vec{b}\);
        !! else, the array should be sized as N-element with the first M
        !! elements containing \(\vec{b}\).  On output, the first N rows will
        !! contain the solution vector \(\vec{x}\).
    integer(int32), intent(inout), target, optional, dimension(:) :: ipvt
        !! An optional input that on input, an N-element array that if 
        !! IPVT(I) .ne. 0, the I-th column of A is permuted to the front
        !! of A * P; if IPVT(I) = 0, the I-th column of A is a free column.  On
        !! output, if IPVT(I) = K, then the I-th column of A * P was the K-th
        !! column of A.  If not supplied, memory is allocated internally, and 
        !! IPVT is set to all zeros such that all columns are treated as free.
    integer(int32), intent(out), optional :: arnk
        !! An optional output, that if provided, will return the rank of \(A\).
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
    integer(int32) :: m, n, maxmn, lwork, istat, flag, rnk
    real(real64), pointer, dimension(:) :: wptr
    real(real64), allocatable, target, dimension(:) :: wrk
    integer(int32), allocatable, target, dimension(:) :: iwrk
    integer(int32), pointer, dimension(:) :: iptr
    real(real64), dimension(1) :: temp
    integer(int32), dimension(1) :: itemp
    real(real64) :: rc
    class(errors), pointer :: errmgr
    type(errors), target :: deferr

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
        call report_array_size_error("solve_least_squares_vec_pvt", errmgr, &
            "b", maxmn, size(b))
        return
    end if

    ! Workspace Query
    call DGELSY(m, n, 1, a, m, b, maxmn, itemp, rc, rnk, temp, -1, flag)
    lwork = int(temp(1), int32)
    if (present(olwork)) then
        olwork = lwork
        return
    end if

    ! Local Memory Allocation
    if (present(ipvt)) then
        if (size(ipvt) < n) then
            call report_array_size_error("solve_least_squares_vec_pvt", &
                errmgr, "ipvt", n, size(ipvt))
            return
        end if
        iptr => ipvt(1:n)
    else
        allocate(iwrk(n), stat = istat)
        if (istat /= 0) then
            call report_memory_error("solve_least_squares_vec_pvt", errmgr, &
                istat)
            return
        end if
        iptr => iwrk
        iptr = 0
    end if

    if (present(work)) then
        if (size(work) < lwork) then
            call report_array_size_error("solve_least_squares_vec_pvt", &
                errmgr, "work", lwork, size(work))
            return
        end if
        wptr => work(1:lwork)
    else
        allocate(wrk(lwork), stat = istat)
        if (istat /= 0) then
            call report_memory_error("solve_least_squares_vec_pvt", errmgr, &
                istat)
            return
        end if
        wptr => wrk
    end if

    ! Process
    call DGELSY(m, n, 1, a, m, b, maxmn, iptr, rc, rnk, wptr, lwork, flag)
    if (present(arnk)) arnk = rnk
end subroutine

! ------------------------------------------------------------------------------
subroutine solve_least_squares_vec_pvt_cmplx(a, b, ipvt, arnk, &
    work, olwork, rwork, err)
    !! Solves the system of equations \(A \vec{x} = \vec{b}\) using a full 
    !! orthogonal factorization of \(A\).
    complex(real64), intent(inout), dimension(:,:) :: a
        !! On input, the M-by-N matrix \(A\).  On output, the matrix is 
        !! overwritten by its orthogonal factorization.
    complex(real64), intent(inout), dimension(:) :: b
        !! If the system is overdetermined, the M-element vector \(\vec{b}\);
        !! else, the array should be sized as N-element with the first M
        !! elements containing \(\vec{b}\).  On output, the first N rows will
        !! contain the solution vector \(\vec{x}\).
    integer(int32), intent(inout), target, optional, dimension(:) :: ipvt
        !! An optional input that on input, an N-element array that if 
        !! IPVT(I) .ne. 0, the I-th column of A is permuted to the front
        !! of A * P; if IPVT(I) = 0, the I-th column of A is a free column.  On
        !! output, if IPVT(I) = K, then the I-th column of A * P was the K-th
        !! column of A.  If not supplied, memory is allocated internally, and 
        !! IPVT is set to all zeros such that all columns are treated as free.
    integer(int32), intent(out), optional :: arnk
        !! An optional output, that if provided, will return the rank of \(A\).
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
        !! must be at least 2 * N.
    class(errors), intent(inout), optional, target :: err
        !! The error object to be updated.

    ! Local Variables
    integer(int32) :: m, n, maxmn, lwork, lrwork, istat, flag, rnk
    complex(real64), pointer, dimension(:) :: wptr
    complex(real64), allocatable, target, dimension(:) :: wrk
    real(real64), pointer, dimension(:) :: rwptr
    real(real64), allocatable, target, dimension(:) :: rwrk
    integer(int32), allocatable, target, dimension(:) :: iwrk
    integer(int32), pointer, dimension(:) :: iptr
    complex(real64), dimension(1) :: temp
    real(real64), dimension(1) :: rtemp
    integer(int32), dimension(1) :: itemp
    real(real64) :: rc
    class(errors), pointer :: errmgr
    type(errors), target :: deferr

    ! Initialization
    m = size(a, 1)
    n = size(a, 2)
    maxmn = max(m, n)
    lrwork = 2 * n
    rc = epsilon(rc)
    if (present(arnk)) arnk = 0
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if

    ! Input Check
    if (size(b, 1) /= maxmn) then
        call report_array_size_error("solve_least_squares_vec_pvt_cmplx", &
            errmgr, "b", maxmn, size(b))
        return
    end if

    ! Workspace Query
    call ZGELSY(m, n, 1, a, m, b, maxmn, itemp, rc, rnk, temp, -1, rtemp, &
        flag)
    lwork = int(temp(1), int32)
    if (present(olwork)) then
        olwork = lwork
        return
    end if

    ! Local Memory Allocation
    if (present(ipvt)) then
        if (size(ipvt) < n) then
            call report_array_size_error("solve_least_squares_vec_pvt_cmplx", &
                errmgr, "ipvt", n, size(ipvt))
            return
        end if
        iptr => ipvt(1:n)
    else
        allocate(iwrk(n), stat = istat)
        if (istat /= 0) then
            call report_memory_error("solve_least_squares_vec_pvt_cmplx", &
                errmgr, istat)
            return
        end if
        iptr => iwrk
        iptr = 0
    end if

    if (present(work)) then
        if (size(work) < lwork) then
            call report_array_size_error("solve_least_squares_vec_pvt_cmplx", &
                errmgr, "work", lwork, size(work))
            return
        end if
        wptr => work(1:lwork)
    else
        allocate(wrk(lwork), stat = istat)
        if (istat /= 0) then
            call report_memory_error("solve_least_squares_vec_pvt_cmplx", &
                errmgr, istat)
            return
        end if
        wptr => wrk
    end if

    if (present(rwork)) then
        if (size(rwork) < lrwork) then
            call report_array_size_error("solve_least_squares_vec_pvt_cmplx", &
                errmgr, "rwork", lrwork, size(rwork))
            return
        end if
        rwptr => rwork(1:lrwork)
    else
        allocate(rwrk(lrwork), stat = istat)
        if (istat /= 0) then
            call report_memory_error("solve_least_squares_vec_pvt_cmplx", &
                errmgr, istat)
            return
        end if
        rwptr => rwrk
    end if

    ! Process
    call ZGELSY(m, n, 1, a, m, b, maxmn, iptr, rc, rnk, wptr, lwork, &
        rwptr, flag)
    if (present(arnk)) arnk = rnk
end subroutine

! ------------------------------------------------------------------------------
subroutine solve_least_squares_mtx_svd(a, b, s, arnk, work, olwork, err)
    !! Solves the system of equations \(A X = B\) using a singular value
    !! decomposition of \(A\).
    real(real64), intent(inout), dimension(:,:) :: a
        !! On input, the M-by-N matrix \(A\).  On output, the matrix is 
        !! overwritten.
    real(real64), intent(inout), dimension(:,:) :: b
        !! If the system is overdetermined, the M-by-NRHS matrix \(B\); else,
        !! the matrix should be sized as N-by-NRHS with the first M rows 
        !! containing \(B\).  On output, the first N rows will contain the
        !! solution matrix \(X\).
    real(real64), intent(out), target, optional, dimension(:) :: s
        !! An optional MIN(M, N)-element array that on output contains the 
        !! singular values of \(A\) in descending order.  Notice, the condition
        !! number of \(A\) can be determined by S(1) / S(MIN(M, N)).
    integer(int32), intent(out), optional :: arnk
        !! An optional output, that if provided, will return the rank of \(A\).
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
    integer(int32) :: m, n, nrhs, mn, maxmn, istat, flag, lwork, rnk
    real(real64), pointer, dimension(:) :: wptr, sptr
    real(real64), allocatable, target, dimension(:) :: wrk, sing
    real(real64), dimension(1) :: temp
    real(real64) :: rcond
    class(errors), pointer :: errmgr
    type(errors), target :: deferr

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
        call report_matrix_size_error("solve_least_squares_mtx_svd", errmgr, &
            "b", maxmn, nrhs, size(b, 1), size(b, 2))
        return
    end if

    ! Workspace Query
    call DGELSS(m, n, nrhs, a, m, b, maxmn, temp, rcond, rnk, temp, -1, &
        flag)
    lwork = int(temp(1), int32)
    if (present(olwork)) then
        olwork = lwork
        return
    end if

    ! Local Memory Allocation
    if (present(s)) then
        if (size(s) < mn) then
            call report_array_size_error("solve_least_squares_mtx_svd", &
                errmgr, "s", mn, size(s))
            return
        end if
        sptr => s(1:mn)
    else
        allocate(sing(mn), stat = istat)
        if (istat /= 0) then
            call report_memory_error("solve_least_squares_mtx_svd", errmgr, &
                istat)
            return
        end if
        sptr => sing
    end if

    if (present(work)) then
        if (size(work) < lwork) then
            call report_array_size_error("solve_least_squares_mtx_svd", &
                errmgr, "work", lwork, size(work))
            return
        end if
        wptr => work(1:lwork)
    else
        allocate(wrk(lwork), stat = istat)
        if (istat /= 0) then
            call report_memory_error("solve_least_squares_mtx_svd", errmgr, &
                istat)
            return
        end if
        wptr => wrk
    end if

    ! Process
    call DGELSS(m, n, nrhs, a, m, b, maxmn, sptr, rcond, rnk, wptr, lwork, &
        flag)
    if (present(arnk)) arnk = rnk
    if (flag > 0) then
        call errmgr%report_error("solve_least_squares_mtx_svd", &
            "The QR iteration process could not converge.", &
            LA_CONVERGENCE_ERROR)
        return
    end if
end subroutine

! ------------------------------------------------------------------------------
subroutine solve_least_squares_mtx_svd_cmplx(a, b, s, arnk, work, &
    olwork, rwork, err)
    !! Solves the system of equations \(A X = B\) using a singular value
    !! decomposition of \(A\).
    complex(real64), intent(inout), dimension(:,:) :: a
        !! On input, the M-by-N matrix \(A\).  On output, the matrix is 
        !! overwritten.
    complex(real64), intent(inout), dimension(:,:) :: b
        !! If the system is overdetermined, the M-by-NRHS matrix \(B\); else,
        !! the matrix should be sized as N-by-NRHS with the first M rows 
        !! containing \(B\).  On output, the first N rows will contain the
        !! solution matrix \(X\).
    real(real64), intent(out), target, optional, dimension(:) :: s
        !! An optional MIN(M, N)-element array that on output contains the 
        !! singular values of \(A\) in descending order.  Notice, the condition
        !! number of \(A\) can be determined by S(1) / S(MIN(M, N)).
    integer(int32), intent(out), optional :: arnk
        !! An optional output, that if provided, will return the rank of \(A\).
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
        !! allocation for real-valued workspaces.  If not provided, the 
        !! memory required is allocated within.  If provided, the length of the 
        !! array must be at least 5 * MIN(M, N).
    class(errors), intent(inout), optional, target :: err
        !! The error object to be updated.

    ! Local Variables
    integer(int32) :: m, n, nrhs, mn, maxmn, istat, flag, lwork, rnk, lrwork
    complex(real64), pointer, dimension(:) :: wptr
    complex(real64), allocatable, target, dimension(:) :: wrk
    real(real64), pointer, dimension(:) :: rwptr, sptr
    real(real64), allocatable, target, dimension(:) :: rwrk, sing
    complex(real64), dimension(1) :: temp
    real(real64), dimension(1) :: rtemp
    real(real64) :: rcond
    class(errors), pointer :: errmgr
    type(errors), target :: deferr

    ! Initialization
    m = size(a, 1)
    n = size(a, 2)
    nrhs = size(b, 2)
    mn = min(m, n)
    lrwork = 5 * mn
    maxmn = max(m, n)
    rcond = epsilon(rcond)
    if (present(arnk)) arnk = 0
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if

    ! Input Check
    if (size(b, 1) /= maxmn) then
        call report_matrix_size_error("solve_least_squares_mtx_svd_cmplx", &
            errmgr, "b", maxmn, nrhs, size(b, 1), size(b, 2))
        return
    end if

    ! Workspace Query
    call ZGELSS(m, n, nrhs, a, m, b, maxmn, rtemp, rcond, rnk, temp, -1, &
        rtemp, flag)
    lwork = int(temp(1), int32)
    if (present(olwork)) then
        olwork = lwork
        return
    end if

    ! Local Memory Allocation
    if (present(s)) then
        if (size(s) < mn) then
            call report_array_size_error("solve_least_squares_mtx_svd_cmplx", &
                errmgr, "s", mn, size(s))
            return
        end if
        sptr => s(1:mn)
    else
        allocate(sing(mn), stat = istat)
        if (istat /= 0) then
            call report_memory_error("solve_least_squares_mtx_svd_cmplx", &
                errmgr, istat)
            return
        end if
        sptr => sing
    end if

    if (present(work)) then
        if (size(work) < lwork) then
            call report_array_size_error("solve_least_squares_mtx_svd_cmplx", &
                errmgr, "work", lwork, size(work))
            return
        end if
        wptr => work(1:lwork)
    else
        allocate(wrk(lwork), stat = istat)
        if (istat /= 0) then
            call report_memory_error("solve_least_squares_mtx_svd_cmplx", &
                errmgr, istat)
            return
        end if
        wptr => wrk
    end if

    if (present(rwork)) then
        if (size(rwork) < lrwork) then
            call report_array_size_error("solve_least_squares_mtx_svd_cmplx", &
                errmgr, "rwork", lrwork, size(rwork))
            return
        end if
        rwptr => rwork(1:lrwork)
    else
        allocate(rwrk(lrwork), stat = istat)
        if (istat /= 0) then
            call report_memory_error("solve_least_squares_mtx_svd_cmplx", &
                errmgr, istat)
            return
        end if
        rwptr => rwrk
    end if

    ! Process
    call ZGELSS(m, n, nrhs, a, m, b, maxmn, sptr, rcond, rnk, wptr, lwork, &
        rwptr, flag)
    if (present(arnk)) arnk = rnk
    if (flag > 0) then
        call errmgr%report_error("solve_least_squares_mtx_svd_cmplx", &
            "The QR iteration process could not converge.", &
            LA_CONVERGENCE_ERROR)
        return
    end if
end subroutine

! ------------------------------------------------------------------------------
subroutine solve_least_squares_vec_svd(a, b, s, arnk, work, olwork, err)
    !! Solves the system of equations \(A \vec{x} = \vec{b}\) using a singular 
    !! value decomposition of \(A\).
    real(real64), intent(inout), dimension(:,:) :: a
        !! On input, the M-by-N matrix \(A\).  On output, the matrix is 
        !! overwritten.
    real(real64), intent(inout), dimension(:) :: b
        !! If the system is overdetermined, the M-element vector \(\vec{b}\);
        !! else, the array should be sized as N-element with the first M
        !! elements containing \(\vec{b}\).  On output, the first N rows will
        !! contain the solution vector \(\vec{x}\).
    real(real64), intent(out), target, optional, dimension(:) :: s
        !! An optional MIN(M, N)-element array that on output contains the 
        !! singular values of \(A\) in descending order.  Notice, the condition
        !! number of \(A\) can be determined by S(1) / S(MIN(M, N)).
    integer(int32), intent(out), optional :: arnk
        !! An optional output, that if provided, will return the rank of \(A\).
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
    integer(int32) :: m, n, mn, maxmn, istat, flag, lwork, rnk
    real(real64), pointer, dimension(:) :: wptr, sptr
    real(real64), allocatable, target, dimension(:) :: wrk, sing
    real(real64), dimension(1) :: temp
    real(real64) :: rcond
    class(errors), pointer :: errmgr
    type(errors), target :: deferr

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
    if (size(b) /= maxmn) then
        call report_array_size_error("solve_least_squares_vec_svd", errmgr, &
            "b", maxmn, size(b))
        return
    end if

    ! Workspace Query
    call DGELSS(m, n, 1, a, m, b, maxmn, temp, rcond, rnk, temp, -1, flag)
    lwork = int(temp(1), int32)
    if (present(olwork)) then
        olwork = lwork
        return
    end if

    ! Local Memory Allocation
    if (present(s)) then
        if (size(s) < mn) then
            call report_array_size_error("solve_least_squares_vec_svd", &
                errmgr, "s", mn, size(s))
            return
        end if
        sptr => s(1:mn)
    else
        allocate(sing(mn), stat = istat)
        if (istat /= 0) then
            call report_memory_error("solve_least_squares_vec_svd", errmgr, &
                istat)
            return
        end if
        sptr => sing
    end if

    if (present(work)) then
        if (size(work) < lwork) then
            call report_array_size_error("solve_least_squares_vec_svd", &
                errmgr, "work", lwork, size(work))
            return
        end if
        wptr => work(1:lwork)
    else
        allocate(wrk(lwork), stat = istat)
        if (istat /= 0) then
            call report_memory_error("solve_least_squares_vec_svd", errmgr, &
                istat)
            return
        end if
        wptr => wrk
    end if

    ! Process
    call DGELSS(m, n, 1, a, m, b, maxmn, sptr, rcond, rnk, wptr, lwork, &
        flag)
    if (present(arnk)) arnk = rnk
    if (flag > 0) then
        call errmgr%report_error("solve_least_squares_vec_svd", &
            "The QR iteration process could not converge.", &
            LA_CONVERGENCE_ERROR)
        return
    end if
end subroutine

! ------------------------------------------------------------------------------
subroutine solve_least_squares_vec_svd_cmplx(a, b, s, arnk, work, &
    olwork, rwork, err)
    !! Solves the system of equations \(A \vec{x} = \vec{b}\) using a singular 
    !! value decomposition of \(A\).
    complex(real64), intent(inout), dimension(:,:) :: a
        !! On input, the M-by-N matrix \(A\).  On output, the matrix is 
        !! overwritten.
    complex(real64), intent(inout), dimension(:) :: b
        !! If the system is overdetermined, the M-element vector \(\vec{b}\);
        !! else, the array should be sized as N-element with the first M
        !! elements containing \(\vec{b}\).  On output, the first N rows will
        !! contain the solution vector \(\vec{x}\).
    real(real64), intent(out), target, optional, dimension(:) :: s
        !! An optional MIN(M, N)-element array that on output contains the 
        !! singular values of \(A\) in descending order.  Notice, the condition
        !! number of \(A\) can be determined by S(1) / S(MIN(M, N)).
    integer(int32), intent(out), optional :: arnk
        !! An optional output, that if provided, will return the rank of \(A\).
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
        !! allocation for real-valued workspaces.  If not provided, the 
        !! memory required is allocated within.  If provided, the length of the 
        !! array must be at least 5 * MIN(M, N).
    class(errors), intent(inout), optional, target :: err
        !! The error object to be updated.

    ! Local Variables
    integer(int32) :: m, n, mn, maxmn, istat, flag, lwork, rnk, lrwork
    real(real64), pointer, dimension(:) :: rwptr, sptr
    real(real64), allocatable, target, dimension(:) :: rwrk, sing
    complex(real64), pointer, dimension(:) :: wptr
    complex(real64), allocatable, target, dimension(:) :: wrk
    complex(real64), dimension(1) :: temp
    real(real64), dimension(1) :: rtemp
    real(real64) :: rcond
    class(errors), pointer :: errmgr
    type(errors), target :: deferr

    ! Initialization
    m = size(a, 1)
    n = size(a, 2)
    mn = min(m, n)
    lrwork = 5 * mn
    maxmn = max(m, n)
    rcond = epsilon(rcond)
    if (present(arnk)) arnk = 0
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if

    ! Input Check
    if (size(b) /= maxmn) then
        call report_array_size_error("solve_least_squares_vec_svd_cmplx", &
            errmgr, "b", maxmn, size(b))
        return
    end if

    ! Workspace Query
    call ZGELSS(m, n, 1, a, m, b, maxmn, rtemp, rcond, rnk, temp, -1, &
        rtemp, flag)
    lwork = int(temp(1), int32)
    if (present(olwork)) then
        olwork = lwork
        return
    end if

    ! Local Memory Allocation
    if (present(s)) then
        if (size(s) < mn) then
            call report_array_size_error("solve_least_squares_vec_svd_cmplx", &
                errmgr, "s", mn, size(s))
            return
        end if
        sptr => s(1:mn)
    else
        allocate(sing(mn), stat = istat)
        if (istat /= 0) then
            call report_memory_error("solve_least_squares_vec_svd_cmplx", &
                errmgr, istat)
            return
        end if
        sptr => sing
    end if

    if (present(work)) then
        if (size(work) < lwork) then
            call report_array_size_error("solve_least_squares_vec_svd_cmplx", &
                errmgr, "work", lwork, size(work))
            return
        end if
        wptr => work(1:lwork)
    else
        allocate(wrk(lwork), stat = istat)
        if (istat /= 0) then
            call report_memory_error("solve_least_squares_vec_svd_cmplx", &
                errmgr, istat)
            return
        end if
        wptr => wrk
    end if

    if (present(rwork)) then
        if (size(rwork) < lrwork) then
            call report_array_size_error("solve_least_squares_vec_svd_cmplx", &
                errmgr, "rwork", lrwork, size(rwork))
            return
        end if
        rwptr => rwork(1:lrwork)
    else
        allocate(rwrk(lrwork), stat = istat)
        if (istat /= 0) then
            call report_memory_error("solve_least_squares_vec_svd_cmplx", &
                errmgr, istat)
            return
        end if
        rwptr => rwrk
    end if

    ! Process
    call ZGELSS(m, n, 1, a, m, b, maxmn, sptr, rcond, rnk, wptr, lwork, &
        rwptr, flag)
    if (present(arnk)) arnk = rnk
    if (flag > 0) then
        call errmgr%report_error("solve_least_squares_vec_svd_cmplx", &
            "The QR iteration process could not converge.", &
            LA_CONVERGENCE_ERROR)
        return
    end if
end subroutine

! ------------------------------------------------------------------------------
end module