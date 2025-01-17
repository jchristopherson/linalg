module linalg_qr
    use iso_fortran_env
    use linalg_errors
    use linalg_rz
    use linalg_tri
    use lapack
    use blas
    use ferror
    use qrupdate
    implicit none
    private
    public :: qr_factor
    public :: form_qr
    public :: mult_qr
    public :: qr_rank1_update
    public :: solve_qr

    interface qr_factor
        module procedure :: qr_factor_no_pivot
        module procedure :: qr_factor_no_pivot_cmplx
        module procedure :: qr_factor_pivot
        module procedure :: qr_factor_pivot_cmplx
    end interface

    interface form_qr
        module procedure :: form_qr_no_pivot
        module procedure :: form_qr_no_pivot_cmplx
        module procedure :: form_qr_pivot
        module procedure :: form_qr_pivot_cmplx
    end interface

    interface mult_qr
        module procedure :: mult_qr_mtx
        module procedure :: mult_qr_mtx_cmplx
        module procedure :: mult_qr_vec
        module procedure :: mult_qr_vec_cmplx
    end interface

    interface qr_rank1_update
        module procedure :: qr_rank1_update_dbl
        module procedure :: qr_rank1_update_cmplx
    end interface

    interface solve_qr
        module procedure :: solve_qr_no_pivot_mtx
        module procedure :: solve_qr_no_pivot_mtx_cmplx
        module procedure :: solve_qr_no_pivot_vec
        module procedure :: solve_qr_no_pivot_vec_cmplx
        module procedure :: solve_qr_pivot_mtx
        module procedure :: solve_qr_pivot_mtx_cmplx
        module procedure :: solve_qr_pivot_vec
        module procedure :: solve_qr_pivot_vec_cmplx
    end interface
contains
! ------------------------------------------------------------------------------
subroutine qr_factor_no_pivot(a, tau, work, olwork, err)
    !! Computes the QR factorization of an M-by-N matrix.
    real(real64), intent(inout), dimension(:,:) :: a
        !! On input, the M-by-N matrix to factor.  On output, the elements on 
        !! and above the diagonal contain the MIN(M, N)-by-N upper trapezoidal 
        !! matrix R (R is upper triangular if M >= N).  The elements below the 
        !! diagonal, along with the array tau, represent the orthogonal matrix
        !! Q as a product of elementary reflectors.
    real(real64), intent(out), dimension(:) :: tau
        !! A MIN(M, N)-element array used to store the scalar factors of the 
        !! elementary reflectors.
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

    ! Local Variables
    integer(int32) :: m, n, mn, istat, lwork, flag
    real(real64), dimension(1) :: temp
    real(real64), pointer, dimension(:) :: wptr
    real(real64), allocatable, target, dimension(:) :: wrk
    class(errors), pointer :: errmgr
    type(errors), target :: deferr

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
    if (size(tau) /= mn) then
        ! ERROR: TAU not sized correctly
        call report_array_size_error("qr_factor_no_pivot", errmgr, "tau", mn, &
            size(tau))
        return
    end if

    ! Workspace Query
    call DGEQRF(m, n, a, m, tau, temp, -1, flag)
    lwork = int(temp(1), int32)
    if (present(olwork)) then
        olwork = lwork
        return
    end if

    ! Local Memory Allocation
    if (present(work)) then
        if (size(work) < lwork) then
            ! ERROR: WORK not sized correctly
            call report_array_size_error("qr_factor_no_pivot", errmgr, "work", &
                lwork, size(work))
            return
        end if
        wptr => work(1:lwork)
    else
        allocate(wrk(lwork), stat = istat)
        if (istat /= 0) then
            call report_memory_error("qr_factor_no_pivot", errmgr, istat)
            return
        end if
        wptr => wrk
    end if

    ! Call DGEQRF
    call DGEQRF(m, n, a, m, tau, wptr, lwork, flag)
end subroutine

! ------------------------------------------------------------------------------
subroutine qr_factor_no_pivot_cmplx(a, tau, work, olwork, err)
    !! Computes the QR factorization of an M-by-N matrix.
    complex(real64), intent(inout), dimension(:,:) :: a
        !! On input, the M-by-N matrix to factor.  On output, the elements on 
        !! and above the diagonal contain the MIN(M, N)-by-N upper trapezoidal 
        !! matrix R (R is upper triangular if M >= N).  The elements below the 
        !! diagonal, along with the array tau, represent the orthogonal matrix
        !! Q as a product of elementary reflectors.
    complex(real64), intent(out), dimension(:) :: tau
        !! A MIN(M, N)-element array used to store the scalar factors of the 
        !! elementary reflectors.
    complex(real64), intent(out), target, dimension(:), optional :: work
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
    integer(int32) :: m, n, mn, istat, lwork, flag
    complex(real64), dimension(1) :: temp
    complex(real64), pointer, dimension(:) :: wptr
    complex(real64), allocatable, target, dimension(:) :: wrk
    class(errors), pointer :: errmgr
    type(errors), target :: deferr

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
    if (size(tau) /= mn) then
        call report_array_size_error("qr_factor_no_pivot_cmplx", errmgr, &
            "tau", mn, size(tau))
        return
    end if

    ! Workspace Query
    call ZGEQRF(m, n, a, m, tau, temp, -1, flag)
    lwork = int(temp(1), int32)
    if (present(olwork)) then
        olwork = lwork
        return
    end if

    ! Local Memory Allocation
    if (present(work)) then
        if (size(work) < lwork) then
            call report_array_size_error("qr_factor_no_pivot_cmplx", errmgr, &
                "tau", lwork, size(work))
            return
        end if
        wptr => work(1:lwork)
    else
        allocate(wrk(lwork), stat = istat)
        if (istat /= 0) then
            call report_memory_error("qr_factor_no_pivot_cmplx", errmgr, istat)
            return
        end if
        wptr => wrk
    end if

    ! Call ZGEQRF
    call ZGEQRF(m, n, a, m, tau, wptr, lwork, flag)
end subroutine

! ------------------------------------------------------------------------------
subroutine qr_factor_pivot(a, tau, jpvt, work, olwork, err)
    !! Computes the QR factorization of an M-by-N matrix using column pivoting
    !! such that \(A P = Q R\).
    real(real64), intent(inout), dimension(:,:) :: a
        !! On input, the M-by-N matrix to factor.  On output, the elements on 
        !! and above the diagonal contain the MIN(M, N)-by-N upper trapezoidal 
        !! matrix R (R is upper triangular if M >= N).  The elements below the 
        !! diagonal, along with the array tau, represent the orthogonal matrix
        !! Q as a product of elementary reflectors.
    real(real64), intent(out), dimension(:) :: tau
        !! A MIN(M, N)-element array used to store the scalar factors of the 
        !! elementary reflectors.
    integer(int32), intent(inout), dimension(:) :: jpvt
        !! On input, an N-element array that if JPVT(I) .ne. 0, the I-th column 
        !! of A is permuted to the front of A * P; if JPVT(I) = 0, the I-th 
        !! column of A is a free column.  On output, if JPVT(I) = K, then the 
        !! I-th column of A * P was the K-th column of A.
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

    ! Local Variables
    integer(int32) :: m, n, mn, istat, lwork, flag
    real(real64), dimension(1) :: temp
    real(real64), pointer, dimension(:) :: wptr
    real(real64), allocatable, target, dimension(:) :: wrk
    class(errors), pointer :: errmgr
    type(errors), target :: deferr

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
    if (size(tau) /= mn) then
        call report_array_size_error("qr_factor_pivot", errmgr, "tau", mn, &
            size(tau))
        return
    else if (size(jpvt) /= n) then
        call report_array_size_error("qr_factor_pivot", errmgr, "jpvt", n, &
            size(jpvt))
        return
    end if

    ! Workspace Query
    call DGEQP3(m, n, a, m, jpvt, tau, temp, -1, flag)
    lwork = int(temp(1), int32)
    if (present(olwork)) then
        olwork = lwork
        return
    end if

    ! Local Memory Allocation
    if (present(work)) then
        if (size(work) < lwork) then
            call report_array_size_error("qr_factor_pivot", errmgr, "work", &
                lwork, size(work))
            return
        end if
        wptr => work(1:lwork)
    else
        allocate(wrk(lwork), stat = istat)
        if (istat /= 0) then
            call report_memory_error("qr_factor_pivot", errmgr, istat)
        end if
        wptr => wrk
    end if

    ! Call DGEQP3
    call DGEQP3(m, n, a, m, jpvt, tau, wptr, lwork, flag)
end subroutine

! ------------------------------------------------------------------------------
subroutine qr_factor_pivot_cmplx(a, tau, jpvt, work, olwork, rwork, err)
    !! Computes the QR factorization of an M-by-N matrix using column pivoting
    !! such that \(A P = Q R\).
    complex(real64), intent(inout), dimension(:,:) :: a
        !! On input, the M-by-N matrix to factor.  On output, the elements on 
        !! and above the diagonal contain the MIN(M, N)-by-N upper trapezoidal 
        !! matrix R (R is upper triangular if M >= N).  The elements below the 
        !! diagonal, along with the array tau, represent the orthogonal matrix
        !! Q as a product of elementary reflectors.
    complex(real64), intent(out), dimension(:) :: tau
        !! A MIN(M, N)-element array used to store the scalar factors of the 
        !! elementary reflectors.
    integer(int32), intent(inout), dimension(:) :: jpvt
        !! On input, an N-element array that if JPVT(I) .ne. 0, the I-th column 
        !! of A is permuted to the front of A * P; if JPVT(I) = 0, the I-th 
        !! column of A is a free column.  On output, if JPVT(I) = K, then the 
        !! I-th column of A * P was the K-th column of A.
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
        !! An optional input, that if provided, prevents any local allocate of 
        !! real-valued memory.  If not provided, the memory required is 
        !! allocated within.  If provided, the length of the array must be at
        !! least 2*N.
    class(errors), intent(inout), optional, target :: err
        !! The error object to be updated.

    ! Local Variables
    integer(int32) :: m, n, mn, istat, lwork, flag
    complex(real64), dimension(1) :: temp
    complex(real64), pointer, dimension(:) :: wptr
    complex(real64), allocatable, target, dimension(:) :: wrk
    real(real64), pointer, dimension(:) :: rptr
    real(real64), allocatable, target, dimension(:) :: rwrk
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    character(len = :), allocatable :: errmsg

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
    if (size(tau) /= mn) then
        call report_array_size_error("qr_factor_pivot_cmplx", errmgr, &
            "tau", mn, size(tau))
        return
    else if (size(jpvt) /= n) then
        call report_array_size_error("qr_factor_pivot_cmplx", errmgr, "jpvt", &
            n, size(jpvt))
        return
    end if
    if (present(rwork)) then
        if (size(rwork) < 2 * n) then
            call report_array_size_error("qr_factor_pivot_cmplx", errmgr, &
                "rwork", 2 * n, size(rwork))
            return
        end if
        rptr => rwork(1:2*n)
    else
        allocate(rwrk(2 * n), stat = flag)
        if (flag /= 0) then
            call report_memory_error("qr_factor_pivot_cmplx", errmgr, flag)
            return
        end if
        rptr => rwrk
    end if

    ! Workspace Query
    call ZGEQP3(m, n, a, m, jpvt, tau, temp, -1, rptr, flag)
    lwork = int(temp(1), int32)
    if (present(olwork)) then
        olwork = lwork
        return
    end if

    ! Local Memory Allocation
    if (present(work)) then
        if (size(work) < lwork) then
            call report_array_size_error("qr_factor_pivot_cmplx", errmgr, &
                "work", lwork, size(work))
            return
        end if
        wptr => work(1:lwork)
    else
        allocate(wrk(lwork), stat = istat)
        if (istat /= 0) then
            call report_memory_error("qr_factor_pivot_cmplx", errmgr, istat)
            return
        end if
        wptr => wrk
    end if

    ! Call ZGEQP3
    call ZGEQP3(m, n, a, m, jpvt, tau, wptr, lwork, rptr, flag)
end subroutine

! ------------------------------------------------------------------------------
subroutine form_qr_no_pivot(r, tau, q, work, olwork, err)
    !! Forms the full M-by-M orthogonal matrix \(Q\) from the elementary
    !! reflectors returned by the base QR factorization algorithm.
    real(real64), intent(inout), dimension(:,:) :: r
        !! On input, an M-by-N matrix where the elements below the diagonal 
        !! contain the elementary reflectors generated from the QR 
        !! factorization.  On and above the diagonal, the matrix contains the
        !! matrix \(R\).  On output, the elements below the diagonal are zeroed 
        !! such that the remaining matrix is simply the M-by-N matrix \(R\).
    real(real64), intent(in), dimension(:) :: tau
        !! A MIN(M, N)-element array containing the scalar factors of each 
        !! elementary reflector defined in \(R\).
    real(real64), intent(out), dimension(:,:) :: q
        !! An M-by-M matrix where the full orthogonal matrix \(Q\) will be
        !! written.  In the event that M > N, \(Q\) may be supplied as M-by-N, 
        !! and therefore only return the useful submatrix \(Q_1\)
        !! \(Q = [Q_1 Q_2]\) as the factorization can be written as 
        !! \(Q R = [Q_1, Q_2] [R1 0]^T\).
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

    ! Local Variables
    integer(int32) :: j, m, n, mn, qcol, istat, flag, lwork
    real(real64), pointer, dimension(:) :: wptr
    real(real64), allocatable, target, dimension(:) :: wrk
    real(real64), dimension(1) :: temp
    class(errors), pointer :: errmgr
    type(errors), target :: deferr

    ! Initialization
    m = size(r, 1)
    n = size(r, 2)
    mn = min(m, n)
    qcol = size(q, 2)
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if

    ! Input Check
    if (size(tau) /= mn) then
        call report_array_size_error("form_qr_no_pivot", errmgr, "tau", &
            mn, size(tau))
        return
    else if (size(q, 1) /= m .or. (qcol /= m .and. qcol /= n)) then
        call report_matrix_size_error("form_qr_no_pivot", errmgr, "q", m, mn, &
            size(q, 1), size(q, 2))
        return
    else if (qcol == n .and. m < n) then
        call report_matrix_size_error("form_qr_no_pivot", errmgr, "q", m, m, &
            size(q, 1), size(q, 2))
        return
    end if

    ! Workspace Query
    call DORGQR(m, qcol, mn, q, m, tau, temp, -1, flag)
    lwork = int(temp(1), int32)
    if (present(olwork)) then
        olwork = lwork
        return
    end if

    ! Local Memory Allocation
    if (present(work)) then
        if (size(work) < lwork) then
            call report_array_size_error("form_qr_no_pivot", errmgr, "work", &
                lwork, size(work))
            return
        end if
        wptr => work(1:lwork)
    else
        allocate(wrk(lwork), stat = istat)
        if (istat /= 0) then
            call report_memory_error("form_qr_no_pivot", errmgr, istat)
            return
        end if
        wptr => wrk
    end if

    ! Copy the sub-diagonal portion of R to Q, and then zero out the
    ! sub-diagonal portion of R
    do j = 1, mn
        q(j+1:m,j) = r(j+1:m,j)
        r(j+1:m,j) = zero
    end do

    ! Build Q - Build M-by-M or M-by-N, but M-by-N only for M >= N
    call DORGQR(m, qcol, mn, q, m, tau, wptr, lwork, flag)
end subroutine

! ------------------------------------------------------------------------------
subroutine form_qr_no_pivot_cmplx(r, tau, q, work, olwork, err)
    !! Forms the full M-by-M orthogonal matrix \(Q\) from the elementary
    !! reflectors returned by the base QR factorization algorithm.
    complex(real64), intent(inout), dimension(:,:) :: r
        !! On input, an M-by-N matrix where the elements below the diagonal 
        !! contain the elementary reflectors generated from the QR 
        !! factorization.  On and above the diagonal, the matrix contains the
        !! matrix \(R\).  On output, the elements below the diagonal are zeroed 
        !! such that the remaining matrix is simply the M-by-N matrix \(R\).
    complex(real64), intent(in), dimension(:) :: tau
        !! A MIN(M, N)-element array containing the scalar factors of each 
        !! elementary reflector defined in \(R\).
    complex(real64), intent(out), dimension(:,:) :: q
        !! An M-by-M matrix where the full orthogonal matrix \(Q\) will be
        !! written.  In the event that M > N, \(Q\) may be supplied as M-by-N, 
        !! and therefore only return the useful submatrix \(Q_1\)
        !! \(Q = [Q_1 Q_2]\) as the factorization can be written as 
        !! \(Q R = [Q_1, Q_2] [R1 0]^T\).
    complex(real64), intent(out), target, dimension(:), optional :: work
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
    complex(real64), parameter :: zero = (0.0d0, 0.0d0)

    ! Local Variables
    integer(int32) :: j, m, n, mn, qcol, istat, flag, lwork
    complex(real64), pointer, dimension(:) :: wptr
    complex(real64), allocatable, target, dimension(:) :: wrk
    complex(real64), dimension(1) :: temp
    class(errors), pointer :: errmgr
    type(errors), target :: deferr

    ! Initialization
    m = size(r, 1)
    n = size(r, 2)
    mn = min(m, n)
    qcol = size(q, 2)
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if

    ! Input Check
    if (size(tau) /= mn) then
        call report_array_size_error("form_qr_no_pivot_cmplx", errmgr, "tau", &
            mn, size(tau))
        return
    else if (size(q, 1) /= m .or. (qcol /= m .and. qcol /= n)) then
        call report_matrix_size_error("form_qr_no_pivot_cmplx", errmgr, "q", &
            m, mn, size(q, 1), size(q, 2))
        return
    else if (qcol == n .and. m < n) then
        call report_matrix_size_error("form_qr_no_pivot_cmplx", errmgr, "q", &
            m, m, size(q, 1), size(q, 2))
        return
    end if

    ! Workspace Query
    call ZUNGQR(m, qcol, mn, q, m, tau, temp, -1, flag)
    lwork = int(temp(1), int32)
    if (present(olwork)) then
        olwork = lwork
        return
    end if

    ! Local Memory Allocation
    if (present(work)) then
        if (size(work) < lwork) then
            call report_array_size_error("form_qr_no_pivot_cmplx", errmgr, &
                "work", lwork, size(work))
            return
        end if
        wptr => work(1:lwork)
    else
        allocate(wrk(lwork), stat = istat)
        if (istat /= 0) then
            call report_memory_error("form_qr_no_pivot_cmplx", errmgr, istat)
            return
        end if
        wptr => wrk
    end if

    ! Copy the sub-diagonal portion of R to Q, and then zero out the
    ! sub-diagonal portion of R
    do j = 1, mn
        q(j+1:m,j) = r(j+1:m,j)
        r(j+1:m,j) = zero
    end do

    ! Build Q - Build M-by-M or M-by-N, but M-by-N only for M >= N
    call ZUNGQR(m, qcol, mn, q, m, tau, wptr, lwork, flag)
end subroutine

! ------------------------------------------------------------------------------
subroutine form_qr_pivot(r, tau, pvt, q, p, work, olwork, err)
    !! Forms the full M-by-M orthogonal matrix \(Q\) from the elementary
    !! reflectors returned by the base QR factorization algorithm.
    real(real64), intent(inout), dimension(:,:) :: r
        !! On input, an M-by-N matrix where the elements below the diagonal 
        !! contain the elementary reflectors generated from the QR 
        !! factorization.  On and above the diagonal, the matrix contains the
        !! matrix \(R\).  On output, the elements below the diagonal are zeroed 
        !! such that the remaining matrix is simply the M-by-N matrix \(R\).
    real(real64), intent(in), dimension(:) :: tau
        !! A MIN(M, N)-element array containing the scalar factors of each 
        !! elementary reflector defined in \(R\).
    integer(int32), intent(in), dimension(:) :: pvt
        !! An N-element column pivot array as returned by the QR factorization.
    real(real64), intent(out), dimension(:,:) :: q
        !! An M-by-M matrix where the full orthogonal matrix \(Q\) will be
        !! written.  In the event that M > N, \(Q\) may be supplied as M-by-N, 
        !! and therefore only return the useful submatrix \(Q_1\)
        !! \(Q = [Q_1 Q_2]\) as the factorization can be written as 
        !! \(Q R = [Q_1, Q_2] [R1 0]^T\).
    real(real64), intent(out), dimension(:,:) :: p
        !! An N-by-N matrix where the pivot matrix will be written.
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
    integer(int32) :: j, jp, m, n, mn
    class(errors), pointer :: errmgr
    type(errors), target :: deferr

    ! Initialization
    m = size(r, 1)
    n = size(r, 2)
    mn = min(m, n)
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if

    ! Input Check
    if (size(p, 1) /= n .or. size(p, 2) /= n) then
        call report_matrix_size_error("form_qr_pivot", errmgr, "p", n, n, &
            size(p, 1), size(p, 2))
        return
    end if

    ! Generate Q and R
    call form_qr_no_pivot(r, tau, q, work = work, olwork = olwork, &
        err = errmgr)
    if (present(olwork)) return ! Just a workspace query
    if (errmgr%has_error_occurred()) return

    ! Form P
    do j = 1, n
        jp = pvt(j)
        p(:,j) = zero
        p(jp,j) = one
    end do
end subroutine

! ------------------------------------------------------------------------------
module subroutine form_qr_pivot_cmplx(r, tau, pvt, q, p, work, olwork, err)
    !! Forms the full M-by-M orthogonal matrix \(Q\) from the elementary
    !! reflectors returned by the base QR factorization algorithm.
    complex(real64), intent(inout), dimension(:,:) :: r
        !! On input, an M-by-N matrix where the elements below the diagonal 
        !! contain the elementary reflectors generated from the QR 
        !! factorization.  On and above the diagonal, the matrix contains the
        !! matrix \(R\).  On output, the elements below the diagonal are zeroed 
        !! such that the remaining matrix is simply the M-by-N matrix \(R\).
    complex(real64), intent(in), dimension(:) :: tau
        !! A MIN(M, N)-element array containing the scalar factors of each 
        !! elementary reflector defined in \(R\).
    integer(int32), intent(in), dimension(:) :: pvt
        !! An N-element column pivot array as returned by the QR factorization.
    complex(real64), intent(out), dimension(:,:) :: q
        !! An M-by-M matrix where the full orthogonal matrix \(Q\) will be
        !! written.  In the event that M > N, \(Q\) may be supplied as M-by-N, 
        !! and therefore only return the useful submatrix \(Q_1\)
        !! \(Q = [Q_1 Q_2]\) as the factorization can be written as 
        !! \(Q R = [Q_1, Q_2] [R1 0]^T\).
    complex(real64), intent(out), dimension(:,:) :: p
        !! An N-by-N matrix where the pivot matrix will be written.
    complex(real64), intent(out), target, dimension(:), optional :: work
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
    complex(real64), parameter :: zero = (0.0d0, 0.0d0)
    complex(real64), parameter :: one = (1.0d0, 0.0d0)

    ! Local Variables
    integer(int32) :: j, jp, m, n, mn, flag
    class(errors), pointer :: errmgr
    type(errors), target :: deferr

    ! Initialization
    m = size(r, 1)
    n = size(r, 2)
    mn = min(m, n)
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if

    ! Input Check
    if (size(p, 1) /= n .or. size(p, 2) /= n) then
        call report_matrix_size_error("form_qr_pivot_cmplx", errmgr, "p", &
            n, n, size(p, 1), size(p, 2))
        return
    end if

    ! Generate Q and R
    call form_qr_no_pivot_cmplx(r, tau, q, work = work, olwork = olwork, &
        err = errmgr)
    if (present(olwork)) return ! Just a workspace query
    if (errmgr%has_error_occurred()) return

    ! Form P
    do j = 1, n
        jp = pvt(j)
        p(:,j) = zero
        p(jp,j) = one
    end do
end subroutine

! ------------------------------------------------------------------------------
subroutine mult_qr_mtx(lside, trans, a, tau, c, work, olwork, err)
    !! Multiplies a general matrix by the orthogonal matrix \(Q\) from a QR
    !! factorization such that \(C = op(Q) C\) or \(C = C op(Q)\).
    logical, intent(in) :: lside
        !! Set to true to apply \(Q\) or \(Q^T\) from the left; else, set to 
        !! false to apply \(Q\) or \(Q^T\) from the right.
    logical, intent(in) :: trans
        !! Set to true to apply \(Q^T\); else, set to false to apply \(Q\).
    real(real64), intent(inout), dimension(:,:) :: a
        !! On input, an LDA-by-K matrix containing the elementary reflectors 
        !! output from the QR factorization.  If lside is set to true, LDA = M, 
        !! and M >= K >= 0; else, if lside is set to false, LDA = N, and 
        !! N >= K >= 0.  Notice, the contents of this matrix are
        !! restored on exit.
    real(real64), intent(in), dimension(:) :: tau
        !! A K-element array containing the scalar factors of each elementary 
        !! reflector defined in\(A\).
    real(real64), intent(inout), dimension(:,:) :: c
        !! On input, the M-by-N matrix \(C\).  On output, the product of the 
        !! orthogonal matrix \(Q\) and the original matrix \(C\).
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
    real(real64), parameter :: one = 1.0d0

    ! Local Variables
    character :: side, t
    integer(int32) :: m, n, k, nrowa, istat, flag, lwork
    real(real64), pointer, dimension(:) :: wptr
    real(real64), allocatable, target, dimension(:) :: wrk
    real(real64), dimension(1) :: temp
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    character(len = :), allocatable :: errmsg

    ! Initialization
    m = size(c, 1)
    n = size(c, 2)
    k = size(tau)
    if (lside) then
        side = 'L'
        nrowa = m
    else
        side = 'R'
        nrowa = n
    end if
    if (trans) then
        t = 'T'
    else
        t = 'N'
    end if
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if

    ! Input Check
    if (lside) then
        ! A is M-by-K, M >= K >= 0
        if (size(a, 1) /= m .or. size(a, 2) < k) then
            call report_matrix_size_error("mult_qr_mtx", errmgr, "a", m, k, &
                size(a, 1), size(a, 2))
            return
        end if
    else
        ! A is N-by-K, N >= K >= 0
        if (size(a, 1) /= n .or. size(a, 2) < k) then
            call report_matrix_size_error("mult_qr_mtx", errmgr, "a", n, k, &
                size(a, 1), size(a, 2))
            return
        end if
    end if

    ! Workspace Query
    call DORMQR(side, t, m, n, k, a, nrowa, tau, c, m, temp, -1, flag)
    lwork = int(temp(1), int32)
    if (present(olwork)) then
        olwork = lwork
        return
    end if

    ! Local Memory Allocation
    if (present(work)) then
        if (size(work) < lwork) then
            call report_array_size_error("mult_qr_mtx", errmgr, "work", lwork, &
                size(work))
            return
        end if
        wptr => work(1:lwork)
    else
        allocate(wrk(lwork), stat = istat)
        if (istat /= 0) then
            call report_memory_error("mult_qr_mtx", errmgr, istat)
            return
        end if
        wptr => wrk
    end if

    ! Call DORMQR
    call DORMQR(side, t, m, n, k, a, nrowa, tau, c, m, wptr, lwork, flag)
end subroutine

! ------------------------------------------------------------------------------
subroutine mult_qr_mtx_cmplx(lside, trans, a, tau, c, work, olwork, err)
    !! Multiplies a general matrix by the orthogonal matrix \(Q\) from a QR
    !! factorization such that \(C = op(Q) C\) or \(C = C op(Q)\).
    logical, intent(in) :: lside
        !! Set to true to apply \(Q\) or \(Q^H\) from the left; else, set to 
        !! false to apply \(Q\) or \(Q^H\) from the right.
    logical, intent(in) :: trans
        !! Set to true to apply \(Q^H\); else, set to false to apply \(Q\).
    complex(real64), intent(inout), dimension(:,:) :: a
        !! On input, an LDA-by-K matrix containing the elementary reflectors 
        !! output from the QR factorization.  If lside is set to true, LDA = M, 
        !! and M >= K >= 0; else, if lside is set to false, LDA = N, and 
        !! N >= K >= 0.  Notice, the contents of this matrix are
        !! restored on exit.
    complex(real64), intent(in), dimension(:) :: tau
        !! A K-element array containing the scalar factors of each elementary 
        !! reflector defined in\(A\).
    complex(real64), intent(inout), dimension(:,:) :: c
        !! On input, the M-by-N matrix \(C\).  On output, the product of the 
        !! orthogonal matrix \(Q\) and the original matrix \(C\).
    complex(real64), intent(out), target, dimension(:), optional :: work
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
    complex(real64), parameter :: one = (1.0d0, 0.0d0)

    ! Local Variables
    character :: side, t
    integer(int32) :: m, n, k, nrowa, istat, flag, lwork
    complex(real64), pointer, dimension(:) :: wptr
    complex(real64), allocatable, target, dimension(:) :: wrk
    complex(real64), dimension(1) :: temp
    class(errors), pointer :: errmgr
    type(errors), target :: deferr

    ! Initialization
    m = size(c, 1)
    n = size(c, 2)
    k = size(tau)
    if (lside) then
        side = 'L'
        nrowa = m
    else
        side = 'R'
        nrowa = n
    end if
    if (trans) then
        t = 'C'
    else
        t = 'N'
    end if
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if

    ! Input Check
    if (lside) then
        ! A is M-by-K, M >= K >= 0
        if (size(a, 1) /= m .or. size(a, 2) < k) then
            call report_matrix_size_error("mult_qr_mtx_cmplx", errmgr, "a", &
                m, k, size(a, 1), size(a, 2))
            return
        end if
    else
        ! A is N-by-K, N >= K >= 0
        if (size(a, 1) /= n .or. size(a, 2) < k) then
            call report_matrix_size_error("mult_qr_mtx_cmplx", errmgr, "a", &
                n, k, size(a, 1), size(a, 2))
            return
        end if
    end if

    ! Workspace Query
    call ZUNMQR(side, t, m, n, k, a, nrowa, tau, c, m, temp, -1, flag)
    lwork = int(temp(1), int32)
    if (present(olwork)) then
        olwork = lwork
        return
    end if

    ! Local Memory Allocation
    if (present(work)) then
        if (size(work) < lwork) then
            ! ERROR: WORK not sized correctly
            call report_array_size_error("mult_qr_mtx_cmplx", errmgr, "work", &
                lwork, size(work))
            return
        end if
        wptr => work(1:lwork)
    else
        allocate(wrk(lwork), stat = istat)
        if (istat /= 0) then
            call report_memory_error("mult_qr_mtx_cmplx", errmgr, istat)
            return
        end if
        wptr => wrk
    end if

    ! Call ZUNMQR
    call ZUNMQR(side, t, m, n, k, a, nrowa, tau, c, m, wptr, lwork, flag)
end subroutine

! ------------------------------------------------------------------------------
subroutine mult_qr_vec(trans, a, tau, c, work, olwork, err)
    !! Multiplies a vector by the orthogonal matrix \(Q\) from a QR 
    !! factorization such that \(\vec{c} = op(Q) \vec{c}\).
    logical, intent(in) :: trans
        !! Set to true to apply \(Q^T\); else, set to false to apply \(Q\).
    real(real64), intent(inout), dimension(:,:) :: a
        !! On input, an M-by-K matrix containing the elementary reflectors 
        !! output from the QR factorization. Notice, the contents of this matrix
        !! are restored on exit.
    real(real64), intent(in), dimension(:) :: tau
        !! A K-element array containing the scalar factors of each elementary 
        !! reflector defined in\(A\).
    real(real64), intent(inout), dimension(:) :: c
        !! On input, the M-element vector \(\vec{c}\).  On output, the
        !! product of the orthogonal matrix \(Q\) and the original vector 
        !! \(\vec{c}\).
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
    real(real64), parameter :: one = 1.0d0

    ! Local Variables
    character :: side, t
    integer(int32) :: m, k, nrowa, istat, flag, lwork
    real(real64), pointer, dimension(:) :: wptr
    real(real64), allocatable, target, dimension(:) :: wrk
    real(real64), dimension(1) :: temp
    class(errors), pointer :: errmgr
    type(errors), target :: deferr

    ! Initialization
    m = size(c)
    k = size(tau)
    side = 'L'
    nrowa = m
    if (trans) then
        t = 'T'
    else
        t = 'N'
    end if
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if

    ! Input Check
    flag = 0
    if (size(a, 1) /= m .or. size(a, 2) < k) then
        call report_matrix_size_error("mult_qr_vec", errmgr, "a", m, k, &
            size(a, 1), size(a, 2))
        return
    end if

    ! Workspace Query
    call DORMQR(side, t, m, 1, k, a, nrowa, tau, c, m, temp, -1, flag)
    lwork = int(temp(1), int32)
    if (present(olwork)) then
        olwork = lwork
        return
    end if

    ! Local Memory Allocation
    if (present(work)) then
        if (size(work) < lwork) then
            call report_array_size_error("mult_qr_vec", errmgr, "work", &
                lwork, size(work))
            return
        end if
        wptr => work(1:lwork)
    else
        allocate(wrk(lwork), stat = istat)
        if (istat /= 0) then
            call report_memory_error("mult_qr_vec", errmgr, istat)
            return
        end if
        wptr => wrk
    end if

    ! Call DORMQR
    call DORMQR(side, t, m, 1, k, a, nrowa, tau, c, m, wptr, lwork, flag)
end subroutine

! ------------------------------------------------------------------------------
subroutine mult_qr_vec_cmplx(trans, a, tau, c, work, olwork, err)
    !! Multiplies a vector by the orthogonal matrix \(Q\) from a QR 
    !! factorization such that \(\vec{c} = op(Q) \vec{c}\).
    logical, intent(in) :: trans
        !! Set to true to apply \(Q^H\); else, set to false to apply \(Q\).
    complex(real64), intent(inout), dimension(:,:) :: a
        !! On input, an M-by-K matrix containing the elementary reflectors 
        !! output from the QR factorization. Notice, the contents of this matrix
        !! are restored on exit.
    complex(real64), intent(in), dimension(:) :: tau
        !! A K-element array containing the scalar factors of each elementary 
        !! reflector defined in\(A\).
    complex(real64), intent(inout), dimension(:) :: c
        !! On input, the M-element vector \(\vec{c}\).  On output, the
        !! product of the orthogonal matrix \(Q\) and the original vector 
        !! \(\vec{c}\).
    complex(real64), intent(out), target, dimension(:), optional :: work
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
    complex(real64), parameter :: one = (1.0d0, 0.0d0)

    ! Local Variables
    character :: side, t
    integer(int32) :: m, k, nrowa, istat, flag, lwork
    complex(real64), pointer, dimension(:) :: wptr
    complex(real64), allocatable, target, dimension(:) :: wrk
    complex(real64), dimension(1) :: temp
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    character(len = :), allocatable :: errmsg

    ! Initialization
    m = size(c)
    k = size(tau)
    side = 'L'
    nrowa = m
    if (trans) then
        t = 'C'
    else
        t = 'N'
    end if
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if

    ! Input Check
    flag = 0
    if (size(a, 1) /= m .or. size(a, 2) < k) then
        call report_matrix_size_error("mult_qr_vec_cmplx", errmgr, "a", m, k, &
            size(a, 1), size(a, 2))
        return
    end if

    ! Workspace Query
    call ZUNMQR(side, t, m, 1, k, a, nrowa, tau, c, m, temp, -1, flag)
    lwork = int(temp(1), int32)
    if (present(olwork)) then
        olwork = lwork
        return
    end if

    ! Local Memory Allocation
    if (present(work)) then
        if (size(work) < lwork) then
            call report_array_size_error("mult_qr_vec_cmplx", errmgr, "work", &
                lwork, size(work))
            return
        end if
        wptr => work(1:lwork)
    else
        allocate(wrk(lwork), stat = istat)
        if (istat /= 0) then
            call report_memory_error("mult_qr_vec_cmplx", errmgr, istat)
            return
        end if
        wptr => wrk
    end if

    ! Call ZUNMQR
    call ZUNMQR(side, t, m, 1, k, a, nrowa, tau, c, m, wptr, lwork, flag)
end subroutine

! ------------------------------------------------------------------------------
subroutine qr_rank1_update_dbl(q, r, u, v, work, err)
    !! Computes the rank-1 update to an M-by-N QR factored matrix \(A\) where
    !! \(M \ge N\), \(A = Q R\), and \(A_1 = A + \vec{u} \vec{v}^T\) such that
    !! \(A_1 = Q_1 R_1\).
    real(real64), intent(inout), dimension(:,:) :: q
        !! On input, the original M-by-K orthogonal matrix \(Q\).  On output, 
        !! the updated matrix \(Q_1\).
    real(real64), intent(inout), dimension(:,:) :: r
        !! On input, the M-by-N matrix \(R\).  On output, the updated matrix 
        !! \(R_1\).
    real(real64), intent(inout), dimension(:) :: u
        !! On input, the M-element \(\vec{u}\) update vector.  On output, the 
        !! original content of the array is overwritten.
    real(real64), intent(inout), dimension(:) :: v
        !! On input, the N-element \(\vec{v}\) update vector.  On output, the 
        !! original content of the array is overwritten.
    real(real64), intent(out), target, optional, dimension(:) :: work
        !! An optional argument that if supplied prevents local memory 
        !! allocation.  If provided, the array must have at least K elements.
    class(errors), intent(inout), optional, target :: err
        !! The error object to be updated.

    ! Local Variables
    logical :: full
    integer(int32) :: m, n, k, lwork, istat
    real(real64), pointer, dimension(:) :: wptr
    real(real64), allocatable, target, dimension(:) :: wrk
    class(errors), pointer :: errmgr
    type(errors), target :: deferr

    ! Initialization
    m = size(u, 1)
    n = size(r, 2)
    k = min(m, n)
    full = size(q, 2) == m
    lwork = 2 * k
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if

    ! Input Check
    if (m < n) then
        call errmgr%report_error("qr_rank1_update_dbl", &
            "The problem must not be underdetermined.", LA_INVALID_INPUT_ERROR)
        return
    else if (.not.full .and. size(q, 2) /= k) then
        call report_matrix_size_error("qr_rank1_update_dbl", errmgr, "q", &
            m, m, size(q, 1), size(q, 2))
        return
    else if (size(r, 1) /= m) then
        call report_inner_matrix_dimension_error("qr_rank1_update_dbl", &
            errmgr, "q", "r", m, size(r, 1))
        return
    else if (size(u) /= m) then
        call report_array_size_error("qr_rank1_update_dbl", errmgr, "u", m, &
            size(u))
        return
    else if (size(v) /= n) then
        call report_array_size_error("qr_rank1_update_dbl", errmgr, "v", n, &
            size(v))
        return
    end if

    ! Local Memory Allocation
    if (present(work)) then
        if (size(work) < lwork) then
            call report_array_size_error("qr_rank1_update_dbl", errmgr, &
                "work", lwork, size(work))
            return
        end if
        wptr => work(1:lwork)
    else
        allocate(wrk(lwork), stat = istat)
        if (istat /= 0) then
            call report_memory_error("qr_rank1_update_dbl", errmgr, istat)
            return
        end if
        wptr => wrk
    end if

    ! Process
    call DQR1UP(m, n, k, q, m, r, m, u, v, wptr)
end subroutine

! ------------------------------------------------------------------------------
module subroutine qr_rank1_update_cmplx(q, r, u, v, work, rwork, err)
    !! Computes the rank-1 update to an M-by-N QR factored matrix \(A\) where
    !! \(M \ge N\), \(A = Q R\), and \(A_1 = A + \vec{u} \vec{v}^H\) such that
    !! \(A_1 = Q_1 R_1\).
    complex(real64), intent(inout), dimension(:,:) :: q
        !! On input, the original M-by-K orthogonal matrix \(Q\).  On output, 
        !! the updated matrix \(Q_1\).
    complex(real64), intent(inout), dimension(:,:) :: r
        !! On input, the M-by-N matrix \(R\).  On output, the updated matrix 
        !! \(R_1\).
    complex(real64), intent(inout), dimension(:) :: u
        !! On input, the M-element \(\vec{u}\) update vector.  On output, the 
        !! original content of the array is overwritten.
    complex(real64), intent(inout), dimension(:) :: v
        !! On input, the N-element \(\vec{v}\) update vector.  On output, the 
        !! original content of the array is overwritten.
    complex(real64), intent(out), target, optional, dimension(:) :: work
        !! An optional argument that if supplied prevents local memory 
        !! allocation.  If provided, the array must have at least K elements.
    real(real64), intent(out), target, optional, dimension(:) :: rwork
        !! An optional argument that if supplied prevents local memory 
        !! allocation.  If provided, the array must have at least K elements.
    class(errors), intent(inout), optional, target :: err
        !! The error object to be updated.

    ! Local Variables
    logical :: full
    integer(int32) :: m, n, k, lwork, istat, lrwork
    complex(real64), pointer, dimension(:) :: wptr
    complex(real64), allocatable, target, dimension(:) :: wrk
    real(real64), pointer, dimension(:) :: rwptr
    real(real64), allocatable, target, dimension(:) :: rwrk
    class(errors), pointer :: errmgr
    type(errors), target :: deferr

    ! Initialization
    m = size(u, 1)
    n = size(r, 2)
    k = min(m, n)
    full = size(q, 2) == m
    lwork = k
    lrwork = k
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if

    ! Input Check
    if (m < n) then
        call errmgr%report_error("qr_rank1_update_cmplx", &
            "The problem must not be underdetermined.", LA_INVALID_INPUT_ERROR)
        return
    else if (.not.full .and. size(q, 2) /= k) then
        call report_matrix_size_error("qr_rank1_update_cmplx", errmgr, "q", &
            m, m, size(q, 1), size(q, 2))
        return
    else if (size(r, 1) /= m) then
        call report_inner_matrix_dimension_error("qr_rank1_update_cmplx", &
            errmgr, "q", "r", m, size(r, 1))
        return
    else if (size(u) /= m) then
        call report_array_size_error("qr_rank1_update_cmplx", errmgr, "u", m, &
            size(u))
        return
    else if (size(v) /= n) then
        call report_array_size_error("qr_rank1_update_cmplx", errmgr, "v", n, &
            size(v))
        return
    end if

    ! Local Memory Allocation
    if (present(work)) then
        if (size(work) < lwork) then
            call report_array_size_error("qr_rank1_update_cmplx", errmgr, &
                "work", lwork, size(work))
            return
        end if
        wptr => work(1:lwork)
    else
        allocate(wrk(lwork), stat = istat)
        if (istat /= 0) then
            call report_memory_error("qr_rank1_update_cmplx", errmgr, istat)
            return
        end if
        wptr => wrk
    end if

    if (present(rwork)) then
        if (size(rwork) < lrwork) then
            call report_array_size_error("qr_rank1_update_cmplx", errmgr, &
                "rwork", lrwork, size(rwork))
            return
        end if
        wptr => work(1:lrwork)
    else
        allocate(rwrk(lrwork), stat = istat)
        if (istat /= 0) then
            call report_memory_error("qr_rank1_update_cmplx", errmgr, istat)
            return
        end if
        rwptr => rwrk
    end if

    ! Process
    call ZQR1UP(m, n, k, q, m, r, m, u, v, wptr, rwptr)
end subroutine

! ------------------------------------------------------------------------------
subroutine solve_qr_no_pivot_mtx(a, tau, b, work, olwork, err)
    !! Solves a system of M QR-factored equations of N unknowns.  M must be
    !! greater than or equal to N.
    real(real64), intent(inout), dimension(:,:) :: a
        !! On input, the M-by-N QR factored matrix as returned by qr_factor.  
        !! On output, the contents of this matrix are restored.  Notice, M must
        !! be greater than or equal to N.
    real(real64), intent(in), dimension(:) :: tau
        !! A MIN(M, N)-element array containing the scalar factors of the 
        !! elementary reflectors as returned by qr_factor.
    real(real64), intent(inout), dimension(:,:) :: b
        !! On input, the M-by-NRHS right-hand-side matrix.  On output, the 
        !! first N rows are overwritten by the solution matrix.
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

    ! Parameters
    real(real64), parameter :: one = 1.0d0

    ! Local Variables
    integer(int32) :: m, n, nrhs, k, lwork, istat
    real(real64), pointer, dimension(:) :: wptr
    real(real64), allocatable, target, dimension(:) :: wrk
    class(errors), pointer :: errmgr
    type(errors), target :: deferr

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
    if (m < n) then
        call errmgr%report_error("solve_qr_no_pivot_mtx", &
            "The problem must not be underdetermined.", LA_INVALID_INPUT_ERROR)
        return
    else if (size(tau) /= k) then
        call report_array_size_error("solve_qr_no_pivot_mtx", errmgr, "tau", &
            k, size(tau))
        return
    else if (size(b, 1) /= m) then
        call report_matrix_size_error("solve_qr_no_pivot_mtx", errmgr, "b", &
            m, nrhs, size(b, 1), size(b, 2))
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
subroutine solve_qr_no_pivot_mtx_cmplx(a, tau, b, work, olwork, err)
    !! Solves a system of M QR-factored equations of N unknowns.  M must be
    !! greater than or equal to N.
    complex(real64), intent(inout), dimension(:,:) :: a
        !! On input, the M-by-N QR factored matrix as returned by qr_factor.  
        !! On output, the contents of this matrix are restored.  Notice, M must
        !! be greater than or equal to N.
    complex(real64), intent(in), dimension(:) :: tau
        !! A MIN(M, N)-element array containing the scalar factors of the 
        !! elementary reflectors as returned by qr_factor.
    complex(real64), intent(inout), dimension(:,:) :: b
        !! On input, the M-by-NRHS right-hand-side matrix.  On output, the 
        !! first N rows are overwritten by the solution matrix.
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

    ! Parameters
    complex(real64), parameter :: one = (1.0d0, 0.0d0)

    ! Local Variables
    integer(int32) :: m, n, nrhs, k, lwork, istat
    complex(real64), pointer, dimension(:) :: wptr
    complex(real64), allocatable, target, dimension(:) :: wrk
    class(errors), pointer :: errmgr
    type(errors), target :: deferr

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
    if (m < n) then
        call errmgr%report_error("solve_qr_no_pivot_mtx_cmplx", &
            "The problem must not be underdetermined.", LA_INVALID_INPUT_ERROR)
        return
    else if (size(tau) /= k) then
        call report_array_size_error("solve_qr_no_pivot_mtx_cmplx", errmgr, &
            "tau", k, size(tau))
        return
    else if (size(b, 1) /= m) then
        call report_matrix_size_error("solve_qr_no_pivot_mtx_cmplx", errmgr, &
            "b", m, nrhs, size(b, 1), size(b, 2))
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
            call report_array_size_error("solve_qr_no_pivot_mtx_cmplx", &
                errmgr, "work", lwork, size(work))
            return
        end if
        wptr => work(1:lwork)
    else
        allocate(wrk(lwork), stat = istat)
        if (istat /= 0) then
            call report_memory_error("solve_qr_no_pivot_mtx_cmplx", &
                errmgr, istat)
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
subroutine solve_qr_no_pivot_vec(a, tau, b, work, olwork, err)
    !! Solves a system of M QR-factored equations of N unknowns.  M must be
    !! greater than or equal to N.
    real(real64), intent(inout), dimension(:,:) :: a
        !! On input, the M-by-N QR factored matrix as returned by qr_factor.  
        !! On output, the contents of this matrix are restored.  Notice, M must
        !! be greater than or equal to N.
    real(real64), intent(in), dimension(:) :: tau
        !! A MIN(M, N)-element array containing the scalar factors of the 
        !! elementary reflectors as returned by qr_factor.
    real(real64), intent(inout), dimension(:) :: b
        !! On input, the M-element right-hand-side vector.  On output, the first
        !! N elements are overwritten with the solution vector.
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
    integer(int32) :: m, n, k, flag, lwork, istat
    real(real64), pointer, dimension(:) :: wptr
    real(real64), allocatable, target, dimension(:) :: wrk
    class(errors), pointer :: errmgr
    type(errors), target :: deferr

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
    if (m < n) then
        call errmgr%report_error("solve_qr_no_pivot_vec", &
            "The problem must not be underdetermined.", LA_INVALID_INPUT_ERROR)
        return
    else if (size(tau) /= k) then
        call report_array_size_error("solve_qr_no_pivot_vec", errmgr, "tau", &
            k, size(tau))
        return
    else if (size(b) /= m) then
        call report_array_size_error("solve_qr_no_pivot_vec", errmgr, "b", &
            m, size(b))
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
            call report_array_size_error("solve_qr_no_pivot_vec", errmgr, &
                "work", lwork, size(work))
            return
        end if
        wptr => work(1:lwork)
    else
        allocate(wrk(lwork), stat = istat)
        if (istat /= 0) then
            call report_memory_error("solve_qr_no_pivot_vec", errmgr, istat)
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
module subroutine solve_qr_no_pivot_vec_cmplx(a, tau, b, work, olwork, err)
    !! Solves a system of M QR-factored equations of N unknowns.  M must be
    !! greater than or equal to N.
    complex(real64), intent(inout), dimension(:,:) :: a
        !! On input, the M-by-N QR factored matrix as returned by qr_factor.  
        !! On output, the contents of this matrix are restored.  Notice, M must
        !! be greater than or equal to N.
    complex(real64), intent(in), dimension(:) :: tau
        !! A MIN(M, N)-element array containing the scalar factors of the 
        !! elementary reflectors as returned by qr_factor.
    complex(real64), intent(inout), dimension(:) :: b
        !! On input, the M-element right-hand-side vector.  On output, the first
        !! N elements are overwritten with the solution vector.
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
    integer(int32) :: m, n, k, lwork, istat
    complex(real64), pointer, dimension(:) :: wptr
    complex(real64), allocatable, target, dimension(:) :: wrk
    class(errors), pointer :: errmgr
    type(errors), target :: deferr

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
    if (m < n) then
        call errmgr%report_error("solve_qr_no_pivot_vec_cmplx", &
            "The problem must not be underdetermined.", LA_INVALID_INPUT_ERROR)
        return
    else if (size(tau) /= k) then
        call report_array_size_error("solve_qr_no_pivot_vec_cmplx", errmgr, &
            "tau", k, size(tau))
        return
    else if (size(b) /= m) then
        call report_array_size_error("solve_qr_no_pivot_vec_cmplx", errmgr, &
            "b", m, size(b))
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
            call report_array_size_error("solve_qr_no_pivot_vec_cmplx", &
                errmgr, "work", lwork, size(work))
            return
        end if
        wptr => work(1:lwork)
    else
        allocate(wrk(lwork), stat = istat)
        if (istat /= 0) then
            call report_memory_error("solve_qr_no_pivot_vec_cmplx", &
                errmgr, istat)
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
subroutine solve_qr_pivot_mtx(a, tau, jpvt, b, work, olwork, err)
    !! Solves a system of M QR-factored equations of N unknowns.
    real(real64), intent(inout), dimension(:,:) :: a
        !! On input, the M-by-N QR factored matrix as returned by qr_factor.  
        !! On output, the contents of this matrix are restored.
    real(real64), intent(in), dimension(:) :: tau
        !! A MIN(M, N)-element array containing the scalar factors of the 
        !! elementary reflectors as returned by qr_factor.
    integer(int32), intent(in), dimension(:) :: jpvt
        !! An N-element array, as output by qr_factor, used to track the 
        !! column pivots.
    real(real64), intent(inout), dimension(:,:) :: b
        !! On input, the MAX(M, N)-by-NRHS right-hand-side matrix.  On output,
        !! the first N rows are overwritten by the solution matrix.
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

    ! Parameters
    integer(int32), parameter :: imin = 2
    integer(int32), parameter :: imax = 1
    real(real64), parameter :: zero = 0.0d0
    real(real64), parameter :: one = 1.0d0

    ! Local Variables
    integer(int32) :: i, j, m, n, mn, nrhs, lwork, ismin, ismax, &
        rnk, maxmn, istat, lwork1, lwork2, lwork3
    real(real64) :: rcond, smax, smin, smaxpr, sminpr, s1, c1, s2, c2
    real(real64), pointer, dimension(:) :: wptr, w, tau2
    real(real64), allocatable, target, dimension(:) :: wrk
    class(errors), pointer :: errmgr
    type(errors), target :: deferr

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
    if (size(tau) /= mn) then
        call report_array_size_error("solve_qr_pivot_mtx", errmgr, "tau", &
            mn, size(tau))
        return
    else if (size(jpvt) /= n) then
        call report_array_size_error("solve_qr_pivot_mtx", errmgr, "jpvt", &
            n, size(jpvt))
        return
    else if (size(b, 1) /= maxmn) then
        call report_matrix_size_error("solve_qr_pivot_mtx", errmgr, "b", &
            maxmn, nrhs, size(b, 1), size(b, 2))
        return
    end if

    ! Workspace Query
    call rz_factor(a(1:mn,:), a(1:mn,1), olwork = lwork1)
    call mult_qr(.true., .true., a, tau, b(1:m,:), olwork = lwork2)
    call mult_rz(.true., .true., n, a(1:mn,:), a(1:mn,1), b(1:n,:), &
        olwork = lwork3)
    lwork = max(lwork1, lwork2, lwork3, 2 * mn + 1) + mn
    if (present(olwork)) then
        olwork = lwork
        return
    end if

    ! Local Memory Allocation
    if (present(work)) then
        if (size(work) < lwork) then
            call report_array_size_error("solve_qr_pivot_mtx", errmgr, &
                "work", lwork, size(work))
            return
        end if
        wptr => work(1:lwork)
    else
        allocate(wrk(lwork), stat = istat)
        if (istat /= 0) then
            call report_memory_error("solve_qr_pivot_mtx", errmgr, istat)
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
        b(1:maxmn,:) = zero
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
        call mult_rz(.true., .true., n - rnk, a(1:rnk,:), tau2, b(1:n,:), w)
    end if

    ! Apply the pivoting: B(1:N,1:NRHS) = P * B(1:N,1:NRHS)
    do j = 1, nrhs
        do i = 1, n
            wptr(jpvt(i)) = b(i,j)
        end do
        b(1:n,j) = wptr(1:n)
    end do
end subroutine

! ------------------------------------------------------------------------------
subroutine solve_qr_pivot_mtx_cmplx(a, tau, jpvt, b, work, olwork, err)
    !! Solves a system of M QR-factored equations of N unknowns.
    complex(real64), intent(inout), dimension(:,:) :: a
        !! On input, the M-by-N QR factored matrix as returned by qr_factor.  
        !! On output, the contents of this matrix are restored.
    complex(real64), intent(in), dimension(:) :: tau
        !! A MIN(M, N)-element array containing the scalar factors of the 
        !! elementary reflectors as returned by qr_factor.
    integer(int32), intent(in), dimension(:) :: jpvt
        !! An N-element array, as output by qr_factor, used to track the 
        !! column pivots.
    complex(real64), intent(inout), dimension(:,:) :: b
        !! On input, the MAX(M, N)-by-NRHS right-hand-side matrix.  On output,
        !! the first N rows are overwritten by the solution matrix.
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

    ! Parameters
    integer(int32), parameter :: imin = 2
    integer(int32), parameter :: imax = 1
    complex(real64), parameter :: zero = (0.0d0, 0.0d0)
    complex(real64), parameter :: one = (1.0d0, 0.0d0)

    ! Local Variables
    integer(int32) :: i, j, m, n, mn, nrhs, lwork, ismin, ismax, &
        rnk, maxmn, istat, lwork1, lwork2, lwork3
    real(real64) :: rcond, smax, smin, smaxpr, sminpr
    complex(real64) :: s1, c1, s2, c2
    complex(real64), pointer, dimension(:) :: wptr, w, tau2
    complex(real64), allocatable, target, dimension(:) :: wrk
    class(errors), pointer :: errmgr
    type(errors), target :: deferr

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
    if (size(tau) /= mn) then
        call report_array_size_error("solve_qr_pivot_mtx_cmplx", errmgr, &
            "tau", mn, size(tau))
        return
    else if (size(jpvt) /= n) then
        call report_array_size_error("solve_qr_pivot_mtx_cmplx", errmgr, &
            "jpvt", n, size(jpvt))
        return
    else if (size(b, 1) /= maxmn) then
        call report_matrix_size_error("solve_qr_pivot_mtx_cmplx", errmgr, "b", &
            maxmn, nrhs, size(b, 1), size(b, 2))
        return
    end if

    ! Workspace Query
    call rz_factor(a(1:mn,:), a(1:mn,1), olwork = lwork1)
    call mult_qr(.true., .true., a, tau, b(1:m,:), olwork = lwork2)
    call mult_rz(.true., .true., n, a(1:mn,:), a(1:mn,1), b(1:n,:), &
        olwork = lwork3)
    lwork = max(lwork1, lwork2, lwork3, 2 * mn + 1) + mn
    if (present(olwork)) then
        olwork = lwork
        return
    end if

    ! Local Memory Allocation
    if (present(work)) then
        if (size(work) < lwork) then
            call report_array_size_error("solve_qr_pivot_mtx_cmplx", errmgr, &
                "work", lwork, size(work))
            return
        end if
        wptr => work(1:lwork)
    else
        allocate(wrk(lwork), stat = istat)
        if (istat /= 0) then
            call report_memory_error("solve_qr_pivot_mtx_cmplx", errmgr, istat)
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
        b(1:maxmn,:) = zero
        return
    else
        rnk = 1
    end if
    do
        if (rnk < mn) then
            i = rnk + 1
            call ZLAIC1(imin, rnk, wptr(ismin:ismin+rnk-1), smin, &
                a(1:rnk-1,i), a(i,i), sminpr, s1, c1)
            call ZLAIC1(imax, rnk, wptr(ismax:ismax+rnk-1), smax, &
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
        call mult_rz(.true., .true., n - rnk, a(1:rnk,:), tau2, b(1:n,:), w)
    end if

    ! Apply the pivoting: B(1:N,1:NRHS) = P * B(1:N,1:NRHS)
    do j = 1, nrhs
        do i = 1, n
            wptr(jpvt(i)) = b(i,j)
        end do
        b(1:n,j) = wptr(1:n)
    end do
end subroutine

! ------------------------------------------------------------------------------
module subroutine solve_qr_pivot_vec(a, tau, jpvt, b, work, olwork, err)
    !! Solves a system of M QR-factored equations of N unknowns.
    real(real64), intent(inout), dimension(:,:) :: a
        !! On input, the M-by-N QR factored matrix as returned by qr_factor.  
        !! On output, the contents of this matrix are restored.
    real(real64), intent(in), dimension(:) :: tau
        !! A MIN(M, N)-element array containing the scalar factors of the 
        !! elementary reflectors as returned by qr_factor.
    integer(int32), intent(in), dimension(:) :: jpvt
        !! An N-element array, as output by qr_factor, used to track the 
        !! column pivots.
    real(real64), intent(inout), dimension(:) :: b
        !! On input, the MAX(M, N)-by-NRHS right-hand-side vector.  On output,
        !! the first N rows are overwritten by the solution vector.
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

    ! Parameters
    integer(int32), parameter :: imin = 2
    integer(int32), parameter :: imax = 1
    real(real64), parameter :: zero = 0.0d0
    real(real64), parameter :: one = 1.0d0

    ! Local Variables
    integer(int32) :: i, m, n, mn, lwork, ismin, ismax, rnk, maxmn, &
        istat, lwork1, lwork2
    real(real64) :: rcond, smax, smin, smaxpr, sminpr, s1, c1, s2, c2
    real(real64), pointer, dimension(:) :: wptr, w, tau2
    real(real64), allocatable, target, dimension(:) :: wrk
    class(errors), pointer :: errmgr
    type(errors), target :: deferr

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
    if (size(tau) /= mn) then
        call report_array_size_error("solve_qr_pivot_vec", errmgr, "tau", &
            mn, size(tau))
        return
    else if (size(jpvt) /= n) then
        call report_array_size_error("solve_qr_pivot_vec", errmgr, "jpvt", &
            n, size(jpvt))
        return
    else if (size(b) /= maxmn) then
        call report_array_size_error("solve_qr_pivot_vec", errmgr, "b", &
            maxmn, size(b))
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
            call report_array_size_error("solve_qr_pivot_vec", errmgr, &
                "work", lwork, size(work))
            return
        end if
        wptr => work(1:lwork)
    else
        allocate(wrk(lwork), stat = istat)
        if (istat /= 0) then
            call report_memory_error("solve_qr_pivot_vec", errmgr, istat)
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

! ------------------------------------------------------------------------------
module subroutine solve_qr_pivot_vec_cmplx(a, tau, jpvt, b, work, olwork, err)
    !! Solves a system of M QR-factored equations of N unknowns.
    complex(real64), intent(inout), dimension(:,:) :: a
        !! On input, the M-by-N QR factored matrix as returned by qr_factor.  
        !! On output, the contents of this matrix are restored.
    complex(real64), intent(in), dimension(:) :: tau
        !! A MIN(M, N)-element array containing the scalar factors of the 
        !! elementary reflectors as returned by qr_factor.
    integer(int32), intent(in), dimension(:) :: jpvt
        !! An N-element array, as output by qr_factor, used to track the 
        !! column pivots.
    complex(real64), intent(inout), dimension(:) :: b
        !! On input, the MAX(M, N)-by-NRHS right-hand-side vector.  On output,
        !! the first N rows are overwritten by the solution vector.
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

    ! Parameters
    integer(int32), parameter :: imin = 2
    integer(int32), parameter :: imax = 1
    complex(real64), parameter :: zero = (0.0d0, 0.0d0)
    complex(real64), parameter :: one = (1.0d0, 0.0d0)

    ! Local Variables
    integer(int32) :: i, m, n, mn, lwork, ismin, ismax, rnk, maxmn, &
        istat, lwork1, lwork2
    real(real64) :: rcond, smax, smin, smaxpr, sminpr
    complex(real64) :: s1, c1, s2, c2
    complex(real64), pointer, dimension(:) :: wptr, w, tau2
    complex(real64), allocatable, target, dimension(:) :: wrk
    class(errors), pointer :: errmgr
    type(errors), target :: deferr

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
    if (size(tau) /= mn) then
        call report_array_size_error("solve_qr_pivot_vec_cmplx", errmgr, &
            "tau", mn, size(tau))
        return
    else if (size(jpvt) /= n) then
        call report_array_size_error("solve_qr_pivot_vec_cmplx", errmgr, &
            "jpvt", n, size(jpvt))
        return
    else if (size(b) /= maxmn) then
        call report_array_size_error("solve_qr_pivot_vec_cmplx", errmgr, "b", &
            maxmn, size(b))
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
            call report_array_size_error("solve_qr_pivot_vec_cmplx", errmgr, &
                "work", lwork, size(work))
            return
        end if
        wptr => work(1:lwork)
    else
        allocate(wrk(lwork), stat = istat)
        if (istat /= 0) then
            call report_memory_error("solve_qr_pivot_vec_cmplx", errmgr, istat)
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
            call ZLAIC1(imin, rnk, wptr(ismin:ismin+rnk-1), smin, &
                a(1:rnk-1,i), a(i,i), sminpr, s1, c1)
            call ZLAIC1(imax, rnk, wptr(ismax:ismax+rnk-1), smax, &
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

! ------------------------------------------------------------------------------
end module