module linalg_lq
    use iso_fortran_env, only : int32, real64
    use lapack
    use linalg_errors
    use linalg_tri
    implicit none
    private
    public :: lq_factor
    public :: form_lq
    public :: mult_lq
    public :: solve_lq

    interface lq_factor
        module procedure :: lq_factor_no_pivot
        module procedure :: lq_factor_no_pivot_cmplx
    end interface

    interface form_lq
        module procedure :: form_lq_no_pivot
        module procedure :: form_lq_no_pivot_cmplx
    end interface

    interface mult_lq
        module procedure :: mult_lq_mtx
        module procedure :: mult_lq_mtx_cmplx
        module procedure :: mult_lq_vec
        module procedure :: mult_lq_vec_cmplx
    end interface

    interface solve_lq
        module procedure :: solve_lq_mtx
        module procedure :: solve_lq_mtx_cmplx
        module procedure :: solve_lq_vec
        module procedure :: solve_lq_vec_cmplx
    end interface
contains
! ------------------------------------------------------------------------------
subroutine lq_factor_no_pivot(a, tau, work, olwork, err)
    !! Computes the LQ factorization of an M-by-N matrix \(A = L Q\) where
    !! \(L\) is a lower triangular (or lower trapezoidal) matrix and \(Q\) is
    !! a orthogonal matrix.
    real(real64), intent(inout), dimension(:,:) :: a
        !! On input, the M-by-N matrix to factor.  On output, the elements on 
        !! and below the diagonal contain the MIN(M, N)-by-N lower trapezoidal 
        !! matrix \(L\) (\(L\) is lower triangular if M >= N).  The elements
        !! above the diagonal, along with the array tau, represent the 
        !! orthogonal matrix \(Q\) as a product of elementary reflectors.
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
        call report_array_size_error("lq_factor_no_pivot", errmgr, "tau", &
            mn, size(tau))
        return
    end if

    ! Workspace Query
    call DGELQF(m, n, a, m, tau, temp, -1, flag)
    lwork = int(temp(1), int32)
    if (present(olwork)) then
        olwork = lwork
        return
    end if

    ! Local Memory Allocation
    if (present(work)) then
        if (size(work) < lwork) then
            call report_array_size_error("lq_factor_no_pivot", errmgr, &
                "work", lwork, size(work))
            return
        end if
        wptr => work(1:lwork)
    else
        allocate(wrk(lwork), stat = istat)
        if (istat /= 0) then
            call report_memory_error("lq_factor_no_pivot", errmgr, istat)
            return
        end if
        wptr => wrk
    end if

    ! Call DGELQF
    call DGELQF(m, n, a, m, tau, wptr, lwork, flag)
end subroutine

! ------------------------------------------------------------------------------
subroutine lq_factor_no_pivot_cmplx(a, tau, work, olwork, err)
    !! Computes the LQ factorization of an M-by-N matrix \(A = L Q\) where
    !! \(L\) is a lower triangular (or lower trapezoidal) matrix and \(Q\) is
    !! a orthogonal matrix.
    complex(real64), intent(inout), dimension(:,:) :: a
        !! On input, the M-by-N matrix to factor.  On output, the elements on 
        !! and below the diagonal contain the MIN(M, N)-by-N lower trapezoidal 
        !! matrix \(L\) (\(L\) is lower triangular if M >= N).  The elements
        !! above the diagonal, along with the array tau, represent the 
        !! orthogonal matrix \(Q\) as a product of elementary reflectors.
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
        call report_array_size_error("lq_factor_no_pivot_cmplx", errmgr, &
            "tau", mn, size(tau))
        return
    end if

    ! Workspace Query
    call ZGELQF(m, n, a, m, tau, temp, -1, flag)
    lwork = int(temp(1), int32)
    if (present(olwork)) then
        olwork = lwork
        return
    end if

    ! Local Memory Allocation
    if (present(work)) then
        if (size(work) < lwork) then
            call report_array_size_error("lq_factor_no_pivot_cmplx", errmgr, &
                "work", lwork, size(work))
            return
        end if
        wptr => work(1:lwork)
    else
        allocate(wrk(lwork), stat = istat)
        if (istat /= 0) then
            call report_memory_error("lq_factor_no_pivot_cmplx", errmgr, istat)
            return
        end if
        wptr => wrk
    end if

    ! Call ZGELQF
    call ZGELQF(m, n, a, m, tau, wptr, lwork, flag)
end subroutine

! ------------------------------------------------------------------------------
subroutine form_lq_no_pivot(l, tau, q, work, olwork, err)
    !! Forms the orthogonal matrix \(Q\) from the elementary reflectors returned 
    !! by the LQ factorization algorithm.
    real(real64), intent(inout), dimension(:,:) :: l
        !! On input, an M-by-N matrix where the elements above the diagonal 
        !! contain the elementary reflectors generated from the LQ factorization
        !! performed by lq_factor.  On and below the diagonal the matrix 
        !! contains the matrix \(L\).  On output, the elements above the 
        !! diagonal are zeroed sucht hat the remaining matrix is the M-by-N 
        !! lower trapezoidal matrix \(L\) where only the M-by-M submatrix is 
        !! the lower triangular matrix \(L\).  Notice, M must be less than or 
        !! equal to N for this routine.
    real(real64), intent(in), dimension(:) :: tau
        !! A MIN(M, N)-element array containing the scalar factors of each 
        !! elementary reflector defined in \(L\).
    real(real64), intent(out), dimension(:,:) :: q
        !! An N-by-N matrix where the orthogonal matrix \(Q\) will be written.
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
    integer(int32) :: i, j, m, n, mn, k, istat, flag, lwork
    real(real64), pointer, dimension(:) :: wptr
    real(real64), allocatable, target, dimension(:) :: wrk
    real(real64), dimension(1) :: temp
    class(errors), pointer :: errmgr
    type(errors), target :: deferr

    ! Initialization
    m = size(l, 1)
    n = size(l, 2)
    mn = min(m, n)
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if

    ! Input Check
    if (m > n) then
        call errmgr%report_error("form_lq_no_pivot", &
            "This routine does not handle the overdetermined case.", &
            LA_INVALID_INPUT_ERROR)
        return
    else if (size(tau) /= mn) then
        call report_array_size_error("form_lq_no_pivot", errmgr, "tau", &
            mn, size(tau))
        return
    else if (size(q, 1) /= n .or. size(q, 2) /= n) then
        call report_matrix_size_error("form_lq_no_pivot", errmgr, "q", &
            n, n, size(q, 1), size(q, 2))
        return
    end if

    ! Workspace Query
    call DORGLQ(n, n, mn, q, n, tau, temp, -1, flag)
    lwork = int(temp(1), int32)
    if (present(olwork)) then
        olwork = lwork
        return
    end if

    ! Local Memory Allocation
    if (present(work)) then
        if (size(work) < lwork) then
            call report_array_size_error("form_lq_no_pivot", errmgr, &
                "work", lwork, size(work))
            return
        end if
        wptr => work(1:lwork)
    else
        allocate(wrk(lwork), stat = istat)
        if (istat /= 0) then
            call report_memory_error("form_lq_no_pivot", errmgr, istat)
            return
        end if
        wptr => wrk
    end if

    ! Copy the upper triangular portion of L to Q, and then zero it out in L
    do j = 2, n
        k = min(j - 1, m)
        q(1:k,j) = l(1:k,j)
        l(1:k,j) = zero
    end do

    ! Build Q
    call DORGLQ(n, n, mn, q, n, tau, wptr, lwork, flag)
end subroutine

! ------------------------------------------------------------------------------
subroutine form_lq_no_pivot_cmplx(l, tau, q, work, olwork, err)
    !! Forms the orthogonal matrix \(Q\) from the elementary reflectors returned 
    !! by the LQ factorization algorithm.
    complex(real64), intent(inout), dimension(:,:) :: l
        !! On input, an M-by-N matrix where the elements above the diagonal 
        !! contain the elementary reflectors generated from the LQ factorization
        !! performed by lq_factor.  On and below the diagonal the matrix 
        !! contains the matrix \(L\).  On output, the elements above the 
        !! diagonal are zeroed sucht hat the remaining matrix is the M-by-N 
        !! lower trapezoidal matrix \(L\) where only the M-by-M submatrix is 
        !! the lower triangular matrix \(L\).  Notice, M must be less than or 
        !! equal to N for this routine.
    complex(real64), intent(in), dimension(:) :: tau
        !! A MIN(M, N)-element array containing the scalar factors of each 
        !! elementary reflector defined in \(L\).
    complex(real64), intent(out), dimension(:,:) :: q
        !! An N-by-N matrix where the orthogonal matrix \(Q\) will be written.
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
    integer(int32) :: i, j, m, n, mn, k, istat, flag, lwork
    complex(real64), pointer, dimension(:) :: wptr
    complex(real64), allocatable, target, dimension(:) :: wrk
    complex(real64), dimension(1) :: temp
    class(errors), pointer :: errmgr
    type(errors), target :: deferr

    ! Initialization
    m = size(l, 1)
    n = size(l, 2)
    mn = min(m, n)
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if

    ! Input Check
    if (m > n) then
        call errmgr%report_error("form_lq_no_pivot_cmplx", &
            "This routine does not handle the overdetermined case.", &
            LA_INVALID_INPUT_ERROR)
        return
    else if (size(tau) /= mn) then
        call report_array_size_error("form_lq_no_pivot_cmplx", errmgr, "tau", &
            mn, size(tau))
        return
    else if (size(q, 1) /= n .or. size(q, 2) /= n) then
        call report_matrix_size_error("form_lq_no_pivot_cmplx", errmgr, "q", &
            n, n, size(q, 1), size(q, 2))
        return
    end if

    ! Workspace Query
    call ZUNGLQ(n, n, mn, q, n, tau, temp, -1, flag)
    lwork = int(temp(1), int32)
    if (present(olwork)) then
        olwork = lwork
        return
    end if

    ! Local Memory Allocation
    if (present(work)) then
        if (size(work) < lwork) then
            call report_array_size_error("form_lq_no_pivot_cmplx", errmgr, &
                "work", lwork, size(work))
            return
        end if
        wptr => work(1:lwork)
    else
        allocate(wrk(lwork), stat = istat)
        if (istat /= 0) then
            call report_memory_error("form_lq_no_pivot_cmplx", errmgr, istat)
            return
        end if
        wptr => wrk
    end if

    ! Copy the upper triangular portion of L to Q, and then zero it out in L
    do j = 2, n
        k = min(j - 1, m)
        q(1:k,j) = l(1:k,j)
        l(1:k,j) = zero
    end do

    ! Build Q
    call ZUNGLQ(n, n, mn, q, n, tau, wptr, lwork, flag)
end subroutine

! ------------------------------------------------------------------------------
subroutine mult_lq_mtx(lside, trans, a, tau, c, work, olwork, err)
    !! Multiplies a matrix by the orthogonal matrix \(Q\) from an LQ
    !! factorization.
    logical, intent(in) :: lside
        !! Set to true to compute \(C = op(Q) C\); else, set to false to
        !! compute \(C = C op(Q)\).
    logical, intent(in) :: trans
        !! Set to true to compute \(op(Q) = Q^T\); else, set to false to 
        !! compute \(op(Q) = Q\).
    real(real64), intent(in), dimension(:,:) :: a
        !! On input, an K-by-P matrix containing the elementary reflectors 
        !! output from the LQ factorization.  If lside is set to true, P = M; 
        !! else, if lside is set to false, P = N.
    real(real64), intent(in), dimension(:) :: tau
        !! A K-element array containing the scalar factors of each elementary 
        !! reflector defined in a.
    real(real64), intent(inout), dimension(:,:) :: c
        !! On input, the M-by-N matrix C.  On output, the product of the 
        !! orthogonal matrix Q and the original matrix C.
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
    character :: side, t
    integer(int32) :: m, n, k, ncola, istat, flag, lwork
    real(real64), pointer, dimension(:) :: wptr
    real(real64), allocatable, target, dimension(:) :: wrk
    real(real64), dimension(1) :: temp
    class(errors), pointer :: errmgr
    type(errors), target :: deferr

    ! Initialization
    m = size(c, 1)
    n = size(c, 2)
    k = size(tau)
    if (lside) then
        side = 'L'
        ncola = m
    else
        side = 'R'
        ncola = n
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
    if (size(a, 1) /= k .or. size(a, 2) /= ncola) then
        call report_matrix_size_error("mult_lq_mtx", errmgr, "a", k, ncola, &
            size(a, 1), size(a, 2))
        return
    end if

    ! Workspace Query
    call DORMLQ(side, t, m, n, k, a, k, tau, c, m, temp, -1, flag)
    lwork = int(temp(1), int32)
    if (present(olwork)) then
        olwork = lwork
        return
    end if

    ! Local Memory Allocation
    if (present(work)) then
        if (size(work) < lwork) then
            call report_array_size_error("mult_lq_mtx", errmgr, "work", &
                lwork, size(work))
            return
        end if
        wptr => work(1:lwork)
    else
        allocate(wrk(lwork), stat = istat)
        if (istat /= 0) then
            call report_memory_error("mult_lq_mtx", errmgr, istat)
            return
        end if
        wptr => wrk
    end if

    ! Call DORMLQ
    call DORMLQ(side, t, m, n, k, a, k, tau, c, m, wptr, lwork, flag)
end subroutine

! ------------------------------------------------------------------------------
subroutine mult_lq_mtx_cmplx(lside, trans, a, tau, c, work, olwork, err)
    !! Multiplies a matrix by the orthogonal matrix \(Q\) from an LQ
    !! factorization.
    logical, intent(in) :: lside
        !! Set to true to compute \(C = op(Q) C\); else, set to false to
        !! compute \(C = C op(Q)\).
    logical, intent(in) :: trans
        !! Set to true to compute \(op(Q) = Q^H\); else, set to false to 
        !! compute \(op(Q) = Q\).
    complex(real64), intent(in), dimension(:,:) :: a
        !! On input, an K-by-P matrix containing the elementary reflectors 
        !! output from the LQ factorization.  If lside is set to true, P = M; 
        !! else, if lside is set to false, P = N.
    complex(real64), intent(in), dimension(:) :: tau
        !! A K-element array containing the scalar factors of each elementary 
        !! reflector defined in a.
    complex(real64), intent(inout), dimension(:,:) :: c
        !! On input, the M-by-N matrix C.  On output, the product of the 
        !! orthogonal matrix Q and the original matrix C.
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
    character :: side, t
    integer(int32) :: m, n, k, ncola, istat, flag, lwork
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
        ncola = m
    else
        side = 'R'
        ncola = n
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
    if (size(a, 1) /= k .or. size(a, 2) /= ncola) then
        call report_matrix_size_error("mult_lq_mtx_cmplx", errmgr, "a", k, &
            ncola, size(a, 1), size(a, 2))
        return
    end if

    ! Workspace Query
    call ZUNMLQ(side, t, m, n, k, a, k, tau, c, m, temp, -1, flag)
    lwork = int(temp(1), int32)
    if (present(olwork)) then
        olwork = lwork
        return
    end if

    ! Local Memory Allocation
    if (present(work)) then
        if (size(work) < lwork) then
            call report_array_size_error("mult_lq_mtx_cmplx", errmgr, "work", &
                lwork, size(work))
            return
        end if
        wptr => work(1:lwork)
    else
        allocate(wrk(lwork), stat = istat)
        if (istat /= 0) then
            call report_memory_error("mult_lq_mtx_cmplx", errmgr, istat)
            return
        end if
        wptr => wrk
    end if

    ! Call ZUNMLQ
    call ZUNMLQ(side, t, m, n, k, a, k, tau, c, m, wptr, lwork, flag)
end subroutine

! ------------------------------------------------------------------------------
subroutine mult_lq_vec(trans, a, tau, c, work, olwork, err)
    !! Multiplies a vector with the orthogonal matrix \(Q\) from an LQ 
    !! factorization such that \(\vec{c} = op(Q) \vec{c}\).
    logical, intent(in) :: trans
        !! Set to true to compute \(op(Q) = Q^T\); else, set to false to 
        !! compute \(op(Q) = Q\).
    real(real64), intent(in), dimension(:,:) :: a
        !! On input, an K-by-M matrix containing the elementary reflectors 
        !! output from the LQ factorization.  Notice, the contents of this 
        !! matrix are restored on exit.
    real(real64), intent(in), dimension(:) :: tau
        !! A K-element array containing the scalar factors of each elementary 
        !! reflector defined in a.
    real(real64), intent(inout), dimension(:) :: c
        !! On input, the M-element vector \(\vec{c}\).  On output, the product 
        !! of the orthogonal matrix \(Q\) and the original vector \(\vec{c}\).
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
    character :: side, t
    integer(int32) :: m, n, k, istat, flag, lwork
    real(real64), pointer, dimension(:) :: wptr
    real(real64), allocatable, target, dimension(:) :: wrk
    real(real64), dimension(1) :: temp
    class(errors), pointer :: errmgr
    type(errors), target :: deferr

    ! Initialization
    m = size(c)
    n = 1
    k = size(tau)
    side = 'L'
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
    if (size(a, 1) /= k .or. size(a, 2) /= m) then
        call report_matrix_size_error("mult_lq_vec", errmgr, "a", k, m, &
            size(a, 1), size(a, 2))
        return
    end if

    ! Workspace Query
    call DORMLQ(side, t, m, n, k, a, k, tau, c, m, temp, -1, flag)
    lwork = int(temp(1), int32)
    if (present(olwork)) then
        olwork = lwork
        return
    end if

    ! Local Memory Allocation
    if (present(work)) then
        if (size(work) < lwork) then
            call report_array_size_error("mult_lq_vec", errmgr, "work", &
                lwork, size(work))
            return
        end if
        wptr => work(1:lwork)
    else
        allocate(wrk(lwork), stat = istat)
        if (istat /= 0) then
            call report_memory_error("mult_lq_vec", errmgr, istat)
            return
        end if
        wptr => wrk
    end if

    ! Call DORMLQ
    call DORMLQ(side, t, m, n, k, a, k, tau, c, m, wptr, lwork, flag)
end subroutine

! ------------------------------------------------------------------------------
subroutine mult_lq_vec_cmplx(trans, a, tau, c, work, olwork, err)
    !! Multiplies a vector with the orthogonal matrix \(Q\) from an LQ 
    !! factorization such that \(\vec{c} = op(Q) \vec{c}\).
    logical, intent(in) :: trans
        !! Set to true to compute \(op(Q) = Q^H\); else, set to false to 
        !! compute \(op(Q) = Q\).
    complex(real64), intent(in), dimension(:,:) :: a
        !! On input, an K-by-M matrix containing the elementary reflectors 
        !! output from the LQ factorization.  Notice, the contents of this 
        !! matrix are restored on exit.
    complex(real64), intent(in), dimension(:) :: tau
        !! A K-element array containing the scalar factors of each elementary 
        !! reflector defined in a.
    complex(real64), intent(inout), dimension(:) :: c
        !! On input, the M-element vector \(\vec{c}\).  On output, the product 
        !! of the orthogonal matrix \(Q\) and the original vector \(\vec{c}\).
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
    character :: side, t
    integer(int32) :: m, n, k, istat, flag, lwork
    complex(real64), pointer, dimension(:) :: wptr
    complex(real64), allocatable, target, dimension(:) :: wrk
    complex(real64), dimension(1) :: temp
    class(errors), pointer :: errmgr
    type(errors), target :: deferr

    ! Initialization
    m = size(c)
    n = 1
    k = size(tau)
    side = 'L'
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
    if (size(a, 1) /= k .or. size(a, 2) /= m) then
        call report_matrix_size_error("mult_lq_vec_cmplx", errmgr, "a", k, m, &
            size(a, 1), size(a, 2))
        return
    end if

    ! Workspace Query
    call ZUNMLQ(side, t, m, n, k, a, k, tau, c, m, temp, -1, flag)
    lwork = int(temp(1), int32)
    if (present(olwork)) then
        olwork = lwork
        return
    end if

    ! Local Memory Allocation
    if (present(work)) then
        if (size(work) < lwork) then
            call report_array_size_error("mult_lq_vec_cmplx", errmgr, "work", &
                lwork, size(work))
            return
        end if
        wptr => work(1:lwork)
    else
        allocate(wrk(lwork), stat = istat)
        if (istat /= 0) then
            call report_memory_error("mult_lq_vec_cmplx", errmgr, istat)
            return
        end if
        wptr => wrk
    end if

    ! Call ZUNMLQ
    call ZUNMLQ(side, t, m, n, k, a, k, tau, c, m, wptr, lwork, flag)
end subroutine

! ------------------------------------------------------------------------------
subroutine solve_lq_mtx(a, tau, b, work, olwork, err)
    !! Solves a system of LQ factored equations of the form \(A X = L Q X = B\).
    real(real64), intent(in), dimension(:,:) :: a
        !! On input, the M-by-N LQ factored matrix as returned by lq_factor.  
        !! On output, the contents of this matrix are restored.  Notice, N must
        !! be greater than or equal to M.
    real(real64), intent(in), dimension(:) :: tau
        !! A MIN(M, N)-element array containing the scalar factors of the 
        !! elementary reflectors as returned by lq_factor.
    real(real64), intent(inout), dimension(:,:) :: b
        !! On input, an N-by-NRHS matrix where the first M rows contain 
        !! the right-hand-side matrix \(B\).  On output, the N-by-NRHS solution
        !! matrix \(X\).
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
    if (m > n) then
        call errmgr%report_error("solve_lq_mtx", &
            "This routine does not handle the overdetermined case.", &
            LA_INVALID_INPUT_ERROR)
        return
    else if (size(tau) /= k) then
        call report_array_size_error("solve_lq_mtx", errmgr, "tau", k, &
            size(tau))
        return
    else if (size(b, 1) /= n) then
        call report_matrix_size_error("solve_lq_mtx", errmgr, "b", n, nrhs, &
            size(b, 1), size(b, 2))
        return
    end if

    ! Workspace Query
    call mult_lq(.true., .true., a, tau, b, olwork = lwork)

    if (present(olwork)) then
        olwork = lwork
        return
    end if

    ! Local Memory Allocation
    if (present(work)) then
        if (size(work) < lwork) then
            call report_array_size_error("solve_lq_mtx", errmgr, "work", &
                lwork, size(work))
            return
        end if
        wptr => work(1:lwork)
    else
        allocate(wrk(lwork), stat = istat)
        if (istat /= 0) then
            call report_memory_error("solve_lq_mtx", errmgr, istat)
            return
        end if
        wptr => wrk
    end if

    ! Solve the lower triangular system L * Y = B for Y, where Y = Q * X.
    ! The lower triangular system is M-by-M and Y is M-by-NHRS.
    call solve_triangular_system(.true., .false., .false., .true., one, &
        a(1:m,1:m), b(1:m,:), errmgr)
    if (errmgr%has_error_occurred()) return

    ! Compute Q**T * Y = X
    call mult_lq(.true., .true., a, tau, b, work = wptr, err = errmgr)
    if (errmgr%has_error_occurred()) return
end subroutine

! ------------------------------------------------------------------------------
module subroutine solve_lq_mtx_cmplx(a, tau, b, work, olwork, err)
    !! Solves a system of LQ factored equations of the form \(A X = L Q X = B\).
    complex(real64), intent(in), dimension(:,:) :: a
        !! On input, the M-by-N LQ factored matrix as returned by lq_factor.  
        !! On output, the contents of this matrix are restored.  Notice, N must
        !! be greater than or equal to M.
    complex(real64), intent(in), dimension(:) :: tau
        !! A MIN(M, N)-element array containing the scalar factors of the 
        !! elementary reflectors as returned by lq_factor.
    complex(real64), intent(inout), dimension(:,:) :: b
        !! On input, an N-by-NRHS matrix where the first M rows contain 
        !! the right-hand-side matrix \(B\).  On output, the N-by-NRHS solution
        !! matrix \(X\).
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
    if (m > n) then
        call errmgr%report_error("solve_lq_mtx_cmplx", &
            "This routine does not handle the overdetermined case.", &
            LA_INVALID_INPUT_ERROR)
        return
    else if (size(tau) /= k) then
        call report_array_size_error("solve_lq_mtx_cmplx", errmgr, "tau", k, &
            size(tau))
        return
    else if (size(b, 1) /= n) then
        call report_matrix_size_error("solve_lq_mtx_cmplx", errmgr, "b", n, &
            nrhs, size(b, 1), size(b, 2))
        return
    end if

    ! Workspace Query
    call mult_lq(.true., .true., a, tau, b, olwork = lwork)

    if (present(olwork)) then
        olwork = lwork
        return
    end if

    ! Local Memory Allocation
    if (present(work)) then
        if (size(work) < lwork) then
            call report_array_size_error("solve_lq_mtx_cmplx", errmgr, "work", &
                lwork, size(work))
            return
        end if
        wptr => work(1:lwork)
    else
        allocate(wrk(lwork), stat = istat)
        if (istat /= 0) then
            call report_memory_error("solve_lq_mtx_cmplx", errmgr, istat)
            return
        end if
        wptr => wrk
    end if

    ! Solve the lower triangular system L * Y = B for Y, where Y = Q * X.
    ! The lower triangular system is M-by-M and Y is M-by-NHRS.
    call solve_triangular_system(.true., .false., .false., .true., one, &
        a(1:m,1:m), b(1:m,:), errmgr)
    if (errmgr%has_error_occurred()) return

    ! Compute Q**T * Y = X
    call mult_lq(.true., .true., a, tau, b, work = wptr, err = errmgr)
    if (errmgr%has_error_occurred()) return
end subroutine

! ------------------------------------------------------------------------------
subroutine solve_lq_vec(a, tau, b, work, olwork, err)
    !! Solves a system of LQ factored equations of the form 
    !! \(A \vec{x} = L Q \vec{x} = \vec{b}\).
    real(real64), intent(in), dimension(:,:) :: a
        !! !! On input, the M-by-N LQ factored matrix as returned by lq_factor.  
        !! On output, the contents of this matrix are restored.  Notice, N must
        !! be greater than or equal to M.
    real(real64), intent(in), dimension(:) :: tau
        !! A MIN(M, N)-element array containing the scalar factors of the 
        !! elementary reflectors as returned by lq_factor.
    real(real64), intent(inout), dimension(:) :: b
        !! On input, an N-element vector where the first M rows contain the 
        !! right-hand-side vector \(\vec{b}\).  On output, the N-element vector
        !! \(\vec{x}\).
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
    integer(int32) :: m, n, k, lwork, istat
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
    if (m > n) then
        call errmgr%report_error("solve_lq_vec", &
            "This routine does not handle the overdetermined case.", &
            LA_INVALID_INPUT_ERROR)
        return
    else if (size(tau) /= k) then
        call report_array_size_error("solve_lq_vec", errmgr, "tau", k, &
            size(tau))
        return
    else if (size(b) /= n) then
        call report_memory_error("solve_lq_vec", errmgr, istat)
        return
    end if

    ! Workspace Query
    call mult_lq(.true., a, tau, b, olwork = lwork)

    if (present(olwork)) then
        olwork = lwork
        return
    end if

    ! Local Memory Allocation
    if (present(work)) then
        if (size(work) < lwork) then
            call report_array_size_error("solve_lq_vec", errmgr, "work", &
                lwork, size(work))
            return
        end if
        wptr => work(1:lwork)
    else
        allocate(wrk(lwork), stat = istat)
        if (istat /= 0) then
            call report_memory_error("solve_lq_vec", errmgr, istat)
            return
        end if
        wptr => wrk
    end if

    ! Solve the lower triangular system L * Y = B for Y, where Y = Q * X.
    ! The lower triangular system is M-by-M and Y is M-by-NHRS.
    call solve_triangular_system(.false., .false., .true., a(1:m,1:m), &
        b(1:m), errmgr)
    if (errmgr%has_error_occurred()) return

    ! Compute Q**T * Y = X
    call mult_lq(.true., a, tau, b, work = wptr, err = errmgr)
    if (errmgr%has_error_occurred()) return
end subroutine

! ------------------------------------------------------------------------------
subroutine solve_lq_vec_cmplx(a, tau, b, work, olwork, err)
    !! Solves a system of LQ factored equations of the form 
    !! \(A \vec{x} = L Q \vec{x} = \vec{b}\).
    complex(real64), intent(in), dimension(:,:) :: a
        !! !! On input, the M-by-N LQ factored matrix as returned by lq_factor.  
        !! On output, the contents of this matrix are restored.  Notice, N must
        !! be greater than or equal to M.
    complex(real64), intent(in), dimension(:) :: tau
        !! A MIN(M, N)-element array containing the scalar factors of the 
        !! elementary reflectors as returned by lq_factor.
    complex(real64), intent(inout), dimension(:) :: b
        !! On input, an N-element vector where the first M rows contain the 
        !! right-hand-side vector \(\vec{b}\).  On output, the N-element vector
        !! \(\vec{x}\).
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
    if (m > n) then
        call errmgr%report_error("solve_lq_vec_cmplx", &
            "This routine does not handle the overdetermined case.", &
            LA_INVALID_INPUT_ERROR)
        return
    else if (size(tau) /= k) then
        call report_array_size_error("solve_lq_vec_cmplx", errmgr, "tau", k, &
            size(tau))
        return
    else if (size(b) /= n) then
        call report_memory_error("solve_lq_vec_cmplx", errmgr, istat)
        return
    end if

    ! Workspace Query
    call mult_lq(.true., a, tau, b, olwork = lwork)

    if (present(olwork)) then
        olwork = lwork
        return
    end if

    ! Local Memory Allocation
    if (present(work)) then
        if (size(work) < lwork) then
            call report_array_size_error("solve_lq_vec_cmplx", errmgr, "work", &
                lwork, size(work))
            return
        end if
        wptr => work(1:lwork)
    else
        allocate(wrk(lwork), stat = istat)
        if (istat /= 0) then
            call report_memory_error("solve_lq_vec_cmplx", errmgr, istat)
            return
        end if
        wptr => wrk
    end if

    ! Solve the lower triangular system L * Y = B for Y, where Y = Q * X.
    ! The lower triangular system is M-by-M and Y is M-by-NHRS.
    call solve_triangular_system(.false., .false., .true., a(1:m,1:m), &
        b(1:m), errmgr)
    if (errmgr%has_error_occurred()) return

    ! Compute Q**T * Y = X
    call mult_lq(.true., a, tau, b, work = wptr, err = errmgr)
    if (errmgr%has_error_occurred()) return
end subroutine

! ------------------------------------------------------------------------------
end module