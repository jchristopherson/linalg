module linalg_rz
    use iso_fortran_env, only : int32, real64
    use linalg_errors
    use lapack
    use ferror
    implicit none
    private
    public :: rz_factor
    public :: mult_rz

    interface rz_factor
        module procedure :: rz_factor_dbl
        module procedure :: rz_factor_cmplx
    end interface

    interface mult_rz
        module procedure :: mult_rz_mtx
        module procedure :: mult_rz_mtx_cmplx
        module procedure :: mult_rz_vec
        module procedure :: mult_rz_vec_cmplx
    end interface

contains
! ------------------------------------------------------------------------------
subroutine rz_factor_dbl(a, tau, work, olwork, err)
    !! Factors an upper trapezoidal matrix by means of orthogonal 
    !! transformations such that \(A = R Z = (R 0) Z \). \(Z\) is an orthogonal
    !! matrix of dimension N-by-N, and \(R\) is an M-by-M upper triangular
    !! matrix.
    real(real64), intent(inout), dimension(:,:) :: a
        !! On input, the M-by-N upper trapezoidal matrix to factor.  On output,
        !! the leading M-by-M upper triangular part of the matrix contains the 
        !! upper triangular matrix \(R\), and elements N-L+1 to N of the
        !! first M rows of \(A\), with the array tau, represent the orthogonal
        !! matrix \(Z\) as a product of M elementary reflectors.
    real(real64), intent(out), dimension(:) :: tau
        !! An M-element array used to store the scalar factors of the 
        !! elementary reflectors.
    real(real64), intent(out), target, optional, dimension(:) :: work
        !! An optional input, that if provided, prevents any local memory 
        !! allocation.  If not provided, the memory required is allocated 
        !! within.  If provided, the length of the array must be at least
        !! olwork.
    integer(int32), intent(out), optional :: olwork
        !! An optional output used to determine workspace size.  If supplied, 
        !! the routine determines the optimal size for @p work, and returns 
        !! without performing any actual calculations.
    class(errors), intent(inout), optional, target :: err
        !! The error object to be updated.

    ! Local Variables
    integer(int32) :: m, n, lwork, flag, istat
    real(real64), pointer, dimension(:) :: wptr
    real(real64), allocatable, target, dimension(:) :: wrk
    real(real64), dimension(1) :: temp
    class(errors), pointer :: errmgr
    type(errors), target :: deferr

    ! Initialization
    m = size(a, 1)
    n = size(a, 2)
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if

    ! Input Check
    if (size(tau) /= m) then
        call report_array_size_error("rz_factor_dbl", errmgr, "tau", m, &
            size(tau))
        return
    end if

    ! Workspace Query
    call DTZRZF(m, n, a, m, tau, temp, -1, flag)
    lwork = int(temp(1), int32)
    if (present(olwork)) then
        olwork = lwork
        return
    end if

    ! Local Memory Allocation
    if (present(work)) then
        if (size(work) < lwork) then
            call report_array_size_error("rz_factor_dbl", errmgr, "lwork", &
                lwork, size(work))
            return
        end if
        wptr => work(1:lwork)
    else
        allocate(wrk(lwork), stat = istat)
        if (istat /= 0) then
            call report_memory_error("rz_factor_dbl", errmgr, istat)
            return
        end if
        wptr => wrk
    end if

    ! Call DTZRZF
    call DTZRZF(m, n, a, m, tau, wptr, lwork, flag)
end subroutine


! ------------------------------------------------------------------------------
subroutine rz_factor_cmplx(a, tau, work, olwork, err)
    !! Factors an upper trapezoidal matrix by means of orthogonal 
    !! transformations such that \(A = R Z = (R 0) Z \). \(Z\) is an orthogonal
    !! matrix of dimension N-by-N, and \(R\) is an M-by-M upper triangular
    !! matrix.
    complex(real64), intent(inout), dimension(:,:) :: a
        !! On input, the M-by-N upper trapezoidal matrix to factor.  On output,
        !! the leading M-by-M upper triangular part of the matrix contains the 
        !! upper triangular matrix \(R\), and elements N-L+1 to N of the
        !! first M rows of \(A\), with the array tau, represent the orthogonal
        !! matrix \(Z\) as a product of M elementary reflectors.
    complex(real64), intent(out), dimension(:) :: tau
        !! An M-element array used to store the scalar factors of the 
        !! elementary reflectors.
    complex(real64), intent(out), target, optional, dimension(:) :: work
        !! An optional input, that if provided, prevents any local memory 
        !! allocation.  If not provided, the memory required is allocated 
        !! within.  If provided, the length of the array must be at least
        !! olwork.
    integer(int32), intent(out), optional :: olwork
        !! An optional output used to determine workspace size.  If supplied, 
        !! the routine determines the optimal size for @p work, and returns 
        !! without performing any actual calculations.
    class(errors), intent(inout), optional, target :: err
        !! The error object to be updated.

    ! Local Variables
    integer(int32) :: m, n, lwork, flag, istat
    complex(real64), pointer, dimension(:) :: wptr
    complex(real64), allocatable, target, dimension(:) :: wrk
    complex(real64), dimension(1) :: temp
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    character(len = :), allocatable :: errmsg

    ! Initialization
    m = size(a, 1)
    n = size(a, 2)
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if

    ! Input Check
    flag = 0
    if (size(tau) /= m) then
        call report_array_size_error("rz_factor_cmplx", errmgr, "tau", m, &
            size(tau))
        return
    end if

    ! Workspace Query
    call ZTZRZF(m, n, a, m, tau, temp, -1, flag)
    lwork = int(temp(1), int32)
    if (present(olwork)) then
        olwork = lwork
        return
    end if

    ! Local Memory Allocation
    if (present(work)) then
        if (size(work) < lwork) then
            call report_array_size_error("rz_factor_cmplx", errmgr, "work", &
                lwork, size(work))
            return
        end if
        wptr => work(1:lwork)
    else
        allocate(wrk(lwork), stat = istat)
        if (istat /= 0) then
            call report_memory_error("rz_factor_cmplx", errmgr, istat)
            return
        end if
        wptr => wrk
    end if

    ! Call ZTZRZF
    call ZTZRZF(m, n, a, m, tau, wptr, lwork, flag)
end subroutine

! ------------------------------------------------------------------------------
subroutine mult_rz_mtx(lside, trans, l, a, tau, c, work, olwork, err)
    !! Multiplies a general matrix by the orthogonal matrix Z from an 
    !! RZ factorization such that \(C = op(Z) C\) or \(C = C op(Z)\)
    logical, intent(in) :: lside
        !! Set to true to compute \(C = op(Z) C\); else, set to false to 
        !! compute \(C = C op(Z)\).
    logical, intent(in) :: trans
        !! Set to true if \(op(Z) = Z^{T}\); else, set to false if 
        !! \(op(Z) = Z\).
    integer(int32), intent(in) :: l
        !! The number of columns in matrix \(A\) containing the meaningful part 
        !! of the Householder vectors.  If lside is true, \(M \ge L \ge 0\); 
        !! else, if lside is false, \(N \ge L \ge 0\).
    real(real64), intent(inout), dimension(:,:) :: a
        !! On input the \(K\)-by-\(LTA\) matrix \(Z\), where \(LTA = M\) if 
        !! lside is true; else, \(LTA = N\) if lside is false.  The I-th row 
        !! must contain the Householder vector in the last \(k\) rows. Notice, 
        !! the contents of this matrix are restored on exit.
    real(real64), intent(inout), dimension(:,:) :: c
        !! On input, the \(M\)-by-\(N\) matrix \(C\).  On output, the product 
        !! of the orthogonal matrix \(Z\) and the original matrix \(C\).
    real(real64), intent(in), dimension(:) :: tau
        !! A \(K\)-element array containing the scalar factors of the elementary 
        !! reflectors, where \(M \ge K \ge 0\) if lside is true; else,
        !! \(N \ge K \ge 0\) if lside is false.
    real(real64), intent(out), target, optional, dimension(:) :: work
        !! An optional input, that if provided, prevents any local memory 
        !! allocation.  If not provided, the memory required is allocated 
        !! within.  If provided, the length of the array must be at least
        !! olwork.
    integer(int32), intent(out), optional :: olwork
        !! An optional output used to determine workspace size.  If supplied, 
        !! the routine determines the optimal size for @p work, and returns 
        !! without performing any actual calculations.
    class(errors), intent(inout), optional, target :: err
        !! The error object to be updated.

    ! Local Variables
    character :: side, t
    integer(int32) :: m, n, k, lwork, flag, istat, lda
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
    lda = size(a, 1)
    if (lside) then
        side = 'L'
    else
        side = 'R'
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
    flag = 0
    if (lside) then
        if (l > m .or. l < 0) then
           flag = 3
        else if (k > m) then
            flag = 5
        else if (size(a, 1) < k .or. size(a, 2) /= m) then
            flag = 4
        end if
    else
        if (l > n .or. l < 0) then
            flag = 3
        else if (k > n) then
            flag = 5
        else if (size(a, 1) < k .or. size(a, 2) /= n) then
            flag = 4
        end if
    end if
    if (flag /= 0) then
        ! ERROR: One of the input arrays is not sized correctly
        allocate(character(len = 256) :: errmsg)
        write(errmsg, 100) "Input number ", flag, &
            " is not sized correctly."
        call errmgr%report_error("mult_rz_mtx", trim(errmsg), &
            LA_ARRAY_SIZE_ERROR)
        return
    end if

    ! Workspace Query
    call DORMRZ(side, t, m, n, k, l, a, lda, tau, c, m, temp, -1, flag)
    lwork = int(temp(1), int32)
    if (present(olwork)) then
        olwork = lwork
        return
    end if

    ! Local Memory Allocation
    if (present(work)) then
        if (size(work) < lwork) then
            call report_array_size_error("mult_rz_mtx", errmgr, "work", lwork, &
                size(work))
            return
        end if
        wptr => work(1:lwork)
    else
        allocate(wrk(lwork), stat = istat)
        if (istat /= 0) then
            call report_memory_error("mult_rz_mtx", errmgr, istat)
            return
        end if
        wptr => wrk
    end if

    ! Call DORMRZ
    call DORMRZ(side, t, m, n, k, l, a, lda, tau, c, m, wptr, lwork, flag)

    ! Formatting
100 format(A, I0, A)
end subroutine

! ------------------------------------------------------------------------------
subroutine mult_rz_mtx_cmplx(lside, trans, l, a, tau, c, work, olwork, err)
    !! Multiplies a general matrix by the orthogonal matrix Z from an 
    !! RZ factorization such that \(C = op(Z) C\) or \(C = C op(Z)\).
    logical, intent(in) :: lside
        !! Set to true to compute \(C = op(Z) C\); else, set to false to 
        !! compute \(C = C op(Z)\).
    logical, intent(in) :: trans
        !! Set to true if \(op(Z) = Z^{T}\); else, set to false if 
        !! \(op(Z) = Z\).
    integer(int32), intent(in) :: l
        !! The number of columns in matrix \(A\) containing the meaningful part 
        !! of the Householder vectors.  If lside is true, \(M \ge L \ge 0\); 
        !! else, if lside is false, \(N \ge L \ge 0\).
    complex(real64), intent(inout), dimension(:,:) :: a
        !! On input the \(K\)-by-\(LTA\) matrix \(Z\), where \(LTA = M\) if 
        !! lside is true; else, \(LTA = N\) if lside is false.  The I-th row 
        !! must contain the Householder vector in the last \(k\) rows. Notice, 
        !! the contents of this matrix are restored on exit.
    complex(real64), intent(inout), dimension(:,:) :: c
        !! On input, the \(M\)-by-\(N\) matrix \(C\).  On output, the product 
        !! of the orthogonal matrix \(Z\) and the original matrix \(C\).
    complex(real64), intent(in), dimension(:) :: tau
        !! A \(K\)-element array containing the scalar factors of the elementary 
        !! reflectors, where \(M \ge K \ge 0\) if lside is true; else,
        !! \(N \ge K \ge 0\) if lside is false.
    complex(real64), intent(out), target, optional, dimension(:) :: work
        !! An optional input, that if provided, prevents any local memory 
        !! allocation.  If not provided, the memory required is allocated 
        !! within.  If provided, the length of the array must be at least
        !! olwork.
    integer(int32), intent(out), optional :: olwork
        !! An optional output used to determine workspace size.  If supplied, 
        !! the routine determines the optimal size for @p work, and returns 
        !! without performing any actual calculations.
    class(errors), intent(inout), optional, target :: err
        !! The error object to be updated.

    ! Local Variables
    character :: side, t
    integer(int32) :: m, n, k, lwork, flag, istat, lda
    complex(real64), pointer, dimension(:) :: wptr
    complex(real64), allocatable, target, dimension(:) :: wrk
    complex(real64), dimension(1) :: temp
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    character(len = :), allocatable :: errmsg

    ! Initialization
    m = size(c, 1)
    n = size(c, 2)
    k = size(tau)
    lda = size(a, 1)
    if (lside) then
        side = 'L'
    else
        side = 'R'
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
    flag = 0
    if (lside) then
        if (l > m .or. l < 0) then
           flag = 3
        else if (k > m) then
            flag = 5
        else if (size(a, 1) < k .or. size(a, 2) /= m) then
            flag = 4
        end if
    else
        if (l > n .or. l < 0) then
            flag = 3
        else if (k > n) then
            flag = 5
        else if (size(a, 1) < k .or. size(a, 2) /= n) then
            flag = 4
        end if
    end if
    if (flag /= 0) then
        ! ERROR: One of the input arrays is not sized correctly
        allocate(character(len = 256) :: errmsg)
        write(errmsg, 100) "Input number ", flag, &
            " is not sized correctly."
        call errmgr%report_error("mult_rz_mtx_cmplx", trim(errmsg), &
            LA_ARRAY_SIZE_ERROR)
        return
    end if

    ! Workspace Query
    call ZUNMRZ(side, t, m, n, k, l, a, lda, tau, c, m, temp, -1, flag)
    lwork = int(temp(1), int32)
    if (present(olwork)) then
        olwork = lwork
        return
    end if

    ! Local Memory Allocation
    if (present(work)) then
        if (size(work) < lwork) then
            call report_array_size_error("mult_rz_mtx_cmplx", errmgr, "work", &
                lwork, size(work))
            return
        end if
        wptr => work(1:lwork)
    else
        allocate(wrk(lwork), stat = istat)
        if (istat /= 0) then
            call report_memory_error("mult_rz_mtx_cmplx", errmgr, istat)
            return
        end if
        wptr => wrk
    end if

    ! Call ZUNMRZ
    call ZUNMRZ(side, t, m, n, k, l, a, lda, tau, c, m, wptr, lwork, flag)

    ! Formatting
100 format(A, I0, A)
end subroutine

! ------------------------------------------------------------------------------
subroutine mult_rz_vec(trans, l, a, tau, c, work, olwork, err)
    !! Multiplies a general matrix by the orthogonal matrix Z from an 
    !! RZ factorization such that \(C = op(Z) C\).
    logical, intent(in) :: trans
        !! Set to true if \(op(Z) = Z^{T}\); else, set to false if 
        !! \(op(Z) = Z\).
    integer(int32), intent(in) :: l
        !! The number of columns in matrix \(A\) containing the meaningful part 
        !! of the Householder vectors.
    real(real64), intent(inout), dimension(:,:) :: a
        !! On input the \(M\)-by-\(M\) matrix \(Z\).  The I-th row must contain 
        !! the Householder vector in the last \(k\) rows. Notice, the contents 
        !! of this matrix are restored on exit.
    real(real64), intent(in), dimension(:) :: tau
        !! An \(M\)-element array containing the scalar factors of the
        !! elementary reflectors.
    real(real64), intent(inout), dimension(:) :: c
        !! On input, the \(M\)-element array \(C\).  On output, the product
        !! of \(Z\) and \(C\).
    real(real64), intent(out), target, optional, dimension(:) :: work
        !! An optional input, that if provided, prevents any local memory 
        !! allocation.  If not provided, the memory required is allocated 
        !! within.  If provided, the length of the array must be at least
        !! olwork.
    integer(int32), intent(out), optional :: olwork
        !! An optional output used to determine workspace size.  If supplied, 
        !! the routine determines the optimal size for @p work, and returns 
        !! without performing any actual calculations.
    class(errors), intent(inout), optional, target :: err
        !! The error object to be updated.

    ! Local Variables
    character :: side, t
    integer(int32) :: m, k, lwork, flag, istat, lda
    real(real64), pointer, dimension(:) :: wptr
    real(real64), allocatable, target, dimension(:) :: wrk
    real(real64), dimension(1) :: temp
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    character(len = :), allocatable :: errmsg

    ! Initialization
    m = size(c)
    k = size(tau)
    lda = size(a, 1)
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
    flag = 0
    if (l > m .or. l < 0) then
        flag = 2
    else if (k > m) then
        flag = 4
    else if (size(a, 1) < k .or. size(a, 2) /= m) then
        flag = 3
    end if
    if (flag /= 0) then
        ! ERROR: One of the input arrays is not sized correctly
        allocate(character(len = 256) :: errmsg)
        write(errmsg, 100) "Input number ", flag, &
            " is not sized correctly."
        call errmgr%report_error("mult_rz_vec", trim(errmsg), &
            LA_ARRAY_SIZE_ERROR)
        return
    end if

    ! Workspace Query
    call DORMRZ(side, t, m, 1, k, l, a, lda, tau, c, m, temp, -1, flag)
    lwork = int(temp(1), int32)
    if (present(olwork)) then
        olwork = lwork
        return
    end if

    ! Local Memory Allocation
    if (present(work)) then
        if (size(work) < lwork) then
            call report_array_size_error("mult_rz_vec", errmgr, "work", &
                lwork, size(work))
            return
        end if
        wptr => work(1:lwork)
    else
        allocate(wrk(lwork), stat = istat)
        if (istat /= 0) then
            call report_memory_error("mult_rz_vec", errmgr, istat)
            return
        end if
        wptr => wrk
    end if

    ! Call DORMRZ
    call DORMRZ(side, t, m, 1, k, l, a, lda, tau, c, m, wptr, lwork, flag)

    ! Formatting
100 format(A, I0, A)
end subroutine

! ------------------------------------------------------------------------------
subroutine mult_rz_vec_cmplx(trans, l, a, tau, c, work, olwork, err)
    !! Multiplies a general matrix by the orthogonal matrix Z from an 
    !! RZ factorization such that \(C = op(Z) C\).
    logical, intent(in) :: trans
        !! Set to true if \(op(Z) = Z^{T}\); else, set to false if 
        !! \(op(Z) = Z\).
    integer(int32), intent(in) :: l
        !! The number of columns in matrix \(A\) containing the meaningful part 
        !! of the Householder vectors.
    complex(real64), intent(inout), dimension(:,:) :: a
        !! On input the \(M\)-by-\(M\) matrix \(Z\).  The I-th row must contain 
        !! the Householder vector in the last \(k\) rows. Notice, the contents 
        !! of this matrix are restored on exit.
    complex(real64), intent(in), dimension(:) :: tau
        !! An \(M\)-element array containing the scalar factors of the
        !! elementary reflectors.
    complex(real64), intent(inout), dimension(:) :: c
        !! On input, the \(M\)-element array \(C\).  On output, the product
        !! of \(Z\) and \(C\).
    complex(real64), intent(out), target, optional, dimension(:) :: work
        !! An optional input, that if provided, prevents any local memory 
        !! allocation.  If not provided, the memory required is allocated 
        !! within.  If provided, the length of the array must be at least
        !! olwork.
    integer(int32), intent(out), optional :: olwork
        !! An optional output used to determine workspace size.  If supplied, 
        !! the routine determines the optimal size for @p work, and returns 
        !! without performing any actual calculations.
    class(errors), intent(inout), optional, target :: err
        !! The error object to be updated.

    ! Local Variables
    character :: side, t
    integer(int32) :: m, k, lwork, flag, istat, lda
    complex(real64), pointer, dimension(:) :: wptr
    complex(real64), allocatable, target, dimension(:) :: wrk
    complex(real64), dimension(1) :: temp
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    character(len = :), allocatable :: errmsg

    ! Initialization
    m = size(c)
    k = size(tau)
    lda = size(a, 1)
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
    flag = 0
    if (l > m .or. l < 0) then
        flag = 2
    else if (k > m) then
        flag = 4
    else if (size(a, 1) < k .or. size(a, 2) /= m) then
        flag = 3
    end if
    if (flag /= 0) then
        ! ERROR: One of the input arrays is not sized correctly
        allocate(character(len = 256) :: errmsg)
        write(errmsg, 100) "Input number ", flag, &
            " is not sized correctly."
        call errmgr%report_error("mult_rz_vec_cmplx", trim(errmsg), &
            LA_ARRAY_SIZE_ERROR)
        return
    end if

    ! Workspace Query
    call ZUNMRZ(side, t, m, 1, k, l, a, lda, tau, c, m, temp, -1, flag)
    lwork = int(temp(1), int32)
    if (present(olwork)) then
        olwork = lwork
        return
    end if

    ! Local Memory Allocation
    if (present(work)) then
        if (size(work) < lwork) then
            call report_array_size_error("mult_rz_vec_cmplx", errmgr, "work", &
                lwork, size(work))
            return
        end if
        wptr => work(1:lwork)
    else
        allocate(wrk(lwork), stat = istat)
        if (istat /= 0) then
            call report_memory_error("mult_rz_vec_cmplx", errmgr, istat)
            return
        end if
        wptr => wrk
    end if

    ! Call ZUNMRZ
    call ZUNMRZ(side, t, m, 1, k, l, a, lda, tau, c, m, wptr, lwork, flag)

    ! Formatting
100     format(A, I0, A)
end subroutine

! ------------------------------------------------------------------------------
end module