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
pure subroutine lq_factor_no_pivot(a, tau, lq, l, q)
    !! Computes the LQ factorization of an M-by-N matrix \(A = L Q\) where
    !! \(L\) is a lower triangular (or lower trapezoidal) matrix and \(Q\) is
    !! a orthogonal matrix.
    real(real64), intent(in), dimension(:,:) :: a
        !! The M-by-N matrix to factor.
    real(real64), intent(out), allocatable, optional, target, dimension(:) :: tau
        !! A MIN(M, N)-element array used to store the scalar factors of the 
        !! elementary reflectors.
    real(real64), intent(out), allocatable, optional, target, dimension(:,:) :: lq
        !! An M-by-N matrix with the elements below the diagonal containing the
        !! MIN(M,N)-by-N lower trapezoidal matrix \(L\) (\(L\) is lower 
        !! triangluar if M >= N).  The elements above the diagonal, along with
        !! the array tau, represent the orthogonal matrix \(Q\) as a product
        !! of elementary reflectors.
    real(real64), intent(out), allocatable, optional, target, dimension(:,:) :: l
        !! The M-by-N lower trapezoidal matrix \(L\).
    real(real64), intent(out), allocatable, optional, target, dimension(:,:) :: q
        !! The N-by-N orthogonal matrix \(Q\).

    ! Local Variables
    logical :: buildlq
    integer(int32) :: m, n, mn, lwork, flag
    real(real64), dimension(1) :: temp
    real(real64), allocatable, target, dimension(:) :: w, tc
    real(real64), allocatable, target, dimension(:,:) :: ac, lc, qc
    real(real64), pointer, dimension(:) :: tptr
    real(real64), pointer, dimension(:,:) :: aptr, lptr, qptr

    ! Initialization
    m = size(a, 1)
    n = size(a, 2)
    mn = min(m, n)
    if (present(tau)) then
        allocate(tau(mn))
        tptr => tau
    else
        allocate(tc(mn))
        tptr => tc
    end if
    if (present(lq)) then
        allocate(lq(m, n), source = a)
        aptr => lq
    else
        allocate(ac(m, n), source = a)
        aptr => ac
    end if
    buildlq = present(l) .or. present(q)

    ! Workspace Query
    call DGELQF(m, n, temp, m, temp, temp, -1, flag)
    lwork = int(temp(1), int32)
    allocate(w(lwork))

    ! Call DGELQF
    call DGELQF(m, n, aptr, m, tptr, w, lwork, flag)

    ! Build L & Q?
    if (buildlq) then
        ! L
        if (present(l)) then
            allocate(l(m, n), source = aptr)
            lptr => l
        else
            if (allocated(ac)) then
                lptr => ac
            else
                allocate(lc(m, n), source = aptr)
                lptr => lc
            end if
        end if

        ! Q
        if (present(q)) then
            allocate(q(n, n))
            qptr => q
        else
            allocate(qc(n, n))
            qptr => q
        end if

        ! Build L & Q
        call form_lq(lptr, tptr, qptr)
    end if
end subroutine

! ------------------------------------------------------------------------------
pure subroutine lq_factor_no_pivot_cmplx(a, tau, lq, l, q)
    !! Computes the LQ factorization of an M-by-N matrix \(A = L Q\) where
    !! \(L\) is a lower triangular (or lower trapezoidal) matrix and \(Q\) is
    !! a orthogonal matrix.
    complex(real64), intent(in), dimension(:,:) :: a
        !! The M-by-N matrix to factor.
    complex(real64), intent(out), allocatable, optional, target, dimension(:) :: tau
        !! A MIN(M, N)-element array used to store the scalar factors of the 
        !! elementary reflectors.
    complex(real64), intent(out), allocatable, optional, target, dimension(:,:) :: lq
        !! An M-by-N matrix with the elements below the diagonal containing the
        !! MIN(M,N)-by-N lower trapezoidal matrix \(L\) (\(L\) is lower 
        !! triangluar if M >= N).  The elements above the diagonal, along with
        !! the array tau, represent the orthogonal matrix \(Q\) as a product
        !! of elementary reflectors.
    complex(real64), intent(out), allocatable, optional, target, dimension(:,:) :: l
        !! The M-by-N lower trapezoidal matrix \(L\).
    complex(real64), intent(out), allocatable, optional, target, dimension(:,:) :: q
        !! The N-by-N orthogonal matrix \(Q\).

    ! Local Variables
    logical :: buildlq
    integer(int32) :: m, n, mn, lwork, flag
    complex(real64), dimension(1) :: temp
    complex(real64), allocatable, target, dimension(:) :: w, tc
    complex(real64), allocatable, target, dimension(:,:) :: ac, lc, qc
    complex(real64), pointer, dimension(:) :: tptr
    complex(real64), pointer, dimension(:,:) :: aptr, lptr, qptr

    ! Initialization
    m = size(a, 1)
    n = size(a, 2)
    mn = min(m, n)
    if (present(tau)) then
        allocate(tau(mn))
        tptr => tau
    else
        allocate(tc(mn))
        tptr => tc
    end if
    if (present(lq)) then
        allocate(lq(m, n), source = a)
        aptr => lq
    else
        allocate(ac(m, n), source = a)
        aptr => ac
    end if
    buildlq = present(l) .or. present(q)

    ! Workspace Query
    call ZGELQF(m, n, temp, m, temp, temp, -1, flag)
    lwork = int(temp(1), int32)
    allocate(w(lwork))

    ! Call ZGELQF
    call ZGELQF(m, n, aptr, m, tptr, w, lwork, flag)

    ! Build L & Q?
    if (buildlq) then
        ! L
        if (present(l)) then
            allocate(l(m, n), source = aptr)
            lptr => l
        else
            if (allocated(ac)) then
                lptr => ac
            else
                allocate(lc(m, n), source = aptr)
                lptr => lc
            end if
        end if

        ! Q
        if (present(q)) then
            allocate(q(n, n))
            qptr => q
        else
            allocate(qc(n, n))
            qptr => q
        end if

        ! Build L & Q
        call form_lq(lptr, tptr, qptr)
    end if
end subroutine

! ------------------------------------------------------------------------------
pure subroutine form_lq_no_pivot(l, tau, q)
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

    ! Parameters
    real(real64), parameter :: zero = 0.0d0

    ! Local Variables
    integer(int32) :: i, j, m, n, mn, k, flag, lwork
    real(real64), allocatable, dimension(:) :: w
    real(real64), dimension(1) :: temp
    
    ! Initialization
    m = size(l, 1)
    n = size(l, 2)
    mn = min(m, n)

    ! Input Check
    if (m > n) then
        error stop LA_INVALID_INPUT_ERROR
    end if
    if (size(tau) /= mn) then
        error stop 2
    end if
    if (size(q, 1) /= n .or. size(q, 2) /= n) then
        error stop 3
    end if

    ! Workspace Query
    call DORGLQ(n, n, mn, temp, n, tau, temp, -1, flag)
    lwork = int(temp(1), int32)
    allocate(w(lwork))

    ! Copy the upper triangular portion of L to Q, and then zero it out in L
    do j = 2, n
        k = min(j - 1, m)
        q(1:k,j) = l(1:k,j)
        l(1:k,j) = zero
    end do

    ! Build Q
    call DORGLQ(n, n, mn, q, n, tau, w, lwork, flag)
end subroutine

! ------------------------------------------------------------------------------
pure subroutine form_lq_no_pivot_cmplx(l, tau, q)
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

    ! Parameters
    complex(real64), parameter :: zero = (0.0d0, 0.0d0)

    ! Local Variables
    integer(int32) :: i, j, m, n, mn, k, flag, lwork
    complex(real64), allocatable, dimension(:) :: w
    complex(real64), dimension(1) :: temp

    ! Initialization
    m = size(l, 1)
    n = size(l, 2)
    mn = min(m, n)

    ! Input Check
    if (m > n) then
        error stop LA_INVALID_INPUT_ERROR
    end if
    if (size(tau) /= mn) then
        error stop 2
    end if
    if (size(q, 1) /= n .or. size(q, 2) /= n) then
        error stop 3
    end if

    ! Workspace Query
    call ZUNGLQ(n, n, mn, temp, n, tau, temp, -1, flag)
    lwork = int(temp(1), int32)
    allocate(w(lwork))

    ! Copy the upper triangular portion of L to Q, and then zero it out in L
    do j = 2, n
        k = min(j - 1, m)
        q(1:k,j) = l(1:k,j)
        l(1:k,j) = zero
    end do

    ! Build Q
    call ZUNGLQ(n, n, mn, q, n, tau, w, lwork, flag)
end subroutine

! ------------------------------------------------------------------------------
pure function mult_lq_mtx(lside, trans, a, tau, c) result(qc)
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
    real(real64), intent(in), dimension(:,:) :: c
        !! The M-by-N matrix \(C\).  
    real(real64), allocatable, dimension(:,:) :: qc
        !! The M-by-N product of the orthogonal \(Q\) and the \(C\).

    ! Local Variables
    character :: side, t
    integer(int32) :: m, n, k, ncola, flag, lwork
    real(real64), allocatable, dimension(:) :: w
    real(real64), dimension(1) :: temp

    ! Initialization
    m = size(c, 1)
    n = size(c, 2)
    k = size(tau)
    allocate(qc(m,n), source = c)
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

    ! Input Check
    if (size(a, 1) /= k .or. size(a, 2) /= ncola) then
        error stop 3
    end if

    ! Workspace Query
    call DORMLQ(side, t, m, n, k, a, k, tau, qc, m, temp, -1, flag)
    lwork = int(temp(1), int32)
    allocate(w(lwork))

    ! Call DORMLQ
    call DORMLQ(side, t, m, n, k, a, k, tau, qc, m, w, lwork, flag)
end function

! ------------------------------------------------------------------------------
pure function mult_lq_mtx_cmplx(lside, trans, a, tau, c) result(qc)
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
    complex(real64), intent(in), dimension(:,:) :: c
        !! The M-by-N matrix \(C\).
    complex(real64), allocatable, dimension(:,:) :: qc
        !! The M-by-N product of the orthogonal \(Q\) and the \(C\).

    ! Local Variables
    character :: side, t
    integer(int32) :: m, n, k, ncola, flag, lwork
    complex(real64), allocatable, dimension(:) :: w
    complex(real64), dimension(1) :: temp

    ! Initialization
    m = size(c, 1)
    n = size(c, 2)
    k = size(tau)
    allocate(qc(m, n), source = c)
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

    ! Input Check
    if (size(a, 1) /= k .or. size(a, 2) /= ncola) then
        error stop 2
    end if

    ! Workspace Query
    call ZUNMLQ(side, t, m, n, k, a, k, tau, qc, m, temp, -1, flag)
    lwork = int(temp(1), int32)
    allocate(w(lwork))

    ! Call ZUNMLQ
    call ZUNMLQ(side, t, m, n, k, a, k, tau, qc, m, w, lwork, flag)
end function

! ------------------------------------------------------------------------------
pure function mult_lq_vec(trans, a, tau, c) result(qc)
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
    real(real64), intent(in), dimension(:) :: c
        !! On input, the M-element vector \(\vec{c}\).
    real(real64), allocatable, dimension(:) :: qc
        !! The M-element product of the orthogonal matrix \(Q\) and the 
        !! vector \(\vec{c}\).

    ! Local Variables
    character :: side, t
    integer(int32) :: m, n, k, istat, flag, lwork
    real(real64), allocatable, dimension(:) :: w
    real(real64), dimension(1) :: temp

    ! Initialization
    m = size(c)
    n = 1
    k = size(tau)
    allocate(qc(m), source = c)
    side = 'L'
    if (trans) then
        t = 'T'
    else
        t = 'N'
    end if

    ! Input Check
    if (size(a, 1) /= k .or. size(a, 2) /= m) then
        error stop 2
    end if

    ! Workspace Query
    call DORMLQ(side, t, m, n, k, a, k, tau, qc, m, temp, -1, flag)
    lwork = int(temp(1), int32)
    allocate(w(lwork))

    ! Call DORMLQ
    call DORMLQ(side, t, m, n, k, a, k, tau, qc, m, w, lwork, flag)
end function

! ------------------------------------------------------------------------------
pure function mult_lq_vec_cmplx(trans, a, tau, c) result(qc)
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
    complex(real64), intent(in), dimension(:) :: c
        !! The M-element vector \(\vec{c}\).
    complex(real64), allocatable, dimension(:) :: qc
        !! The M-element product of the orthogonal matrix \(Q\) and the 
        !! vector \(\vec{c}\).

    ! Local Variables
    character :: side, t
    integer(int32) :: m, n, k, flag, lwork
    complex(real64), allocatable, dimension(:) :: w
    complex(real64), dimension(1) :: temp

    ! Initialization
    m = size(c)
    n = 1
    k = size(tau)
    allocate(qc(m), source = c)
    side = 'L'
    if (trans) then
        t = 'C'
    else
        t = 'N'
    end if

    ! Input Check
    if (size(a, 1) /= k .or. size(a, 2) /= m) then
        error stop 2
    end if

    ! Workspace Query
    call ZUNMLQ(side, t, m, n, k, a, k, tau, qc, m, temp, -1, flag)
    lwork = int(temp(1), int32)
    allocate(w(lwork))

    ! Call ZUNMLQ
    call ZUNMLQ(side, t, m, n, k, a, k, tau, qc, m, w, lwork, flag)
end function

! ------------------------------------------------------------------------------
pure function solve_lq_mtx(a, tau, b) result(x)
    !! Solves a system of LQ factored equations of the form \(A X = L Q X = B\).
    real(real64), intent(in), dimension(:,:) :: a
        !! The M-by-N LQ factored matrix as returned by lq_factor.  Notice, N 
        !! must be greater than or equal to M.
    real(real64), intent(in), dimension(:) :: tau
        !! A MIN(M, N)-element array containing the scalar factors of the 
        !! elementary reflectors as returned by lq_factor.
    real(real64), intent(in), dimension(:,:) :: b
        !! The M-by-NRHS matrix \(B\).
    real(real64), allocatable, dimension(:,:) :: x
        !! The N-by-NRHS matrix \(X\).

    ! Parameters
    real(real64), parameter :: zero = 0.0d0
    real(real64), parameter :: one = 1.0d0

    ! Local Variables
    integer(int32) :: m, n, nrhs, k

    ! Initialization
    m = size(a, 1)
    n = size(a, 2)
    nrhs = size(b, 2)
    k = min(m, n)
    allocate(x(n, nrhs), source = zero)

    ! Input Check
    if (m > n) then
        error stop 1
    else if (size(tau) /= k) then
        error stop 2
    else if (size(b, 1) /= m) then
        error stop 3
    end if

    ! Solve the lower triangular system L * Y = B for Y, where Y = Q * X.
    ! The lower triangular system is M-by-M and Y is M-by-NHRS.
    x(1:m,:) = solve_triangular_system(.true., .false., .false., .true., one, &
        a(1:m,1:m), b)

    ! Compute Q**T * Y = X
    x = mult_lq(.true., .true., a, tau, x)
end function

! ------------------------------------------------------------------------------
pure function solve_lq_mtx_cmplx(a, tau, b) result(x)
    !! Solves a system of LQ factored equations of the form \(A X = L Q X = B\).
    complex(real64), intent(in), dimension(:,:) :: a
        !! The M-by-N LQ factored matrix as returned by lq_factor.  Notice, N 
        !! must be greater than or equal to M.
    complex(real64), intent(in), dimension(:) :: tau
        !! A MIN(M, N)-element array containing the scalar factors of the 
        !! elementary reflectors as returned by lq_factor.
    complex(real64), intent(in), dimension(:,:) :: b
        !! The M-by-NRHS matrix \(B\).
    complex(real64), allocatable, dimension(:,:) :: x
        !! The N-by-NRHS matrix \(X\).

    ! Parameters
    complex(real64), parameter :: zero = (0.0d0, 0.0d0)
    complex(real64), parameter :: one = (1.0d0, 0.0d0)

    ! Local Variables
    integer(int32) :: m, n, nrhs, k

    ! Initialization
    m = size(a, 1)
    n = size(a, 2)
    nrhs = size(b, 2)
    k = min(m, n)
    allocate(x(n, nrhs), source = zero)

    ! Input Check
    if (m > n) then
        error stop 1
    end if
    if (size(tau) /= k) then
        error stop 2
    end if
    if (size(b, 1) /= m) then
        error stop 3
    end if

    ! Solve the lower triangular system L * Y = B for Y, where Y = Q * X.
    ! The lower triangular system is M-by-M and Y is M-by-NHRS.
    x(1:m,:) = solve_triangular_system(.true., .false., .false., .true., one, &
        a(1:m,1:m), b)

    ! Compute Q**T * Y = X
    x = mult_lq(.true., .true., a, tau, x)
end function

! ------------------------------------------------------------------------------
pure function solve_lq_vec(a, tau, b) result(x)
    !! Solves a system of LQ factored equations of the form 
    !! \(A \vec{x} = L Q \vec{x} = \vec{b}\).
    real(real64), intent(in), dimension(:,:) :: a
        !! The M-by-N LQ factored matrix as returned by lq_factor.  Notice, N 
        !! must be greater than or equal to M.
    real(real64), intent(in), dimension(:) :: tau
        !! A MIN(M, N)-element array containing the scalar factors of the 
        !! elementary reflectors as returned by lq_factor.
    real(real64), intent(in), dimension(:) :: b
        !! The M-element vector \(\vec{b}\).
    real(real64), allocatable, dimension(:) :: x
        !! The N-element vector \(\vec{x}\).

    ! Local Variables
    integer(int32) :: m, n, k

    ! Initialization
    m = size(a, 1)
    n = size(a, 2)
    k = min(m, n)
    allocate(x(n), source = 0.0d0)

    ! Solve the lower triangular system L * Y = B for Y, where Y = Q * X.
    ! The lower triangular system is M-by-M and Y is M-by-NHRS.
    x(1:m) = solve_triangular_system(.false., .false., .true., a(1:m,1:m), b)

    ! Compute Q**T * Y = X
    x = mult_lq(.true., a, tau, x)
end function

! ------------------------------------------------------------------------------
pure function solve_lq_vec_cmplx(a, tau, b) result(x)
    !! Solves a system of LQ factored equations of the form 
    !! \(A \vec{x} = L Q \vec{x} = \vec{b}\).
    complex(real64), intent(in), dimension(:,:) :: a
        !! The M-by-N LQ factored matrix as returned by lq_factor.  Notice, N 
        !! must be greater than or equal to M.
    complex(real64), intent(in), dimension(:) :: tau
        !! A MIN(M, N)-element array containing the scalar factors of the 
        !! elementary reflectors as returned by lq_factor.
    complex(real64), intent(in), dimension(:) :: b
        !! The M-element vector \(\vec{b}\).
    complex(real64), allocatable, dimension(:) :: x
        !! The N-element vector \(\vec{x}\).

    ! Local Variables
    integer(int32) :: m, n, k

    ! Initialization
    m = size(a, 1)
    n = size(a, 2)
    k = min(m, n)

    ! Input Check
    if (m > n) then
        error stop 1
    end if
    if (size(tau) /= k) then
        error stop 2
    end if
    if (size(b) /= m) then
        error stop 3
    end if

    ! Solve the lower triangular system L * Y = B for Y, where Y = Q * X.
    ! The lower triangular system is M-by-M and Y is M-by-NHRS.
    x(1:m) = solve_triangular_system(.false., .false., .true., a(1:m,1:m), b)

    ! Compute Q**T * Y = X
    x = mult_lq(.true., a, tau, x)
end function

! ------------------------------------------------------------------------------
end module