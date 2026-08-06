module linalg_qr
    use iso_fortran_env
    use linalg_errors
    use linalg_rz
    use linalg_tri
    use lapack
    use blas
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
pure subroutine qr_factor_no_pivot(a, tau, qr, q, r)
    !! Computes the QR factorization of an M-by-N matrix.
    real(real64), intent(in), dimension(:,:) :: a
        !! The M-by-N matrix to factor.
    real(real64), intent(out), allocatable, optional, target, dimension(:) :: tau
        !! A MIN(M, N)-element array used to store the scalar factors of the 
        !! elementary reflectors.
    real(real64), intent(out), allocatable, optional, target, dimension(:,:) :: qr
        !! The M-by-N factored matrix stored such that the elements on and above
        !! the diagonal contain the MIN(M, N)-by-N upper trapezoidal matrix 
        !! \(R\) (\(R\) is upper triangular if M >= N).  The elements below the
        !! diagonal, along with the array tau, represent the orthogonal matrix
        !! \(Q\) as a product of elementary reflectors.
    real(real64), intent(out), allocatable, optional, target, dimension(:,:) :: q
        !! The M-by-M orthogonal matrix \(Q\).
    real(real64), intent(out), allocatable, optional, target, dimension(:,:) :: r
        !! The M-by-N upper trapezoidal matrix \(R\).

    ! Local Variables
    logical :: buildqr
    integer(int32) :: m, n, mn, lwork, flag
    real(real64), dimension(1) :: temp
    real(real64), allocatable, target, dimension(:) :: w, tc
    real(real64), allocatable, target, dimension(:,:) :: ac, qc, rc
    real(real64), pointer, dimension(:) :: tptr
    real(real64), pointer, dimension(:,:) :: aptr, qptr, rptr

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
    if (present(qr)) then
        allocate(qr(m, n), source = a)
        aptr => qr
    else
        allocate(ac(m, n), source = a)
        aptr => ac
    end if
    buildqr = present(q) .or. present(r)

    ! Workspace Query
    call DGEQRF(m, n, temp, m, temp, temp, -1, flag)
    lwork = int(temp(1), int32)
    allocate(w(lwork))

    ! Call DGEQRF
    call DGEQRF(m, n, aptr, m, tptr, w, lwork, flag)

    ! Build Q & R
    if (buildqr) then
        ! Q
        if (present(q)) then
            allocate(q(m, m))
            qptr => q
        else
            allocate(qc(m, m))
            qptr => qc
        end if

        ! R
        if (present(r)) then
            allocate(r(m, n), source = aptr)
            rptr => r
        else
            if (allocated(ac)) then
                rptr => ac
            else
                allocate(rc(m, n), source = aptr)
                rptr => rc
            end if
        end if
        call form_qr(rptr, tptr, qptr)
    end if
end subroutine

! ------------------------------------------------------------------------------
pure subroutine qr_factor_no_pivot_cmplx(a, tau, qr, q, r)
    !! Computes the QR factorization of an M-by-N matrix.
    complex(real64), intent(in), dimension(:,:) :: a
        !! The M-by-N matrix to factor.
    complex(real64), intent(out), allocatable, optional, target, dimension(:) :: tau
        !! A MIN(M, N)-element array used to store the scalar factors of the 
        !! elementary reflectors.
    complex(real64), intent(out), allocatable, optional, target, dimension(:,:) :: qr
        !! The M-by-N factored matrix stored such that the elements on and above
        !! the diagonal contain the MIN(M, N)-by-N upper trapezoidal matrix 
        !! \(R\) (\(R\) is upper triangular if M >= N).  The elements below the
        !! diagonal, along with the array tau, represent the orthogonal matrix
        !! \(Q\) as a product of elementary reflectors.
    complex(real64), intent(out), allocatable, optional, target, dimension(:,:) :: q
        !! The M-by-M orthogonal matrix \(Q\).
    complex(real64), intent(out), allocatable, optional, target, dimension(:,:) :: r
        !! The M-by-N upper trapezoidal matrix \(R\).

    ! Local Variables
    logical :: buildqr
    integer(int32) :: m, n, mn, lwork, flag
    complex(real64), dimension(1) :: temp
    complex(real64), allocatable, target, dimension(:) :: w, tc
    complex(real64), allocatable, target, dimension(:,:) :: ac, qc, rc
    complex(real64), pointer, dimension(:) :: tptr
    complex(real64), pointer, dimension(:,:) :: aptr, qptr, rptr

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
    if (present(qr)) then
        allocate(qr(m, n), source = a)
        aptr => qr
    else
        allocate(ac(m, n), source = a)
        aptr => ac
    end if
    buildqr = present(q) .or. present(r)
    
    ! Workspace Query
    call ZGEQRF(m, n, temp, m, temp, temp, -1, flag)
    lwork = int(temp(1), int32)
    allocate(w(lwork))

    ! Call ZGEQRF
    call ZGEQRF(m, n, aptr, m, tptr, w, lwork, flag)

    ! Build Q & R
    if (buildqr) then
        ! Q
        if (present(q)) then
            allocate(q(m, m))
            qptr => q
        else
            allocate(qc(m, m))
            qptr => qc
        end if

        ! R
        if (present(r)) then
            allocate(r(m, n), source = aptr)
            rptr => r
        else
            if (allocated(ac)) then
                rptr => ac
            else
                allocate(rc(m, n), source = aptr)
                rptr => rc
            end if
        end if
        call form_qr(rptr, tptr, qptr)
    end if
end subroutine

! ------------------------------------------------------------------------------
pure subroutine qr_factor_pivot(a, jpvt, tau, qr, q, r, p)
    !! Computes the QR factorization of an M-by-N matrix using column pivoting
    !! such that \(A P = Q R\).
    real(real64), intent(in), dimension(:,:) :: a
        !! The M-by-N matrix to factor.
    integer(int32), intent(inout), dimension(:) :: jpvt
        !! On input, an N-element array that if JPVT(I) .ne. 0, the I-th column 
        !! of A is permuted to the front of A * P; if JPVT(I) = 0, the I-th 
        !! column of A is a free column.  On output, if JPVT(I) = K, then the 
        !! I-th column of A * P was the K-th column of A.
    real(real64), intent(out), allocatable, optional, target, dimension(:) :: tau
        !! A MIN(M, N)-element array used to store the scalar factors of the 
        !! elementary reflectors.
    real(real64), intent(out), allocatable, optional, target, dimension(:,:) :: qr
        !! The M-by-N factored matrix stored such that the elements on and above
        !! the diagonal contain the MIN(M, N)-by-N upper trapezoidal matrix 
        !! \(R\) (\(R\) is upper triangular if M >= N).  The elements below the
        !! diagonal, along with the array tau, represent the orthogonal matrix
        !! \(Q\) as a product of elementary reflectors.
    real(real64), intent(out), allocatable, optional, target, dimension(:,:) :: q
        !! The M-by-M orthogonal matrix \(Q\).
    real(real64), intent(out), allocatable, optional, target, dimension(:,:) :: r
        !! The M-by-N upper trapezoidal matrix \(R\).
    real(real64), intent(out), allocatable, optional, target, dimension(:,:) :: p
        !! The N-by-N column-pivot tracking matrix \(P\) where \(A P = Q R\).

    ! Local Variables
    logical :: buildqr
    integer(int32) :: m, n, mn, lwork, flag
    real(real64), dimension(1) :: temp
    real(real64), allocatable, target, dimension(:) :: w, tc
    real(real64), allocatable, target, dimension(:,:) :: ac, qc, rc, pc
    real(real64), pointer, dimension(:) :: tptr
    real(real64), pointer, dimension(:,:) :: aptr, qptr, rptr, pptr

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
    if (present(qr)) then
        allocate(qr(m, n), source = a)
        aptr => qr
    else
        allocate(ac(m, n), source = a)
        aptr => ac
    end if
    buildqr = present(q) .or. present(r) .or. present(p)

    ! Input Check
    if (size(jpvt) /= n) then
        error stop 2
    end if

    ! Workspace Query
    call DGEQP3(m, n, temp, m, jpvt, temp, temp, -1, flag)
    lwork = int(temp(1), int32)
    allocate(w(lwork))

    ! Call DGEQP3
    call DGEQP3(m, n, aptr, m, jpvt, tptr, w, lwork, flag)

    ! Build Q & R
    if (buildqr) then
        ! Q
        if (present(q)) then
            allocate(q(m, m))
            qptr => q
        else
            allocate(qc(m, m))
            qptr => qc
        end if

        ! R
        if (present(r)) then
            allocate(r(m, n), source = aptr)
            rptr => r
        else
            if (allocated(ac)) then
                rptr => ac
            else
                allocate(rc(m, n), source = aptr)
                rptr => rc
            end if
        end if

        ! P
        if (present(p)) then
            allocate(p(n, n))
            pptr => p
        else
            allocate(pc(n, n))
            pptr => pc
        end if

        ! Form the matrices
        call form_qr(rptr, tptr, jpvt, qptr, pptr)
    end if
end subroutine

! ------------------------------------------------------------------------------
pure subroutine qr_factor_pivot_cmplx(a, jpvt, tau, qr, q, r, p)
    !! Computes the QR factorization of an M-by-N matrix using column pivoting
    !! such that \(A P = Q R\).
    complex(real64), intent(in), dimension(:,:) :: a
        !! The M-by-N matrix to factor.
    integer(int32), intent(inout), dimension(:) :: jpvt
        !! On input, an N-element array that if JPVT(I) .ne. 0, the I-th column 
        !! of A is permuted to the front of A * P; if JPVT(I) = 0, the I-th 
        !! column of A is a free column.  On output, if JPVT(I) = K, then the 
        !! I-th column of A * P was the K-th column of A.
    complex(real64), intent(out), allocatable, optional, target, dimension(:) :: tau
        !! A MIN(M, N)-element array used to store the scalar factors of the 
        !! elementary reflectors.
    complex(real64), intent(out), allocatable, optional, target, dimension(:,:) :: qr
        !! The M-by-N factored matrix stored such that the elements on and above
        !! the diagonal contain the MIN(M, N)-by-N upper trapezoidal matrix 
        !! \(R\) (\(R\) is upper triangular if M >= N).  The elements below the
        !! diagonal, along with the array tau, represent the orthogonal matrix
        !! \(Q\) as a product of elementary reflectors.
    complex(real64), intent(out), allocatable, optional, target, dimension(:,:) :: q
        !! The M-by-M orthogonal matrix \(Q\).
    complex(real64), intent(out), allocatable, optional, target, dimension(:,:) :: r
        !! The M-by-N upper trapezoidal matrix \(R\).
    complex(real64), intent(out), allocatable, optional, target, dimension(:,:) :: p
        !! The N-by-N column-pivot tracking matrix \(P\) where \(A P = Q R\).

    ! Local Variables
    logical :: buildqr
    integer(int32) :: m, n, mn, lwork, flag
    complex(real64), dimension(1) :: temp
    real(real64), dimension(1) :: rtemp
    complex(real64), allocatable, target, dimension(:) :: w, tc
    real(real64), allocatable, dimension(:) :: rw
    complex(real64), allocatable, target, dimension(:,:) :: ac, qc, rc, pc
    complex(real64), pointer, dimension(:) :: tptr
    complex(real64), pointer, dimension(:,:) :: aptr, qptr, rptr, pptr

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
    if (present(qr)) then
        allocate(qr(m, n), source = a)
        aptr => qr
    else
        allocate(ac(m, n), source = a)
        aptr => ac
    end if
    buildqr = present(q) .or. present(r) .or. present(p)

    ! Input Check
    if (size(jpvt) /= n) then
        error stop 2
    end if

    ! Workspace Query
    call ZGEQP3(m, n, temp, m, jpvt, temp, temp, -1, rtemp, flag)
    lwork = int(temp(1), int32)
    allocate(w(lwork), rw(2 * n))

    ! Call ZGEQP3
    call ZGEQP3(m, n, aptr, m, jpvt, tptr, w, lwork, rw, flag)

    ! Build Q & R
    if (buildqr) then
        ! Q
        if (present(q)) then
            allocate(q(m, m))
            qptr => q
        else
            allocate(qc(m, m))
            qptr => qc
        end if

        ! R
        if (present(r)) then
            allocate(r(m, n), source = aptr)
            rptr => r
        else
            if (allocated(ac)) then
                rptr => ac
            else
                allocate(rc(m, n), source = aptr)
                rptr => rc
            end if
        end if

        ! P
        if (present(p)) then
            allocate(p(n, n))
            pptr => p
        else
            allocate(pc(n, n))
            pptr => pc
        end if

        ! Form the matrices
        call form_qr(rptr, tptr, jpvt, qptr, pptr)
    end if
end subroutine

! ------------------------------------------------------------------------------
pure subroutine form_qr_no_pivot(r, tau, q)
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

    ! Parameters
    real(real64), parameter :: zero = 0.0d0

    ! Local Variables
    integer(int32) :: j, m, n, mn, qcol, flag, lwork
    real(real64), allocatable, dimension(:) :: w
    real(real64), dimension(1) :: temp

    ! Initialization
    m = size(r, 1)
    n = size(r, 2)
    mn = min(m, n)
    qcol = size(q, 2)

    ! Input Check
    if (size(tau) /= mn) then
        error stop 2
    end if
    if (size(q, 1) /= m .or. (qcol /= m .and. qcol /= n)) then
        error stop 3
    end if
    if (qcol == n .and. m < n) then
        error stop 3
    end if

    ! Workspace Query
    call DORGQR(m, qcol, mn, q, m, tau, temp, -1, flag)
    lwork = int(temp(1), int32)
    allocate(w(lwork))

    ! Copy the sub-diagonal portion of R to Q, and then zero out the
    ! sub-diagonal portion of R
    do j = 1, mn
        q(j+1:m,j) = r(j+1:m,j)
        r(j+1:m,j) = zero
    end do

    ! Build Q - Build M-by-M or M-by-N, but M-by-N only for M >= N
    call DORGQR(m, qcol, mn, q, m, tau, w, lwork, flag)
end subroutine

! ------------------------------------------------------------------------------
pure subroutine form_qr_no_pivot_cmplx(r, tau, q)
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

    ! Parameters
    complex(real64), parameter :: zero = (0.0d0, 0.0d0)

    ! Local Variables
    integer(int32) :: j, m, n, mn, qcol, flag, lwork
    complex(real64), allocatable, dimension(:) :: w
    complex(real64), dimension(1) :: temp

    ! Initialization
    m = size(r, 1)
    n = size(r, 2)
    mn = min(m, n)
    qcol = size(q, 2)

    ! Input Check
    if (size(tau) /= mn) then
        error stop 2
    end if
    if (size(q, 1) /= m .or. (qcol /= m .and. qcol /= n)) then
        error stop 3
    end if
    if (qcol == n .and. m < n) then
        error stop 3
    end if

    ! Workspace Query
    call ZUNGQR(m, qcol, mn, q, m, tau, temp, -1, flag)
    lwork = int(temp(1), int32)
    allocate(w(lwork))

    ! Copy the sub-diagonal portion of R to Q, and then zero out the
    ! sub-diagonal portion of R
    do j = 1, mn
        q(j+1:m,j) = r(j+1:m,j)
        r(j+1:m,j) = zero
    end do

    ! Build Q - Build M-by-M or M-by-N, but M-by-N only for M >= N
    call ZUNGQR(m, qcol, mn, q, m, tau, w, lwork, flag)
end subroutine

! ------------------------------------------------------------------------------
pure subroutine form_qr_pivot(r, tau, pvt, q, p)
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

    ! Parameters
    real(real64), parameter :: zero = 0.0d0
    real(real64), parameter :: one = 1.0d0

    ! Local Variables
    integer(int32) :: j, jp, m, n, mn

    ! Initialization
    m = size(r, 1)
    n = size(r, 2)
    mn = min(m, n)

    ! Input Check
    if (size(p, 1) /= n .or. size(p, 2) /= n) then
        error stop 5
    end if

    ! Generate Q and R
    call form_qr_no_pivot(r, tau, q)

    ! Form P
    do j = 1, n
        jp = pvt(j)
        p(:,j) = zero
        p(jp,j) = one
    end do
end subroutine

! ------------------------------------------------------------------------------
pure subroutine form_qr_pivot_cmplx(r, tau, pvt, q, p)
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

    ! Parameters
    complex(real64), parameter :: zero = (0.0d0, 0.0d0)
    complex(real64), parameter :: one = (1.0d0, 0.0d0)

    ! Local Variables
    integer(int32) :: j, jp, m, n, mn, flag

    ! Initialization
    m = size(r, 1)
    n = size(r, 2)
    mn = min(m, n)

    ! Input Check
    if (size(p, 1) /= n .or. size(p, 2) /= n) then
        error stop 5
    end if

    ! Generate Q and R
    call form_qr_no_pivot_cmplx(r, tau, q)

    ! Form P
    do j = 1, n
        jp = pvt(j)
        p(:,j) = zero
        p(jp,j) = one
    end do
end subroutine

! ------------------------------------------------------------------------------
pure function mult_qr_mtx(lside, trans, a, tau, c) result(qc)
    !! Multiplies a general matrix by the orthogonal matrix \(Q\) from a QR
    !! factorization such that \(C = op(Q) C\) or \(C = C op(Q)\).
    logical, intent(in) :: lside
        !! Set to true to apply \(Q\) or \(Q^T\) from the left; else, set to 
        !! false to apply \(Q\) or \(Q^T\) from the right.
    logical, intent(in) :: trans
        !! Set to true to apply \(Q^T\); else, set to false to apply \(Q\).
    real(real64), intent(in), dimension(:,:) :: a
        !! On input, an LDA-by-K matrix containing the elementary reflectors 
        !! output from the QR factorization.  If lside is set to true, LDA = M, 
        !! and M >= K >= 0; else, if lside is set to false, LDA = N, and 
        !! N >= K >= 0.
    real(real64), intent(in), dimension(:) :: tau
        !! A K-element array containing the scalar factors of each elementary 
        !! reflector defined in\(A\).
    real(real64), intent(in), dimension(:,:) :: c
        !! The M-by-N matrix \(C\).
    real(real64), allocatable, dimension(:,:) :: qc
        !! The M-by-N product of \(C\) and \(Q\).

    ! Parameters
    real(real64), parameter :: one = 1.0d0

    ! Local Variables
    character :: side, t
    integer(int32) :: m, n, k, nrowa, flag, lwork
    real(real64), allocatable, dimension(:) :: w
    real(real64), dimension(1) :: temp

    ! Initialization
    m = size(c, 1)
    n = size(c, 2)
    k = size(tau)
    allocate(qc(m, n), source = c)
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

    ! Input Check
    if (lside) then
        ! A is M-by-K, M >= K >= 0
        if (size(a, 1) /= m .or. size(a, 2) < k) then
            error stop 3
        end if
    else
        ! A is N-by-K, N >= K >= 0
        if (size(a, 1) /= n .or. size(a, 2) < k) then
            error stop 3
        end if
    end if

    ! Workspace Query
    call DORMQR(side, t, m, n, k, a, nrowa, tau, qc, m, temp, -1, flag)
    lwork = int(temp(1), int32)
    allocate(w(lwork))

    ! Call DORMQR
    call DORMQR(side, t, m, n, k, a, nrowa, tau, qc, m, w, lwork, flag)
end function

! ------------------------------------------------------------------------------
pure function mult_qr_mtx_cmplx(lside, trans, a, tau, c) result(qc)
    !! Multiplies a general matrix by the orthogonal matrix \(Q\) from a QR
    !! factorization such that \(C = op(Q) C\) or \(C = C op(Q)\).
    logical, intent(in) :: lside
        !! Set to true to apply \(Q\) or \(Q^H\) from the left; else, set to 
        !! false to apply \(Q\) or \(Q^H\) from the right.
    logical, intent(in) :: trans
        !! Set to true to apply \(Q^H\); else, set to false to apply \(Q\).
    complex(real64), intent(in), dimension(:,:) :: a
        !! On input, an LDA-by-K matrix containing the elementary reflectors 
        !! output from the QR factorization.  If lside is set to true, LDA = M, 
        !! and M >= K >= 0; else, if lside is set to false, LDA = N, and 
        !! N >= K >= 0.
    complex(real64), intent(in), dimension(:) :: tau
        !! A K-element array containing the scalar factors of each elementary 
        !! reflector defined in\(A\).
    complex(real64), intent(in), dimension(:,:) :: c
        !! The M-by-N matrix \(C\).
    complex(real64), allocatable, dimension(:,:) :: qc
        !! The M-by-N product of \(C\) and \(Q\).

    ! Parameters
    complex(real64), parameter :: one = (1.0d0, 0.0d0)

    ! Local Variables
    character :: side, t
    integer(int32) :: m, n, k, nrowa, flag, lwork
    complex(real64), allocatable, dimension(:) :: w
    complex(real64), dimension(1) :: temp

    ! Initialization
    m = size(c, 1)
    n = size(c, 2)
    k = size(tau)
    allocate(qc(m, n), source = c)
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

    ! Input Check
    if (lside) then
        ! A is M-by-K, M >= K >= 0
        if (size(a, 1) /= m .or. size(a, 2) < k) then
            error stop 3
        end if
    else
        ! A is N-by-K, N >= K >= 0
        if (size(a, 1) /= n .or. size(a, 2) < k) then
            error stop 3
        end if
    end if

    ! Workspace Query
    call ZUNMQR(side, t, m, n, k, a, nrowa, tau, qc, m, temp, -1, flag)
    lwork = int(temp(1), int32)
    allocate(w(lwork))

    ! Call ZUNMQR
    call ZUNMQR(side, t, m, n, k, a, nrowa, tau, qc, m, w, lwork, flag)
end function

! ------------------------------------------------------------------------------
pure function mult_qr_vec(trans, a, tau, c) result(qc)
    !! Multiplies a vector by the orthogonal matrix \(Q\) from a QR 
    !! factorization such that \(\vec{c} = op(Q) \vec{c}\).
    logical, intent(in) :: trans
        !! Set to true to apply \(Q^T\); else, set to false to apply \(Q\).
    real(real64), intent(in), dimension(:,:) :: a
        !! On input, an M-by-K matrix containing the elementary reflectors 
        !! output from the QR factorization.
    real(real64), intent(in), dimension(:) :: tau
        !! A K-element array containing the scalar factors of each elementary 
        !! reflector defined in\(A\).
    real(real64), intent(in), dimension(:) :: c
        !! The M-element vector \(\vec{c}\).
    real(real64), allocatable, dimension(:) :: qc
        !! The product of the orthogonal matrix \(Q\) and the original vector 
        !! \(\vec{c}\).

    ! Parameters
    real(real64), parameter :: one = 1.0d0

    ! Local Variables
    character :: side, t
    integer(int32) :: m, k, nrowa, flag, lwork
    real(real64), allocatable, dimension(:) :: w
    real(real64), dimension(1) :: temp

    ! Initialization
    m = size(c)
    k = size(tau)
    side = 'L'
    nrowa = m
    allocate(qc(m), source = c)
    if (trans) then
        t = 'T'
    else
        t = 'N'
    end if

    ! Input Check
    flag = 0
    if (size(a, 1) /= m .or. size(a, 2) < k) then
        error stop 2
    end if

    ! Workspace Query
    call DORMQR(side, t, m, 1, k, a, nrowa, tau, qc, m, temp, -1, flag)
    lwork = int(temp(1), int32)
    allocate(w(lwork))

    ! Call DORMQR
    call DORMQR(side, t, m, 1, k, a, nrowa, tau, qc, m, w, lwork, flag)
end function

! ------------------------------------------------------------------------------
pure function mult_qr_vec_cmplx(trans, a, tau, c) result(qc)
    !! Multiplies a vector by the orthogonal matrix \(Q\) from a QR 
    !! factorization such that \(\vec{c} = op(Q) \vec{c}\).
    logical, intent(in) :: trans
        !! Set to true to apply \(Q^H\); else, set to false to apply \(Q\).
    complex(real64), intent(in), dimension(:,:) :: a
        !! On input, an M-by-K matrix containing the elementary reflectors 
        !! output from the QR factorization. Notice, the contents of this matrix
        !! are restored on exit.
    complex(real64), intent(in), dimension(:) :: tau
        !! A K-element array containing the scalar factors of each elementary 
        !! reflector defined in\(A\).
    complex(real64), intent(in), dimension(:) :: c
        !! The M-element vector \(\vec{c}\).
    complex(real64), allocatable, dimension(:) :: qc
        !! The product of the orthogonal matrix \(Q\) and the original vector 
        !! \(\vec{c}\).

    ! Parameters
    complex(real64), parameter :: one = (1.0d0, 0.0d0)

    ! Local Variables
    character :: side, t
    integer(int32) :: m, k, nrowa, flag, lwork
    complex(real64), allocatable, dimension(:) :: w
    complex(real64), dimension(1) :: temp

    ! Initialization
    m = size(c)
    k = size(tau)
    side = 'L'
    nrowa = m
    allocate(qc(m), source = c)
    if (trans) then
        t = 'C'
    else
        t = 'N'
    end if

    ! Input Check
    flag = 0
    if (size(a, 1) /= m .or. size(a, 2) < k) then
        error stop 2
    end if

    ! Workspace Query
    call ZUNMQR(side, t, m, 1, k, a, nrowa, tau, qc, m, temp, -1, flag)
    lwork = int(temp(1), int32)
    allocate(w(lwork))

    ! Call ZUNMQR
    call ZUNMQR(side, t, m, 1, k, a, nrowa, tau, qc, m, w, lwork, flag)
end function

! ------------------------------------------------------------------------------
pure subroutine qr_rank1_update_dbl(q, r, u, v)
    !! Computes the rank-1 update to an M-by-N QR factored matrix \(A\) where
    !! \(M \ge N\), \(A = Q R\), and \(A_1 = A + \vec{u} \vec{v}^T\) such that
    !! \(A_1 = Q_1 R_1\).
    real(real64), intent(inout), dimension(:,:) :: q
        !! On input, the original M-by-K orthogonal matrix \(Q\).  On output, 
        !! the updated matrix \(Q_1\).
    real(real64), intent(inout), dimension(:,:) :: r
        !! On input, the M-by-N matrix \(R\).  On output, the updated matrix 
        !! \(R_1\).
    real(real64), intent(in), dimension(:) :: u
        !! The M-element \(\vec{u}\) update vector.
    real(real64), intent(in), dimension(:) :: v
        !! The N-element \(\vec{v}\) update vector.

    ! Local Variables
    logical :: full
    integer(int32) :: m, n, k, lwork
    real(real64), allocatable, dimension(:) :: w, uc, vc

    ! Initialization
    m = size(u, 1)
    n = size(r, 2)
    k = min(m, n)
    full = size(q, 2) == m
    lwork = 2 * k
    allocate(w(lwork))

    ! Input Check
    if (m < n) then
        error stop 2
    end if
    if (.not.full .and. size(q, 2) /= k) then
        error stop 1
    end if
    if (size(r, 1) /= m) then
        error stop 2
    end if
    if (size(u) /= m) then
        error stop 3
    end if
    if (size(v) /= n) then
        error stop 4
    end if

    ! Local Memory Allocation
    allocate(uc(m), source = u)
    allocate(vc(n), source = v)

    ! Process
    call DQR1UP(m, n, k, q, m, r, m, uc, vc, w)
end subroutine

! ------------------------------------------------------------------------------
pure subroutine qr_rank1_update_cmplx(q, r, u, v)
    !! Computes the rank-1 update to an M-by-N QR factored matrix \(A\) where
    !! \(M \ge N\), \(A = Q R\), and \(A_1 = A + \vec{u} \vec{v}^H\) such that
    !! \(A_1 = Q_1 R_1\).
    complex(real64), intent(inout), dimension(:,:) :: q
        !! On input, the original M-by-K orthogonal matrix \(Q\).  On output, 
        !! the updated matrix \(Q_1\).
    complex(real64), intent(inout), dimension(:,:) :: r
        !! On input, the M-by-N matrix \(R\).  On output, the updated matrix 
        !! \(R_1\).
    complex(real64), intent(in), dimension(:) :: u
        !! On input, the M-element \(\vec{u}\) update vector.  On output, the 
        !! original content of the array is overwritten.
    complex(real64), intent(in), dimension(:) :: v
        !! On input, the N-element \(\vec{v}\) update vector.  On output, the 
        !! original content of the array is overwritten.

    ! Local Variables
    logical :: full
    integer(int32) :: m, n, k, lwork, lrwork
    complex(real64), allocatable, dimension(:) :: w, uc, vc
    real(real64), allocatable, dimension(:) :: rw

    ! Initialization
    m = size(u, 1)
    n = size(r, 2)
    k = min(m, n)
    full = size(q, 2) == m
    lwork = k
    lrwork = k
    allocate(w(lwork), rw(lrwork))

    ! Input Check
    if (m < n) then
        error stop 2
    end if
    if (.not.full .and. size(q, 2) /= k) then
        error stop 1
    end if
    if (size(r, 1) /= m) then
        error stop 2
    end if
    if (size(u) /= m) then
        error stop 3
    end if
    if (size(v) /= n) then
        error stop 4
    end if

    ! Local Memory Allocation
    allocate(uc(m), source = u)
    allocate(vc(n), source = v)

    ! Process
    call ZQR1UP(m, n, k, q, m, r, m, uc, vc, w, rw)
end subroutine

! ------------------------------------------------------------------------------
pure function solve_qr_no_pivot_mtx(a, tau, b) result(x)
    !! Solves a system of M QR-factored equations of N unknowns.  M must be
    !! greater than or equal to N.
    real(real64), intent(in), dimension(:,:) :: a
        !! The M-by-N QR factored matrix as returned by qr_factor.
    real(real64), intent(in), dimension(:) :: tau
        !! A MIN(M, N)-element array containing the scalar factors of the 
        !! elementary reflectors as returned by qr_factor.
    real(real64), intent(in), dimension(:,:) :: b
        !! The M-by-NRHS right-hand-side matrix.
    real(real64), allocatable, dimension(:,:) :: x
        !! The N-by-NRHS solution matrix.

    ! Parameters
    real(real64), parameter :: one = 1.0d0

    ! Local Variables
    integer(int32) :: m, n, nrhs, k
    real(real64), allocatable, dimension(:,:) :: qtb

    ! Initialization
    m = size(a, 1)
    n = size(a, 2)
    nrhs = size(b, 2)
    k = min(m, n)

    ! Input Check
    if (m < n) then
        error stop 1
    end if
    if (size(tau) /= k) then
        error stop 2
    end if
    if (size(b, 1) /= m) then
        error stop 3
    end if

    ! Compute Q**T * B
    qtb = mult_qr(.true., .true., a, tau, b)

    ! Solve the triangular system: A(1:N,1:N)*X = B(1:N,:)
    x = solve_triangular_system(.true., .true., .false., .true., one, &
        a(1:n,1:n), qtb(1:n,:))
end function

! ------------------------------------------------------------------------------
pure function solve_qr_no_pivot_mtx_cmplx(a, tau, b) result(x)
    !! Solves a system of M QR-factored equations of N unknowns.  M must be
    !! greater than or equal to N.
    complex(real64), intent(in), dimension(:,:) :: a
        !! The M-by-N QR factored matrix as returned by qr_factor.
    complex(real64), intent(in), dimension(:) :: tau
        !! A MIN(M, N)-element array containing the scalar factors of the 
        !! elementary reflectors as returned by qr_factor.
    complex(real64), intent(in), dimension(:,:) :: b
        !! The M-by-NRHS right-hand-side matrix.
    complex(real64), allocatable, dimension(:,:) :: x
        !! The N-by-NRHS right-hand-side matrix.

    ! Parameters
    complex(real64), parameter :: one = (1.0d0, 0.0d0)

    ! Local Variables
    integer(int32) :: m, n, nrhs, k
    complex(real64), allocatable, dimension(:,:) :: qtb

    ! Initialization
    m = size(a, 1)
    n = size(a, 2)
    nrhs = size(b, 2)
    k = min(m, n)

    ! Input Check
    if (m < n) then
        error stop 1
    end if
    if (size(tau) /= k) then
        error stop 2
    end if
    if (size(b, 1) /= m) then
        error stop 3
    end if

    ! Compute Q**T * B
    qtb = mult_qr(.true., .true., a, tau, b)

    ! Solve the triangular system: A(1:N,1:N)*X = B(1:N,:)
    x = solve_triangular_system(.true., .true., .false., .true., one, &
        a(1:n,1:n), qtb(1:n,:))
end function

! ------------------------------------------------------------------------------
pure function solve_qr_no_pivot_vec(a, tau, b) result(x)
    !! Solves a system of M QR-factored equations of N unknowns.  M must be
    !! greater than or equal to N.
    real(real64), intent(in), dimension(:,:) :: a
        !! The M-by-N QR factored matrix as returned by qr_factor.
    real(real64), intent(in), dimension(:) :: tau
        !! A MIN(M, N)-element array containing the scalar factors of the 
        !! elementary reflectors as returned by qr_factor.
    real(real64), intent(in), dimension(:) :: b
        !! The M-element right-hand-side vector.
    real(real64), allocatable, dimension(:) :: x
        !! The N-element solution vector.

    ! Local Variables
    integer(int32) :: m, n, k
    real(real64), allocatable, dimension(:) :: qtb

    ! Initialization
    m = size(a, 1)
    n = size(a, 2)
    k = min(m, n)

    ! Input Check
    if (m < n) then
        error stop 1
    end if
    if (size(tau) /= k) then
        error stop 2
    end if
    if (size(b) /= m) then
        error stop 3
    end if

    ! Compute Q**T * B
    qtb = mult_qr(.true., a, tau, b)

    ! Solve the triangular system: A(1:N,1:N)*X = B(1:N)
    x = solve_triangular_system(.true., .false., .true., a(1:n,1:n), qtb(1:n))
end function

! ------------------------------------------------------------------------------
pure function solve_qr_no_pivot_vec_cmplx(a, tau, b) result(x)
    !! Solves a system of M QR-factored equations of N unknowns.  M must be
    !! greater than or equal to N.
    complex(real64), intent(in), dimension(:,:) :: a
        !! The M-by-N QR factored matrix as returned by qr_factor.
    complex(real64), intent(in), dimension(:) :: tau
        !! A MIN(M, N)-element array containing the scalar factors of the 
        !! elementary reflectors as returned by qr_factor.
    complex(real64), intent(in), dimension(:) :: b
        !! The M-element right-hand-side vector.
    complex(real64), allocatable, dimension(:) :: x
        !! The N-element solution vector.

    ! Local Variables
    integer(int32) :: m, n, k
    complex(real64), allocatable, dimension(:) :: qtb

    ! Initialization
    m = size(a, 1)
    n = size(a, 2)
    k = min(m, n)

    ! Input Check
    if (m < n) then
        error stop 1
    end if
    if (size(tau) /= k) then
        error stop 2
    end if
    if (size(b) /= m) then
        error stop 3
    end if

    ! Compute Q**T * B
    qtb = mult_qr(.true., a, tau, b)

    ! Solve the triangular system: A(1:N,1:N)*X = B(1:N)
    x = solve_triangular_system(.true., .false., .true., a(1:n,1:n), qtb(1:n))
end function

! ------------------------------------------------------------------------------
pure function solve_qr_pivot_mtx(a, tau, jpvt, b) result(x)
    !! Solves a system of M QR-factored equations of N unknowns.
    real(real64), intent(in), dimension(:,:) :: a
        !! The M-by-N QR factored matrix as returned by qr_factor.
    real(real64), intent(in), dimension(:) :: tau
        !! A MIN(M, N)-element array containing the scalar factors of the 
        !! elementary reflectors as returned by qr_factor.
    integer(int32), intent(in), dimension(:) :: jpvt
        !! An N-element array, as output by qr_factor, used to track the 
        !! column pivots.
    real(real64), intent(in), dimension(:,:) :: b
        !! The MAX(M, N)-by-NRHS right-hand-side matrix.
    real(real64), allocatable, dimension(:,:) :: x
        !! The N-by-NRHS solution matrix.

    ! Parameters
    integer(int32), parameter :: imin = 2
    integer(int32), parameter :: imax = 1
    real(real64), parameter :: zero = 0.0d0
    real(real64), parameter :: one = 1.0d0

    ! Local Variables
    integer(int32) :: i, j, m, n, mn, nrhs, lwork, ismin, ismax, &
        rnk, maxmn
    real(real64) :: rcond, smax, smin, smaxpr, sminpr, s1, c1, s2, c2
    real(real64), allocatable, dimension(:) :: w, t2
    real(real64), allocatable, dimension(:,:) :: rz, qtb

    ! Initialization
    m = size(a, 1)
    n = size(a, 2)
    mn = min(m, n)
    maxmn = max(m, n)
    nrhs = size(b, 2)
    ismin = 1
    ismax = mn + 1
    rcond = epsilon(rcond)
    lwork = 2 * maxmn + 1
    allocate(w(lwork), source = one)
    allocate(x(n, nrhs))

    ! Input Check
    if (size(tau) /= mn) then
        error stop 2
    end if
    if (size(jpvt) /= n) then
        error stop 3
    end if
    if (size(b, 1) /= m) then
        error stop 4
    end if

    ! Determine the rank of R11 using an incremental condition estimation
    smax = abs(a(1,1))
    smin = smax
    if (abs(a(1,1)) == zero) then
        rnk = 0
        x = zero
        return
    else
        rnk = 1
    end if
    do
        if (rnk < mn) then
            i = rnk + 1
            call DLAIC1(imin, rnk, w(ismin:ismin+rnk-1), smin, &
                a(1:rnk-1,i), a(i,i), sminpr, s1, c1)
            call DLAIC1(imax, rnk, w(ismax:ismax+rnk-1), smax, &
                a(1:rnk-1,i), a(i,i), smaxpr, s2, c2)
            if (smaxpr * rcond <= sminpr) then
                do i = 1, rnk
                    w(ismin+i-1) = s1 * w(ismin+i-1)
                    w(ismax+i-1) = s2 * w(ismax+i-1)
                end do
                w(ismin+rnk) = c1
                w(ismax+rnk) = c2
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
    ! [R11 R12] = [T11 0] * Y
    if (rnk < n) then
        call rz_factor(a(1:rnk,:), t2, rz)
    end if

    ! Compute B(1:m,1:NRHS) = Q**T * B(1:M,1:NRHS)
    qtb = mult_qr(.true., .true., a, tau, b)

    ! Solve the triangular system T11 * X(1:rnk,1:nrhs) = B(1:rnk,1:nrhs)
    if (rnk < n) then
        x(1:rnk,:) = solve_triangular_system(.true., .true., .false., .true., &
            one, rz(1:rnk,1:rnk), qtb(1:rnk,:))
        x(rnk+1:n,:) = zero

        ! Compute X(1:n,1:nrhs) = Y**T * X(1:n,1:nrhs)
        x = mult_rz(.true., .true., n - rnk, rz, t2, x)
    else
        x = solve_triangular_system(.true., .true., .false., .true., one, &
            a(1:rnk,1:rnk), qtb(1:rnk,:))
    end if

    ! Apply the pivoting: X(1:N,1:NRHS) = P * X(1:N,1:NRHS)
    do j = 1, nrhs
        do i = 1, n
            w(jpvt(i)) = x(i,j)
        end do
        x(1:n,j) = w(1:n)
    end do
end function

! ------------------------------------------------------------------------------
pure function solve_qr_pivot_mtx_cmplx(a, tau, jpvt, b) result(x)
    !! Solves a system of M QR-factored equations of N unknowns.
    complex(real64), intent(in), dimension(:,:) :: a
        !! The M-by-N QR factored matrix as returned by qr_factor.
    complex(real64), intent(in), dimension(:) :: tau
        !! A MIN(M, N)-element array containing the scalar factors of the 
        !! elementary reflectors as returned by qr_factor.
    integer(int32), intent(in), dimension(:) :: jpvt
        !! An N-element array, as output by qr_factor, used to track the 
        !! column pivots.
    complex(real64), intent(in), dimension(:,:) :: b
        !! The M-by-NRHS right-hand-side matrix.
    complex(real64), allocatable, dimension(:,:) :: x
        !! The N-by-NRHS solution matrix.

    ! Parameters
    integer(int32), parameter :: imin = 2
    integer(int32), parameter :: imax = 1
    complex(real64), parameter :: zero = (0.0d0, 0.0d0)
    complex(real64), parameter :: one = (1.0d0, 0.0d0)

    ! Local Variables
    integer(int32) :: i, j, m, n, mn, nrhs, lwork, ismin, ismax, &
        rnk, maxmn
    real(real64) :: rcond, smax, smin, smaxpr, sminpr
    complex(real64) :: s1, c1, s2, c2
    complex(real64), allocatable, dimension(:) :: w, t2
    complex(real64), allocatable, dimension(:,:) :: rz, qtb

    ! Initialization
    m = size(a, 1)
    n = size(a, 2)
    mn = min(m, n)
    maxmn = max(m, n)
    nrhs = size(b, 2)
    ismin = 1
    ismax = mn + 1
    rcond = epsilon(rcond)
    lwork = 2 * maxmn + 1
    allocate(w(lwork), source = one)
    allocate(x(n, nrhs))

    ! Input Check
    if (size(tau) /= mn) then
        error stop 2
    end if
    if (size(jpvt) /= n) then
        error stop 3
    end if
    if (size(b, 1) /= m) then
        error stop 4
    end if

    ! Determine the rank of R11 using an incremental condition estimation
    smax = abs(a(1,1))
    smin = smax
    if (abs(a(1,1)) == zero) then
        rnk = 0
        x = zero
        return
    else
        rnk = 1
    end if
    do
        if (rnk < mn) then
            i = rnk + 1
            call ZLAIC1(imin, rnk, w(ismin:ismin+rnk-1), smin, &
                a(1:rnk-1,i), a(i,i), sminpr, s1, c1)
            call ZLAIC1(imax, rnk, w(ismax:ismax+rnk-1), smax, &
                a(1:rnk-1,i), a(i,i), smaxpr, s2, c2)
            if (smaxpr * rcond <= sminpr) then
                do i = 1, rnk
                    w(ismin+i-1) = s1 * w(ismin+i-1)
                    w(ismax+i-1) = s2 * w(ismax+i-1)
                end do
                w(ismin+rnk) = c1
                w(ismax+rnk) = c2
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
    ! [R11 R12] = [T11 0] * Y
    if (rnk < n) then
        call rz_factor(a(1:rnk,:), t2, rz)
    end if

    ! Compute B(1:m,1:NRHS) = Q**T * B(1:M,1:NRHS)
    qtb = mult_qr(.true., .true., a, tau, b)

    ! Solve the triangular system T11 * B(1:rnk,1:nrhs) = B(1:rnk,1:nrhs)
    if (rnk < n) then
        x(1:rnk,:) = solve_triangular_system(.true., .true., .false., .true., &
            one, rz(1:rnk,1:rnk), qtb(1:rnk,:))
        x(rnk+1:n,:) = zero

        ! Compute X(1:n,1:nrhs) = Y**T * X(1:n,1:nrhs)
        x = mult_rz(.true., .true., n - rnk, rz, t2, x)
    else
        x = solve_triangular_system(.true., .true., .false., .true., one, &
            a(1:rnk,1:rnk), qtb(1:rnk,:))
    end if

    ! Apply the pivoting: X(1:N,1:NRHS) = P * X(1:N,1:NRHS)
    do j = 1, nrhs
        do i = 1, n
            w(jpvt(i)) = x(i,j)
        end do
        x(1:n,j) = w(1:n)
    end do
end function

! ------------------------------------------------------------------------------
pure function solve_qr_pivot_vec(a, tau, jpvt, b) result(x)
    !! Solves a system of M QR-factored equations of N unknowns.
    real(real64), intent(in), dimension(:,:) :: a
        !! The M-by-N QR factored matrix as returned by qr_factor.
    real(real64), intent(in), dimension(:) :: tau
        !! A MIN(M, N)-element array containing the scalar factors of the 
        !! elementary reflectors as returned by qr_factor.
    integer(int32), intent(in), dimension(:) :: jpvt
        !! An N-element array, as output by qr_factor, used to track the 
        !! column pivots.
    real(real64), intent(in), dimension(:) :: b
        !! The M-element right-hand-side vector.
    real(real64), allocatable, dimension(:) :: x
        !! The N-element solution vector.

    ! Parameters
    integer(int32), parameter :: imin = 2
    integer(int32), parameter :: imax = 1
    real(real64), parameter :: zero = 0.0d0
    real(real64), parameter :: one = 1.0d0

    ! Local Variables
    integer(int32) :: i, m, n, mn, lwork, ismin, ismax, rnk, maxmn
    real(real64) :: rcond, smax, smin, smaxpr, sminpr, s1, c1, s2, c2
    real(real64), allocatable, dimension(:) :: w, t2, qtb
    real(real64), allocatable, dimension(:,:) :: rz

    ! Initialization
    m = size(a, 1)
    n = size(a, 2)
    mn = min(m, n)
    maxmn = max(m, n)
    ismin = 1
    ismax = mn + 1
    rcond = epsilon(rcond)
    lwork = 2 * maxmn + 1
    allocate(w(lwork), source = one)
    allocate(x(n))

    ! Input Check
    if (size(tau) /= mn) then
        error stop 2
    end if
    if (size(jpvt) /= n) then
        error stop 3
    end if
    if (size(b) /= m) then
        error stop 4
    end if

    ! Determine the rank of R11 using an incremental condition estimation
    smax = abs(a(1,1))
    smin = smax
    if (abs(a(1,1)) == zero) then
        rnk = 0
        x = zero
        return
    else
        rnk = 1
    end if
    do
        if (rnk < mn) then
            i = rnk + 1
            call DLAIC1(imin, rnk, w(ismin:ismin+rnk-1), smin, &
                a(1:rnk-1,i), a(i,i), sminpr, s1, c1)
            call DLAIC1(imax, rnk, w(ismax:ismax+rnk-1), smax, &
                a(1:rnk-1,i), a(i,i), smaxpr, s2, c2)
            if (smaxpr * rcond <= sminpr) then
                do i = 1, rnk
                    w(ismin+i-1) = s1 * w(ismin+i-1)
                    w(ismax+i-1) = s2 * w(ismax+i-1)
                end do
                w(ismin+rnk) = c1
                w(ismax+rnk) = c2
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
    ! [R11 R12] = [T11 0] * Y
    if (rnk < n) then
        call rz_factor(a(1:rnk,:), t2, rz)
    end if

    ! Compute B(1:m,1:NRHS) = Q**T * B(1:M,1:NRHS)
    qtb = mult_qr(.true., a, tau, b)

    ! Solve the triangular system T11 * B(1:rnk) = B(1:rnk)
    if (rnk < n) then
        x(1:rnk) = solve_triangular_system(.true., .false., .true., &
            rz(1:rnk,1:rnk), qtb(1:rnk))
        x(rnk+1:n) = zero

        ! Compute X(1:n) = Y**T * X(1:n)
        x = mult_rz(.true., n - rnk, rz, t2, x)
    else
        x = solve_triangular_system(.true., .false., .true., a(1:rnk,1:rnk), &
            qtb(1:rnk))
    end if

    ! Apply the pivoting: X(1:N) = P * X(1:N)
    do i = 1, n
        w(jpvt(i)) = x(i)
    end do
    x(1:n) = w(1:n)
end function

! ------------------------------------------------------------------------------
pure function solve_qr_pivot_vec_cmplx(a, tau, jpvt, b) result(x)
    !! Solves a system of M QR-factored equations of N unknowns.
    complex(real64), intent(in), dimension(:,:) :: a
        !! The M-by-N QR factored matrix as returned by qr_factor.
    complex(real64), intent(in), dimension(:) :: tau
        !! A MIN(M, N)-element array containing the scalar factors of the 
        !! elementary reflectors as returned by qr_factor.
    integer(int32), intent(in), dimension(:) :: jpvt
        !! An N-element array, as output by qr_factor, used to track the 
        !! column pivots.
    complex(real64), intent(in), dimension(:) :: b
        !! The M-element right-hand-side vector.
    complex(real64), allocatable, dimension(:) :: x
        !! The N-element solution vector.

    ! Parameters
    integer(int32), parameter :: imin = 2
    integer(int32), parameter :: imax = 1
    complex(real64), parameter :: zero = (0.0d0, 0.0d0)
    complex(real64), parameter :: one = (1.0d0, 0.0d0)

    ! Local Variables
    integer(int32) :: i, m, n, mn, lwork, ismin, ismax, rnk, maxmn
    real(real64) :: rcond, smax, smin, smaxpr, sminpr
    complex(real64) :: s1, c1, s2, c2
    complex(real64), allocatable, dimension(:) :: w, t2, qtb
    complex(real64), allocatable, dimension(:,:) :: rz

    ! Initialization
    m = size(a, 1)
    n = size(a, 2)
    mn = min(m, n)
    maxmn = max(m, n)
    ismin = 1
    ismax = mn + 1
    rcond = epsilon(rcond)
    lwork = 2 * maxmn + 1
    allocate(w(lwork), source = one)
    allocate(x(n))

    ! Input Check
    if (size(tau) /= mn) then
        error stop 2
    end if
    if (size(jpvt) /= n) then
        error stop 3
    end if
    if (size(b, 1) /= m) then
        error stop 4
    end if

    ! Determine the rank of R11 using an incremental condition estimation
    smax = abs(a(1,1))
    smin = smax
    if (abs(a(1,1)) == zero) then
        rnk = 0
        x = zero
        return
    else
        rnk = 1
    end if
    do
        if (rnk < mn) then
            i = rnk + 1
            call ZLAIC1(imin, rnk, w(ismin:ismin+rnk-1), smin, &
                a(1:rnk-1,i), a(i,i), sminpr, s1, c1)
            call ZLAIC1(imax, rnk, w(ismax:ismax+rnk-1), smax, &
                a(1:rnk-1,i), a(i,i), smaxpr, s2, c2)
            if (smaxpr * rcond <= sminpr) then
                do i = 1, rnk
                    w(ismin+i-1) = s1 * w(ismin+i-1)
                    w(ismax+i-1) = s2 * w(ismax+i-1)
                end do
                w(ismin+rnk) = c1
                w(ismax+rnk) = c2
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
    ! [R11 R12] = [T11 0] * Y
    if (rnk < n) then
        call rz_factor(a(1:rnk,:), t2, rz)
    end if

    ! Compute B(1:m,1:NRHS) = Q**T * B(1:M,1:NRHS)
    qtb = mult_qr(.true., a, tau, b)

    ! Solve the triangular system T11 * B(1:rnk) = B(1:rnk)
    if (rnk < n) then
        x(1:rnk) = solve_triangular_system(.true., .false., .true., &
            rz(1:rnk,1:rnk), qtb(1:rnk))
        x(rnk+1:n) = zero

        ! Compute X(1:n) = Y**T * X(1:n)
        x = mult_rz(.true., n - rnk, rz, t2, x)
    else
        x = solve_triangular_system(.true., .false., .true., a(1:rnk,1:rnk), &
            qtb(1:rnk))
    end if

    ! Apply the pivoting: X(1:N) = P * X(1:N)
    do i = 1, n
        w(jpvt(i)) = x(i)
    end do
    x(1:n) = w(1:n)
end function

! ------------------------------------------------------------------------------
end module