module linear_algebra
    !! This module contains pure functions for basic linear algebra operations.
    use iso_fortran_env
    use ieee_arithmetic
    use lapack
    use blas
    implicit none
    private
    public :: swap_arrays
    public :: identity
    public :: lu_factors
    public :: qr_factors
    public :: svd_factors
    public :: lu_factor
    public :: qr_factor
    public :: cholesky_factor
    public :: svd
    public :: solve_triangular_system
    public :: solve_linear_system
    public :: solve_least_squares
    public :: inverse
    public :: pinverse

    type :: lu_factors
        !! A container for the results of a LU factorization.
        real(real64), allocatable, dimension(:,:) :: L
            !! The lower-triangular factorization.
        real(real64), allocatable, dimension(:,:) :: U
            !! The upper-triangular factorization.
        real(real64), allocatable, dimension(:,:) :: P
            !! The pivot tracking matrix.
    end type

    type :: qr_factors
        !! A container for the results of a QR factorization of an M-by-N 
        !! matrix.
        real(real64), allocatable, dimension(:,:) :: Q
            !! The M-by-M orthogonal matrix, \(Q\).
        real(real64), allocatable, dimension(:,:) :: R
            !! The M-by-N upper trapezoidal matrix, \(R\).
        real(real64), allocatable, dimension(:,:) :: P
            !! The N-by-N pivot tracking matrix.
    end type

    type :: svd_factors
        !! A container for the results of a singular value decomposition of
        !! an M-by-N matrix.
        real(real64), allocatable, dimension(:,:) :: U
            !! The M-by-M orthogonal matrix \( U \).
        real(real64), allocatable, dimension(:,:) :: S
            !! The M-by-N diagonal matrix \(S\) containing the singular values
            !! on the diagonal.
        real(real64), allocatable, dimension(:,:) :: Vt
            !! The transpose of the N-by-N right singular vector matrix \(V\).
    end type

    interface solve_triangular_system
        module procedure :: solve_triangular_system_mtx
        module procedure :: solve_triangular_system_vec
    end interface

    interface solve_linear_system
        module procedure :: solve_linear_system_mtx
        module procedure :: solve_linear_system_vec
    end interface

    interface solve_least_squares
        module procedure :: solve_least_squares_mtx
        module procedure :: solve_least_squares_vec
    end interface
contains
! ******************************************************************************
! COMMON OPERATIONS
! ------------------------------------------------------------------------------
pure subroutine swap_arrays(x, y)
    !! Swaps the contents of two arrays.
    real(real64), intent(inout), dimension(:) :: x
        !! The first array.
    real(real64), intent(inout), dimension(:) :: y
        !! The second array.

    ! Local Variables
    integer(int32) :: i, m, mp1, nx, ny, n
    real(real64) :: temp

    ! Initialization
    nx = size(x)
    ny = size(y)
    n = min(nx, ny)
    m = mod(n, 3)
    mp1 = m + 1

    ! Process
    if (m /= 0) then
        do i = 1, m
            temp = x(i)
            x(i) = y(i)
            y(i) = temp
        end do
        if (n < 3) return
    end if
    do i = mp1, n, 3
        temp = x(i)
        x(i) = y(i)
        y(i) = temp

        temp = x(i+1)
        x(i+1) = y(i+1)
        y(i+1) = temp

        temp = x(i+2)
        x(i+2) = y(i+2)
        y(i+2) = temp
    end do
end subroutine

! ------------------------------------------------------------------------------
pure function identity(n) result(rst)
    !! Constructs an N-by-N identity matrix.
    integer(int32), intent(in) :: n
        !! The size of the matrix.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The resulting matrix.

    ! Local Variables
    integer(int32) :: i

    ! Process
    if (n < 1) return
    allocate(rst(n, n), source = 0.0d0)
    do i = 1, n
        rst(i,i) = 1.0d0
    end do
end function

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ******************************************************************************
! FACTORIZATIONS
! ------------------------------------------------------------------------------
pure function lu_factor(a) result(rst)
    use linalg_lu, only : form_lu
    !! Computes the LU factorization of a square matrix such that
    !! \( P A = L U \).
    real(real64), intent(in), dimension(:,:) :: a
        !! The N-by-N matrix to factor.
    type(lu_factors) :: rst
        !! The factored form of the matrix.

    ! Local Variables
    integer(int32) :: i, ip, n, info
    real(real64), allocatable, dimension(:,:) :: u
    integer(int32), allocatable, dimension(:) :: p

    ! Initialization
    n = size(a, 1)
    if (size(a, 2) /= n) return
    allocate(rst%L(n, n), source = a)
    allocate(rst%U(n, n), source = 0.0d0)
    allocate(p(n))

    ! Process
    call DGETRF(n, n, rst%L, n, p, info)
    rst%P = identity(n)
    do i = 1, n
        ! Build the pivot matrix
        ip = p(i)
        if (i /= ip) call swap_arrays(rst%P(i,:), rst%P(ip,:))

        ! Build L & U
        rst%U(1:i,i) = rst%L(1:i,i)
        if (i > 1) rst%L(1:i-1,i) = 0.0d0
        rst%L(i,i) = 1.0d0
    end do
end function

! ------------------------------------------------------------------------------
pure function qr_factor(a, pivot) result(rst)
    !! Computes the QR factorization of an M-by-N matrix such that either
    !! \(A = Q R \) (no pivoting), or \(A P = Q R\) (with pivoting).
    real(real64), intent(in), dimension(:,:) :: a
        !! The M-by-N matrix to factor.
    logical, intent(in), optional :: pivot
        !! An optional parameter used to specifiy if pivoting should be used
        !! (true); else, false if no pivoting is used.  The default is false
        !! such that no pivoting is performed.
    type(qr_factors) :: rst
        !! The factored form of the matrix.

    ! Local Variables
    logical :: pvt
    integer(int32) :: j, jp, m, n, mn, lwork, info
    integer(int32), allocatable, dimension(:) :: jpvt
    real(real64) :: temp1(1), temp2(1)
    real(real64), allocatable, dimension(:) :: tau, work

    ! Initialization
    m = size(a, 1)
    n = size(a, 2)
    mn = min(m, n)
    allocate(tau(mn))
    pvt = .false.
    if (present(pivot)) pvt = pivot
    allocate(rst%Q(m, m), source = 0.0d0)
    allocate(rst%R(m, n), source = a)
    rst%P = identity(n)

    ! Determine the workspace requirements
    if (pvt) then
        allocate(jpvt(n), source = 0)
        call DGEQP3(m, n, rst%R, m, jpvt, tau, temp1, -1, info)
    else
        call DGEQRF(m, n, rst%R, m, tau, temp1, -1, info)
    end if
    call DORGQR(m, m, mn, rst%R, m, tau, temp2, -1, info)
    lwork = int(max(temp1(1), temp2(1)), int32)
    allocate(work(lwork))

    ! Compute the factorization
    if (pvt) then
        call DGEQP3(m, n, rst%R, m, jpvt, tau, work, lwork, info)
    else
        call DGEQRF(m, n, rst%R, m, tau, work, lwork, info)
    end if

    ! Build the matrices
    do j = 1, mn
        rst%Q(j+1:m,j) = rst%R(j+1:m,j)
        rst%R(j+1:m,j) = 0.0d0
    end do
    call DORGQR(m, m, mn, rst%Q, m, tau, work, lwork, info)

    ! Construct the pivot matrix, if necessary
    if (pvt) then
        do j = 1, n
            jp = jpvt(j)
            rst%P(:,j) = 0.0d0
            rst%P(jp,j) = 1.0d0
        end do
    end if
end function

! ------------------------------------------------------------------------------
pure function cholesky_factor(a, upper) result(rst)
    !! Computes the Cholesky factorization of a positive-definite matrix.
    real(real64), intent(in), dimension(:,:) :: a
        !! The matrix to factor.
    logical, intent(in), optional :: upper
        !! An optional parameter to specifiy if the upper factorization
        !! \(A = R^{T} R \) should be computed (true); else, false for the lower
        !! factorization \( A = L L^{T} \).  The default is to compute the upper
        !! factorization.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The factored matrix, either \(R\) or \(L\).

    ! Local Variables
    integer(int32) :: i, n, info
    character :: uplo

    ! Initialization
    n = size(a, 1)
    if (size(a, 2) /= n) return
    uplo = 'U'
    if (present(upper)) then
        if (.not.upper) uplo = 'L'
    end if
    allocate(rst(n, n), source = a)

    ! Factor
    call DPOTRF(uplo, n, rst, n, info)

    ! Zero out the non-used upper or lower diagonal portions of the matrix
    if (uplo == 'U') then
        do i = 1, n - 1
            rst(i+1:n,i) = 0.0d0
        end do
    else
        do i = 2, n
            rst(1:i-1,i) = 0.0d0
        end do
    end if
end function

! ------------------------------------------------------------------------------
pure function svd(a) result(rst)
    !! Computes the singular value decomposition of an M-by-N matrix such that
    !! \( A = U S V^{T} \).
    real(real64), intent(in), dimension(:,:) :: a
        !! The M-by-N matrix to factor.
    type(svd_factors) :: rst
        !! The factored form of the matrix.

    ! Local Variables
    integer(int32) :: i, m, n, mn, lwork, info
    real(real64) :: temp(1)
    real(real64), allocatable, dimension(:) :: sigma, work
    real(real64), allocatable, dimension(:,:) :: ac

    ! Initialization
    m = size(a, 1)
    n = size(a, 2)
    mn = min(m, n)
    ac = a
    allocate(rst%U(m, m), rst%S(m, n), rst%Vt(n, n), sigma(mn), source = 0.0d0)

    ! Determine the workspace
    call DGESVD('A', 'A', m, n, ac, m, sigma, rst%U, m, rst%Vt, n, temp, -1, info)
    lwork = int(temp(1), int32)
    allocate(work(lwork))

    ! Compute the factorization
    call DGESVD('A', 'A', m, n, ac, m, sigma, rst%U, m, rst%Vt, n, work, lwork, info)

    ! Populate the diagonal singular-value matrix
    do i = 1, mn
        rst%S(i,i) = sigma(i)
    end do
end function

! ******************************************************************************
! SOLVERS
! ------------------------------------------------------------------------------
pure function solve_triangular_system_mtx(a, b, upper) result(rst)
    !! Solves a triangular system of the form \(A X = B\) where \(A\) is a
    !! triangular matrix, either upper or lower, for equation \(X\).
    real(real64), intent(in), dimension(:,:) :: a
        !! The N-by-N triangular \(A\) matrix.
    real(real64), intent(in), dimension(:,:) :: b
        !! The N-by-NRHS \(B\) matrix.
    logical, intent(in), optional :: upper
        !! An optional argument specifying if the \(A\) matrix is upper 
        !! triangular (true), or lower triangular (false).  The default 
        !! assumption is that \(A\) is an upper triangular matrix.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The N-by-NRHS solution matrix, \(X\).

    ! Local Variables
    character :: uplo, side, transa, diag
    integer(int32) :: n, nrhs

    ! Initialization
    n = size(a, 1)
    nrhs = size(b, 2)
    side = 'L'
    transa = 'N'
    diag = 'N'
    uplo = 'U'
    if (present(upper)) then
        if (upper) then
            uplo = 'U'
        else
            uplo = 'L'
        end if
    end if
    if (size(a, 2) /= n .or. size(b, 1) /= n) return
    allocate(rst(n, nrhs), source = b)

    ! Process
    call DTRSM(side, uplo, transa, diag, n, nrhs, 1.0d0, a, n, rst, n)
end function

! ------------------------------------------------------------------------------
pure function solve_triangular_system_vec(a, b, upper) result(rst)
    !! Solves a triangular system of the form \(A X = B\) where \(A\) is a
    !! triangular matrix, either upper or lower, for equation \(X\).
    real(real64), intent(in), dimension(:,:) :: a
        !! The N-by-N triangular \(A\) matrix.
    real(real64), intent(in), dimension(:) :: b
        !! The N-element \(B\) array.
    logical, intent(in), optional :: upper
        !! An optional argument specifying if the \(A\) matrix is upper 
        !! triangular (true), or lower triangular (false).  The default 
        !! assumption is that \(A\) is an upper triangular matrix.
    real(real64), allocatable, dimension(:) :: rst
        !! The N-element solution array, \(X\).

    ! Local Variables
    character :: uplo, side, transa, diag
    integer(int32) :: n, nrhs

    ! Initialization
    n = size(a, 1)
    nrhs = 1
    side = 'L'
    transa = 'N'
    diag = 'N'
    uplo = 'U'
    if (present(upper)) then
        if (upper) then
            uplo = 'U'
        else
            uplo = 'L'
        end if
    end if
    if (size(a, 2) /= n .or. size(b, 1) /= n) return
    allocate(rst(n), source = b)

    ! Process
    call DTRSM(side, uplo, transa, diag, n, nrhs, 1.0d0, a, n, rst, n)
end function

! ------------------------------------------------------------------------------
pure function solve_linear_system_mtx(a, b) result(rst)
    !! Solves the M-by-N linear system \(A X = B\) for \(X\).
    real(real64), intent(in), dimension(:,:) :: a
        !! The M-by-N matrix \(A\).
    real(real64), intent(in), dimension(:,:) :: b
        !! The M-by-NRHS matrix \(B\).
    real(real64), allocatable, dimension(:,:) :: rst
        !! The resulting N-by-NRHS matrix \(X\).

    ! Local Variables
    logical :: usedgels
    integer(int32) :: m, n, maxmn, nrhs, lwork, info
    integer(int32), allocatable, dimension(:) :: ipiv
    real(real64) :: temp(1)
    real(real64), allocatable, dimension(:) :: work
    real(real64), allocatable, dimension(:,:) :: ac, xc

    ! Initialization
    m = size(a, 1)
    n = size(a, 2)
    nrhs = size(b, 2)
    if (size(b, 1) /= m) return
    usedgels = .false.
    if (m /= n) usedgels = .true.
    allocate(ac(m, n), source = a)
    maxmn = max(m, n)

    ! Process
    if (usedgels) then
        allocate(xc(maxmn,nrhs))
        if (m >= n) then
            xc = b
        else
            xc(1:m,:) = b
        end if
        call DGELS('N', m, n, nrhs, ac, m, xc, maxmn, temp, -1, info)
        lwork = int(temp(1), int32)
        allocate(work(lwork))
        call DGELS('N', m, n, nrhs, ac, m, xc, maxmn, work, lwork, info)
        if (m >= n) then
            allocate(rst(n, nrhs), source = xc(1:n,:))
        else
            allocate(rst(n, nrhs), source = xc)
        end if
    else
        ! Use LU factorization to solve
        allocate(rst(n, nrhs), source = b)
        allocate(ipiv(n))
        call DGESV(n, nrhs, ac, n, ipiv, rst, n, info)
    end if
end function

! ------------------------------------------------------------------------------
pure function solve_linear_system_vec(a, b) result(rst)
    !! Solves the M-by-N linear system \(A X = B\) for \(X\).
    real(real64), intent(in), dimension(:,:) :: a
        !! The M-by-N matrix \(A\).
    real(real64), intent(in), dimension(:) :: b
        !! The M-element array \(B\).
    real(real64), allocatable, dimension(:) :: rst
        !! The resulting N-element array \(X\).

    ! Local Variables
    logical :: usedgels
    integer(int32) :: m, n, maxmn, nrhs, lwork, info
    integer(int32), allocatable, dimension(:) :: ipiv
    real(real64) :: temp(1)
    real(real64), allocatable, dimension(:) :: work, xc
    real(real64), allocatable, dimension(:,:) :: ac

    ! Initialization
    m = size(a, 1)
    n = size(a, 2)
    nrhs = 1
    if (size(b) /= m) return
    usedgels = .false.
    if (m /= n) usedgels = .true.
    allocate(ac(m, n), source = a)
    maxmn = max(m, n)

    ! Process
    if (usedgels) then
        allocate(xc(maxmn))
        if (m >= n) then
            xc = b
        else
            xc(1:m) = b
        end if
        call DGELS('N', m, n, nrhs, ac, m, xc, maxmn, temp, -1, info)
        lwork = int(temp(1), int32)
        allocate(work(lwork))
        call DGELS('N', m, n, nrhs, ac, m, xc, maxmn, work, lwork, info)
        if (m >= n) then
            allocate(rst(n), source = xc(1:n))
        else
            allocate(rst(n), source = xc)
        end if
    else
        ! Use LU factorization to solve
        allocate(rst(n), source = b)
        allocate(ipiv(n))
        call DGESV(n, nrhs, ac, n, ipiv, rst, n, info)
    end if
end function

! ------------------------------------------------------------------------------
pure function solve_least_squares_mtx(a, b) result(rst)
    !! Solves the least squares problem by minimizing \(|| A X - B ||\) using
    !! a complete orthogonal factorization of \(A\).
    real(real64), intent(in), dimension(:,:) :: a
        !! The M-by-N matrix \(A\).
    real(real64), intent(in), dimension(:,:) :: b
        !! The M-by-NRHS matrix \(B\).
    real(real64), allocatable, dimension(:,:) :: rst
        !! The resulting N-by-NRHS matrix \(X\).

    ! Local Variables
    integer(int32) :: m, n, maxmn, nrhs, lwork, info, rnk
    integer(int32), allocatable, dimension(:) :: jpvt
    real(real64) :: rcond, temp(1)
    real(real64), allocatable, dimension(:) :: work
    real(real64), allocatable, dimension(:,:) :: ac, x

    ! Initialization
    m = size(a, 1)
    n = size(a, 2)
    nrhs = size(b, 2)
    maxmn = max(m, n)
    if (size(b, 1) /= m) return
    allocate(jpvt(n), source = 0)
    allocate(ac(m, n), source = a)
    allocate(x(maxmn, nrhs))
    if (m >= n) then
        x = b
    else
        x(1:m,:) = b
    end if
    call DGELSY(m, n, nrhs, ac, m, x, maxmn, jpvt, rcond, rnk, temp, -1, info)
    lwork = int(temp(1), int32)
    allocate(work(lwork))

    ! Process
    call DGELSY(m, n, nrhs, ac, m, x, maxmn, jpvt, rcond, rnk, work, lwork, info)
    allocate(rst(n, nrhs), source = x(1:n,:))
end function

! ------------------------------------------------------------------------------
pure function solve_least_squares_vec(a, b) result(rst)
    !! Solves the least squares problem by minimizing \(|| A X - B ||\) using
    !! a complete orthogonal factorization of \(A\).
    real(real64), intent(in), dimension(:,:) :: a
        !! The M-by-N matrix \(A\).
    real(real64), intent(in), dimension(:) :: b
        !! The M-element array \(B\).
    real(real64), allocatable, dimension(:) :: rst
        !! The resulting N-element array \(X\).

    ! Local Variables
    integer(int32) :: m, n, maxmn, nrhs, lwork, info, rnk
    integer(int32), allocatable, dimension(:) :: jpvt
    real(real64) :: rcond, temp(1)
    real(real64), allocatable, dimension(:) :: work, x
    real(real64), allocatable, dimension(:,:) :: ac

    ! Initialization
    m = size(a, 1)
    n = size(a, 2)
    nrhs = 1
    maxmn = max(m, n)
    if (size(b) /= m) return
    allocate(jpvt(n), source = 0)
    allocate(ac(m, n), source = a)
    allocate(x(maxmn))
    if (m >= n) then
        x = b
    else
        x(1:m) = b
    end if
    call DGELSY(m, n, nrhs, ac, m, x, maxmn, jpvt, rcond, rnk, temp, -1, info)
    lwork = int(temp(1), int32)
    allocate(work(lwork))

    ! Process
    call DGELSY(m, n, nrhs, ac, m, x, maxmn, jpvt, rcond, rnk, work, lwork, info)
    allocate(rst(n), source = x(1:n))
end function

! ******************************************************************************
! INVERSE OPERATIONS
! ------------------------------------------------------------------------------
pure function inverse(a) result(rst)
    !! Computes the inverse of a square matrix.
    real(real64), intent(in), dimension(:,:) :: a
        !! The N-by-N matrix to invert.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The N-by-N inverted matrix.

    ! Local Variables
    integer(int32) :: n, lwork, info
    integer(int32), allocatable, dimension(:) :: ipvt
    real(real64) :: temp(1), nan
    real(real64), allocatable, dimension(:) :: work

    ! Process
    n = size(a, 1)
    if (size(a, 2) /= n) return
    allocate(rst(n, n), source = a)
    allocate(ipvt(n))
    call DGETRI(n, rst, n, ipvt, temp, -1, info)
    lwork = int(temp(1), int32)
    allocate(work(lwork))

    ! Compute the LU factorization of A
    call DGETRF(n, n, rst, n, ipvt, info)
    if (info > 0) then
        nan = ieee_value(nan, IEEE_QUIET_NAN)
        rst = nan
        return
    end if

    ! Compute the inverse of the LU factored matrix
    call DGETRI(n, rst, n, ipvt, work, lwork, info)
end function

! ------------------------------------------------------------------------------
pure function pinverse(a) result(rst)
    !! Computes the Moore-Penrose pseudoinverse of an M-by-N matrix.
    real(real64), intent(in), dimension(:,:) :: a
        !! The M-by-N matrix to invert.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The N-by-M inverted matrix.

    ! Local Variables
    integer(int32) :: i, m, n, mn, info, lwork
    real(real64) :: temp(1), t, ss
    real(real64), allocatable, dimension(:) :: s, work
    real(real64), allocatable, dimension(:,:) :: ac, u, vt

    ! Initialization
    m = size(a, 1)
    n = size(a, 2)
    mn = min(m, n)
    allocate(ac(m, n), source = a)
    allocate(s(mn), u(m,mn), vt(mn,n))
    call DGESVD('S', 'S', m, n, ac, m, s, u, m, vt, mn, temp, -1, info)
    lwork = int(temp(1), int32)
    allocate(work(lwork))
    allocate(rst(n, m))

    ! Compute the SVD of A
    call DGESVD('S', 'S', m, n, ac, m, s, u, m, vt, mn, work, lwork, info)

    ! Compute u = U * inv(S)
    t = max(m, n) * epsilon(t) * s(1)
    do i = 1, mn
        if (s(i) > t) then
            ss = 1.0d0 / s(i)
        else
            ss = s(i)
        end if
        u(:,i) = ss * u(:,i)
    end do

    ! Compute pinv(A) = (V**T)**T * u**T
    call DGEMM('T', 'T', n, m, mn, 1.0d0, vt, mn, u, m, 0.0d0, rst, n)
end function

! ******************************************************************************
! EIGEN ANALYSIS
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end module