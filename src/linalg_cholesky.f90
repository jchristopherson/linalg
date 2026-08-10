module linalg_cholesky
    use iso_fortran_env, only : int32, real64
    use linalg_errors
    use lapack
    use qrupdate
    implicit none
    private
    public :: cholesky_factor
    public :: cholesky_rank1_update
    public :: cholesky_rank1_downdate
    public :: solve_cholesky

    interface cholesky_factor
        module procedure :: cholesky_factor_dbl
        module procedure :: cholesky_factor_cmplx
    end interface

    interface cholesky_rank1_update
        module procedure :: cholesky_rank1_update_dbl
        module procedure :: cholesky_rank1_update_cmplx
    end interface

    interface cholesky_rank1_downdate
        module procedure :: cholesky_rank1_downdate_dbl
        module procedure :: cholesky_rank1_downdate_cmplx
    end interface

    interface solve_cholesky
        module procedure :: solve_cholesky_mtx
        module procedure :: solve_cholesky_mtx_cmplx
        module procedure :: solve_cholesky_vec
        module procedure :: solve_cholesky_vec_cmplx
    end interface
contains
! ------------------------------------------------------------------------------
pure function cholesky_factor_dbl(a, upper) result(rst)
    !! Computes the Cholesky factorization of a symmetric, positive definite
    !! matrix.
    real(real64), intent(in), dimension(:,:) :: a
        !! The N-by-N matrix to factor.
    logical, intent(in), optional :: upper
        !! An optional input that, if specified, provides control over whether
        !! the factorization is computed as \(A = U^T U\) (set to true), or
        !! as \(A = L L^T\) (set to false).  The default is true such that
        !! \(A = U^T U\).
    real(real64), allocatable, dimension(:,:) :: rst
        !! The factored matrix.

    ! Parameters
    real(real64), parameter :: zero = 0.0d0

    ! Local Variables
    character :: uplo
    integer(int32) :: i, n, flag

    ! Initialization
    n = size(a, 1)
    if (present(upper)) then
        if (upper) then
            uplo = 'U'
        else
            uplo = 'L'
        end if
    else
        uplo = 'U'
    end if

    ! Input Check
    if (size(a, 2) /= n) then
        error stop 1
    end if

    ! Process
    allocate(rst(n, n), source = a)
    call DPOTRF(uplo, n, rst, n, flag)
    if (flag > 0) then
        error stop LA_MATRIX_FORMAT_ERROR
    end if

    ! Zero out the non-used upper or lower diagonal
    if (uplo == 'U') then
        ! Zero out the lower
        do i = 1, n - 1
            rst(i+1:n,i) = zero
        end do
    else
        ! Zero out the upper
        do i = 2, n
            rst(1:i-1,i) = zero
        end do
    end if
end function

! ------------------------------------------------------------------------------
pure function cholesky_factor_cmplx(a, upper) result(rst)
    !! Computes the Cholesky factorization of a symmetric, positive definite
    !! matrix.
    complex(real64), intent(in), dimension(:,:) :: a
        !! The N-by-N matrix to factor.
    logical, intent(in), optional :: upper
        !! An optional input that, if specified, provides control over whether
        !! the factorization is computed as \(A = U^H U\) (set to true), or
        !! as \(A = L L^H\) (set to false).  The default is true such that
        !! \(A = U^H U\).
    complex(real64), allocatable, dimension(:,:) :: rst
        !! The factored matrix.

    ! Parameters
    complex(real64), parameter :: zero = (0.0d0, 0.0d0)

    ! Local Variables
    character :: uplo
    integer(int32) :: i, n, flag

    ! Initialization
    n = size(a, 1)
    if (present(upper)) then
        if (upper) then
            uplo = 'U'
        else
            uplo = 'L'
        end if
    else
        uplo = 'U'
    end if

    ! Input Check
    if (size(a, 2) /= n) then
        error stop 1
    end if

    ! Process
    allocate(rst(n, n), source = a)
    call ZPOTRF(uplo, n, rst, n, flag)
    if (flag > 0) then
        ! ERROR: Matrix is not positive definite
        error stop LA_MATRIX_FORMAT_ERROR
    end if

    ! Zero out the non-used upper or lower diagonal
    if (uplo == 'U') then
        ! Zero out the lower
        do i = 1, n - 1
            rst(i+1:n,i) = zero
        end do
    else
        ! Zero out the upper
        do i = 2, n
            rst(1:i-1,i) = zero
        end do
    end if
end function

! ------------------------------------------------------------------------------
pure subroutine cholesky_rank1_update_dbl(r, u)
    !! Computes the rank 1 update to a Cholesky factored matrix \(A = R^T R\) 
    !! such that \(A_1 = A + \vec{u} \vec{u}^T\).
    real(real64), intent(inout), dimension(:,:) :: r
        !! On input, the N-by-N upper triangular matrix \(R\).  On output, the 
        !! updated matrix \(R_1\).
    real(real64), intent(inout), dimension(:) :: u
        !! On input, the N-element vector \(\vec{u}\).  On output, the rotation
        !! sines used to transform \(R\) to \(R_1\).

    ! Local Variables
    integer(int32) :: n, lwork
    real(real64), allocatable, dimension(:) :: w

    ! Initialization
    n = size(r, 1)
    lwork = n

    ! Input Check
    if (size(r, 2) /= n) then
        error stop 1
    else if (size(u) /= n) then
        error stop 2
    end if

    ! Local Memory Allocation
    allocate(w(lwork))

    ! Process
    call DCH1UP(n, r, n, u, w)
end subroutine

! ------------------------------------------------------------------------------
pure subroutine cholesky_rank1_update_cmplx(r, u)
    !! Computes the rank 1 update to a Cholesky factored matrix \(A = R^H R\) 
    !! such that \(A_1 = A + \vec{u} \vec{u}^H\).
    complex(real64), intent(inout), dimension(:,:) :: r
        !! On input, the N-by-N upper triangular matrix \(R\).  On output, the 
        !! updated matrix \(R_1\).
    complex(real64), intent(inout), dimension(:) :: u
        !! On input, the N-element vector \(\vec{u}\).  On output, the rotation
        !! sines used to transform \(R\) to \(R_1\).

    ! Local Variables
    integer(int32) :: n, lwork
    real(real64), allocatable, dimension(:) :: w

    ! Initialization
    n = size(r, 1)
    lwork = n

    ! Input Check
    if (size(r, 2) /= n) then
        error stop 1
    else if (size(u) /= n) then
        error stop 2
    end if

    ! Local Memory Allocation
    allocate(w(lwork))

    ! Process
    call ZCH1UP(n, r, n, u, w)
end subroutine

! ------------------------------------------------------------------------------
pure subroutine cholesky_rank1_downdate_dbl(r, u)
    !! Computes the rank 1 downdate to a Cholesky factored matrix \(A = R^T R\) 
    !! such that \(A_1 = A - \vec{u} \vec{u}^T\).  This operation only works if
    !! the new matrix \(A_1\) is positive definite.
    real(real64), intent(inout), dimension(:,:) :: r
        !! On input, the N-by-N upper triangular matrix \(R\).  On output, the 
        !! updated matrix \(R_1\).
    real(real64), intent(inout), dimension(:) :: u
        !! On input, the N-element vector \(\vec{u}\).  On output, the rotation
        !! sines used to transform \(R\) to \(R_1\).

    ! Local Variables
    integer(int32) :: n, lwork, flag
    real(real64), allocatable, dimension(:) :: w

    ! Initialization
    n = size(r, 1)
    lwork = n

    ! Input Check
    if (size(r, 2) /= n) then
        error stop 1
    else if (size(u) /= n) then
        error stop 2
    end if

    ! Local Memory Allocation
    allocate(w(lwork))

    ! Process
    call DCH1DN(n, r, n, u, w, flag)
    if (flag == 1) then
        ! ERROR: The matrix is not positive definite
        error stop LA_MATRIX_FORMAT_ERROR
    else if (flag == 2) then
        ! ERROR: The matrix is singular
        error stop LA_SINGULAR_MATRIX_ERROR
    end if
end subroutine

! ------------------------------------------------------------------------------
pure subroutine cholesky_rank1_downdate_cmplx(r, u)
    !! Computes the rank 1 downdate to a Cholesky factored matrix \(A = R^H R\) 
    !! such that \(A_1 = A - \vec{u} \vec{u}^H\).  This operation only works if
    !! the new matrix \(A_1\) is positive definite.
    complex(real64), intent(inout), dimension(:,:) :: r
        !! On input, the N-by-N upper triangular matrix \(R\).  On output, the 
        !! updated matrix \(R_1\).
    complex(real64), intent(inout), dimension(:) :: u
        !! On input, the N-element vector \(\vec{u}\).  On output, the rotation
        !! sines used to transform \(R\) to \(R_1\).

    ! Local Variables
    integer(int32) :: n, lwork, flag
    real(real64), allocatable, dimension(:) :: w

    ! Initialization
    n = size(r, 1)
    lwork = n

    ! Input Check
    if (size(r, 2) /= n) then
        error stop 1
    else if (size(u) /= n) then
        error stop 2
    end if

    ! Local Memory Allocation
    allocate(w(lwork))

    ! Process
    call ZCH1DN(n, r, n, u, w, flag)
    if (flag == 1) then
        ! ERROR: The matrix is not positive definite
        error stop LA_MATRIX_FORMAT_ERROR
    else if (flag == 2) then
        ! ERROR: The matrix is singular
        error stop LA_SINGULAR_MATRIX_ERROR
    end if
end subroutine

! ------------------------------------------------------------------------------
pure function solve_cholesky_mtx(upper, a, b) result(x)
    !! Solves the system of Cholesky factored equations \(A X = R^T R X = B\) or
    !! \(A X = L L^T X = B\).
    logical, intent(in) :: upper
        !! Set to true if \(A\) is factored such that \(A = R^T R\); else, set
        !! to false if \(A\) is factored such that \(A = L L^T\).
    real(real64), intent(in), dimension(:,:) :: a
        !! The N-by-N Cholesky factored matrix as returned by cholesky_factor.
    real(real64), intent(in), dimension(:,:) :: b
        !! The N-by-NRHS matrix \(B\).
    real(real64), allocatable, dimension(:,:) :: x 
        !! The resulting N-by-NRHS matrix \(X\).

    ! Local Variables
    character :: uplo
    integer(int32) :: n, nrhs, flag

    ! Initialization
    n = size(a, 1)
    nrhs = size(b, 2)
    if (upper) then
        uplo = 'U'
    else
        uplo = 'L'
    end if

    ! Input Check
    if (size(a, 2) /= n) then
        error stop 1
    else if (size(b, 1) /= n) then
        error stop 2
    end if

    ! Process
    allocate(x(n,nrhs), source = b)
    call DPOTRS(uplo, n, nrhs, a, n, x, n, flag)
end function

! ------------------------------------------------------------------------------
pure function solve_cholesky_mtx_cmplx(upper, a, b) result(x)
    !! Solves the system of Cholesky factored equations \(A X = R^H R X = B\) or
    !! \(A X = L L^H X = B\).
    logical, intent(in) :: upper
        !! Set to true if \(A\) is factored such that \(A = R^H R\); else, set
        !! to false if \(A\) is factored such that \(A = L L^H\).
    complex(real64), intent(in), dimension(:,:) :: a
        !! The N-by-N Cholesky factored matrix as returned by cholesky_factor.
    complex(real64), intent(in), dimension(:,:) :: b
        !! The N-by-NRHS matrix \(B\).
    complex(real64), allocatable, dimension(:,:) :: x
        !! The resulting N-by-NRHS matrix \(X\).

    ! Local Variables
    character :: uplo
    integer(int32) :: n, nrhs, flag

    ! Initialization
    n = size(a, 1)
    nrhs = size(b, 2)
    if (upper) then
        uplo = 'U'
    else
        uplo = 'L'
    end if

    ! Input Check
    if (size(a, 2) /= n) then
        error stop 1
    else if (size(b, 1) /= n) then
        error stop 2
    end if

    ! Process
    allocate(x(n,nrhs), source = b)
    call ZPOTRS(uplo, n, nrhs, a, n, x, n, flag)
end function

! ------------------------------------------------------------------------------
pure function solve_cholesky_vec(upper, a, b) result(x)
    !! Solves the system of Cholesky factored equations 
    !! \(A \vec{x} = R^T R \vec{x} = \vec{b}\) or
    !! \(A \vec{x} = L L^T \vec{x} = \vec{b}\).
    logical, intent(in) :: upper
        !! Set to true if \(A\) is factored such that \(A = R^T R\); else, set
        !! to false if \(A\) is factored such that \(A = L L^T\).
    real(real64), intent(in), dimension(:,:) :: a
        !! The N-by-N Cholesky factored matrix as returned by cholesky_factor.
    real(real64), intent(in), dimension(:) :: b
        !! The N-element vector \(\vec{b}\).
    real(real64), allocatable, dimension(:) :: x
        !! The resulting N-element vector \(\vec{x}\).

    ! Local Variables
    character :: uplo
    integer(int32) :: n, flag

    ! Initialization
    n = size(a, 1)
    if (upper) then
        uplo = 'U'
    else
        uplo = 'L'
    end if

    ! Input Check
    if (size(a, 2) /= n) then
        error stop 1
    else if (size(b) /= n) then
        error stop 2
    end if

    ! Process
    allocate(x(n), source = b)
    call DPOTRS(uplo, n, 1, a, n, x, n, flag)
end function

! ------------------------------------------------------------------------------
pure function solve_cholesky_vec_cmplx(upper, a, b) result(x)
    !! Solves the system of Cholesky factored equations 
    !! \(A \vec{x} = R^H R \vec{x} = \vec{b}\) or
    !! \(A \vec{x} = L L^H \vec{x} = \vec{b}\).
    logical, intent(in) :: upper
        !! Set to true if \(A\) is factored such that \(A = R^H R\); else, set
        !! to false if \(A\) is factored such that \(A = L L^H\).
    complex(real64), intent(in), dimension(:,:) :: a
        !! The N-by-N Cholesky factored matrix as returned by cholesky_factor.
    complex(real64), intent(in), dimension(:) :: b
        !! The N-element vector \(\vec{b}\).
    complex(real64), allocatable, dimension(:) :: x
        !! The resulting N-element vector \(\vec{x}\).

    ! Local Variables
    character :: uplo
    integer(int32) :: n, flag

    ! Initialization
    n = size(a, 1)
    if (upper) then
        uplo = 'U'
    else
        uplo = 'L'
    end if

    ! Input Check
    if (size(a, 2) /= n) then
        error stop 1
    else if (size(b) /= n) then
        error stop 2
    end if

    ! Process
    allocate(x(n), source = b)
    call ZPOTRS(uplo, n, 1, a, n, x, n, flag)
end function

! ------------------------------------------------------------------------------
end module