module linalg_tri
    use iso_fortran_env
    use blas
    use linalg_errors
    implicit none
    private
    public :: solve_triangular_system

    interface solve_triangular_system
        module procedure :: solve_tri_mtx
        module procedure :: solve_tri_mtx_cmplx
        module procedure :: solve_tri_vec
        module procedure :: solve_tri_vec_cmplx
    end interface
contains
! ------------------------------------------------------------------------------
pure function solve_tri_mtx(lside, upper, trans, nounit, alpha, a, b) result(x)
    !! Solves a triangular system of equations of the form 
    !! \(op(A) X = \alpha B\) or \(X op(A) = \alpha B\) where \(A\) is a 
    !! triangular matrix (either upper or lower) for the unknown \(X\).
    logical, intent(in) :: lside
        !! Set to true to solve \(op(A) X = \alpha B\); else, set to false to
        !! solve \(X op(A) = \alpha B\).
    logical, intent(in) :: upper
        !! Set to true if \(A\) is upper triangular; else, set to false if
        !! \(A\) is lower triangular.
    logical, intent(in) :: trans
        !! Set to true if \(op(A) = A^T\); else, set to false if \(op(A) = A\).
    logical, intent(in) :: nounit
        !! Set to true if \(A\) is unit-triangular (ones on the diagonal); else,
        !! false if \(A\) is not unit-triangular.
    real(real64), intent(in) :: alpha
        !! The scalar multiplier \(\alpha\).
    real(real64), intent(in), dimension(:,:) :: a
        !! If lside is true, the M-by-M triangular matrix \(A\); else, \(A\) is
        !! N-by-N if lside is false.
    real(real64), intent(in), dimension(:,:) :: b
        !! The M-by-N matrix \(B\).
    real(real64), allocatable, dimension(:,:) :: x
        !! The M-by-N matrix \(X\).

    ! Local Variables
    character :: side, uplo, transa, diag
    integer(int32) :: m, n, nrowa

    ! Initialization
    m = size(b, 1)
    n = size(b, 2)
    if (lside) then
        nrowa = m
        side = 'L'
    else
        nrowa = n
        side = 'R'
    end if
    if (upper) then
        uplo = 'U'
    else
        uplo = 'L'
    end if
    if (trans) then
        transa = 'T'
    else
        transa = 'N'
    end if
    if (nounit) then
        diag = 'N'
    else
        diag = 'U'
    end if

    ! Input Check - matrix A must be square
    if (size(a, 1) /= nrowa .or. size(a, 2) /= nrowa) then
        error stop 6
    end if

    ! Call DTRSM
    allocate(x(m, n), source = b)
    call DTRSM(side, uplo, transa, diag, m, n, alpha, a, nrowa, x, m)
end function

! ------------------------------------------------------------------------------
pure function solve_tri_mtx_cmplx(lside, upper, trans, nounit, alpha, a, b) result(x)
    !! Solves a triangular system of equations of the form 
    !! \(op(A) X = \alpha B\) or \(X op(A) = \alpha B\) where \(A\) is a 
    !! triangular matrix (either upper or lower) for the unknown \(X\).
    logical, intent(in) :: lside
        !! Set to true to solve \(op(A) X = \alpha B\); else, set to false to
        !! solve \(X op(A) = \alpha B\).
    logical, intent(in) :: upper
        !! Set to true if \(A\) is upper triangular; else, set to false if
        !! \(A\) is lower triangular.
    logical, intent(in) :: trans
        !! Set to true if \(op(A) = A^H\); else, set to false if \(op(A) = A\).
    logical, intent(in) :: nounit
        !! Set to true if \(A\) is unit-triangular (ones on the diagonal); else,
        !! false if \(A\) is not unit-triangular.
    complex(real64), intent(in) :: alpha
        !! The scalar multiplier \(\alpha\).
    complex(real64), intent(in), dimension(:,:) :: a
        !! If lside is true, the M-by-M triangular matrix \(A\); else, \(A\) is
        !! N-by-N if lside is false.
    complex(real64), intent(in), dimension(:,:) :: b
        !! On input, The M-by-N matrix \(B\).
    complex(real64), allocatable, dimension(:,:) :: x
        !! The M-by-N matrix \(X\).

    ! Local Variables
    character :: side, uplo, transa, diag
    integer(int32) :: m, n, nrowa

    ! Initialization
    m = size(b, 1)
    n = size(b, 2)
    if (lside) then
        nrowa = m
        side = 'L'
    else
        nrowa = n
        side = 'R'
    end if
    if (upper) then
        uplo = 'U'
    else
        uplo = 'L'
    end if
    if (trans) then
        transa = 'C'
    else
        transa = 'N'
    end if
    if (nounit) then
        diag = 'N'
    else
        diag = 'U'
    end if

    ! Input Check - matrix A must be square
    if (size(a, 1) /= nrowa .or. size(a, 2) /= nrowa) then
        error stop 6
    end if

    ! Call ZTRSM
    allocate(x(m, n), source = b)
    call ZTRSM(side, uplo, transa, diag, m, n, alpha, a, nrowa, x, m)
end function

! ------------------------------------------------------------------------------
pure function solve_tri_vec(upper, trans, nounit, a, b) result(x)
    !! Solves the triangular system \(op(A) \vec{x} = \vec{b}\) where \(A\)
    !! is a triangular matrix.
    logical, intent(in) :: upper
        !! Set to true if \(A\) is upper triangular; else, set to false if \(A\)
        !! is lower triangular.
    logical, intent(in) :: trans
        !! Set to true if \(op(A) = A^T\); else, set to false if \(op(A) = A\).
    logical, intent(in) :: nounit
        !! Set to true if \(A\) is unit-triangular (ones on the diagonal); else,
        !! false if \(A\) is not unit-triangular.
    real(real64), intent(in), dimension(:,:) :: a
        !! The N-by-N triangular matrix \(A\).
    real(real64), intent(in), dimension(:) :: b
        !! The N-element vector \(\vec{b}\).
    real(real64), allocatable, dimension(:) :: x  
        !! The N-element vector \(\vec{x}\).

    ! Parameters
    real(real64), parameter :: zero = 0.0d0

    ! Local Variables
    character :: uplo, t, diag
    integer(int32) :: n

    ! Initialization
    n = size(a, 1)
    if (upper) then
        uplo = 'U'
    else
        uplo = 'L'
    end if
    if (trans) then
        t = 'T'
    else
        t = 'N'
    end if
    if (nounit) then
        diag = 'N'
    else
        diag = 'U'
    end if

    ! Input Check
    if (size(a, 2) /= n) then
        error stop 4
    end if
    if (size(b) /= n) then
        error stop 5
    end if

    ! Call DTRSV
    allocate(x(n), source = b)
    call DTRSV(uplo, t, diag, n, a, n, x, 1)
end function

! ------------------------------------------------------------------------------
pure function solve_tri_vec_cmplx(upper, trans, nounit, a, b) result(x)
    !! Solves the triangular system \(op(A) \vec{x} = \vec{b}\) where \(A\)
    !! is a triangular matrix.
    logical, intent(in) :: upper
        !! Set to true if \(A\) is upper triangular; else, set to false if \(A\)
        !! is lower triangular.
    logical, intent(in) :: trans
        !! Set to true if \(op(A) = A^T\); else, set to false if \(op(A) = A\).
    logical, intent(in) :: nounit
        !! Set to true if \(A\) is unit-triangular (ones on the diagonal); else,
        !! false if \(A\) is not unit-triangular.
    complex(real64), intent(in), dimension(:,:) :: a
        !! The N-by-N triangular matrix \(A\).
    complex(real64), intent(in), dimension(:) :: b
        !! The N-element vector \(\vec{b}\).
    complex(real64), allocatable, dimension(:) :: x  
        !! The N-element vector \(\vec{x}\).

    ! Parameters
    real(real64), parameter :: zero = 0.0d0

    ! Local Variables
    character :: uplo, t, diag
    integer(int32) :: n

    ! Initialization
    n = size(a, 1)
    if (upper) then
        uplo = 'U'
    else
        uplo = 'L'
    end if
    if (trans) then
        t = 'C'
    else
        t = 'N'
    end if
    if (nounit) then
        diag = 'N'
    else
        diag = 'U'
    end if

    ! Input Check
    if (size(a, 2) /= n) then
        error stop 4
    else if (size(b) /= n) then
        error stop 5
    end if

    ! Call ZTRSV
    allocate(x(n), source = b)
    call ZTRSV(uplo, t, diag, n, a, n, x, 1)
end function

! ------------------------------------------------------------------------------
end module