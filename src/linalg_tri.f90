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
subroutine solve_tri_mtx(lside, upper, trans, nounit, alpha, a, b, err)
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
    real(real64), intent(inout), dimension(:,:) :: b
        !! On input, the M-by-N matrix \(B\).  On output, the M-by-N solution 
        !! matrix \(X\).
    class(errors), intent(inout), optional, target :: err
        !! The error object to be updated.

    ! Parameters
    character :: side, uplo, transa, diag

    ! Local Variables
    integer(int32) :: m, n, nrowa
    class(errors), pointer :: errmgr
    type(errors), target :: deferr

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
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if

    ! Input Check - matrix A must be square
    if (size(a, 1) /= nrowa .or. size(a, 2) /= nrowa) then
        call report_matrix_size_error("solve_tri_mtx", errmgr, "a", &
            nrowa, nrowa, size(a, 1), size(a, 2))
        return
    end if

    ! Call DTRSM
    call DTRSM(side, uplo, transa, diag, m, n, alpha, a, nrowa, b, m)
end subroutine

! ------------------------------------------------------------------------------
subroutine solve_tri_mtx_cmplx(lside, upper, trans, nounit, alpha, a, b, err)
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
    complex(real64), intent(inout), dimension(:,:) :: b
        !! On input, the M-by-N matrix \(B\).  On output, the M-by-N solution 
        !! matrix \(X\).
    class(errors), intent(inout), optional, target :: err
        !! The error object to be updated.

    ! Parameters
    character :: side, uplo, transa, diag

    ! Local Variables
    integer(int32) :: m, n, nrowa
    class(errors), pointer :: errmgr
    type(errors), target :: deferr

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
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if

    ! Input Check - matrix A must be square
    if (size(a, 1) /= nrowa .or. size(a, 2) /= nrowa) then
        call report_matrix_size_error("solve_tri_mtx_cmplx", errmgr, "a", &
            nrowa, nrowa, size(a, 1), size(a, 2))
        return
    end if

    ! Call ZTRSM
    call ZTRSM(side, uplo, transa, diag, m, n, alpha, a, nrowa, b, m)
end subroutine

! ------------------------------------------------------------------------------
subroutine solve_tri_vec(upper, trans, nounit, a, x, err)
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
    real(real64), intent(inout), dimension(:) :: x
        !! On input, the N-element vector \(\vec{b}\).  On output, the 
        !! N-element solution vector \(\vec{x}\).
    class(errors), intent(inout), optional, target :: err
        !! The error object to be updated.

    ! Parameters
    real(real64), parameter :: zero = 0.0d0

    ! Local Variables
    character :: uplo, t, diag
    integer(int32) :: n
    class(errors), pointer :: errmgr
    type(errors), target :: deferr

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
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if

    ! Input Check
    if (size(a, 2) /= n) then
        call report_square_matrix_error("solve_tri_vec", errmgr, "a", &
            n, size(a, 1), size(a, 2))
        return
    else if (size(x) /= n) then
        call report_inner_matrix_dimension_error("solve_tri_vec", errmgr, &
            "a", "x", n, size(x))
        return
    end if

    ! Call DTRSV
    call DTRSV(uplo, t, diag, n, a, n, x, 1)
end subroutine

! ------------------------------------------------------------------------------
subroutine solve_tri_vec_cmplx(upper, trans, nounit, a, x, err)
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
    complex(real64), intent(inout), dimension(:) :: x
        !! On input, the N-element vector \(\vec{b}\).  On output, the 
        !! N-element solution vector \(\vec{x}\).
    class(errors), intent(inout), optional, target :: err
        !! The error object to be updated.

    ! Parameters
    real(real64), parameter :: zero = 0.0d0

    ! Local Variables
    character :: uplo, t, diag
    integer(int32) :: n
    class(errors), pointer :: errmgr
    type(errors), target :: deferr

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
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if

    ! Input Check
    if (size(a, 2) /= n) then
        call report_square_matrix_error("solve_tri_vec_cmplx", errmgr, "a", &
            n, size(a, 1), size(a, 2))
        return
    else if (size(x) /= n) then
        call report_inner_matrix_dimension_error("solve_tri_vec_cmplx", &
            errmgr, "a", "x", n, size(x))
        return
    end if

    ! Call ZTRSV
    call ZTRSV(uplo, t, diag, n, a, n, x, 1)
end subroutine

! ------------------------------------------------------------------------------
end module