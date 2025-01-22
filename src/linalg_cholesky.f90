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
subroutine cholesky_factor_dbl(a, upper, err)
    !! Computes the Cholesky factorization of a symmetric, positive definite
    !! matrix.
    real(real64), intent(inout), dimension(:,:) :: a
        !! On input, the N-by-N matrix to factor.  On output, the factored 
        !! matrix is returned in either the upper or lower triangular portion 
        !! of the matrix, dependent upon the value of upper.
    logical, intent(in), optional :: upper
        !! An optional input that, if specified, provides control over whether
        !! the factorization is computed as \(A = U^T U\) (set to true), or
        !! as \(A = L L^T\) (set to false).  The default is true such that
        !! \(A = U^T U\).
    class(errors), intent(inout), optional, target :: err
        !! The error object to be updated.

    ! Parameters
    real(real64), parameter :: zero = 0.0d0

    ! Local Variables
    character :: uplo
    integer(int32) :: i, n, flag
    class(errors), pointer :: errmgr
    type(errors), target :: deferr

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
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if

    ! Input Check
    if (size(a, 2) /= n) then
        call report_square_matrix_error("cholesky_factor", errmgr, "a", &
            n, size(a, 1), size(a, 2))
        return
    end if

    ! Process
    call DPOTRF(uplo, n, a, n, flag)
    if (flag > 0) then
        call errmgr%report_error("cholesky_factor", &
            "The matrix is not positive-definite.", LA_MATRIX_FORMAT_ERROR)
        return
    end if

    ! Zero out the non-used upper or lower diagonal
    if (uplo == 'U') then
        ! Zero out the lower
        do i = 1, n - 1
            a(i+1:n,i) = zero
        end do
    else
        ! Zero out the upper
        do i = 2, n
            a(1:i-1,i) = zero
        end do
    end if
end subroutine

! ------------------------------------------------------------------------------
subroutine cholesky_factor_cmplx(a, upper, err)
    !! Computes the Cholesky factorization of a symmetric, positive definite
    !! matrix.
    complex(real64), intent(inout), dimension(:,:) :: a
        !! On input, the N-by-N matrix to factor.  On output, the factored 
        !! matrix is returned in either the upper or lower triangular portion 
        !! of the matrix, dependent upon the value of upper.
    logical, intent(in), optional :: upper
        !! An optional input that, if specified, provides control over whether
        !! the factorization is computed as \(A = U^H U\) (set to true), or
        !! as \(A = L L^H\) (set to false).  The default is true such that
        !! \(A = U^H U\).
    class(errors), intent(inout), optional, target :: err
        !! The error object to be updated.

    ! Parameters
    complex(real64), parameter :: zero = (0.0d0, 0.0d0)

    ! Local Variables
    character :: uplo
    integer(int32) :: i, n, flag
    class(errors), pointer :: errmgr
    type(errors), target :: deferr

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
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if

    ! Input Check
    if (size(a, 2) /= n) then
        call report_square_matrix_error("cholesky_factor_cmplx", errmgr, "a", &
            n, size(a, 1), size(a, 2))
        return
    end if

    ! Process
    call ZPOTRF(uplo, n, a, n, flag)
    if (flag > 0) then
        ! ERROR: Matrix is not positive definite
        call errmgr%report_error("cholesky_factor_cmplx", &
            "The matrix is not positive-definite.", LA_MATRIX_FORMAT_ERROR)
    end if

    ! Zero out the non-used upper or lower diagonal
    if (uplo == 'U') then
        ! Zero out the lower
        do i = 1, n - 1
            a(i+1:n,i) = zero
        end do
    else
        ! Zero out the upper
        do i = 2, n
            a(1:i-1,i) = zero
        end do
    end if
end subroutine

! ------------------------------------------------------------------------------
subroutine cholesky_rank1_update_dbl(r, u, work, err)
    !! Computes the rank 1 update to a Cholesky factored matrix \(A = R^T R\) 
    !! such that \(A_1 = A + \vec{u} \vec{u}^T\).
    real(real64), intent(inout), dimension(:,:) :: r
        !! On input, the N-by-N upper triangular matrix \(R\).  On output, the 
        !! updated matrix \(R_1\).
    real(real64), intent(inout), dimension(:) :: u
        !! On input, the N-element vector \(\vec{u}\).  On output, the rotation
        !! sines used to transform \(R\) to \(R_1\).
    real(real64), intent(out), target, optional, dimension(:) :: work
        !! An optional argument that if supplied prevents local memory 
        !! allocation.  If provided, the array must have at least N elements.  
        !! Additionally, this workspace array is used to contain the rotation 
        !! cosines used to transform \(R\) to \(R_1\).
    class(errors), intent(inout), optional, target :: err
        !! An error object to report any errors that occur.

    ! Local Variables
    integer(int32) :: n, lwork, istat
    real(real64), pointer, dimension(:) :: wptr
    real(real64), allocatable, target, dimension(:) :: wrk
    class(errors), pointer :: errmgr
    type(errors), target :: deferr

    ! Initialization
    n = size(r, 1)
    lwork = n
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if

    ! Input Check
    if (size(r, 2) /= n) then
        call report_square_matrix_error("cholesky_rank1_update_dbl", errmgr, &
            "r", n, size(r, 1), size(r, 2))
        return
    else if (size(u) /= n) then
        call report_array_size_error("cholesky_rank1_update_dbl", errmgr, &
            "u", n, size(u))
        return
    end if

    ! Local Memory Allocation
    if (present(work)) then
        if (size(work) < lwork) then
            call report_array_size_error("cholesky_rank1_update_dbl", errmgr, &
                "work", lwork, size(work))
            return
        end if
        wptr => work(1:lwork)
    else
        allocate(wrk(lwork), stat = istat)
        if (istat /= 0) then
            call report_memory_error("cholesky_rank1_update_dbl", errmgr, istat)
            return
        end if
        wptr => wrk
    end if

    ! Process
    call DCH1UP(n, r, n, u, wptr)
end subroutine

! ------------------------------------------------------------------------------
subroutine cholesky_rank1_update_cmplx(r, u, work, err)
    !! Computes the rank 1 update to a Cholesky factored matrix \(A = R^H R\) 
    !! such that \(A_1 = A + \vec{u} \vec{u}^H\).
    complex(real64), intent(inout), dimension(:,:) :: r
        !! On input, the N-by-N upper triangular matrix \(R\).  On output, the 
        !! updated matrix \(R_1\).
    complex(real64), intent(inout), dimension(:) :: u
        !! On input, the N-element vector \(\vec{u}\).  On output, the rotation
        !! sines used to transform \(R\) to \(R_1\).
    real(real64), intent(out), target, optional, dimension(:) :: work
        !! An optional argument that if supplied prevents local memory 
        !! allocation.  If provided, the array must have at least N elements.  
        !! Additionally, this workspace array is used to contain the rotation 
        !! cosines used to transform \(R\) to \(R_1\).
    class(errors), intent(inout), optional, target :: err
        !! An error object to report any errors that occur.

    ! Local Variables
    integer(int32) :: n, lwork, istat
    real(real64), pointer, dimension(:) :: wptr
    real(real64), allocatable, target, dimension(:) :: wrk
    class(errors), pointer :: errmgr
    type(errors), target :: deferr

    ! Initialization
    n = size(r, 1)
    lwork = n
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if

    ! Input Check
    if (size(r, 2) /= n) then
        call report_square_matrix_error("cholesky_rank1_update_cmplx", errmgr, &
            "r", n, size(r, 1), size(r, 2))
        return
    else if (size(u) /= n) then
        call report_array_size_error("cholesky_rank1_update_cmplx", errmgr, &
            "u", n, size(u))
        return
    end if

    ! Local Memory Allocation
    if (present(work)) then
        if (size(work) < lwork) then
            call report_array_size_error("cholesky_rank1_update_cmplx", &
                errmgr, "work", lwork, size(work))
            return
        end if
        wptr => work(1:lwork)
    else
        allocate(wrk(lwork), stat = istat)
        if (istat /= 0) then
            call report_memory_error("cholesky_rank1_update_cmplx", errmgr, &
                istat)
            return
        end if
        wptr => wrk
    end if

    ! Process
    call ZCH1UP(n, r, n, u, wptr)
end subroutine

! ------------------------------------------------------------------------------
subroutine cholesky_rank1_downdate_dbl(r, u, work, err)
    !! Computes the rank 1 downdate to a Cholesky factored matrix \(A = R^T R\) 
    !! such that \(A_1 = A - \vec{u} \vec{u}^T\).  This operation only works if
    !! the new matrix \(A_1\) is positive definite.
    real(real64), intent(inout), dimension(:,:) :: r
        !! On input, the N-by-N upper triangular matrix \(R\).  On output, the 
        !! updated matrix \(R_1\).
    real(real64), intent(inout), dimension(:) :: u
        !! On input, the N-element vector \(\vec{u}\).  On output, the rotation
        !! sines used to transform \(R\) to \(R_1\).
    real(real64), intent(out), target, optional, dimension(:) :: work
        !! An optional argument that if supplied prevents local memory 
        !! allocation.  If provided, the array must have at least N elements.  
        !! Additionally, this workspace array is used to contain the rotation 
        !! cosines used to transform \(R\) to \(R_1\).
    class(errors), intent(inout), optional, target :: err
        !! An error object to report any errors that occur.

    ! Local Variables
    integer(int32) :: n, lwork, istat, flag
    real(real64), pointer, dimension(:) :: wptr
    real(real64), allocatable, target, dimension(:) :: wrk
    class(errors), pointer :: errmgr
    type(errors), target :: deferr

    ! Initialization
    n = size(r, 1)
    lwork = n
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if

    ! Input Check
    flag = 0
    if (size(r, 2) /= n) then
        call report_square_matrix_error("cholesky_rank1_downdate_dbl", &
            errmgr, "r", n, size(r, 1), size(r, 2))
        return
    else if (size(u) /= n) then
        call report_array_size_error("cholesky_rank1_downdate_dbl", errmgr, &
            "u", n, size(u))
    end if

    ! Local Memory Allocation
    if (present(work)) then
        if (size(work) < lwork) then
            call report_array_size_error("cholesky_rank1_downdate_dbl", &
                errmgr, "work", lwork, size(work))
            return
        end if
        wptr => work(1:lwork)
    else
        allocate(wrk(lwork), stat = istat)
        if (istat /= 0) then
            call report_memory_error("cholesky_rank1_downdate_dbl", errmgr, &
                istat)
            return
        end if
        wptr => wrk
    end if

    ! Process
    call DCH1DN(n, r, n, u, wptr, flag)
    if (flag == 1) then
        ! ERROR: The matrix is not positive definite
        call errmgr%report_error("cholesky_rank1_downdate", &
            "The downdated matrix is not positive definite.", &
            LA_MATRIX_FORMAT_ERROR)
    else if (flag == 2) then
        ! ERROR: The matrix is singular
        call errmgr%report_error("cholesky_rank1_downdate", &
            "The input matrix is singular.", LA_SINGULAR_MATRIX_ERROR)
    end if
end subroutine

! ------------------------------------------------------------------------------
subroutine cholesky_rank1_downdate_cmplx(r, u, work, err)
    !! Computes the rank 1 downdate to a Cholesky factored matrix \(A = R^H R\) 
    !! such that \(A_1 = A - \vec{u} \vec{u}^H\).  This operation only works if
    !! the new matrix \(A_1\) is positive definite.
    complex(real64), intent(inout), dimension(:,:) :: r
        !! On input, the N-by-N upper triangular matrix \(R\).  On output, the 
        !! updated matrix \(R_1\).
    complex(real64), intent(inout), dimension(:) :: u
        !! On input, the N-element vector \(\vec{u}\).  On output, the rotation
        !! sines used to transform \(R\) to \(R_1\).
    real(real64), intent(out), target, optional, dimension(:) :: work
        !! An optional argument that if supplied prevents local memory 
        !! allocation.  If provided, the array must have at least N elements.  
        !! Additionally, this workspace array is used to contain the rotation 
        !! cosines used to transform \(R\) to \(R_1\).
    class(errors), intent(inout), optional, target :: err
        !! An error object to report any errors that occur.

    ! Local Variables
    integer(int32) :: n, lwork, istat, flag
    real(real64), pointer, dimension(:) :: wptr
    real(real64), allocatable, target, dimension(:) :: wrk
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    character(len = :), allocatable :: errmsg

    ! Initialization
    n = size(r, 1)
    lwork = n
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if

    ! Input Check
    flag = 0
    if (size(r, 2) /= n) then
        call report_square_matrix_error("cholesky_rank1_downdate_cmplx", &
            errmgr, "r", n, size(r, 1), size(r, 2))
        return
    else if (size(u) /= n) then
        call report_array_size_error("cholesky_rank1_downdate_cmplx", errmgr, &
            "u", n, size(u))
        return
    end if

    ! Local Memory Allocation
    if (present(work)) then
        if (size(work) < lwork) then
            call report_array_size_error("cholesky_rank1_downdate_cmplx", &
                errmgr, "work", lwork, size(work))
            return
        end if
        wptr => work(1:lwork)
    else
        allocate(wrk(lwork), stat = istat)
        if (istat /= 0) then
            call report_memory_error("cholesky_rank1_downdate_cmplx", errmgr, &
                istat)
            return
        end if
        wptr => wrk
    end if

    ! Process
    call ZCH1DN(n, r, n, u, wptr, flag)
    if (flag == 1) then
        ! ERROR: The matrix is not positive definite
        call errmgr%report_error("cholesky_rank1_downdate_cmplx", &
            "The downdated matrix is not positive definite.", &
            LA_MATRIX_FORMAT_ERROR)
    else if (flag == 2) then
        ! ERROR: The matrix is singular
        call errmgr%report_error("cholesky_rank1_downdate_cmplx", &
            "The input matrix is singular.", LA_SINGULAR_MATRIX_ERROR)
    end if
end subroutine

! ------------------------------------------------------------------------------
subroutine solve_cholesky_mtx(upper, a, b, err)
    !! Solves the system of Cholesky factored equations \(A X = R^T R X = B\) or
    !! \(A X = L L^T X = B\).
    logical, intent(in) :: upper
        !! Set to true if \(A\) is factored such that \(A = R^T R\); else, set
        !! to false if \(A\) is factored such that \(A = L L^T\).
    real(real64), intent(in), dimension(:,:) :: a
        !! The N-by-N Cholesky factored matrix as returned by cholesky_factor.
    real(real64), intent(inout), dimension(:,:) :: b
        !! On input, the N-by-NRHS matrix \(B\).  On output, the resulting
        !! N-by-NRHS matrix \(X\).
    class(errors), intent(inout), optional, target :: err
        !! An error object to report any errors that occur.

    ! Local Variables
    character :: uplo
    integer(int32) :: n, nrhs, flag
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    character(len = :), allocatable :: errmsg

    ! Initialization
    n = size(a, 1)
    nrhs = size(b, 2)
    if (upper) then
        uplo = 'U'
    else
        uplo = 'L'
    end if
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if

    ! Input Check
    if (size(a, 2) /= n) then
        call report_square_matrix_error("solve_cholesky_mtx", errmgr, &
            "a", n, size(a, 1), size(a, 2))
        return
    else if (size(b, 1) /= n) then
        call report_matrix_size_error("solve_cholesky_mtx", errmgr, "b", &
            n, nrhs, size(b, 1), size(b, 2))
        return
    end if

    ! Process
    call DPOTRS(uplo, n, nrhs, a, n, b, n, flag)
end subroutine

! ------------------------------------------------------------------------------
subroutine solve_cholesky_mtx_cmplx(upper, a, b, err)
    !! Solves the system of Cholesky factored equations \(A X = R^H R X = B\) or
    !! \(A X = L L^H X = B\).
    logical, intent(in) :: upper
        !! Set to true if \(A\) is factored such that \(A = R^H R\); else, set
        !! to false if \(A\) is factored such that \(A = L L^H\).
    complex(real64), intent(in), dimension(:,:) :: a
        !! The N-by-N Cholesky factored matrix as returned by cholesky_factor.
    complex(real64), intent(inout), dimension(:,:) :: b
        !! On input, the N-by-NRHS matrix \(B\).  On output, the resulting
        !! N-by-NRHS matrix \(X\).
    class(errors), intent(inout), optional, target :: err
        !! An error object to report any errors that occur.

    ! Local Variables
    character :: uplo
    integer(int32) :: n, nrhs, flag
    class(errors), pointer :: errmgr
    type(errors), target :: deferr

    ! Initialization
    n = size(a, 1)
    nrhs = size(b, 2)
    if (upper) then
        uplo = 'U'
    else
        uplo = 'L'
    end if
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if

    ! Input Check
    if (size(a, 2) /= n) then
        call report_square_matrix_error("solve_cholesky_mtx_cmplx", errmgr, &
            "a", n, size(a, 1), size(a, 2))
        return
    else if (size(b, 1) /= n) then
        call report_matrix_size_error("solve_cholesky_mtx_cmplx", errmgr, "b", &
            n, nrhs, size(b, 1), size(b, 2))
        return
    end if

    ! Process
    call ZPOTRS(uplo, n, nrhs, a, n, b, n, flag)
end subroutine

! ------------------------------------------------------------------------------
subroutine solve_cholesky_vec(upper, a, b, err)
    !! Solves the system of Cholesky factored equations 
    !! \(A \vec{x} = R^T R \vec{x} = \vec{b}\) or
    !! \(A \vec{x} = L L^T \vec{x} = \vec{b}\).
    logical, intent(in) :: upper
        !! Set to true if \(A\) is factored such that \(A = R^T R\); else, set
        !! to false if \(A\) is factored such that \(A = L L^T\).
    real(real64), intent(in), dimension(:,:) :: a
        !! The N-by-N Cholesky factored matrix as returned by cholesky_factor.
    real(real64), intent(inout), dimension(:) :: b
        !! On input, the N-element vector \(\vec{b}\).  On output, the resulting
        !! N-element vector \(\vec{x}\).
    class(errors), intent(inout), optional, target :: err
        !! An error object to report any errors that occur.

    ! Local Variables
    character :: uplo
    integer(int32) :: n, flag
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    character(len = :), allocatable :: errmsg

    ! Initialization
    n = size(a, 1)
    if (upper) then
        uplo = 'U'
    else
        uplo = 'L'
    end if
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if

    ! Input Check
    if (size(a, 2) /= n) then
        call report_square_matrix_error("solve_cholesky_vec", errmgr, &
            "a", n, size(a, 1), size(a, 2))
        return
    else if (size(b) /= n) then
        call report_array_size_error("solve_cholesky_vec", errmgr, "b", &
            n, size(b))
        return
    end if

    ! Process
    call DPOTRS(uplo, n, 1, a, n, b, n, flag)
end subroutine

! ------------------------------------------------------------------------------
module subroutine solve_cholesky_vec_cmplx(upper, a, b, err)
    !! Solves the system of Cholesky factored equations 
    !! \(A \vec{x} = R^H R \vec{x} = \vec{b}\) or
    !! \(A \vec{x} = L L^H \vec{x} = \vec{b}\).
    logical, intent(in) :: upper
        !! Set to true if \(A\) is factored such that \(A = R^H R\); else, set
        !! to false if \(A\) is factored such that \(A = L L^H\).
    complex(real64), intent(in), dimension(:,:) :: a
        !! The N-by-N Cholesky factored matrix as returned by cholesky_factor.
    complex(real64), intent(inout), dimension(:) :: b
        !! On input, the N-element vector \(\vec{b}\).  On output, the resulting
        !! N-element vector \(\vec{x}\).
    class(errors), intent(inout), optional, target :: err
        !! An error object to report any errors that occur.

    ! Local Variables
    character :: uplo
    integer(int32) :: n, flag
    class(errors), pointer :: errmgr
    type(errors), target :: deferr

    ! Initialization
    n = size(a, 1)
    if (upper) then
        uplo = 'U'
    else
        uplo = 'L'
    end if
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if

    ! Input Check
    if (size(a, 2) /= n) then
        call report_square_matrix_error("solve_cholesky_vec_cmplx", errmgr, &
            "a", n, size(a, 1), size(a, 2))
        return
    else if (size(b) /= n) then
        call report_array_size_error("solve_cholesky_vec_cmplx", errmgr, "b", &
            n, size(b))
        return
    end if

    ! Process
    call ZPOTRS(uplo, n, 1, a, n, b, n, flag)
end subroutine

! ------------------------------------------------------------------------------
end module