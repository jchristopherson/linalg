module linalg_lu
    use iso_fortran_env
    use linalg_errors
    use linalg_sparse
    use linalg_basic
    use ferror
    use lapack
    use sparskit
    implicit none
    private
    public :: lu_factor
    public :: form_lu
    public :: solve_lu

    interface lu_factor
        module procedure :: lu_factor_dbl
        module procedure :: lu_factor_cmplx
        module procedure :: csr_lu_factor
    end interface

    interface form_lu
        module procedure :: form_lu_all
        module procedure :: form_lu_all_cmplx
        module procedure :: form_lu_only
        module procedure :: form_lu_only_cmplx
    end interface

    interface solve_lu
        module procedure :: solve_lu_mtx
        module procedure :: solve_lu_mtx_cmplx
        module procedure :: solve_lu_vec
        module procedure :: solve_lu_vec_cmplx
        module procedure :: csr_lu_solve
    end interface
contains
! ------------------------------------------------------------------------------
subroutine lu_factor_dbl(a, ipvt, err)
    !! Computes the LU factorization of an M-by-N matrix.
    real(real64), intent(inout), dimension(:,:) :: a
        !! On input, the M-by-N matrix on which to operate.  On output, the 
        !! LU factored matrix in the form [L\\U] where the unit diagonal
        !! elements of L are not stored.
    integer(int32), intent(out), dimension(:) :: ipvt
        !! An MIN(M, N)-element array used to track row-pivot operations.  The 
        !! array stored pivot information such that row I is interchanged with 
        !! row IPVT(I).
    class(errors), intent(inout), optional, target :: err
        !! An error object to report any errors that occur.

    ! Local Variables
    integer(int32) :: m, n, mn, flag
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
    flag = 0
    if (size(ipvt) /= mn) then
        ! ERROR: IPVT not sized correctly
        call report_array_size_error("lu_factor_dbl", errmgr, "ipvt", mn, &
            size(ipvt))
        return
    end if

    ! Compute the LU factorization by calling the LAPACK routine DGETRF
    call DGETRF(m, n, a, m, ipvt, flag)

    ! If flag > 0, the matrix is singular.  Notice, flag should not be
    ! able to be < 0 as we've already verrified inputs prior to making the
    ! call to LAPACK
    if (flag > 0) then
        ! WARNING: Singular matrix
        call report_singular_matrix_warning("lu_factor_dbl", errmgr, flag)
    end if
end subroutine

! ------------------------------------------------------------------------------
subroutine lu_factor_cmplx(a, ipvt, err)
    !! Computes the LU factorization of an M-by-N matrix.
    complex(real64), intent(inout), dimension(:,:) :: a
        !! On input, the M-by-N matrix on which to operate.  On output, the 
        !! LU factored matrix in the form [L\\U] where the unit diagonal
        !! elements of L are not stored.
    integer(int32), intent(out), dimension(:) :: ipvt
        !! An MIN(M, N)-element array used to track row-pivot operations.  The 
        !! array stored pivot information such that row I is interchanged with 
        !! row IPVT(I).
    class(errors), intent(inout), optional, target :: err
        !! An error object to report any errors that occur.

    ! Local Variables
    integer(int32) :: m, n, mn, flag
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
    flag = 0
    if (size(ipvt) /= mn) then
        ! ERROR: IPVT not sized correctly
        call errmgr%report_error("lu_factor_cmplx", &
            "Incorrectly sized input array IPVT, argument 2.", &
            LA_ARRAY_SIZE_ERROR)
        return
    end if

    ! Compute the LU factorization by calling the LAPACK routine ZGETRF
    call ZGETRF(m, n, a, m, ipvt, flag)

    ! If flag > 0, the matrix is singular.  Notice, flag should not be
    ! able to be < 0 as we've already verrified inputs prior to making the
    ! call to LAPACK
    if (flag > 0) then
        ! WARNING: Singular matrix
        call report_singular_matrix_warning("lu_factor_cmplx", errmgr, flag)
    end if
end subroutine

! ------------------------------------------------------------------------------
subroutine csr_lu_factor(a, lu, ju, droptol, err)
    !! Factors a matrix using an LU decomposition.
    class(csr_matrix), intent(in) :: a
        !! The matrix to factor.
    type(msr_matrix), intent(out) :: lu
        !! The LU matrix.
    integer(int32), intent(out), dimension(:) :: ju
        !! The row tracking array.
    real(real64), intent(in), optional :: droptol
        !! The drop tolerance for the ILUT factorization.
    class(errors), intent(inout), optional, target :: err
        !! The error object to be updated.

    ! Local Variables
    integer(int32) :: i, m, n, nn, nnz, lfil, iwk, ierr, flag
    integer(int32), allocatable, dimension(:) :: jlu, jw
    real(real64), allocatable, dimension(:) :: alu, w
    real(real64) :: dt
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    if (present(droptol)) then
        dt = droptol
    else
        dt = sqrt(epsilon(dt))
    end if
    m = size(a, 1)
    n = size(a, 2)
    nnz = nonzero_count(a)

    ! Input Check
    if (size(ju) /= m) then
        call report_array_size_error("csr_lu_factor", errmgr, "ju", m, size(ju))
        return
    end if

    ! Parameter Determination
    lfil = 1
    do i = 1, m
        lfil = max(lfil, a%row_indices(i+1) - a%row_indices(i))
    end do
    iwk = max(lfil * m, nnz)  ! somewhat arbitrary - can be adjusted

    ! Local Memory Allocation
    allocate(alu(iwk), w(n+1), jlu(iwk), jw(2 * n), stat = flag)
    if (flag /= 0) go to 10

    ! Factorization
    do
        ! Factor the matrix
        call ilut(n, a%values, a%column_indices, a%row_indices, lfil, dt, &
            alu, jlu, ju, iwk, w, jw, ierr)

        ! Check the error flag
        if (ierr == 0) then
            ! Success
            exit
        else if (ierr > 0) then
            ! Zero pivot
        else if (ierr == -1) then
            ! The input matrix is not formatted correctly
            go to 20
        else if (ierr == -2 .or. ierr == -3) then
            ! ALU and JLU are too small - try something larger
            ! This is the main reason for the loop - to offload worrying about
            ! workspace size from the user
            iwk = min(iwk + m + n, m * n)
            deallocate(alu)
            deallocate(jlu)
            allocate(alu(iwk), jlu(iwk), stat = flag)
            if (flag /= 0) go to 10
        else if (ierr == -4) then
            ! Illegal value for LFIL - reset and try again
            lfil = n
        else if (ierr == -5) then
            ! Zero row encountered
            go to 30
        else
            ! We should never get here, but just in case
            go to 40
        end if
    end do

    ! Determine the actual number of non-zero elements
    nnz = jlu(m+1) - 1

    ! Copy the contents to the output arrays
    lu%m = m
    lu%n = n
    lu%nnz = nnz
    nn = m + 1 + nnz - min(m, n)
    allocate(lu%values(nn), source = alu(:nn), stat = flag)
    if (flag /= 0) go to 10
    allocate(lu%indices(nn), source = jlu(:nn), stat = flag)

    ! End
    return

    ! Memory Error
10  continue
    call report_memory_error("csr_lu_factor", errmgr, flag) 
    return

    ! Matrix Format Error
20  continue
    call errmgr%report_error("csr_lu_factor", &
        "The input matrix was incorrectly formatted.  A row with more " // &
        "than N entries was found.", LA_MATRIX_FORMAT_ERROR)
    return

    ! Zero Row Error
30  continue
    call errmgr%report_error("csr_lu_factor", &
        "A row with all zeros was encountered in the matrix.", &
        LA_SINGULAR_MATRIX_ERROR)
    return

    ! Unknown Error
40  continue
    call errmgr%report_error("csr_solve_sparse_direct", "ILUT encountered " // &
        "an unknown error.  The error code from the ILUT routine is " // &
        "provided in the output.", ierr)
    return

    ! Zero Pivot Error
50  continue
    call errmgr%report_error("csr_lu_factor", &
        "A zero pivot was encountered.", LA_SINGULAR_MATRIX_ERROR)
    return
end subroutine

! ------------------------------------------------------------------------------
subroutine form_lu_all(lu, ipvt, u, p, err)
    !! Extracts the L and U matrices from the condensed [L\\U] storage format 
    !! used by the lu_factor.
    real(real64), intent(inout), dimension(:,:) :: lu
        !! On input, the N-by-N matrix as output by lu_factor.  On output, the 
        !! N-by-N lower triangular matrix L.
    integer(int32), intent(in), dimension(:) :: ipvt
        !! The N-element pivot array as output by lu_factor.
    real(real64), intent(out), dimension(:,:) :: u
        !! An N-by-N matrix where the U matrix will be written.
    real(real64), intent(out), dimension(:,:) :: p
        !! An N-by-N matrix where the row permutation matrix will be written.
    class(errors), intent(inout), optional, target :: err
        !! The error object to be updated.

    ! Local Variables
    integer(int32) :: j, jp, n
    class(errors), pointer :: errmgr
    type(errors), target :: deferr

    ! Parameters
    real(real64), parameter :: zero = 0.0d0
    real(real64), parameter :: one = 1.0d0

    ! Initialization
    n = size(lu, 1)
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if

    ! Input Check
    if (size(lu, 2) /= n) then
        call report_square_matrix_error("form_lu_all", errmgr, "lu", n, &
            size(lu, 1), size(lu, 2))
        return
    else if (size(ipvt) /= n) then
        call report_array_size_error("form_lu_all", errmgr, "ipvt", n, &
            size(ipvt))
        return
    else if (size(u, 1) /= n .or. size(u, 2) /= n) then
        call report_matrix_size_error("form_lu_all", errmgr, "u", n, n, &
            size(u, 1), size(u, 2))
        return
    else if (size(p, 1) /= n .or. size(p, 2) /= n) then
        call report_matrix_size_error("form_lu_all", errmgr, "p", n, n, &
            size(p, 1), size(p, 2))
        return
    end if

    ! Ensure P starts off as an identity matrix
    call DLASET('A', n, n, zero, one, p, n)

    ! Process
    do j = 1, n
        ! Define the pivot matrix
        jp = ipvt(j)
        if (j /= jp) call swap(p(j,1:n), p(jp,1:n))

        ! Build L and U
        u(1:j,j) = lu(1:j,j)
        u(j+1:n,j) = zero

        if (j > 1) lu(1:j-1,j) = zero
        lu(j,j) = one
    end do
end subroutine

! ------------------------------------------------------------------------------
subroutine form_lu_all_cmplx(lu, ipvt, u, p, err)
    !! Extracts the L and U matrices from the condensed [L\\U] storage format 
    !! used by the lu_factor.
    complex(real64), intent(inout), dimension(:,:) :: lu
        !! On input, the N-by-N matrix as output by lu_factor.  On output, the 
        !! N-by-N lower triangular matrix L.
    integer(int32), intent(in), dimension(:) :: ipvt
        !! The N-element pivot array as output by lu_factor.
    complex(real64), intent(out), dimension(:,:) :: u
        !! An N-by-N matrix where the U matrix will be written.
    real(real64), intent(out), dimension(:,:) :: p
        !! An N-by-N matrix where the row permutation matrix will be written.
    class(errors), intent(inout), optional, target :: err
        !! The error object to be updated.

    ! Local Variables
    integer(int32) :: j, jp, n
    class(errors), pointer :: errmgr
    type(errors), target :: deferr

    ! Parameters
    real(real64), parameter :: zero = 0.0d0
    real(real64), parameter :: one = 1.0d0
    complex(real64), parameter :: c_zero = (0.0d0, 0.0d0)
    complex(real64), parameter :: c_one = (1.0d0, 0.0d0)

    ! Initialization
    n = size(lu, 1)
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if

    ! Input Check
    if (size(lu, 2) /= n) then
        call report_square_matrix_error("form_lu_all_cmplx", errmgr, "lu", n, &
            size(lu, 1), size(lu, 2))
        return
    else if (size(ipvt) /= n) then
        call report_array_size_error("form_lu_all_cmplx", errmgr, "ipvt", n, &
            size(ipvt))
        return
    else if (size(u, 1) /= n .or. size(u, 2) /= n) then
        call report_matrix_size_error("form_lu_all_cmplx", errmgr, "u", n, n, &
            size(u, 1), size(u, 2))
        return
    else if (size(p, 1) /= n .or. size(p, 2) /= n) then
        call report_matrix_size_error("form_lu_all_cmplx", errmgr, "p", n, n, &
            size(p, 1), size(p, 2))
        return
    end if

    ! Ensure P starts off as an identity matrix
    call DLASET('A', n, n, zero, one, p, n)

    ! Process
    do j = 1, n
        ! Define the pivot matrix
        jp = ipvt(j)
        if (j /= jp) call swap(p(j,1:n), p(jp,1:n))

        ! Build L and U
        u(1:j,j) = lu(1:j,j)
        u(j+1:n,j) = c_zero

        if (j > 1) lu(1:j-1,j) = c_zero
        lu(j,j) = c_one
    end do
end subroutine
! ------------------------------------------------------------------------------
subroutine form_lu_only(lu, u, err)
    !! Extracts the L and U matrices from the condensed [L\\U] storage format 
    !! used by the lu_factor.
    real(real64), intent(inout), dimension(:,:) :: lu
        !! On input, the N-by-N matrix as output by lu_factor.  On output, the 
        !! N-by-N lower triangular matrix L.
    real(real64), intent(out), dimension(:,:) :: u
        !! An N-by-N matrix where the U matrix will be written.
    class(errors), intent(inout), optional, target :: err
        !! The error object to be updated.

    ! Local Variables
    integer(int32) :: j, n
    class(errors), pointer :: errmgr
    type(errors), target :: deferr

    ! Parameters
    real(real64), parameter :: zero = 0.0d0
    real(real64), parameter :: one = 1.0d0

    ! Initialization
    n = size(lu, 1)
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if

    ! Input Check
    if (size(lu, 2) /= n) then
        call report_square_matrix_error("form_lu_only", errmgr, "lu", n, &
            size(lu, 1), size(lu, 2))
        return
    else if (size(u, 1) /= n .or. size(u, 2) /= n) then
        call report_matrix_size_error("form_lu_only", errmgr, "u", n, n, &
            size(u, 1), size(u, 2))
        return
    end if

    ! Process
    do j = 1, n
        ! Build L and U
        u(1:j,j) = lu(1:j,j)
        u(j+1:n,j) = zero

        if (j > 1) lu(1:j-1,j) = zero
        lu(j,j) = one
    end do
end subroutine

! ------------------------------------------------------------------------------
subroutine form_lu_only_cmplx(lu, u, err)
    !! Extracts the L and U matrices from the condensed [L\\U] storage format 
    !! used by the lu_factor.
    complex(real64), intent(inout), dimension(:,:) :: lu
        !! On input, the N-by-N matrix as output by lu_factor.  On output, the 
        !! N-by-N lower triangular matrix L.
    complex(real64), intent(out), dimension(:,:) :: u
        !! An N-by-N matrix where the U matrix will be written.
    class(errors), intent(inout), optional, target :: err
        !! The error object to be updated.

    ! Local Variables
    integer(int32) :: j, n
    class(errors), pointer :: errmgr
    type(errors), target :: deferr

    ! Parameters
    complex(real64), parameter :: zero = (0.0d0, 0.0d0)
    complex(real64), parameter :: one = (1.0d0, 0.0d0)

    ! Initialization
    n = size(lu, 1)
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if

    ! Input Check
    if (size(lu, 2) /= n) then
        call report_square_matrix_error("form_lu_only_cmplx", errmgr, "lu", n, &
            size(lu, 1), size(lu, 2))
        return
    else if (size(u, 1) /= n .or. size(u, 2) /= n) then
        call report_matrix_size_error("form_lu_only_cmplx", errmgr, "u", n, n, &
            size(u, 1), size(u, 2))
        return
    end if

    ! Process
    do j = 1, n
        ! Build L and U
        u(1:j,j) = lu(1:j,j)
        u(j+1:n,j) = zero

        if (j > 1) lu(1:j-1,j) = zero
        lu(j,j) = one
    end do
end subroutine

! ------------------------------------------------------------------------------
subroutine solve_lu_mtx(a, ipvt, b, err)
    !! Solves a system of LU-factored equations.
    real(real64), intent(in), dimension(:,:) :: a
        !! The N-by-N LU factored matrix as output by lu_factor.
    integer(int32), intent(in), dimension(:) :: ipvt
        !! The N-element pivot array as output by lu_factor.
    real(real64), intent(inout), dimension(:,:) :: b
        !! On input, the N-by-NRHS right-hand-side matrix.  On output, the 
        !! N-by-NRHS solution matrix.
    class(errors), intent(inout), optional, target :: err
        !! The error object to be updated.

    ! Local Variables
    integer(int32) :: n, nrhs, flag
    class(errors), pointer :: errmgr
    type(errors), target :: deferr

    ! Initialization
    n = size(a, 1)
    nrhs = size(b, 2)
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if

    ! Input Check
    if (size(a, 2) /= n) then
        call report_square_matrix_error("solve_lu_mtx", errmgr, "a", n, &
            size(a, 1), size(a, 2))
        return
    else if (size(ipvt) /= n) then
        call report_array_size_error("solve_lu_mtx", errmgr, "ipvt", n, &
            size(ipvt))
        return
    else if (size(b, 1) /= n) then
        call report_matrix_size_error("solve_lu_mtx", errmgr, "b", n, &
            size(b, 2), size(b, 1), size(b, 2))
        return
    end if

    ! Call DGETRS
    call DGETRS("N", n, nrhs, a, n, ipvt, b, n, flag)
end subroutine

! ------------------------------------------------------------------------------
subroutine solve_lu_mtx_cmplx(a, ipvt, b, err)
    !! Solves a system of LU-factored equations.
    complex(real64), intent(in), dimension(:,:) :: a
        !! The N-by-N LU factored matrix as output by lu_factor.
    integer(int32), intent(in), dimension(:) :: ipvt
        !! The N-element pivot array as output by lu_factor.
    complex(real64), intent(inout), dimension(:,:) :: b
        !! On input, the N-by-NRHS right-hand-side matrix.  On output, the 
        !! N-by-NRHS solution matrix.
    class(errors), intent(inout), optional, target :: err
        !! The error object to be updated.

    ! Local Variables
    integer(int32) :: n, nrhs, flag
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    character(len = :), allocatable :: errmsg

    ! Initialization
    n = size(a, 1)
    nrhs = size(b, 2)
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if

    ! Input Check
    if (size(a, 2) /= n) then
        call report_square_matrix_error("solve_lu_mtx_cmplx", errmgr, "a", n, &
            size(a, 1), size(a, 2))
        return
    else if (size(ipvt) /= n) then
        call report_array_size_error("solve_lu_mtx_cmplx", errmgr, "ipvt", n, &
            size(ipvt))
        return
    else if (size(b, 1) /= n) then
        call report_matrix_size_error("solve_lu_mtx_cmplx", errmgr, "b", n, &
            size(b, 2), size(b, 1), size(b, 2))
        return
    end if

    ! Call ZGETRS
    call ZGETRS("N", n, nrhs, a, n, ipvt, b, n, flag)
end subroutine

! ------------------------------------------------------------------------------
subroutine solve_lu_vec(a, ipvt, b, err)
    !! Solves a system of LU-factored equations.
    real(real64), intent(in), dimension(:,:) :: a
        !! The N-by-N LU factored matrix as output by lu_factor.
    integer(int32), intent(in), dimension(:) :: ipvt
        !! The N-element pivot array as output by lu_factor.
    real(real64), intent(inout), dimension(:) :: b
        !! The N-element right-hand-side array.
    class(errors), intent(inout), optional, target :: err
        !! The error object to be updated.

    ! Local Variables
    integer(int32) :: n, flag
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    character(len = :), allocatable :: errmsg

    ! Initialization
    n = size(a, 1)
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if

    ! Input Check
    if (size(a, 2) /= n) then
        call report_square_matrix_error("solve_lu_vec", errmgr, "a", n, &
            size(a, 1), size(a, 2))
        return
    else if (size(ipvt) /= n) then
        call report_array_size_error("solve_lu_vec", errmgr, "ipvt", n, &
            size(ipvt))
        return
    else if (size(b) /= n) then
        call report_array_size_error("solve_lu_vec", errmgr, "b", n, &
            size(b))
        return
    end if

    ! Call DGETRS
    call DGETRS("N", n, 1, a, n, ipvt, b, n, flag)
end subroutine

! ------------------------------------------------------------------------------
subroutine solve_lu_vec_cmplx(a, ipvt, b, err)
    !! Solves a system of LU-factored equations.
    complex(real64), intent(in), dimension(:,:) :: a
        !! The N-by-N LU factored matrix as output by lu_factor.
    integer(int32), intent(in), dimension(:) :: ipvt
        !! The N-element pivot array as output by lu_factor.
    complex(real64), intent(inout), dimension(:) :: b
        !! The N-element right-hand-side array.
    class(errors), intent(inout), optional, target :: err
        !! The error object to be updated.

    ! Local Variables
    integer(int32) :: n, flag
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    character(len = :), allocatable :: errmsg

    ! Initialization
    n = size(a, 1)
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if

    ! Input Check
    if (size(a, 2) /= n) then
        call report_square_matrix_error("solve_lu_vec_cmplx", errmgr, "a", n, &
            size(a, 1), size(a, 2))
        return
    else if (size(ipvt) /= n) then
        call report_array_size_error("solve_lu_vec_cmplx", errmgr, "ipvt", n, &
            size(ipvt))
        return
    else if (size(b) /= n) then
        call report_array_size_error("solve_lu_vec_cmplx", errmgr, "b", n, &
            size(b))
        return
    end if

    ! Call ZGETRS
    call ZGETRS("N", n, 1, a, n, ipvt, b, n, flag)
end subroutine

! ------------------------------------------------------------------------------
subroutine csr_lu_solve(lu, ju, b, x, err)
    !! Solves a linear system using an LU decomposition.
    class(msr_matrix), intent(in) :: lu
        !! The LU matrix.
    integer(int32), intent(in), dimension(:) :: ju
        !! The row tracking array.
    real(real64), intent(in), dimension(:) :: b
        !! The right-hand side.
    real(real64), intent(out), dimension(:) :: x
        !! The solution.
    class(errors), intent(inout), optional, target :: err
        !! The error object to be updated.

    ! Local Variables
    integer(int32) :: m, n
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    m = size(lu, 1)
    n = size(lu, 2)

    ! Input Check
    if (m /= n) then
        call report_square_matrix_error("csr_lu_solve", errmgr, "lu", m, m, n)
        return
    end if
    if (size(x) /= m) then
        call report_inner_matrix_dimension_error("csr_lu_solve", errmgr, &
            "lu", "x", m, size(x))
        return
    end if
    if (size(b) /= m) then
        call report_array_size_error("csr_lu_solve", errmgr, "b", m, size(b))
        return
    end if
    if (size(ju) /= m) then
        call report_array_size_error("csr_lu_solve", errmgr, "ju", m, size(ju))
        return
    end if

    ! Process
    call lusol(m, b, x, lu%values, lu%indices, ju)
end subroutine

! ------------------------------------------------------------------------------
end module