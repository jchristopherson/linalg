module linalg_lu
    use iso_fortran_env
    use linalg_errors
    use linalg_sparse
    use linalg_basic
    use lapack
    use sparskit
    use ieee_arithmetic, only : ieee_value, ieee_quiet_nan
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
pure subroutine lu_factor_dbl(a, ipvt, lu, l, u, p)
    !! Computes the LU factorization of an M-by-N matrix.  In the event of a 
    !! singular matrix, the output matrices are populated with NaN's.
    real(real64), intent(inout), dimension(:,:) :: a
        !! The N-by-N matrix to factor.
    integer(int32), intent(out), allocatable, optional, target, dimension(:) :: ipvt
        !! An N-element array used to track row-pivot operations.  The 
        !! array stored pivot information such that row I is interchanged with 
        !! row IPVT(I).
    real(real64), intent(out), allocatable, optional, target, dimension(:,:) :: lu
        !! The N-by-N factored matrix in the form [\(L\)\\ \(U\)] where the unit
        !! diagonal elements of \(L\) are not stored.
    real(real64), intent(out), allocatable, optional, target, dimension(:,:) :: l
        !! The N-by-N lower triangular matrix \(L\).
    real(real64), intent(out), allocatable, optional, target, dimension(:,:) :: u
        !! The N-by-N upper triangular matrix \(U\).
    real(real64), intent(out), allocatable, optional, target, dimension(:,:) :: p
        !! The N-by-N row permutation matrix.

    ! Local Variables
    logical :: buildlu
    integer(int32) :: n, flag
    integer(int32), allocatable, target, dimension(:) :: ic
    integer(int32), pointer, dimension(:) :: iptr
    real(real64) :: nan
    real(real64), allocatable, target, dimension(:,:) :: ac, lc, uc, pc
    real(real64), pointer, dimension(:,:) :: aptr, lptr, uptr, pptr

    ! Initialization
    n = size(a, 1)
    if (size(a, 2) /= n) error stop 1
    if (present(ipvt)) then
        allocate(ipvt(n))
        iptr => ipvt
    else
        allocate(ic(n))
        iptr => ic
    end if
    if (present(lu)) then
        allocate(lu(n, n), source = a)
        aptr => lu
    else
        allocate(ac(n, n), source = a)
        aptr => ac
    end if
    buildlu = present(l) .or. present(u) .or. present(p)

    ! Compute the LU factorization by calling the LAPACK routine DGETRF
    call DGETRF(n, n, aptr, n, iptr, flag)

    ! If flag > 0, the matrix is singular.  Notice, flag should not be
    ! able to be < 0 as we've already verrified inputs prior to making the
    ! call to LAPACK
    if (flag > 0) then
        ! WARNING: Singular matrix
        nan = ieee_value(nan, ieee_quiet_nan)
        if (present(lu)) lu = nan
        if (present(l)) allocate(l(n,n), source = nan)
        if (present(u)) allocate(u(n,n), source = nan)
        if (present(p)) allocate(p(n,n), source = nan)
        if (present(ipvt)) ipvt = 0
        return
    end if

    ! Build L & U?
    if (buildlu) then
        ! L
        if (present(l)) then
            allocate(l(n, n), source = aptr)
            lptr => l
        else
            if (allocated(ac)) then
                lptr => ac
            else
                allocate(lc(n, n), source = aptr)
                lptr => lc
            end if
        end if

        ! U
        if (present(u)) then
            allocate(u(n, n))
            uptr => u
        else
            allocate(uc(n, n))
            uptr => uc
        end if

        ! P
        if (present(p)) then
            allocate(p(n, n))
            pptr => p
        else
            allocate(pc(n, n))
            pptr => pc
        end if

        ! Build the matrices
        call form_lu(lptr, iptr, uptr, pptr)
    end if
end subroutine

! ------------------------------------------------------------------------------
pure subroutine lu_factor_cmplx(a, ipvt, lu, l, u, p)
    !! Computes the LU factorization of an M-by-N matrix.  In the event of a 
    !! singular matrix, the output matrices are populated with NaN's.
    complex(real64), intent(inout), dimension(:,:) :: a
        !! The N-by-N matrix to factor.
    integer(int32), intent(out), allocatable, optional, target, dimension(:) :: ipvt
        !! An N-element array used to track row-pivot operations.  The 
        !! array stored pivot information such that row I is interchanged with 
        !! row IPVT(I).
    complex(real64), intent(out), allocatable, optional, target, dimension(:,:) :: lu
        !! The N-by-N factored matrix in the form [\(L\)\\ \(U\)] where the unit
        !! diagonal elements of \(L\) are not stored.
    complex(real64), intent(out), allocatable, optional, target, dimension(:,:) :: l
        !! The N-by-N lower triangular matrix \(L\).
    complex(real64), intent(out), allocatable, optional, target, dimension(:,:) :: u
        !! The N-by-N upper triangular matrix \(U\).
    real(real64), intent(out), allocatable, optional, target, dimension(:,:) :: p
        !! The N-by-N row permutation matrix.

    ! Local Variables
    logical :: buildlu
    integer(int32) :: n, flag
    integer(int32), allocatable, target, dimension(:) :: ic
    integer(int32), pointer, dimension(:) :: iptr
    real(real64), allocatable, target, dimension(:,:) :: pc
    real(real64), pointer, dimension(:,:) :: pptr
    real(real64) :: nan
    complex(real64) :: cnan
    complex(real64), allocatable, target, dimension(:,:) :: ac, lc, uc
    complex(real64), pointer, dimension(:,:) :: aptr, lptr, uptr

    ! Initialization
    n = size(a, 1)
    if (size(a, 2) /= n) error stop 1
    if (present(ipvt)) then
        allocate(ipvt(n))
        iptr => ipvt
    else
        allocate(ic(n))
        iptr => ic
    end if
    if (present(lu)) then
        allocate(lu(n, n), source = a)
        aptr => lu
    else
        allocate(ac(n, n), source = a)
        aptr => ac
    end if
    buildlu = present(l) .or. present(u) .or. present(p)

    ! Compute the LU factorization by calling the LAPACK routine ZGETRF
    call ZGETRF(n, n, aptr, n, iptr, flag)

    ! If flag > 0, the matrix is singular.  Notice, flag should not be
    ! able to be < 0 as we've already verrified inputs prior to making the
    ! call to LAPACK
    if (flag > 0) then
        ! WARNING: Singular matrix
        nan = ieee_value(nan, ieee_quiet_nan)
        cnan = cmplx(nan, nan)
        if (present(lu)) lu = cnan
        if (present(l)) allocate(l(n,n), source = cnan)
        if (present(u)) allocate(u(n,n), source = cnan)
        if (present(p)) allocate(p(n,n), source = nan)
        if (present(ipvt)) ipvt = 0
        return
    end if

    ! Build L & U?
    if (buildlu) then
        ! L
        if (present(l)) then
            allocate(l(n, n), source = aptr)
            lptr => l
        else
            if (allocated(ac)) then
                lptr => ac
            else
                allocate(lc(n, n), source = aptr)
                lptr => lc
            end if
        end if

        ! U
        if (present(u)) then
            allocate(u(n, n))
            uptr => u
        else
            allocate(uc(n, n))
            uptr => uc
        end if

        ! P
        if (present(p)) then
            allocate(p(n, n))
            pptr => p
        else
            allocate(pc(n, n))
            pptr => pc
        end if

        ! Build the matrices
        call form_lu(lptr, iptr, uptr, pptr)
    end if
end subroutine

! ------------------------------------------------------------------------------
pure subroutine csr_lu_factor(a, lu, ju, droptol)
    !! Factors a matrix using an LU decomposition.
    class(csr_matrix), intent(in) :: a
        !! The matrix to factor.
    type(msr_matrix), intent(out) :: lu
        !! The LU matrix.
    integer(int32), intent(out), dimension(:) :: ju
        !! The row tracking array.
    real(real64), intent(in), optional :: droptol
        !! The drop tolerance for the ILUT factorization.

    ! Local Variables
    integer(int32) :: i, m, n, nn, nnz, lfil, iwk, ierr
    integer(int32), allocatable, dimension(:) :: jlu, jw
    real(real64), allocatable, dimension(:) :: alu, w
    real(real64) :: dt
    
    ! Initialization
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
        error stop 3
    end if

    ! Parameter Determination
    lfil = 1
    do i = 1, m
        lfil = max(lfil, a%row_indices(i+1) - a%row_indices(i))
    end do
    iwk = max(lfil * m, nnz)  ! somewhat arbitrary - can be adjusted

    ! Local Memory Allocation
    allocate(alu(iwk), w(n+1), jlu(iwk), jw(2 * n))

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
            error stop LA_MATRIX_FORMAT_ERROR
        else if (ierr == -2 .or. ierr == -3) then
            ! ALU and JLU are too small - try something larger
            ! This is the main reason for the loop - to offload worrying about
            ! workspace size from the user
            iwk = min(iwk + m + n, m * n)
            deallocate(alu)
            deallocate(jlu)
            allocate(alu(iwk), jlu(iwk))
        else if (ierr == -4) then
            ! Illegal value for LFIL - reset and try again
            lfil = n
        else if (ierr == -5) then
            ! Zero row encountered
            error stop LA_MATRIX_FORMAT_ERROR
        else
            ! We should never get here, but just in case
            error stop LA_INVALID_OPERATION_ERROR
        end if
    end do

    ! Determine the actual number of non-zero elements
    nnz = jlu(m+1) - 1

    ! Copy the contents to the output arrays
    lu%m = m
    lu%n = n
    lu%nnz = nnz
    nn = m + 1 + nnz - min(m, n)
    allocate(lu%values(nn), source = alu(:nn))
    allocate(lu%indices(nn), source = jlu(:nn))
end subroutine

! ------------------------------------------------------------------------------
pure subroutine form_lu_all(lu, ipvt, u, p)
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

    ! Local Variables
    integer(int32) :: j, jp, n

    ! Parameters
    real(real64), parameter :: zero = 0.0d0
    real(real64), parameter :: one = 1.0d0

    ! Initialization
    n = size(lu, 1)

    ! Input Check
    if (size(lu, 2) /= n) then
        error stop 1
    end if
    if (size(ipvt) /= n) then
        error stop 2
    end if
    if (size(u, 1) /= n .or. size(u, 2) /= n) then
        error stop 3
    end if
    if (size(p, 1) /= n .or. size(p, 2) /= n) then
        error stop 4
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
pure subroutine form_lu_all_cmplx(lu, ipvt, u, p)
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

    ! Local Variables
    integer(int32) :: j, jp, n

    ! Parameters
    real(real64), parameter :: zero = 0.0d0
    real(real64), parameter :: one = 1.0d0
    complex(real64), parameter :: c_zero = (0.0d0, 0.0d0)
    complex(real64), parameter :: c_one = (1.0d0, 0.0d0)

    ! Initialization
    n = size(lu, 1)

    ! Input Check
    if (size(lu, 2) /= n) then
        error stop 1
    end if
    if (size(ipvt) /= n) then
        error stop 2
    end if
    if (size(u, 1) /= n .or. size(u, 2) /= n) then
        error stop 3
    end if
    if (size(p, 1) /= n .or. size(p, 2) /= n) then
        error stop 4
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
pure subroutine form_lu_only(lu, u)
    !! Extracts the L and U matrices from the condensed [L\\U] storage format 
    !! used by the lu_factor.
    real(real64), intent(inout), dimension(:,:) :: lu
        !! On input, the N-by-N matrix as output by lu_factor.  On output, the 
        !! N-by-N lower triangular matrix L.
    real(real64), intent(out), dimension(:,:) :: u
        !! An N-by-N matrix where the U matrix will be written.

    ! Local Variables
    integer(int32) :: j, n

    ! Parameters
    real(real64), parameter :: zero = 0.0d0
    real(real64), parameter :: one = 1.0d0

    ! Initialization
    n = size(lu, 1)

    ! Input Check
    if (size(lu, 2) /= n) then
        error stop 1
    end if
    if (size(u, 1) /= n .or. size(u, 2) /= n) then
        error stop 2
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
pure subroutine form_lu_only_cmplx(lu, u)
    !! Extracts the L and U matrices from the condensed [L\\U] storage format 
    !! used by the lu_factor.
    complex(real64), intent(inout), dimension(:,:) :: lu
        !! On input, the N-by-N matrix as output by lu_factor.  On output, the 
        !! N-by-N lower triangular matrix L.
    complex(real64), intent(out), dimension(:,:) :: u
        !! An N-by-N matrix where the U matrix will be written.

    ! Local Variables
    integer(int32) :: j, n

    ! Parameters
    complex(real64), parameter :: zero = (0.0d0, 0.0d0)
    complex(real64), parameter :: one = (1.0d0, 0.0d0)

    ! Initialization
    n = size(lu, 1)

    ! Input Check
    if (size(lu, 2) /= n) then
        error stop 1
    end if
    if (size(u, 1) /= n .or. size(u, 2) /= n) then
        error stop 2
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
pure function solve_lu_mtx(a, ipvt, b) result(x)
    !! Solves a system of LU-factored equations.
    real(real64), intent(in), dimension(:,:) :: a
        !! The N-by-N LU factored matrix as output by lu_factor.
    integer(int32), intent(in), dimension(:) :: ipvt
        !! The N-element pivot array as output by lu_factor.
    real(real64), intent(in), dimension(:,:) :: b
        !! The N-by-NRHS right-hand-side matrix.
    real(real64), allocatable, dimension(:,:) :: x
        !! The N-by-NRHS solution matrix.

    ! Local Variables
    integer(int32) :: n, nrhs, flag

    ! Initialization
    n = size(a, 1)
    nrhs = size(b, 2)

    ! Input Check
    if (size(a, 2) /= n) then
        error stop 1
    end if
    if (size(ipvt) /= n) then
        error stop 2
    end if
    if (size(b, 1) /= n) then
        error stop 3
    end if

    ! Call DGETRS
    allocate(x(n, nrhs), source = b)
    call DGETRS("N", n, nrhs, a, n, ipvt, x, n, flag)
end function

! ------------------------------------------------------------------------------
pure function solve_lu_mtx_cmplx(a, ipvt, b) result(x)
    !! Solves a system of LU-factored equations.
    complex(real64), intent(in), dimension(:,:) :: a
        !! The N-by-N LU factored matrix as output by lu_factor.
    integer(int32), intent(in), dimension(:) :: ipvt
        !! The N-element pivot array as output by lu_factor.
    complex(real64), intent(in), dimension(:,:) :: b
        !! The N-by-NRHS right-hand-side matrix.
    complex(real64), allocatable, dimension(:,:) :: x
        !! The N-by-NRHS solution matrix.

    ! Local Variables
    integer(int32) :: n, nrhs, flag

    ! Initialization
    n = size(a, 1)
    nrhs = size(b, 2)

    ! Input Check
    if (size(a, 2) /= n) then
        error stop 1
    end if
    if (size(ipvt) /= n) then
        error stop 2
    end if
    if (size(b, 1) /= n) then
        error stop 3
    end if

    ! Call ZGETRS
    allocate(x(n, nrhs), source = b)
    call ZGETRS("N", n, nrhs, a, n, ipvt, x, n, flag)
end function

! ------------------------------------------------------------------------------
pure function solve_lu_vec(a, ipvt, b) result(x)
    !! Solves a system of LU-factored equations.
    real(real64), intent(in), dimension(:,:) :: a
        !! The N-by-N LU factored matrix as output by lu_factor.
    integer(int32), intent(in), dimension(:) :: ipvt
        !! The N-element pivot array as output by lu_factor.
    real(real64), intent(in), dimension(:) :: b
        !! The N-element right-hand-side array.
    real(real64), allocatable, dimension(:) :: x
        !! The N-element solution array.

    ! Local Variables
    integer(int32) :: n, flag

    ! Initialization
    n = size(a, 1)

    ! Input Check
    if (size(a, 2) /= n) then
        error stop 1
    end if
    if (size(ipvt) /= n) then
        error stop 2
    end if
    if (size(b) /= n) then
        error stop 3
    end if

    ! Call DGETRS
    allocate(x(n), source = b)
    call DGETRS("N", n, 1, a, n, ipvt, x, n, flag)
end function

! ------------------------------------------------------------------------------
pure function solve_lu_vec_cmplx(a, ipvt, b) result(x)
    !! Solves a system of LU-factored equations.
    complex(real64), intent(in), dimension(:,:) :: a
        !! The N-by-N LU factored matrix as output by lu_factor.
    integer(int32), intent(in), dimension(:) :: ipvt
        !! The N-element pivot array as output by lu_factor.
    complex(real64), intent(in), dimension(:) :: b
        !! The N-element right-hand-side array.
    complex(real64), allocatable, dimension(:) :: x
        !! The N-element solution array.

    ! Local Variables
    integer(int32) :: n, flag

    ! Initialization
    n = size(a, 1)

    ! Input Check
    if (size(a, 2) /= n) then
        error stop 1
    end if
    if (size(ipvt) /= n) then
        error stop 2
    end if
    if (size(b) /= n) then
        error stop 3
    end if

    ! Call ZGETRS
    allocate(x(n), source = b)
    call ZGETRS("N", n, 1, a, n, ipvt, x, n, flag)
end function

! ------------------------------------------------------------------------------
pure function csr_lu_solve(lu, ju, b) result(x)
    !! Solves a linear system using an LU decomposition.
    class(msr_matrix), intent(in) :: lu
        !! The LU matrix.
    integer(int32), intent(in), dimension(:) :: ju
        !! The row tracking array.
    real(real64), intent(in), dimension(:) :: b
        !! The right-hand side.
    real(real64), allocatable, dimension(:) :: x
        !! The solution.

    ! Local Variables
    integer(int32) :: m, n
    
    ! Initialization
    m = size(lu, 1)
    n = size(lu, 2)

    ! Input Check
    if (m /= n) then
        error stop 1
    end if
    if (size(ju) /= m) then
        error stop 2
    end if
    if (size(b) /= m) then
        error stop 3
    end if

    ! Process
    allocate(x(m), source = b)
    call lusol(m, b, x, lu%values, lu%indices, ju)
end function

! ------------------------------------------------------------------------------
end module