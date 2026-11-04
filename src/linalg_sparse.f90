! linalg_sparse.f90
! SPDX-License-Identifier: MIT
module linalg_sparse
    !! Provides sparse matrix types, conversions, arithmetic, and iterative or direct solvers.
    use iso_fortran_env, only : real64
    use sparskit
    use blas
    use linalg_errors
    use ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_is_finite
    implicit none
    private
    public :: csr_matrix
    public :: msr_matrix
    public :: size
    public :: create_empty_csr_matrix
    public :: create_empty_msr_matrix
    public :: nonzero_count
    public :: dense_to_csr
    public :: diag_to_csr
    public :: banded_to_csr
    public :: csr_to_dense
    public :: csr_to_msr
    public :: msr_to_csr
    public :: dense_to_msr
    public :: msr_to_dense
    public :: create_csr_matrix
    public :: matmul
    public :: operator(+)
    public :: operator(-)
    public :: operator(*)
    public :: operator(/)
    public :: assignment(=)
    public :: transpose
    public :: sparse_direct_solve
    public :: pgmres_solver
    public :: pgmres_solve
    public :: validate_csr_storage
    public :: validate_msr_lu

type :: csr_matrix
    !! A sparse matrix stored in compressed sparse row (CSR) format.
    integer, allocatable, dimension(:) :: row_indices
        !! An M+1 element array containing the indices in V an JA at which the
        !! requested row starts.
    integer, allocatable, dimension(:) :: column_indices
        !! An NNZ-element array, where NNZ is the number of non-zero values,
        !! containing the column indices of each value.
    real(real64), allocatable, dimension(:) :: values
        !! An NNZ-element array, where NNZ is the number of non-zero values,
        !! containing the non-zero values of the matrix.
    integer, private :: n = 0
        !! The number of columns in the matrix.
contains
    procedure, public :: get => csr_get_element
    procedure, public :: extract_diagonal => csr_extract_diagonal
end type

! ------------------------------------------------------------------------------
type :: msr_matrix
    !! A sparse matrix stored in modified sparse row format.
    integer, allocatable, dimension(:) :: indices
        !! An NNZ-element array containing the index information.
    real(real64), allocatable, dimension(:) :: values
        !! An NNZ-element array containing the non-zero values from the
        !! matrix.  The first MIN(M,N) elements contain the diagonal.
    integer :: m = 0
        !! The number of rows in the matrix.
    integer :: n = 0
        !! The number of columns in the matrix.
    integer :: nnz = 0
        !! The number of nonzero values in the matrix.
end type

! ------------------------------------------------------------------------------

interface nonzero_count
    module procedure :: nonzero_count_csr
    module procedure :: nonzero_count_msr
end interface

interface size
    module procedure :: csr_size
    module procedure :: msr_size
end interface

interface matmul
    module procedure :: csr_mtx_mtx_mult
    module procedure :: csr_mtx_vec_mult
end interface

interface operator(+)
    module procedure :: csr_mtx_add
end interface

interface operator(-)
    module procedure :: csr_mtx_sub
end interface

interface operator(*)
    module procedure :: csr_mtx_mult_scalar_1
    module procedure :: csr_mtx_mult_scalar_2
end interface

interface operator(/)
    module procedure :: csr_mtx_divide_scalar_1
end interface

interface assignment(=)
    module procedure :: csr_assign_to_dense
    module procedure :: dense_assign_to_csr
    module procedure :: msr_assign_to_dense
    module procedure :: dense_assign_to_msr
    module procedure :: csr_assign_to_msr
    module procedure :: msr_assign_to_csr
end interface

interface transpose
    module procedure :: csr_transpose
end interface

interface sparse_direct_solve
    module procedure :: csr_solve_sparse_direct
end interface

interface pgmres_solver
    module procedure :: csr_pgmres_solver
end interface

interface pgmres_solve
    !! Solves a CSR system with optional nonterminating status reporting.
    module procedure :: csr_pgmres_solve
end interface
contains
! ******************************************************************************
! CSR ROUTINES
! ------------------------------------------------------------------------------
function csr_get_element(this, i, j) result(rst)
    !! Retrieves the element at the specified row and column.
    class(csr_matrix), intent(in) :: this
        !! The CSR matrix object.
    integer, intent(in) :: i
        !! The row index.
    integer, intent(in) :: j
        !! The column index.
    real(real64) :: rst
        !! The value at the specified row and column.

    ! Local Variables
    integer :: iadd
    logical :: sorted

    ! Initialization
    sorted = .false.

    ! Input Checking
    if (i < 1 .or. i > size(this, 1)) error stop 2
    if (j < 1 .or. j > size(this, 2)) error stop 3

    ! Process
    if (.not.allocated(this%row_indices) .or. &
        .not.allocated(this%column_indices) .or. &
        .not.allocated(this%values)) &
    then
        rst = 0.0d0
        return
    end if
    rst = getelm(i, j, this%values, this%column_indices, this%row_indices, &
        iadd, sorted)
end function

! ------------------------------------------------------------------------------
pure function csr_size(x, dim) result(rst)
    !! Returns the size of the matrix along the specified dimension.
    class(csr_matrix), intent(in) :: x
        !! The CSR matrix object.
    integer, intent(in) :: dim
        !! The dimension to return the size of.
    integer :: rst
        !! The size of the matrix along the specified dimension.

    ! Process
    select case (dim)
    case (1)
        if (allocated(x%row_indices)) then
            rst = size(x%row_indices) - 1
        else
            rst = 0
        end if
    case (2)
        rst = x%n
    case default
        rst = 0
    end select
end function

! ------------------------------------------------------------------------------
pure function nonzero_count_csr(x) result(rst)
    !! Returns the number of non-zero values in the matrix.
    class(csr_matrix), intent(in) :: x
        !! The CSR matrix object.
    integer :: rst
        !! The number of non-zero values in the matrix.

    ! Process
    if (allocated(x%values)) then
        rst = size(x%values)
    else
        rst = 0
    end if
end function

! ------------------------------------------------------------------------------
pure function create_empty_csr_matrix(m, n, nnz) result(rst)
    !! Creates an empty CSR matrix.
    integer, intent(in) :: m
        !! The number of rows in the matrix.
    integer, intent(in) :: n
        !! The number of columns in the matrix.
    integer, intent(in) :: nnz
        !! The number of non-zero values in the matrix.
    type(csr_matrix) :: rst
        !! The empty CSR matrix.

    ! Local Variables
    integer :: m1
    
    ! Initialization
    m1 = m + 1

    ! Input Checking
    if (m < 0) then
        error stop 1
    end if
    if (n < 0) then
        error stop 2
    end if
    if (nnz < 0) then
        error stop 3
    end if

    ! Allocation
    rst%n = n
    allocate(rst%row_indices(m1), rst%column_indices(nnz), source = 0)
    allocate(rst%values(nnz), source = 0.0d0)
end function

! ------------------------------------------------------------------------------
pure function dense_to_csr(a) result(rst)
    !! Converts a dense matrix to a CSR matrix.
    real(real64), intent(in), dimension(:,:) :: a
        !! The dense matrix to convert.
    type(csr_matrix) :: rst
        !! The CSR matrix.

    ! Local Variables
    integer :: i, j, k, m, n, nnz
    real(real64) :: t
    
    ! Initialization
    t = 2.0d0 * epsilon(t)
    m = size(a, 1)
    n = size(a, 2)
    nnz = 0

    ! Determine how many non-zero values exist
    do j = 1, n
        do i = 1, m
            if (abs(a(i,j)) > t) then
                nnz = nnz + 1
            end if
        end do
    end do

    ! Memory Allocation
    rst = create_empty_csr_matrix(m, n, nnz)

    ! Store the non-zero values
    k = 1
    rst%row_indices(1) = 1
    do i = 1, m
        inner_loop : do j = 1, n
            if (abs(a(i,j)) < t) cycle inner_loop
            rst%column_indices(k) = j
            rst%values(k) = a(i,j)
            k = k + 1
        end do inner_loop
        rst%row_indices(i+1) = k
    end do
end function

! ------------------------------------------------------------------------------
pure function banded_to_csr(m, ml, mu, a) result(rst)
    !! Converts a banded matrix to a CSR matrix.
    integer, intent(in) :: m
        !! The number of rows in the banded matrix.
    integer, intent(in) :: ml
        !! The number of lower diagonals in the banded matrix.
    integer, intent(in) :: mu
        !! The number of upper diagonals in the banded matrix.
    real(real64), intent(in), dimension(:,:) :: a
        !! The banded matrix to convert.
    type(csr_matrix) :: rst
        !! The CSR matrix.

    ! Local Variables
    integer :: n, nnz, lowd, lda, flag
    integer, allocatable, dimension(:) :: ia, ja
    real(real64), allocatable, dimension(:) :: v
    
    ! Initialization
    lda = size(a, 1)
    n = size(a, 2)
    nnz = lda * n
    lowd = ml + mu + 1

    ! Input Checking
    if (ml < 0 .or. mu < 0) then
        error stop 2
    end if
    if (lda /= ml + mu + 1) then
        error stop 3
    end if

    ! Allocation
    allocate(ia(m + 1), ja(nnz), v(nnz))

    ! Process
    call bndcsr(m, a, lda, lowd, ml, mu, v, ja, ia, nnz, flag)
    nnz = ia(m + 1) - 1

    ! Put into the sparse matrix structure
    allocate(rst%row_indices(m + 1), source = ia)
    allocate(rst%column_indices(nnz), source = ja(:nnz))
    allocate(rst%values(nnz), source = v(:nnz))
    rst%n = n
end function

! ------------------------------------------------------------------------------
pure function csr_to_dense(a) result(x)
    !! Converts a CSR matrix to a dense matrix.
    class(csr_matrix), intent(in) :: a
        !! The CSR matrix to convert.
    real(real64), allocatable, dimension(:,:) :: x
        !! The dense matrix.

    ! Local Variables
    integer :: i, j, k, m, n, nnz
    
    ! Initialization
    m = size(a, 1)
    n = size(a, 2)
    nnz = nonzero_count(a)
    allocate(x(m, n), source = 0.0d0)
    
    ! Process
    do i = 1, m
        do k = a%row_indices(i), a%row_indices(i+1) - 1
            j = a%column_indices(k)
            x(i,j) = a%values(k)
        end do
    end do
end function

! ------------------------------------------------------------------------------
pure function diag_to_csr(a) result(rst)
    !! Converts a diagonal matrix to a CSR matrix.
    real(real64), intent(in), dimension(:) :: a
        !! The diagonal matrix to convert.
    type(csr_matrix) :: rst
        !! The CSR matrix.

    ! Local Variables
    integer :: i, n, n1
    
    ! Initialization
    n = size(a)
    n1 = n + 1

    ! Allocation
    allocate(rst%row_indices(n1), rst%column_indices(n))
    allocate(rst%values(n), source = a)
    rst%n = n

    ! Populate IA & JA
    do i = 1, n
        rst%column_indices(i) = i
        rst%row_indices(i) = i
    end do
    rst%row_indices(n1) = n1
end function

! ------------------------------------------------------------------------------
pure subroutine csr_assign_to_dense(dense, sparse)
    !! Assigns the values of a CSR matrix to a dense matrix.
    real(real64), intent(out), dimension(:,:) :: dense
        !! The dense matrix.
    class(csr_matrix), intent(in) :: sparse
        !! The CSR matrix.

    ! Process
    dense = csr_to_dense(sparse)
end subroutine

! ------------------------------------------------------------------------------
pure subroutine dense_assign_to_csr(sparse, dense)
    !! Assigns the values of a dense matrix to a CSR matrix.
    type(csr_matrix), intent(out) :: sparse
        !! The CSR matrix.
    real(real64), intent(in), dimension(:,:) :: dense
        !! The dense matrix.

    ! Process
    sparse = dense_to_csr(dense)
end subroutine

! ------------------------------------------------------------------------------
pure function csr_mtx_mtx_mult(a, b) result(rst)
    !! Multiplies two CSR matrices together.
    class(csr_matrix), intent(in) :: a
        !! The first CSR matrix.
    class(csr_matrix), intent(in) :: b
        !! The second CSR matrix.
    type(csr_matrix) :: rst
        !! The resulting CSR matrix.

    ! Local Variables
    integer, parameter :: sym_mult = 0
    integer, parameter :: full_mult = 1
    integer :: m, n, k, nnza, nnzb, nnzc, ierr
    integer, allocatable, dimension(:) :: ic, jc, iw
    real(real64) :: dummy(1)
    
    ! Initialization
    m = size(a, 1)
    n = size(b, 2)
    k = size(a, 2)
    nnza = nonzero_count(a)
    nnzb = nonzero_count(b)
    nnzc = nnza + nnzb

    ! Input Check
    if (size(b, 1) /= k) then
        error stop 2
    end if

    ! Local Memory Allocations
    allocate(ic(m + 1), jc(nnzc), iw(n))

    ! Determine the structure of C
    call amub(m, n, sym_mult, a%values, a%column_indices, a%row_indices, &
        b%values, b%column_indices, b%row_indices, dummy, jc, ic, &
        nnzc, iw, ierr)
    if (ierr /= 0) then
        ! NNZC was too small - try increasing it
        do while (ierr /= 0)
            deallocate(jc)
            nnzc = nnzc + nnza + nnzb
            allocate(jc(nnzc))
            call amub(m, n, sym_mult, a%values, a%column_indices, &
                a%row_indices, b%values, b%column_indices, b%row_indices, &
                dummy, jc, ic, nnzc, iw, ierr)
        end do
    end if

    ! Determine the actual NNZ for C & allocate space for the output
    nnzc = ic(m + 1) - 1
    deallocate(ic)
    deallocate(jc)
    rst = create_empty_csr_matrix(m, n, nnzc)

    ! Compute the actual product
    call amub(m, n, full_mult, a%values, a%column_indices, a%row_indices, &
        b%values, b%column_indices, b%row_indices, rst%values, &
        rst%column_indices, rst%row_indices, nnzc, iw, ierr)
end function

! ------------------------------------------------------------------------------
pure function csr_mtx_vec_mult(a, b) result(rst)
    !! Multiplies a CSR matrix by a vector.
    class(csr_matrix), intent(in) :: a
        !! The CSR matrix.
    real(real64), intent(in), dimension(:) :: b
        !! The vector.
    real(real64), allocatable, dimension(:) :: rst
        !! The resulting vector.

    ! Local Variables
    integer :: i, k, k1, k2, n, p
    real(real64) :: t

    ! Initialization
    n = size(a, 1)
    p = size(a, 2)

    ! Input Check
    if (size(b) /= p) then
        error stop 2
    end if

    ! Memory Allocation
    allocate(rst(n))

    ! Process
    do i = 1, n
        t = 0.0d0
        k1 = a%row_indices(i)
        k2 = a%row_indices(i+1) - 1
        do k = k1, k2
            t = t + a%values(k) * b(a%column_indices(k))
        end do
        rst(i) = t
    end do
end function

! ------------------------------------------------------------------------------
pure function csr_mtx_add(a, b) result(rst)
    !! Adds two CSR matrices.
    class(csr_matrix), intent(in) :: a
        !! The first CSR matrix.
    class(csr_matrix), intent(in) :: b
        !! The second CSR matrix.
    type(csr_matrix) :: rst
        !! The resulting CSR matrix.

    ! Local Variables
    integer, parameter :: sym_add = 0
    integer, parameter :: full_add = 1
    integer :: m, n, nnza, nnzb, nnzc, ierr
    integer, allocatable, dimension(:) :: ic, jc, iw
    real(real64) :: dummy(1)
    
    ! Initialization
    m = size(a, 1)
    n = size(a, 2)
    nnza = nonzero_count(a)
    nnzb = nonzero_count(b)
    nnzc = nnza + nnzb

    ! Input Checking
    if (size(b, 1) /= m .or. size(b, 2) /= n) then
        error stop 2
    end if

    ! Local Memory Allocations
    allocate(ic(m + 1), jc(nnzc), iw(n))

    ! Determine the structure of C
    call aplb(m, n, sym_add, a%values, a%column_indices, a%row_indices, &
        b%values, b%column_indices, b%row_indices, dummy, jc, ic, &
        nnzc, iw, ierr)
    if (ierr /= 0) then
        ! NNZC was too small - try increasing it
        do while (ierr /= 0)
            deallocate(jc)
            nnzc = nnzc + nnza + nnzb
            allocate(jc(nnzc))
            call aplb(m, n, sym_add, a%values, a%column_indices, &
                a%row_indices, b%values, b%column_indices, b%row_indices, &
                dummy, jc, ic, nnzc, iw, ierr)
        end do
    end if

    ! Determine the actuall NNZ for C & allocate space for the output
    nnzc = ic(m + 1) - 1
    deallocate(ic)
    deallocate(jc)
    rst = create_empty_csr_matrix(m, n, nnzc)

    ! Compute the actual sum
    call aplb(m, n, full_add, a%values, a%column_indices, a%row_indices, &
        b%values, b%column_indices, b%row_indices, rst%values, &
        rst%column_indices, rst%row_indices, nnzc, iw, ierr)
end function

! ------------------------------------------------------------------------------
pure function csr_mtx_sub(a, b) result(rst)
    !! Subtracts two CSR matrices.
    class(csr_matrix), intent(in) :: a
        !! The first CSR matrix.
    class(csr_matrix), intent(in) :: b
        !! The second CSR matrix.
    type(csr_matrix) :: rst
        !! The resulting CSR matrix.

    ! Local Variables
    integer, parameter :: sym_add = 0
    integer :: m, n, nnza, nnzb, nnzc, ierr
    integer, allocatable, dimension(:) :: ic, jc, iw
    real(real64) :: dummy(1)
    
    ! Initialization
    m = size(a, 1)
    n = size(a, 2)
    nnza = nonzero_count(a)
    nnzb = nonzero_count(b)
    nnzc = nnza + nnzb

    ! Input Checking
    if (size(b, 1) /= m .or. size(b, 2) /= n) then
        error stop 2
    end if

    ! Local Memory Allocations
    allocate(ic(m + 1), jc(nnzc), iw(n))

    ! Determine the structure of C
    call aplb(m, n, sym_add, a%values, a%column_indices, a%row_indices, &
        b%values, b%column_indices, b%row_indices, dummy, jc, ic, &
        nnzc, iw, ierr)
    if (ierr /= 0) then
        ! NNZC was too small - try increasing it
        do while (ierr /= 0)
            deallocate(jc)
            nnzc = nnzc + nnza + nnzb
            allocate(jc(nnzc))
            call aplb(m, n, sym_add, a%values, a%column_indices, &
                a%row_indices, b%values, b%column_indices, b%row_indices, &
                dummy, jc, ic, nnzc, iw, ierr)
        end do
    end if

    ! Determine the actuall NNZ for C & allocate space for the output
    nnzc = ic(m + 1) - 1
    deallocate(ic)
    deallocate(jc)
    rst = create_empty_csr_matrix(m, n, nnzc)

    ! Compute the actual sum
    call aplsb(m, n, a%values, a%column_indices, a%row_indices, -1.0d0, &
        b%values, b%column_indices, b%row_indices, rst%values, &
        rst%column_indices, rst%row_indices, nnzc, iw, ierr)
end function

! ------------------------------------------------------------------------------
pure function csr_mtx_mult_scalar_1(a, b) result(rst)
    !! Multiplies a CSR matrix by a scalar.
    class(csr_matrix), intent(in) :: a
        !! The CSR matrix.
    real(real64), intent(in) :: b
        !! The scalar.
    type(csr_matrix) :: rst
        !! The resulting CSR matrix.

    ! Local Variables
    integer :: m, n, nnz
    
    ! Initialization
    m = size(a, 1)
    n = size(a, 2)
    nnz = nonzero_count(a)

    ! Process
    rst = create_empty_csr_matrix(m, n, nnz)

    ! Compute the product
    rst%row_indices = a%row_indices
    rst%column_indices = a%column_indices
    rst%values = b * a%values
end function

! ------------------------------------------------------------------------------
pure function csr_mtx_mult_scalar_2(a, b) result(rst)
    !! Multiplies a scalar by a CSR matrix.
    real(real64), intent(in) :: a
        !! The scalar.
    class(csr_matrix), intent(in) :: b
        !! The CSR matrix.
    type(csr_matrix) :: rst
        !! The resulting CSR matrix.

    ! Local Variables
    integer :: m, n, nnz
    
    ! Initialization
    m = size(b, 1)
    n = size(b, 2)
    nnz = nonzero_count(b)

    ! Process
    rst = create_empty_csr_matrix(m, n, nnz)

    ! Compute the product
    rst%row_indices = b%row_indices
    rst%column_indices = b%column_indices
    rst%values = a * b%values
end function

! ------------------------------------------------------------------------------
pure function csr_mtx_divide_scalar_1(a, b) result(rst)
    !! Divides a CSR matrix by a scalar.
    class(csr_matrix), intent(in) :: a
        !! The CSR matrix.
    real(real64), intent(in) :: b
        !! The scalar.
    type(csr_matrix) :: rst

    ! Local Variables
    integer :: m, n, nnz
    
    ! Initialization
    m = size(a, 1)
    n = size(a, 2)
    nnz = nonzero_count(a)

    ! Process
    rst = create_empty_csr_matrix(m, n, nnz)

    ! Compute the product
    rst%row_indices = a%row_indices
    rst%column_indices = a%column_indices
    rst%values = a%values / b
end function

! ------------------------------------------------------------------------------
pure function csr_transpose(a) result(rst)
    !! Transposes a CSR matrix.
    class(csr_matrix), intent(in) :: a
        !! The CSR matrix.
    type(csr_matrix) :: rst
        !! The transposed CSR matrix.

    ! Local Variables
    integer :: m, n, nnz

    ! Initialization
    m = size(a, 1)
    n = size(a, 2)
    nnz = nonzero_count(a)
    rst = create_empty_csr_matrix(n, m, nnz)

    ! Process
    call csrcsc2(m, n, 1, 1, a%values, a%column_indices, a%row_indices, &
        rst%values, rst%column_indices, rst%row_indices)
end function

! ------------------------------------------------------------------------------
pure function csr_extract_diagonal(a) result(diag)
    !! Extracts the diagonal from a CSR matrix.
    class(csr_matrix), intent(in) :: a
        !! The CSR matrix.
    real(real64), allocatable, dimension(:) :: diag
        !! The diagonal values.

    ! Local Variables
    integer :: m, n, mn, len
    integer, allocatable, dimension(:) :: idiag
    
    ! Initialization
    m = size(a, 1)
    n = size(a, 2)
    mn = min(m, n)
    allocate(diag(mn), idiag(mn))

    ! Process
    call getdia(m, n, 0, a%values, a%column_indices, a%row_indices, len, &
        diag, idiag, 0)
end function

! ------------------------------------------------------------------------------
pure function validate_csr_storage(a) result(flag)
    !! Validates nonempty finite CSR storage before raw-array solver calls.
    !! Rectangular matrices and zero rows are accepted by this storage check;
    !! individual solvers impose their own square-system requirements.
    class(csr_matrix), intent(in) :: a
        !! Matrix to inspect without modification. Rows may be unsorted, but
        !! column indices within each row must be unique and in range.
    integer :: flag
        !! LA_NO_ERROR on success; LA_NOT_INITIALIZED_ERROR for missing arrays;
        !! LA_INVALID_INPUT_ERROR for empty dimensions or nonfinite values;
        !! LA_MATRIX_FORMAT_ERROR for inconsistent lengths, pointers, or columns.
    integer :: rows, columns, stored, row, entry
    integer, allocatable :: seen(:)

    flag = LA_NOT_INITIALIZED_ERROR
    if (.not.allocated(a%row_indices) .or. &
        .not.allocated(a%column_indices) .or. .not.allocated(a%values)) return
    rows = size(a, 1)
    columns = size(a, 2)
    flag = LA_INVALID_INPUT_ERROR
    if (rows < 1 .or. columns < 1) return
    if (.not.all(ieee_is_finite(a%values))) return
    stored = size(a%values)
    flag = LA_MATRIX_FORMAT_ERROR
    if (size(a%column_indices) /= stored) return
    if (a%row_indices(1) /= 1 .or. a%row_indices(rows+1) /= stored+1) return
    if (any(a%row_indices < 1) .or. any(a%row_indices > stored+1)) return
    if (any(a%row_indices(2:) < a%row_indices(:rows))) return
    if (any(a%column_indices < 1) .or. any(a%column_indices > columns)) return
    allocate(seen(columns), source = 0)
    do row = 1, rows
        do entry = a%row_indices(row), a%row_indices(row+1)-1
            if (seen(a%column_indices(entry)) == row) return
            seen(a%column_indices(entry)) = row
        end do
    end do
    flag = LA_NO_ERROR
end function

pure function validate_msr_lu(lu, ju) result(flag)
    !! Validates the combined MSR triangular factors consumed by LUSOL.
    !! The diagonal contains inverse U pivots, not the original U diagonal.
    !! The unused value at N+1 and unused trailing capacity are not inspected.
    class(msr_matrix), intent(in) :: lu
        !! Nonempty square factors with row pointers in INDICES(1:N+1),
        !! followed by lower entries and then upper entries within each row.
    integer, intent(in) :: ju(:)
        !! N pointers splitting each row into its strict lower and upper parts.
    integer :: flag
        !! LA_NO_ERROR on success; LA_NOT_INITIALIZED_ERROR for missing arrays;
        !! LA_MATRIX_NOT_SQUARE_ERROR or LA_MATRIX_DIMENSION_ERROR for sizes;
        !! LA_INVALID_INPUT_ERROR for empty factors; LA_MATRIX_FORMAT_ERROR
        !! for invalid pointers or nontriangular columns; or
        !! LA_SINGULAR_MATRIX_ERROR for zero inverse pivots/nonfinite factors.
    integer :: order, stored, row, first, split, last

    flag = LA_MATRIX_NOT_SQUARE_ERROR
    order = size(lu, 1)
    if (size(lu, 2) /= order) return
    flag = LA_MATRIX_DIMENSION_ERROR
    if (size(ju) /= order) return
    flag = LA_NOT_INITIALIZED_ERROR
    if (.not.allocated(lu%values) .or. .not.allocated(lu%indices)) return
    flag = LA_INVALID_INPUT_ERROR
    if (order < 1) return
    stored = size(lu%values)
    flag = LA_MATRIX_FORMAT_ERROR
    if (size(lu%indices) /= stored .or. stored < order+1) return
    if (lu%indices(1) /= order+2) return
    if (any(lu%indices(:order+1) < order+2) .or. &
        any(lu%indices(:order+1) > stored+1)) return
    if (any(lu%indices(2:order+1) < lu%indices(:order))) return
    do row = 1, order
        first = lu%indices(row)
        split = ju(row)
        last = lu%indices(row+1) - 1
        if (split < first .or. split > last+1) return
        if (any(lu%indices(first:split-1) < 1) .or. &
            any(lu%indices(first:split-1) >= row)) return
        if (any(lu%indices(split:last) <= row) .or. &
            any(lu%indices(split:last) > order)) return
    end do
    flag = LA_SINGULAR_MATRIX_ERROR
    if (any(lu%values(:order) == 0.0d0)) return
    if (.not.all(ieee_is_finite(lu%values(:order)))) return
    if (.not.all(ieee_is_finite(lu%values(order+2:lu%indices(order+1)-1)))) return
    flag = LA_NO_ERROR
end function

pure function csr_solve_sparse_direct(a, b, droptol) result(x)
    !! Legacy convenience entry point for solving the original CSR system.
    !! ILUT builds a preconditioner and GMRES refines the full-system solution;
    !! despite the historical name, this is not a pivoted sparse direct solver.
    !! Detected ILUT breakdown yields NaNs; nonconvergence uses ERROR STOP.
    class(csr_matrix), intent(in) :: a
        !! Nonempty N-by-N real coefficient matrix in valid CSR format.
    real(real64), intent(in), dimension(:) :: b
        !! N-element right-hand side, not modified.
    real(real64), intent(in), optional :: droptol
        !! Finite nonnegative ILUT drop threshold; default SQRT(EPSILON).
        !! Controls the preconditioner, not final accuracy. GMRES uses its
        !! default relative residual tolerance SQRT(EPSILON), restart length
        !! MIN(N,50), and iteration budget 100. Use PGMRES_SOLVE for controls.
    real(real64), allocatable, dimension(:) :: x
        !! N-element solution of A * X = B to the GMRES residual tolerance.
        !! Inputs are preserved; an explicitly detected breakdown yields NaNs.

    ! Local Variables
    integer :: i, m, n, nnz, lfil, iwk, ierr
    integer, allocatable, dimension(:) :: jlu, ju, jw
    real(real64), allocatable, dimension(:) :: alu, w
    real(real64) :: dt, nan
    type(msr_matrix) :: factors
    
    ! Initialization
    if (present(droptol)) then
        dt = droptol
    else
        dt = sqrt(epsilon(dt))
    end if
    m = size(a, 1)
    n = size(a, 2)
    nnz = nonzero_count(a)
    allocate(x(n))

    ! Input Checking
    if (m /= n) then
        error stop 1
    end if
    if (size(b) /= n) then
        error stop 2
    end if
    ierr = validate_csr_storage(a)
    if (ierr /= LA_NO_ERROR) error stop ierr
    if (.not.ieee_is_finite(dt)) error stop LA_INVALID_INPUT_ERROR
    if (dt < 0.0d0 .or. .not.all(ieee_is_finite(b))) &
        error stop LA_INVALID_INPUT_ERROR
    do i = 1, n
        if (all(a%values(a%row_indices(i):a%row_indices(i+1)-1) == 0.0d0)) then
            nan = ieee_value(0.0d0, ieee_quiet_nan)
            x = nan
            return
        end if
    end do

    ! Parameter Determination
    lfil = 1
    do i = 1, m
        lfil = max(lfil, a%row_indices(i+1) - a%row_indices(i))
    end do
    iwk = max(n + 1, lfil * m + 1, nnz)

    ! Local Memory Allocation
    allocate(alu(iwk), w(n+1), jlu(iwk), ju(n), jw(2 * n))

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
            ! WARNING: Singular matrix - a zero pivot was found at step IERR
            nan = ieee_value(nan, ieee_quiet_nan)
            x = nan
            return
        else if (ierr == -1) then
            ! The input matrix is not formatted correctly
            error stop LA_MATRIX_FORMAT_ERROR
        else if (ierr == -2 .or. ierr == -3) then
            ! ALU and JLU are too small - try something larger
            iwk = iwk + max(m + n, iwk / 2)
            deallocate(alu)
            deallocate(jlu)
            allocate(alu(iwk), jlu(iwk))
        else if (ierr == -4) then
            ! Illegal value for LFIL - reset and try again
            lfil = n
        else if (ierr == -5) then
            ! Zero row encountered
            x = ieee_value(0.0d0, ieee_quiet_nan)
            return
        else
            ! We should never get here, but just in case
            error stop LA_INVALID_OPERATION_ERROR
        end if
    end do

    ! Incomplete factors precondition, rather than replace, the original system.
    factors%m = n
    factors%n = n
    factors%nnz = jlu(n+1) - 1
    factors%values = alu(:factors%nnz)
    factors%indices = jlu(:factors%nnz)
    factors%values(n+1) = 0.0d0
    x = csr_pgmres_solver(a, factors, ju, b)
end function

! ******************************************************************************
! MSR ROUTINES
! ------------------------------------------------------------------------------
! TO DO: MSR_GET_ELEMENT

! ------------------------------------------------------------------------------
pure function msr_size(x, dim) result(rst)
    !! Returns the size of the specified dimension of an MSR matrix.
    class(msr_matrix), intent(in) :: x
        !! The MSR matrix.
    integer, intent(in) :: dim
        !! The dimension to return the size of.
    integer :: rst
        !! The size of the specified dimension.

    ! Process
    select case (dim)
    case (1)
        rst = x%m
    case (2)
        rst = x%n
    case default
        rst = 0
    end select
end function

! ------------------------------------------------------------------------------
pure function nonzero_count_msr(x) result(rst)
    !! Returns the number of non-zero elements in an MSR matrix.
    class(msr_matrix), intent(in) :: x
        !! The MSR matrix.
    integer :: rst
        !! The number of non-zero elements.

    ! Process
    rst = x%nnz
end function

! ------------------------------------------------------------------------------
pure function create_empty_msr_matrix(m, n, nnz) result(rst)
    !! Creates an empty MSR matrix.
    integer, intent(in) :: m
        !! The number of rows in the matrix.
    integer, intent(in) :: n
        !! The number of columns in the matrix.
    integer, intent(in) :: nnz
        !! The number of non-zero elements in the matrix.
    type(msr_matrix) :: rst
        !! The MSR matrix.

    ! Local Variables
    integer :: nelem, mn
    
    ! Input Checking
    if (m < 0) then
        error stop 1
    end if
    if (n < 0) then
        error stop 2
    end if
    if (nnz < 0) then
        error stop 3
    end if

    ! Allocation
    rst%m = m
    rst%n = n
    rst%nnz = nnz
    mn = min(m, n)
    nelem = m + 1 + nnz - mn
    allocate(rst%indices(nelem), source = 0)
    allocate(rst%values(nelem), source = 0.0d0)
end function

! ------------------------------------------------------------------------------
pure function csr_to_msr(a) result(rst)
    !! Converts a CSR matrix to an MSR matrix.
    class(csr_matrix), intent(in) :: a
        !! The CSR matrix to convert.
    type(msr_matrix) :: rst
        !! The MSR matrix.

    ! Local Variables
    integer :: m, n, nnz
    integer, allocatable, dimension(:) :: iwork, jc, ic
    real(real64), allocatable, dimension(:) :: work, ac
    
    ! Initialization
    m = size(a, 1)
    n = size(a, 2)
    nnz = nonzero_count(a)

    ! Memory Allocation
    rst = create_empty_msr_matrix(m, n, nnz)
    allocate(work(m), iwork(m + 1))
    allocate(ac(nnz), source = a%values)
    allocate(jc(nnz), source = a%column_indices)
    allocate(ic(m+1), source = a%row_indices)

    ! Perform the conversion
    call csrmsr(m, ac, jc, ic, rst%values, rst%indices, work, iwork)
end function

! ------------------------------------------------------------------------------
pure function msr_to_csr(a) result(rst)
    !! Converts an MSR matrix to a CSR matrix.
    class(msr_matrix), intent(in) :: a
        !! The MSR matrix to convert.
    type(csr_matrix) :: rst
        !! The CSR matrix.

    ! Local Variables
    integer :: m, n, nnz, flag
    integer, allocatable, dimension(:) :: iwork
    real(real64), allocatable, dimension(:) :: work
    
    ! Initialization
    m = size(a, 1)
    n = size(a, 2)
    nnz = nonzero_count(a)

    ! Memory Allocation
    rst = create_empty_csr_matrix(m, n, nnz)
    allocate(work(m), iwork(m+1))

    ! Process
    call msrcsr(m, a%values, a%indices, rst%values, rst%column_indices, &
        rst%row_indices, work, iwork)
end function

! ------------------------------------------------------------------------------
pure function dense_to_msr(a) result(rst)
    !! Converts a dense matrix to an MSR matrix.
    real(real64), intent(in), dimension(:,:) :: a
        !! The dense matrix to convert.
    type(msr_matrix) :: rst
        !! The MSR matrix.

    ! Local Variables
    type(csr_matrix) :: csr

    ! Convert to CSR, and then from CSR to MSR
    csr = dense_to_csr(a)

    ! Convert to MSR
    rst = csr_to_msr(csr)
end function

! ------------------------------------------------------------------------------
pure function msr_to_dense(a) result(x)
    !! Converts an MSR matrix to a dense matrix.
    class(msr_matrix), intent(in) :: a
        !! The MSR matrix to convert.
    real(real64), allocatable, dimension(:,:) :: x
        !! The dense matrix.

    ! Local Variables
    integer :: m, n
    type(csr_matrix) :: csr
    
    ! Initialization
    m = size(a, 1)
    n = size(a, 2)

    ! Process
    csr = msr_to_csr(a)
    x = csr_to_dense(csr)
end function

! ------------------------------------------------------------------------------
pure subroutine msr_assign_to_dense(dense, msr)
    !! Assigns an MSR matrix to a dense matrix.
    real(real64), intent(out), dimension(:,:) :: dense
        !! The dense matrix.
    class(msr_matrix), intent(in) :: msr
        !! The MSR matrix.

    ! Process
    dense = msr_to_dense(msr)
end subroutine

! ------------------------------------------------------------------------------
pure subroutine dense_assign_to_msr(msr, dense)
    !! Assigns a dense matrix to an MSR matrix.
    type(msr_matrix), intent(out) :: msr
        !! The MSR matrix.
    real(real64), intent(in), dimension(:,:) :: dense
        !! The dense matrix.

    ! Process
    msr = dense_to_msr(dense)
end subroutine

! ------------------------------------------------------------------------------
pure subroutine csr_assign_to_msr(msr, csr)
    !! Assigns a CSR matrix to an MSR matrix.
    type(msr_matrix), intent(out) :: msr
        !! The MSR matrix.
    class(csr_matrix), intent(in) :: csr
        !! The CSR matrix.

    ! Process
    msr = csr_to_msr(csr)
end subroutine

! ------------------------------------------------------------------------------
pure subroutine msr_assign_to_csr(csr, msr)
    !! Assigns an MSR matrix to a CSR matrix.
    type(csr_matrix), intent(out) :: csr
        !! The CSR matrix.
    class(msr_matrix), intent(in) :: msr
        !! The MSR matrix.

    ! Process
    csr = msr_to_csr(msr)
end subroutine

! ------------------------------------------------------------------------------
pure function create_csr_matrix(m, n, rows, cols, vals) result(rst)
    !! Creates a CSR matrix from the input data.
    integer, intent(in) :: m
        !! The number of rows in the matrix.
    integer, intent(in) :: n
        !! The number of columns in the matrix.
    integer, intent(in), dimension(:) :: rows
        !! The row indices.
    integer, intent(in), dimension(:) :: cols
        !! The column indices.
    real(real64), intent(in), dimension(:) :: vals
        !! The values.
    type(csr_matrix) :: rst
        !! The CSR matrix.

    ! Local Variables
    integer :: i, nnz
    integer, allocatable, dimension(:) :: ir
    
    ! Initialization
    nnz = size(rows)

    ! Input Checking
    if (m < 0) then
        error stop 1
    end if
    if (n < 0) then
        error stop 2
    end if
    if (size(cols) /= nnz .or. size(vals) /= nnz) then
        error stop 5
    end if
    do i = 1, nnz
        if (rows(i) < 1 .or. rows(i) > m) then
            error stop 3
        end if
        if (cols(i) < 1 .or. cols(i) > n) then
            error stop 4
        end if
    end do
    allocate(ir(nnz), source = rows)

    ! Create an empty matrix
    rst = create_empty_csr_matrix(m, n, nnz)

    ! Populate the empty matrix
    call coocsr(m, nnz, vals, ir, cols, rst%values, rst%column_indices, &
        rst%row_indices)
    call csort(m, rst%values, rst%column_indices, rst%row_indices, .true.)
end function

! ******************************************************************************
! ITERATIVE SOLVERS
! ------------------------------------------------------------------------------
! Additional References:
! - https://www.diva-portal.org/smash/get/diva2:360739/FULLTEXT01.pdf
pure function csr_pgmres_solver(a, lu, ju, b, im, tol, maxits, iout) result(x)
    !! Function form of [[pgmres_solve]], retaining the existing pure API.
    !! Invalid inputs and failure of the true-residual test use ERROR STOP.
    class(csr_matrix), intent(in) :: a
        !! Nonempty N-by-N coefficient matrix in valid CSR format.
    class(msr_matrix), intent(in) :: lu
        !! Prepared N-by-N ILUT preconditioner in MSR format.
    integer, intent(in) :: ju(:)
        !! N pointers to upper-factor rows in LU.
    real(real64), intent(in) :: b(:)
        !! N-element right-hand side, not modified.
    integer, intent(in), optional :: im
        !! Restart length in [1, MIN(N, 50)]; default MIN(N, 50).
    real(real64), intent(in), optional :: tol
        !! Finite relative residual tolerance at least EPSILON; default
        !! SQRT(EPSILON). Success requires NORM2(B-A*X) <= TOL*NORM2(B).
    integer, intent(in), optional :: maxits
        !! Positive iteration budget, checked at restart boundaries; default 100.
    integer, intent(in), optional :: iout
        !! Output unit for iteration logs; nonpositive suppresses logs, default 0.
    real(real64), allocatable :: x(:)
        !! Computed N-element solution of the original system A * X = B.

    call csr_pgmres_solve(a, lu, ju, b, x, im, tol, maxits, iout)
end function

pure subroutine csr_pgmres_solve(a, lu, ju, b, x, im, tol, maxits, iout, flag)
    !! Solves A * X = B using right-preconditioned restarted GMRES, starting
    !! from zero. The incomplete factors are only a preconditioner; convergence
    !! is checked using the true residual of the original CSR system.
    class(csr_matrix), intent(in) :: a
        !! The nonempty N-by-N coefficient matrix in valid CSR format.
    class(msr_matrix), intent(in) :: lu
        !! Prepared N-by-N ILUT factors in MSR format, not modified.
    integer, intent(in), dimension(:) :: ju
        !! N pointers to the starts of the upper-factor rows in LU.
    real(real64), intent(in), dimension(:) :: b
        !! N-element right-hand side, not modified. A zero vector yields zero.
    integer, intent(in), optional :: im
        !! Restart length in [1, MIN(N, 50)]; default MIN(N, 50). The upper
        !! bound is imposed by the bundled SPARSKIT implementation.
    integer, intent(in), optional :: maxits
        !! Positive iteration budget, default 100. SPARSKIT checks this budget
        !! at restart boundaries and may exceed it by up to IM-1 iterations.
    integer, intent(in), optional :: iout
        !! Fortran output unit for iteration information. Nonpositive values
        !! suppress output; default zero.
    real(real64), intent(in), optional :: tol
        !! Finite relative residual tolerance, at least machine epsilon;
        !! default SQRT(EPSILON). Success requires NORM2(B-A*X) <= TOL*NORM2(B).
    integer, intent(out), optional :: flag
        !! LA_NO_ERROR on success; LA_CONVERGENCE_ERROR if iteration or the
        !! true-residual check fails. Invalid dimensions return
        !! LA_MATRIX_DIMENSION_ERROR or LA_MATRIX_NOT_SQUARE_ERROR; invalid
        !! controls or nonfinite CSR values/B return LA_INVALID_INPUT_ERROR.
        !! Missing storage returns LA_NOT_INITIALIZED_ERROR; malformed CSR or
        !! MSR pointers/columns return LA_MATRIX_FORMAT_ERROR; singular or
        !! nonfinite factors return LA_SINGULAR_MATRIX_ERROR. If omitted,
        !! errors terminate with ERROR STOP, preserving the existing API.
    real(real64), intent(out), allocatable, dimension(:) :: x
        !! N-element solution on success, or final iterate on nonconvergence.
        !! Unallocated on input errors. Call-local workspace leaves A and LU
        !! unchanged; this does not assert reentrancy of external routines.

    ! Local Variables
    integer :: n, ierr, io, mit, krylov
    real(real64) :: eps
    real(real64), allocatable, dimension(:) :: bc
    real(real64), allocatable, dimension(:,:) :: vv
    
    ! Initialization
    if (present(flag)) flag = LA_NO_ERROR
    n = size(a, 1)
    if (present(im)) then
        krylov = im
    else
        krylov = min(n, 50)
    end if
    if (present(tol)) then
        eps = tol
    else
        eps = sqrt(epsilon(eps))
    end if
    if (present(maxits)) then
        mit = maxits
    else
        mit = 100
    end if
    if (present(iout)) then
        io = iout
    else
        io = 0
    end if

    ! Input Checking
    ierr = LA_NO_ERROR
    if (size(a, 2) /= n) then
        ierr = LA_MATRIX_NOT_SQUARE_ERROR
    else if (size(lu, 1) /= n .or. size(lu, 2) /= n .or. &
        size(b) /= n .or. size(ju) /= n) then
        ierr = LA_MATRIX_DIMENSION_ERROR
    else if (n < 1 .or. .not.ieee_is_finite(eps)) then
        ierr = LA_INVALID_INPUT_ERROR
    else if (eps < epsilon(eps) .or. mit < 1 .or. &
        krylov < 1 .or. krylov > min(n, 50)) then
        ierr = LA_INVALID_INPUT_ERROR
    else if (.not.all(ieee_is_finite(b))) then
        ierr = LA_INVALID_INPUT_ERROR
    end if
    if (ierr == LA_NO_ERROR) ierr = validate_csr_storage(a)
    if (ierr == LA_NO_ERROR) ierr = validate_msr_lu(lu, ju)
    if (ierr /= LA_NO_ERROR) then
        if (present(flag)) then
            flag = ierr
            return
        end if
        error stop ierr
    end if

    ! Memory Allocation
    allocate(x(n), source = 0.0d0)
    if (all(b == 0.0d0)) return
    allocate(vv(n,krylov+1))

    ! Process
    allocate(bc(n), source = b)
    call pgmres(n, krylov, bc, x, vv, eps, mit, io, a%values, a%column_indices, &
        a%row_indices, lu%values, lu%indices, ju, ierr)
    bc = b - matmul(a, x)
    if ((ierr /= 0 .and. ierr /= -1) .or. &
        .not.all(ieee_is_finite(x)) .or. .not.all(ieee_is_finite(bc)) .or. &
        norm2(bc) > eps * norm2(b)) then
        if (present(flag)) then
            flag = LA_CONVERGENCE_ERROR
            return
        end if
        error stop LA_CONVERGENCE_ERROR
    end if
end subroutine

! ------------------------------------------------------------------------------
end module
