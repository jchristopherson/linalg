module linalg_sparse
    use iso_fortran_env, only : int32, real64
    use sparskit
    use blas
    use ferror
    use linalg_errors
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

type :: csr_matrix
    !! A sparse matrix stored in compressed sparse row (CSR) format.
    integer(int32), allocatable, dimension(:) :: row_indices
        !! An M+1 element array containing the indices in V an JA at which the
        !! requested row starts.
    integer(int32), allocatable, dimension(:) :: column_indices
        !! An NNZ-element array, where NNZ is the number of non-zero values,
        !! containing the column indices of each value.
    real(real64), allocatable, dimension(:) :: values
        !! An NNZ-element array, where NNZ is the number of non-zero values,
        !! containing the non-zero values of the matrix.
    integer(int32), private :: n = 0
        !! The number of columns in the matrix.
contains
    procedure, public :: get => csr_get_element
    procedure, public :: extract_diagonal => csr_extract_diagonal
end type

! ------------------------------------------------------------------------------
type :: msr_matrix
    !! A sparse matrix stored in modified sparse row format.
    integer(int32), allocatable, dimension(:) :: indices
        !! An NNZ-element array containing the index information.
    real(real64), allocatable, dimension(:) :: values
        !! An NNZ-element array containing the non-zero values from the
        !! matrix.  The first MIN(M,N) elements contain the diagonal.
    integer(int32) :: m = 0
        !! The number of rows in the matrix.
    integer(int32) :: n = 0
        !! The number of columns in the matrix.
    integer(int32) :: nnz = 0
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
contains
! ******************************************************************************
! CSR ROUTINES
! ------------------------------------------------------------------------------
function csr_get_element(this, i, j) result(rst)
    !! Retrieves the element at the specified row and column.
    class(csr_matrix), intent(in) :: this
        !! The CSR matrix object.
    integer(int32), intent(in) :: i
        !! The row index.
    integer(int32), intent(in) :: j
        !! The column index.
    real(real64) :: rst
        !! The value at the specified row and column.

    ! Local Variables
    integer(int32) :: iadd
    logical :: sorted

    ! Initialization
    sorted = .false.

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
    integer(int32), intent(in) :: dim
        !! The dimension to return the size of.
    integer(int32) :: rst
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
    integer(int32) :: rst
        !! The number of non-zero values in the matrix.

    ! Process
    if (allocated(x%values)) then
        rst = size(x%values)
    else
        rst = 0
    end if
end function

! ------------------------------------------------------------------------------
function create_empty_csr_matrix(m, n, nnz, err) result(rst)
    !! Creates an empty CSR matrix.
    integer(int32), intent(in) :: m
        !! The number of rows in the matrix.
    integer(int32), intent(in) :: n
        !! The number of columns in the matrix.
    integer(int32), intent(in) :: nnz
        !! The number of non-zero values in the matrix.
    class(errors), intent(inout), optional, target :: err
        !! The error object to be updated.
    type(csr_matrix) :: rst
        !! The empty CSR matrix.

    ! Local Variables
    integer(int32) :: flag, m1
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    m1 = m + 1

    ! Input Checking
    if (m < 0) then
        call errmgr%report_error("create_empty_csr_matrix", &
            "The number of rows must be a positive value.", &
            LA_INVALID_INPUT_ERROR)
        return
    end if
    if (n < 0) then
        call errmgr%report_error("create_empty_csr_matrix", &
            "The number of columns must be a positive value.", &
            LA_INVALID_INPUT_ERROR)
        return
    end if
    if (nnz < 0) then
        call errmgr%report_error("create_empty_csr_matrix", &
            "The number of non-zero values must be a positive value.", &
            LA_INVALID_INPUT_ERROR)
        return
    end if

    ! Allocation
    rst%n = n
    allocate(rst%row_indices(m1), rst%column_indices(nnz), source = 0, &
        stat = flag)
    if (flag == 0) allocate(rst%values(nnz), source = 0.0d0, stat = flag)
    if (flag /= 0) then
        call report_memory_error("create_empty_csr_matrix", errmgr, flag)
        return
    end if
end function

! ------------------------------------------------------------------------------
function dense_to_csr(a, err) result(rst)
    !! Converts a dense matrix to a CSR matrix.
    real(real64), intent(in), dimension(:,:) :: a
        !! The dense matrix to convert.
    class(errors), intent(inout), optional, target :: err
        !! The error object to be updated.
    type(csr_matrix) :: rst
        !! The CSR matrix.

    ! Local Variables
    integer(int32) :: i, j, k, m, n, nnz
    real(real64) :: t
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
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
    rst = create_empty_csr_matrix(m, n, nnz, errmgr)
    if (errmgr%has_error_occurred()) return

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
function banded_to_csr(m, ml, mu, a, err) result(rst)
    !! Converts a banded matrix to a CSR matrix.
    integer(int32), intent(in) :: m
        !! The number of rows in the banded matrix.
    integer(int32), intent(in) :: ml
        !! The number of lower diagonals in the banded matrix.
    integer(int32), intent(in) :: mu
        !! The number of upper diagonals in the banded matrix.
    real(real64), intent(in), dimension(:,:) :: a
        !! The banded matrix to convert.
    class(errors), intent(inout), optional, target :: err
        !! The error object to be updated.
    type(csr_matrix) :: rst
        !! The CSR matrix.

    ! Local Variables
    integer(int32) :: n, nnz, flag, lowd, lda
    integer(int32), allocatable, dimension(:) :: ia, ja
    real(real64), allocatable, dimension(:) :: v
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    lda = size(a, 1)
    n = size(a, 2)
    nnz = lda * n
    lowd = ml + mu + 1

    ! Input Checking
    if (ml < 0 .or. mu < 0) then
        call errmgr%report_error("banded_to_csr", "The bandwidth " // &
            "dimensions cannot be negative.", LA_INVALID_INPUT_ERROR)
        return
    end if
    if (lda /= ml + mu + 1) then
        call errmgr%report_error("banded_to_csr", "The number of rows in " // &
            "the banded matrix does not match the supplied bandwidth " // &
            "dimensions.", LA_MATRIX_FORMAT_ERROR)
        return
    end if

    ! Allocation
    allocate(ia(m + 1), ja(nnz), v(nnz), stat = flag)
    if (flag /= 0) go to 10

    ! Process
    call bndcsr(m, a, lda, lowd, ml, mu, v, ja, ia, nnz, flag)
    nnz = ia(m + 1) - 1

    ! Put into the sparse matrix structure
    allocate(rst%row_indices(m + 1), source = ia, stat = flag)
    if (flag == 0) allocate(rst%column_indices(nnz), source = ja(:nnz), &
        stat = flag)
    if (flag == 0) allocate(rst%values(nnz), source = v(:nnz), stat = flag)
    if (flag /= 0) go to 10
    rst%n = n

    ! End
    return

    ! Memory Error
10  continue
    call report_memory_error("banded_to_csr", errmgr, flag)
    return
end function

! ------------------------------------------------------------------------------
subroutine csr_to_dense(a, x, err)
    !! Converts a CSR matrix to a dense matrix.
    class(csr_matrix), intent(in) :: a
        !! The CSR matrix to convert.
    real(real64), intent(out), dimension(:,:) :: x
        !! The dense matrix.
    class(errors), intent(inout), optional, target :: err
        !! The error object to be updated.

    ! Local Variables
    integer(int32) :: i, j, k, m, n, nnz, flag
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    m = size(a, 1)
    n = size(a, 2)
    nnz = nonzero_count(a)
    
    ! Input Check
    if (size(x, 1) /= m .or. size(x, 2) /= n) then
        call report_matrix_size_error("csr_to_dense", errmgr, "x", m, n, &
            size(x, 1), size(x, 2))
        return
    end if

    ! Process
    do i = 1, m
        x(i,:) = 0.0d0
        do k = a%row_indices(i), a%row_indices(i+1) - 1
            j = a%column_indices(k)
            x(i,j) = a%values(k)
        end do
    end do
end subroutine

! ------------------------------------------------------------------------------
function diag_to_csr(a, err) result(rst)
    !! Converts a diagonal matrix to a CSR matrix.
    real(real64), intent(in), dimension(:) :: a
        !! The diagonal matrix to convert.
    class(errors), intent(inout), optional, target :: err
        !! The error object to be updated.
    type(csr_matrix) :: rst
        !! The CSR matrix.

    ! Local Variables
    integer(int32) :: i, n, n1, flag
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    n = size(a)
    n1 = n + 1

    ! Allocation
    allocate(rst%row_indices(n1), rst%column_indices(n), stat = flag)
    if (flag == 0) allocate(rst%values(n), source = a, stat = flag)
    if (flag /= 0) then
        call report_memory_error("diag_to_csr", errmgr, flag)
        return
    end if
    rst%n = n

    ! Populate IA & JA
    do i = 1, n
        rst%column_indices(i) = i
        rst%row_indices(i) = i
    end do
    rst%row_indices(n1) = n1
end function

! ------------------------------------------------------------------------------
subroutine csr_assign_to_dense(dense, sparse)
    !! Assigns the values of a CSR matrix to a dense matrix.
    real(real64), intent(out), dimension(:,:) :: dense
        !! The dense matrix.
    class(csr_matrix), intent(in) :: sparse
        !! The CSR matrix.

    ! Process
    call csr_to_dense(sparse, dense)
end subroutine

! ------------------------------------------------------------------------------
subroutine dense_assign_to_csr(sparse, dense)
    !! Assigns the values of a dense matrix to a CSR matrix.
    type(csr_matrix), intent(out) :: sparse
        !! The CSR matrix.
    real(real64), intent(in), dimension(:,:) :: dense
        !! The dense matrix.

    ! Process
    sparse = dense_to_csr(dense)
end subroutine

! ------------------------------------------------------------------------------
function csr_mtx_mtx_mult(a, b) result(rst)
    !! Multiplies two CSR matrices together.
    class(csr_matrix), intent(in) :: a
        !! The first CSR matrix.
    class(csr_matrix), intent(in) :: b
        !! The second CSR matrix.
    type(csr_matrix) :: rst
        !! The resulting CSR matrix.

    ! Local Variables
    integer(int32), parameter :: sym_mult = 0
    integer(int32), parameter :: full_mult = 1
    integer(int32) :: flag, m, n, k, nnza, nnzb, nnzc, ierr
    integer(int32), allocatable, dimension(:) :: ic, jc, iw
    real(real64) :: dummy(1)
    type(errors) :: errmgr
    
    ! Initialization
    m = size(a, 1)
    n = size(b, 2)
    k = size(a, 2)
    nnza = nonzero_count(a)
    nnzb = nonzero_count(b)
    nnzc = nnza + nnzb

    ! Input Check
    if (size(b, 1) /= k) then
        call report_inner_matrix_dimension_error("csr_mtx_mtx_mult", errmgr, &
            "a", "b", k, size(b, 1))
        return
    end if

    ! Local Memory Allocations
    allocate(ic(m + 1), jc(nnzc), iw(n), stat = flag)
    if (flag /= 0) go to 10

    ! Determine the structure of C
    call amub(m, n, sym_mult, a%values, a%column_indices, a%row_indices, &
        b%values, b%column_indices, b%row_indices, dummy, jc, ic, &
        nnzc, iw, ierr)
    if (ierr /= 0) then
        ! NNZC was too small - try increasing it
        do while (ierr /= 0)
            deallocate(jc)
            nnzc = nnzc + nnza + nnzb
            allocate(jc(nnzc), stat = flag)
            if (flag /= 0) go to 10
            call amub(m, n, sym_mult, a%values, a%column_indices, &
                a%row_indices, b%values, b%column_indices, b%row_indices, &
                dummy, jc, ic, nnzc, iw, ierr)
        end do
    end if

    ! Determine the actual NNZ for C & allocate space for the output
    nnzc = ic(m + 1) - 1
    deallocate(ic)
    deallocate(jc)
    rst = create_empty_csr_matrix(m, n, nnzc, errmgr)
    if (errmgr%has_error_occurred()) return

    ! Compute the actual product
    call amub(m, n, full_mult, a%values, a%column_indices, a%row_indices, &
        b%values, b%column_indices, b%row_indices, rst%values, &
        rst%column_indices, rst%row_indices, nnzc, iw, ierr)

    ! End
    return

    ! Memory Error
10  continue
    call report_memory_error("csr_mtx_mtx_mult", errmgr, flag)
    return
end function

! ------------------------------------------------------------------------------
function csr_mtx_vec_mult(a, b) result(rst)
    !! Multiplies a CSR matrix by a vector.
    class(csr_matrix), intent(in) :: a
        !! The CSR matrix.
    real(real64), intent(in), dimension(:) :: b
        !! The vector.
    real(real64), allocatable, dimension(:) :: rst
        !! The resulting vector.

    ! Local Variables
    integer(int32) :: i, k, k1, k2, n, p, flag
    real(real64) :: t
    type(errors) :: errmgr

    ! Initialization
    n = size(a, 1)
    p = size(a, 2)

    ! Input Check
    if (size(b) /= p) then
        call report_inner_matrix_dimension_error("csr_mtx_vec_mult", errmgr, &
            "a", "b", p, size(b))
        return
    end if

    ! Memory Allocation
    allocate(rst(n), stat = flag)
    if (flag /= 0) then
        call report_memory_error("csr_mtx_vec_mult", errmgr, flag)
        return
    end if

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
function csr_mtx_add(a, b) result(rst)
    !! Adds two CSR matrices.
    class(csr_matrix), intent(in) :: a
        !! The first CSR matrix.
    class(csr_matrix), intent(in) :: b
        !! The second CSR matrix.
    type(csr_matrix) :: rst
        !! The resulting CSR matrix.

    ! Local Variables
    integer(int32), parameter :: sym_add = 0
    integer(int32), parameter :: full_add = 1
    integer(int32) :: m, n, nnza, nnzb, nnzc, ierr, flag
    integer(int32), allocatable, dimension(:) :: ic, jc, iw
    real(real64) :: dummy(1)
    type(errors) :: errmgr
    
    ! Initialization
    m = size(a, 1)
    n = size(a, 2)
    nnza = nonzero_count(a)
    nnzb = nonzero_count(b)
    nnzc = nnza + nnzb

    ! Input Checking
    if (size(b, 1) /= m .or. size(b, 2) /= n) then
        call report_matrix_size_error("csr_mtx_add", errmgr, "b", m, n, &
            size(b, 1), size(b, 2))
        return
    end if

    ! Local Memory Allocations
    allocate(ic(m + 1), jc(nnzc), iw(n), stat = flag)
    if (flag /= 0) go to 10

    ! Determine the structure of C
    call aplb(m, n, sym_add, a%values, a%column_indices, a%row_indices, &
        b%values, b%column_indices, b%row_indices, dummy, jc, ic, &
        nnzc, iw, ierr)
    if (ierr /= 0) then
        ! NNZC was too small - try increasing it
        do while (ierr /= 0)
            deallocate(jc)
            nnzc = nnzc + nnza + nnzb
            allocate(jc(nnzc), stat = flag)
            if (flag /= 0) go to 10
            call aplb(m, n, sym_add, a%values, a%column_indices, &
                a%row_indices, b%values, b%column_indices, b%row_indices, &
                dummy, jc, ic, nnzc, iw, ierr)
        end do
    end if

    ! Determine the actuall NNZ for C & allocate space for the output
    nnzc = ic(m + 1) - 1
    deallocate(ic)
    deallocate(jc)
    rst = create_empty_csr_matrix(m, n, nnzc, errmgr)
    if (errmgr%has_error_occurred()) return

    ! Compute the actual sum
    call aplb(m, n, full_add, a%values, a%column_indices, a%row_indices, &
        b%values, b%column_indices, b%row_indices, rst%values, &
        rst%column_indices, rst%row_indices, nnzc, iw, ierr)

    ! End
    return

    ! Memory Error
10  continue
    call report_memory_error("csr_mtx_add", errmgr, flag)
    return
end function

! ------------------------------------------------------------------------------
function csr_mtx_sub(a, b) result(rst)
    !! Subtracts two CSR matrices.
    class(csr_matrix), intent(in) :: a
        !! The first CSR matrix.
    class(csr_matrix), intent(in) :: b
        !! The second CSR matrix.
    type(csr_matrix) :: rst
        !! The resulting CSR matrix.

    ! Local Variables
    integer(int32), parameter :: sym_add = 0
    integer(int32) :: m, n, nnza, nnzb, nnzc, ierr, flag
    integer(int32), allocatable, dimension(:) :: ic, jc, iw
    real(real64) :: dummy(1)
    type(errors) :: errmgr
    
    ! Initialization
    m = size(a, 1)
    n = size(a, 2)
    nnza = nonzero_count(a)
    nnzb = nonzero_count(b)
    nnzc = nnza + nnzb

    ! Input Checking
    if (size(b, 1) /= m .or. size(b, 2) /= n) then
        call report_matrix_size_error("csr_mtx_sub", errmgr, "b", m, n, &
            size(b, 1), size(b, 2))
        return
    end if

    ! Local Memory Allocations
    allocate(ic(m + 1), jc(nnzc), iw(n), stat = flag)
    if (flag /= 0) go to 10

    ! Determine the structure of C
    call aplb(m, n, sym_add, a%values, a%column_indices, a%row_indices, &
        b%values, b%column_indices, b%row_indices, dummy, jc, ic, &
        nnzc, iw, ierr)
    if (ierr /= 0) then
        ! NNZC was too small - try increasing it
        do while (ierr /= 0)
            deallocate(jc)
            nnzc = nnzc + nnza + nnzb
            allocate(jc(nnzc), stat = flag)
            if (flag /= 0) go to 10
            call aplb(m, n, sym_add, a%values, a%column_indices, &
                a%row_indices, b%values, b%column_indices, b%row_indices, &
                dummy, jc, ic, nnzc, iw, ierr)
        end do
    end if

    ! Determine the actuall NNZ for C & allocate space for the output
    nnzc = ic(m + 1) - 1
    deallocate(ic)
    deallocate(jc)
    rst = create_empty_csr_matrix(m, n, nnzc, errmgr)
    if (errmgr%has_error_occurred()) return

    ! Compute the actual sum
    call aplsb(m, n, a%values, a%column_indices, a%row_indices, -1.0d0, &
        b%values, b%column_indices, b%row_indices, rst%values, &
        rst%column_indices, rst%row_indices, nnzc, iw, ierr)

    ! End
    return

    ! Memory Error
10  continue
    call report_memory_error("csr_mtx_sub", errmgr, flag)
    return
end function

! ------------------------------------------------------------------------------
function csr_mtx_mult_scalar_1(a, b) result(rst)
    !! Multiplies a CSR matrix by a scalar.
    class(csr_matrix), intent(in) :: a
        !! The CSR matrix.
    real(real64), intent(in) :: b
        !! The scalar.
    type(csr_matrix) :: rst
        !! The resulting CSR matrix.

    ! Local Variables
    integer(int32) :: m, n, nnz
    type(errors) :: errmgr
    
    ! Initialization
    m = size(a, 1)
    n = size(a, 2)
    nnz = nonzero_count(a)

    ! Process
    rst = create_empty_csr_matrix(m, n, nnz, errmgr)
    if (errmgr%has_error_occurred()) return

    ! Compute the product
    rst%row_indices = a%row_indices
    rst%column_indices = a%column_indices
    rst%values = b * a%values
end function

! ------------------------------------------------------------------------------
function csr_mtx_mult_scalar_2(a, b) result(rst)
    !! Multiplies a scalar by a CSR matrix.
    real(real64), intent(in) :: a
        !! The scalar.
    class(csr_matrix), intent(in) :: b
        !! The CSR matrix.
    type(csr_matrix) :: rst
        !! The resulting CSR matrix.

    ! Local Variables
    integer(int32) :: m, n, nnz
    type(errors) :: errmgr
    
    ! Initialization
    m = size(b, 1)
    n = size(b, 2)
    nnz = nonzero_count(b)

    ! Process
    rst = create_empty_csr_matrix(m, n, nnz, errmgr)
    if (errmgr%has_error_occurred()) return

    ! Compute the product
    rst%row_indices = b%row_indices
    rst%column_indices = b%column_indices
    rst%values = a * b%values
end function

! ------------------------------------------------------------------------------
function csr_mtx_divide_scalar_1(a, b) result(rst)
    !! Divides a CSR matrix by a scalar.
    class(csr_matrix), intent(in) :: a
        !! The CSR matrix.
    real(real64), intent(in) :: b
        !! The scalar.
    type(csr_matrix) :: rst

    ! Local Variables
    integer(int32) :: m, n, nnz
    type(errors) :: errmgr
    
    ! Initialization
    m = size(a, 1)
    n = size(a, 2)
    nnz = nonzero_count(a)

    ! Process
    rst = create_empty_csr_matrix(m, n, nnz, errmgr)
    if (errmgr%has_error_occurred()) return

    ! Compute the product
    rst%row_indices = a%row_indices
    rst%column_indices = a%column_indices
    rst%values = a%values / b
end function

! ------------------------------------------------------------------------------
function csr_transpose(a) result(rst)
    !! Transposes a CSR matrix.
    class(csr_matrix), intent(in) :: a
        !! The CSR matrix.
    type(csr_matrix) :: rst
        !! The transposed CSR matrix.

    ! Local Variables
    integer(int32) :: m, n, nnz
    type(errors) :: errmgr

    ! Initialization
    m = size(a, 1)
    n = size(a, 2)
    nnz = nonzero_count(a)
    rst = create_empty_csr_matrix(n, m, nnz, errmgr)
    if (errmgr%has_error_occurred()) return

    ! Process
    call csrcsc2(m, n, 1, 1, a%values, a%column_indices, a%row_indices, &
        rst%values, rst%column_indices, rst%row_indices)
end function

! ------------------------------------------------------------------------------
subroutine csr_extract_diagonal(a, diag, err)
    !! Extracts the diagonal from a CSR matrix.
    class(csr_matrix), intent(in) :: a
        !! The CSR matrix.
    real(real64), intent(out), dimension(:) :: diag
        !! The diagonal values.
    class(errors), intent(inout), optional, target :: err
        !! The error object to be updated.

    ! Local Variables
    integer(int32) :: m, n, mn, len, flag
    integer(int32), allocatable, dimension(:) :: idiag
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    m = size(a, 1)
    n = size(a, 2)
    mn = min(m, n)

    ! Input Check
    if (size(diag) /= mn) then
        call report_array_size_error("csr_extract_diagonal", errmgr, "diag", &
            mn, size(diag))
        return
    end if

    ! Memory Allocation
    allocate(idiag(mn), stat = flag)
    if (flag /= 0) then
        call report_memory_error("csr_extract_diagonal", errmgr, flag)
        return
    end if

    ! Process
    call getdia(m, n, 0, a%values, a%column_indices, a%row_indices, len, &
        diag, idiag, 0)
end subroutine

! ------------------------------------------------------------------------------
subroutine csr_solve_sparse_direct(a, b, x, droptol, err)
    !! Solves a linear system using a direct method.
    class(csr_matrix), intent(in) :: a
        !! The matrix.
    real(real64), intent(in), dimension(:) :: b
        !! The right-hand side.
    real(real64), intent(out), dimension(:) :: x
        !! The solution.
    real(real64), intent(in), optional :: droptol
        !! The drop tolerance for the ILUT factorization.
    class(errors), intent(inout), optional, target :: err
        !! The error object to be updated.

    ! Local Variables
    integer(int32) :: i, m, n, nnz, lfil, iwk, ierr, flag
    integer(int32), allocatable, dimension(:) :: jlu, ju, jw
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

    ! Input Checking
    if (m /= n) then
        call report_square_matrix_error("csr_solve_sparse_direct", errmgr, &
            "a", m, m, n)
        return
    end if
    if (size(x) /= n) then
        call report_inner_matrix_dimension_error("csr_solve_sparse_direct", &
            errmgr, "a", "x", n, size(x))
        return
    end if
    if (size(b) /= n) then
        call report_array_size_error("csr_solve_sparse_direct", errmgr, "b", &
            n, size(b))
        return
    end if

    ! Parameter Determination
    lfil = 1
    do i = 1, m
        lfil = max(lfil, a%row_indices(i+1) - a%row_indices(i))
    end do
    iwk = max(lfil * m, nnz)  ! somewhat arbitrary - can be adjusted

    ! Local Memory Allocation
    allocate(alu(iwk), w(n+1), jlu(iwk), ju(n), jw(2 * n), stat = flag)
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

    ! Solution
    call lusol(n, b, x, alu, jlu, ju)

    ! End
    return

    ! Memory Error
10  continue
    call report_memory_error("csr_solve_sparse_direct", errmgr, flag)
    return

    ! Matrix Format Error
20  continue
    call errmgr%report_error("csr_solve_sparse_direct", &
        "The input matrix was incorrectly formatted.  A row with more " // &
        "than N entries was found.", LA_MATRIX_FORMAT_ERROR)
    return

    ! Zero Row Error
30  continue
    call errmgr%report_error("csr_solve_sparse_direct", &
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
    call errmgr%report_error("csr_solve_sparse_direct", &
        "A zero pivot was encountered.", LA_SINGULAR_MATRIX_ERROR)
    return
end subroutine

! ******************************************************************************
! MSR ROUTINES
! ------------------------------------------------------------------------------
! TO DO: MSR_GET_ELEMENT

! ------------------------------------------------------------------------------
pure function msr_size(x, dim) result(rst)
    !! Returns the size of the specified dimension of an MSR matrix.
    class(msr_matrix), intent(in) :: x
        !! The MSR matrix.
    integer(int32), intent(in) :: dim
        !! The dimension to return the size of.
    integer(int32) :: rst
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
    integer(int32) :: rst
        !! The number of non-zero elements.

    ! Process
    rst = x%nnz
end function

! ------------------------------------------------------------------------------
function create_empty_msr_matrix(m, n, nnz, err) result(rst)
    !! Creates an empty MSR matrix.
    integer(int32), intent(in) :: m
        !! The number of rows in the matrix.
    integer(int32), intent(in) :: n
        !! The number of columns in the matrix.
    integer(int32), intent(in) :: nnz
        !! The number of non-zero elements in the matrix.
    class(errors), intent(inout), optional, target :: err
        !! The error object to be updated.
    type(msr_matrix) :: rst
        !! The MSR matrix.

    ! Local Variables
    integer(int32) :: nelem, mn, flag
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if

    ! Input Checking
    if (m < 0) then
        call errmgr%report_error("create_empty_msr_matrix", &
            "The number of rows must be a positive value.", &
            LA_INVALID_INPUT_ERROR)
        return
    end if
    if (n < 0) then
        call errmgr%report_error("create_empty_msr_matrix", &
            "The number of columns must be a positive value.", &
            LA_INVALID_INPUT_ERROR)
        return
    end if
    if (nnz < 0) then
        call errmgr%report_error("create_empty_msr_matrix", &
            "The number of non-zero values must be a positive value.", &
            LA_INVALID_INPUT_ERROR)
        return
    end if

    ! Allocation
    rst%m = m
    rst%n = n
    rst%nnz = nnz
    mn = min(m, n)
    nelem = m + 1 + nnz - mn
    allocate(rst%indices(nelem), source = 0, stat = flag)
    if (flag == 0) allocate(rst%values(nelem), source = 0.0d0, stat = flag)
    if (flag /= 0) then
        call report_memory_error("create_empty_msr_matrix", errmgr, flag)
        return
    end if
end function

! ------------------------------------------------------------------------------
function csr_to_msr(a, err) result(rst)
    !! Converts a CSR matrix to an MSR matrix.
    class(csr_matrix), intent(in) :: a
        !! The CSR matrix to convert.
    class(errors), intent(inout), optional, target :: err
        !! The error object to be updated.
    type(msr_matrix) :: rst
        !! The MSR matrix.

    ! Local Variables
    integer(int32) :: m, n, nnz, flag
    integer(int32), allocatable, dimension(:) :: iwork, jc, ic
    real(real64), allocatable, dimension(:) :: work, ac
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    m = size(a, 1)
    n = size(a, 2)
    nnz = nonzero_count(a)

    ! Memory Allocation
    rst = create_empty_msr_matrix(m, n, nnz, errmgr)
    if (errmgr%has_error_occurred()) return
    allocate(work(m), iwork(m + 1), stat = flag)
    if (flag == 0) allocate(ac(nnz), source = a%values, stat = flag)
    if (flag == 0) allocate(jc(nnz), source = a%column_indices, stat = flag)
    if (flag == 0) allocate(ic(m+1), source = a%row_indices, stat = flag)
    if (flag /= 0) then
        call report_memory_error("csr_to_msr", errmgr, flag)
        return
    end if

    ! Perform the conversion
    call csrmsr(m, ac, jc, ic, rst%values, rst%indices, work, iwork)
end function

! ------------------------------------------------------------------------------
function msr_to_csr(a, err) result(rst)
    !! Converts an MSR matrix to a CSR matrix.
    class(msr_matrix), intent(in) :: a
        !! The MSR matrix to convert.
    class(errors), intent(inout), optional, target :: err
        !! The error object to be updated.
    type(csr_matrix) :: rst
        !! The CSR matrix.

    ! Local Variables
    integer(int32) :: m, n, nnz, flag
    integer(int32), allocatable, dimension(:) :: iwork
    real(real64), allocatable, dimension(:) :: work
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    m = size(a, 1)
    n = size(a, 2)
    nnz = nonzero_count(a)

    ! Memory Allocation
    rst = create_empty_csr_matrix(m, n, nnz, errmgr)
    if (errmgr%has_error_occurred()) return
    allocate(work(m), iwork(m+1), stat = flag)
    if (flag /= 0) then
        call report_memory_error("msr_to_csr", errmgr, flag)
        return
    end if

    ! Process
    call msrcsr(m, a%values, a%indices, rst%values, rst%column_indices, &
        rst%row_indices, work, iwork)
end function

! ------------------------------------------------------------------------------
function dense_to_msr(a, err) result(rst)
    !! Converts a dense matrix to an MSR matrix.
    real(real64), intent(in), dimension(:,:) :: a
        !! The dense matrix to convert.
    class(errors), intent(inout), optional, target :: err
        !! The error object to be updated.
    type(msr_matrix) :: rst
        !! The MSR matrix.

    ! Local Variables
    type(csr_matrix) :: csr
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if

    ! Convert to CSR, and then from CSR to MSR
    csr = dense_to_csr(a, errmgr)

    ! Convert to MSR
    rst = csr_to_msr(csr, errmgr)
end function

! ------------------------------------------------------------------------------
subroutine msr_to_dense(a, x, err)
    !! Converts an MSR matrix to a dense matrix.
    class(msr_matrix), intent(in) :: a
        !! The MSR matrix to convert.
    real(real64), intent(out), dimension(:,:) :: x
        !! The dense matrix.
    class(errors), intent(inout), optional, target :: err
        !! The error object to be updated.

    ! Local Variables
    integer(int32) :: m, n, flag
    type(csr_matrix) :: csr
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    m = size(a, 1)
    n = size(a, 2)

    ! Input Check
    if (size(x, 1) /= m .or. size(x, 2) /= n) then
        call report_matrix_size_error("msr_to_dense", errmgr, "x", m, n, &
            size(x, 1), size(x, 2))
        return
    end if

    ! Process
    csr = msr_to_csr(a, errmgr)
    if (errmgr%has_error_occurred()) return
    call csr_to_dense(csr, x, errmgr)
end subroutine

! ------------------------------------------------------------------------------
subroutine msr_assign_to_dense(dense, msr)
    !! Assigns an MSR matrix to a dense matrix.
    real(real64), intent(out), dimension(:,:) :: dense
        !! The dense matrix.
    class(msr_matrix), intent(in) :: msr
        !! The MSR matrix.

    ! Process
    call msr_to_dense(msr, dense)
end subroutine

! ------------------------------------------------------------------------------
subroutine dense_assign_to_msr(msr, dense)
    !! Assigns a dense matrix to an MSR matrix.
    type(msr_matrix), intent(out) :: msr
        !! The MSR matrix.
    real(real64), intent(in), dimension(:,:) :: dense
        !! The dense matrix.

    ! Process
    msr = dense_to_msr(dense)
end subroutine

! ------------------------------------------------------------------------------
subroutine csr_assign_to_msr(msr, csr)
    !! Assigns a CSR matrix to an MSR matrix.
    type(msr_matrix), intent(out) :: msr
        !! The MSR matrix.
    class(csr_matrix), intent(in) :: csr
        !! The CSR matrix.

    ! Process
    msr = csr_to_msr(csr)
end subroutine

! ------------------------------------------------------------------------------
subroutine msr_assign_to_csr(csr, msr)
    !! Assigns an MSR matrix to a CSR matrix.
    type(csr_matrix), intent(out) :: csr
        !! The CSR matrix.
    class(msr_matrix), intent(in) :: msr
        !! The MSR matrix.

    ! Process
    csr = msr_to_csr(msr)
end subroutine

! ------------------------------------------------------------------------------
function create_csr_matrix(m, n, rows, cols, vals, err) result(rst)
    !! Creates a CSR matrix from the input data.
    integer(int32), intent(in) :: m
        !! The number of rows in the matrix.
    integer(int32), intent(in) :: n
        !! The number of columns in the matrix.
    integer(int32), intent(in), dimension(:) :: rows
        !! The row indices.
    integer(int32), intent(in), dimension(:) :: cols
        !! The column indices.
    real(real64), intent(in), dimension(:) :: vals
        !! The values.
    class(errors), intent(inout), optional, target :: err
        !! The error object to be updated.
    type(csr_matrix) :: rst
        !! The CSR matrix.

    ! Local Variables
    integer(int32) :: i, flag, nnz
    integer(int32), allocatable, dimension(:) :: ir
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    nnz = size(rows)

    ! Input Checking
    if (m < 0) then
        call errmgr%report_error("create_csr_matrix", &
            "The number of rows must be a positive value.", &
            LA_INVALID_INPUT_ERROR)
        return
    end if
    if (n < 0) then
        call errmgr%report_error("create_csr_matrix", &
            "The number of columns must be a positive value.", &
            LA_INVALID_INPUT_ERROR)
        return
    end if
    if (size(cols) /= nnz .or. size(vals) /= nnz) then
        call errmgr%report_error("create_csr_matrix", &
            "The size of the input arrays must be the same.", &
            LA_ARRAY_SIZE_ERROR)
        return
    end if
    do i = 1, nnz
        if (rows(i) < 1 .or. rows(i) > m) then
            call errmgr%report_error("create_csr_matrix", &
                "All row indices must be within the bounds of the matrix.", &
                LA_INVALID_INPUT_ERROR)
            return
        end if
        if (cols(i) < 1 .or. cols(i) > n) then
            call errmgr%report_error("create_csr_matrix", &
                "All column indices must be within the bounds of the matrix.", &
                LA_INVALID_INPUT_ERROR)
            return
        end if
    end do
    allocate(ir(nnz), source = rows, stat = flag)
    if (flag /= 0) then
        call report_memory_error("create_csr_matrix", errmgr, flag)
        return
    end if

    ! Create an empty matrix
    rst = create_empty_csr_matrix(m, n, nnz, errmgr)
    if (errmgr%has_error_occurred()) return

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
subroutine csr_pgmres_solver(a, lu, ju, b, x, im, tol, maxits, iout, err)
    !! Solves a linear system using the PGMRES method.
    class(csr_matrix), intent(in) :: a
        !! The matrix.
    class(msr_matrix), intent(in) :: lu
        !! The LU factored matrix.
    integer(int32), intent(in), dimension(:) :: ju
        !! The row tracking array.
    real(real64), intent(inout), dimension(:) :: b
        !! The right-hand side.
    real(real64), intent(out), dimension(:) :: x
        !! The solution.
    integer(int32), intent(in), optional :: im
        !! The Krylov subspace size.
    integer(int32), intent(in), optional :: maxits
        !! The maximum number of iterations.
    integer(int32), intent(in), optional :: iout
        !! The output level.
    real(real64), intent(in), optional :: tol
        !! The convergence tolerance.
    class(errors), intent(inout), optional, target :: err
        !! The error object to be updated.

    ! Local Variables
    integer(int32) :: n, ierr, flag, io, mit, krylov
    real(real64) :: eps
    real(real64), allocatable, dimension(:,:) :: vv
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    n = size(a, 1)
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
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
    if (size(a, 2) /= n) then
        call report_square_matrix_error("csr_pgmres_solver", errmgr, "a", n, n, &
            size(a, 2))
        return
    end if
    if (size(lu, 1) /= n .or. size(lu, 2) /= n) then
        call report_matrix_size_error("csr_pgmres_solver", errmgr, "lu", n, n, &
            size(lu, 1), size(lu, 2))
        return
    end if
    if (size(b) /= n) then
        call report_array_size_error("csr_pgmres_solver", errmgr, "b", n, size(b))
        return
    end if
    if (size(x) /= n) then
        call report_array_size_error("csr_pgmres_solver", errmgr, "x", n, size(x))
        return
    end if
    if (eps < epsilon(eps)) then
        call errmgr%report_error("csr_pgmres_solver", &
            "The convergence tolerance is too small.", LA_INVALID_INPUT_ERROR)
        return
    end if
    if (mit < 1) then
        call errmgr%report_error("csr_pgmres_solver", &
            "Too few iterations allowed.", LA_INVALID_INPUT_ERROR)
        return
    end if
    if (krylov < 1) then
        call errmgr%report_error("csr_pgmres_solver", &
            "The requested Krylov subspace size is too small.", &
            LA_INVALID_INPUT_ERROR)
        return
    end if

    ! Memory Allocation
    allocate(vv(n,krylov+1), stat = flag)
    if (flag /= 0) then
        call report_memory_error("csr_pgmres_solver", errmgr, flag)
        return
    end if

    ! Process
    call pgmres(n, krylov, b, x, vv, eps, mit, io, a%values, a%column_indices, &
        a%row_indices, lu%values, lu%indices, ju, ierr)
    if (ierr == 1) then
        call errmgr%report_error("csr_pgmres_solver", &
            "Convergence could not be achieved to the requested tolerance " // &
            "in the allowed number of iterations.", LA_CONVERGENCE_ERROR)
        return
    end if
end subroutine

! ------------------------------------------------------------------------------
end module