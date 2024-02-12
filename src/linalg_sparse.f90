submodule (linalg) linalg_sparse
    use sparskit
    use blas
    implicit none
contains
! ******************************************************************************
! CSR ROUTINES
! ------------------------------------------------------------------------------
module function csr_get_element(this, i, j) result(rst)
    ! Arguments
    class(csr_matrix), intent(in) :: this
    integer(int32), intent(in) :: i, j
    real(real64) :: rst

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
pure module function csr_size(x, dim) result(rst)
    ! Arguments
    class(csr_matrix), intent(in) :: x
    integer(int32), intent(in) :: dim
    integer(int32) :: rst

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
pure module function nonzero_count_csr(x) result(rst)
    ! Arguments
    class(csr_matrix), intent(in) :: x
    integer(int32) :: rst

    ! Process
    if (allocated(x%values)) then
        rst = size(x%values)
    else
        rst = 0
    end if
end function

! ------------------------------------------------------------------------------
module function create_empty_csr_matrix(m, n, nnz, err) result(rst)
    ! Arguments
    integer(int32), intent(in) :: m, n, nnz
    class(errors), intent(inout), optional, target :: err
    type(csr_matrix) :: rst

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
        call errmgr%report_error("create_empty_csr_matrix", &
            "Memory allocation error.", LA_OUT_OF_MEMORY_ERROR)
        return
    end if
end function

! ------------------------------------------------------------------------------
module function dense_to_csr(a, err) result(rst)
    ! Arguments
    real(real64), intent(in), dimension(:,:) :: a
    class(errors), intent(inout), optional, target :: err
    type(csr_matrix) :: rst

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
module function banded_to_csr(m, ml, mu, a, err) result(rst)
    ! Arguments
    integer(int32), intent(in) :: m, ml, mu
    real(real64), intent(in), dimension(:,:) :: a
    class(errors), intent(inout), optional, target :: err
    type(csr_matrix) :: rst

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
    call errmgr%report_error("banded_to_csr", "Memory allocation error.", &
        LA_OUT_OF_MEMORY_ERROR)
    return
end function

! ------------------------------------------------------------------------------
module subroutine csr_to_dense(a, x, err)
    ! Arguments
    class(csr_matrix), intent(in) :: a
    real(real64), intent(out), dimension(:,:) :: x
    class(errors), intent(inout), optional, target :: err

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
        call errmgr%report_error("csr_to_dense", &
            "The output matrix dimensions are not correct.", &
            LA_ARRAY_SIZE_ERROR)
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
module function diag_to_csr(a, err) result(rst)
    ! Arguments
    real(real64), intent(in), dimension(:) :: a
    class(errors), intent(inout), optional, target :: err
    type(csr_matrix) :: rst

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
        call errmgr%report_error("diag_to_csr", "Memory allocation error.", &
            LA_OUT_OF_MEMORY_ERROR)
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
module subroutine csr_assign_to_dense(dense, sparse)
    ! Arguments
    real(real64), intent(out), dimension(:,:) :: dense
    class(csr_matrix), intent(in) :: sparse

    ! Process
    call csr_to_dense(sparse, dense)
end subroutine

! ------------------------------------------------------------------------------
module subroutine dense_assign_to_csr(sparse, dense)
    ! Arguments
    type(csr_matrix), intent(out) :: sparse
    real(real64), intent(in), dimension(:,:) :: dense

    ! Process
    sparse = dense_to_csr(dense)
end subroutine

! ------------------------------------------------------------------------------
module function csr_mtx_mtx_mult(a, b) result(rst)
    ! Arguments
    class(csr_matrix), intent(in) :: a, b
    type(csr_matrix) :: rst

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
        call errmgr%report_error("csr_mtx_mtx_mult", &
            "Inner matrix dimension mismatch.", LA_ARRAY_SIZE_ERROR)
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
    call errmgr%report_error("csr_mtx_mtx_mult", &
        "Memory allocation error.", LA_OUT_OF_MEMORY_ERROR)
    return
end function

! ------------------------------------------------------------------------------
module function csr_mtx_vec_mult(a, b) result(rst)
    ! Arguments
    class(csr_matrix), intent(in) :: a
    real(real64), intent(in), dimension(:) :: b
    real(real64), allocatable, dimension(:) :: rst

    ! Local Variables
    integer(int32) :: i, k, k1, k2, n, p, flag
    real(real64) :: t
    type(errors) :: errmgr

    ! Initialization
    n = size(a, 1)
    p = size(a, 2)

    ! Input Check
    if (size(b) /= p) then
        call errmgr%report_error("csr_mtx_vec_mult", &
            "Inner matrix dimension error.", LA_ARRAY_SIZE_ERROR)
        return
    end if

    ! Memory Allocation
    allocate(rst(n), stat = flag)
    if (flag /= 0) then
        call errmgr%report_error("csr_mtx_vec_mult", &
            "Memory allocation error.", LA_OUT_OF_MEMORY_ERROR)
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
module function csr_mtx_add(a, b) result(rst)
    ! Arguments
    class(csr_matrix), intent(in) :: a, b
    type(csr_matrix) :: rst

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
        call errmgr%report_error("csr_mtx_add", &
            "The matrix sizes must match.", LA_ARRAY_SIZE_ERROR)
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
    call errmgr%report_error("csr_mtx_add", "Memory allocation error.", &
        LA_OUT_OF_MEMORY_ERROR)
    return
end function

! ------------------------------------------------------------------------------
module function csr_mtx_sub(a, b) result(rst)
    ! Arguments
    class(csr_matrix), intent(in) :: a, b
    type(csr_matrix) :: rst

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
        call errmgr%report_error("csr_mtx_sub", &
            "The matrix sizes must match.", LA_ARRAY_SIZE_ERROR)
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
    call errmgr%report_error("csr_mtx_sub", "Memory allocation error.", &
        LA_OUT_OF_MEMORY_ERROR)
    return
end function

! ------------------------------------------------------------------------------
module function csr_mtx_mult_scalar_1(a, b) result(rst)
    ! Arguments
    class(csr_matrix), intent(in) :: a
    real(real64), intent(in) :: b
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
    rst%values = b * a%values
end function

! ------------------------------------------------------------------------------
module function csr_mtx_mult_scalar_2(a, b) result(rst)
    ! Arguments
    real(real64), intent(in) :: a
    class(csr_matrix), intent(in) :: b
    type(csr_matrix) :: rst

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
module function csr_mtx_divide_scalar_1(a, b) result(rst)
    ! Arguments
    class(csr_matrix), intent(in) :: a
    real(real64), intent(in) :: b
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
module function csr_transpose(a) result(rst)
    ! Arguments
    class(csr_matrix), intent(in) :: a
    type(csr_matrix) :: rst

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
module subroutine extract_diagonal_csr(a, diag, err)
    ! Arguments
    class(csr_matrix), intent(in) :: a
    real(real64), intent(out), dimension(:) :: diag
    class(errors), intent(inout), optional, target :: err

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
        call errmgr%report_error("extract_diagonal_csr", &
            "The is expected to have MIN(M, N) elements.", &
            LA_ARRAY_SIZE_ERROR)
        return
    end if

    ! Memory Allocation
    allocate(idiag(mn), stat = flag)
    if (flag /= 0) then
        call errmgr%report_error("extract_diagonal_csr", &
            "Memory allocation error.", LA_OUT_OF_MEMORY_ERROR)
        return
    end if

    ! Process
    call getdia(m, n, 0, a%values, a%column_indices, a%row_indices, len, &
        diag, idiag, 0)
end subroutine

! ------------------------------------------------------------------------------
module subroutine csr_solve_sparse_direct(a, b, x, droptol, err)
    ! Arguments
    class(csr_matrix), intent(in) :: a
    real(real64), intent(in), dimension(:) :: b
    real(real64), intent(out), dimension(:) :: x
    real(real64), intent(in), optional :: droptol
    class(errors), intent(inout), optional, target :: err

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
        call errmgr%report_error("csr_solve_sparse_direct", &
            "The matrix must be square.", LA_ARRAY_SIZE_ERROR)
        return
    end if
    if (size(x) /= n) then
        call errmgr%report_error("csr_solve_sparse_direct", &
            "Inner matrix dimension mismatch.", LA_ARRAY_SIZE_ERROR)
        return
    end if
    if (size(b) /= n) then
        call errmgr%report_error("csr_solve_sparse_direct", &
            "The output array dimension does not match the rest of the problem.", &
            LA_ARRAY_SIZE_ERROR)
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
    call errmgr%report_error("csr_solve_sparse_direct", &
        "Memory allocation error.", LA_OUT_OF_MEMORY_ERROR)
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
pure module function msr_size(x, dim) result(rst)
    ! Arguments
    class(msr_matrix), intent(in) :: x
    integer(int32), intent(in) :: dim
    integer(int32) :: rst

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
pure module function nonzero_count_msr(x) result(rst)
    ! Arguments
    class(msr_matrix), intent(in) :: x
    integer(int32) :: rst

    ! Process
    rst = x%nnz
end function

! ------------------------------------------------------------------------------
module function create_empty_msr_matrix(m, n, nnz, err) result(rst)
    ! Arguments
    integer(int32), intent(in) :: m, n, nnz
    class(errors), intent(inout), optional, target :: err
    type(msr_matrix) :: rst

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
        call errmgr%report_error("create_empty_msr_matrix", &
            "Memory allocation error.", LA_OUT_OF_MEMORY_ERROR)
        return
    end if
end function

! ------------------------------------------------------------------------------
module function csr_to_msr(a, err) result(rst)
    ! Arguments
    class(csr_matrix), intent(in) :: a
    class(errors), intent(inout), optional, target :: err
    type(msr_matrix) :: rst

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
        call errmgr%report_error("csr_to_msr", "Memory allocation error.", &
            LA_OUT_OF_MEMORY_ERROR)
        return
    end if

    ! Perform the conversion
    call csrmsr(m, ac, jc, ic, rst%values, rst%indices, work, iwork)
end function

! ------------------------------------------------------------------------------
module function msr_to_csr(a, err) result(rst)
    ! Arguments
    class(msr_matrix), intent(in) :: a
    class(errors), intent(inout), optional, target :: err
    type(csr_matrix) :: rst

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
        call errmgr%report_error("msr_to_csr", "Memory allocation error.", &
            LA_OUT_OF_MEMORY_ERROR)
        return
    end if

    ! Process
    call msrcsr(m, a%values, a%indices, rst%values, rst%column_indices, &
        rst%row_indices, work, iwork)
end function

! ------------------------------------------------------------------------------
module function dense_to_msr(a, err) result(rst)
    ! Arguments
    real(real64), intent(in), dimension(:,:) :: a
    class(errors), intent(inout), optional, target :: err
    type(msr_matrix) :: rst

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
module subroutine msr_to_dense(a, x, err)
    ! Arguments
    class(msr_matrix), intent(in) :: a
    real(real64), intent(out), dimension(:,:) :: x
    class(errors), intent(inout), optional, target :: err

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
        call errmgr%report_error("msr_to_dense", &
            "The output matrix dimensions are not correct.", &
            LA_ARRAY_SIZE_ERROR)
        return
    end if

    ! Process
    csr = msr_to_csr(a, errmgr)
    if (errmgr%has_error_occurred()) return
    call csr_to_dense(csr, x, errmgr)
end subroutine

! ------------------------------------------------------------------------------
module subroutine msr_assign_to_dense(dense, msr)
    ! Arguments
    real(real64), intent(out), dimension(:,:) :: dense
    class(msr_matrix), intent(in) :: msr

    ! Process
    call msr_to_dense(msr, dense)
end subroutine

! ------------------------------------------------------------------------------
module subroutine dense_assign_to_msr(msr, dense)
    ! Arguments
    type(msr_matrix), intent(out) :: msr
    real(real64), intent(in), dimension(:,:) :: dense

    ! Process
    msr = dense_to_msr(dense)
end subroutine

! ------------------------------------------------------------------------------
module subroutine csr_assign_to_msr(msr, csr)
    ! Arguments
    type(msr_matrix), intent(out) :: msr
    class(csr_matrix), intent(in) :: csr

    ! Process
    msr = csr_to_msr(csr)
end subroutine

! ------------------------------------------------------------------------------
module subroutine msr_assign_to_csr(csr, msr)
    ! Arguments
    type(csr_matrix), intent(out) :: csr
    class(msr_matrix), intent(in) :: msr

    ! Process
    csr = msr_to_csr(msr)
end subroutine

! ******************************************************************************
! LU PRECONDITIONER ROUTINES
! ------------------------------------------------------------------------------
module subroutine csr_lu_factor(a, lu, ju, droptol, err)
    ! Arguments
    class(csr_matrix), intent(in) :: a
    type(msr_matrix), intent(out) :: lu
    integer(int32), intent(out), dimension(:) :: ju
    real(real64), intent(in), optional :: droptol
    class(errors), intent(inout), optional, target :: err

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
        call errmgr%report_error("csr_lu_factor", &
            "U row tracking array is not sized correctly.", LA_ARRAY_SIZE_ERROR)
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
    call errmgr%report_error("csr_lu_factor", &
        "Memory allocation error.", LA_OUT_OF_MEMORY_ERROR)
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
module subroutine csr_lu_solve(lu, ju, b, x, err)
    ! Arguments
    class(msr_matrix), intent(in) :: lu
    integer(int32), intent(in), dimension(:) :: ju
    real(real64), intent(in), dimension(:) :: b
    real(real64), intent(out), dimension(:) :: x
    class(errors), intent(inout), optional, target :: err

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
        call errmgr%report_error("csr_lu_solve", &
            "The input matrix is expected to be square.", LA_ARRAY_SIZE_ERROR)
        return
    end if
    if (size(x) /= m) then
        call errmgr%report_error("csr_lu_solve", &
            "Inner matrix dimension mismatch.", LA_ARRAY_SIZE_ERROR)
        return
    end if
    if (size(b) /= m) then
        call errmgr%report_error("csr_lu_solve", &
            "The output array dimension does not match the rest of the problem.", &
            LA_ARRAY_SIZE_ERROR)
        return
    end if
    if (size(ju) /= m) then
        call errmgr%report_error("csr_lu_solve", &
            "The U row tracking array is not sized correctly.", &
            LA_ARRAY_SIZE_ERROR)
        return
    end if

    ! Process
    call lusol(m, b, x, lu%values, lu%indices, ju)
end subroutine

! ------------------------------------------------------------------------------
end submodule