submodule (linalg) linalg_sparse
    use sparskit
    implicit none
contains
! ------------------------------------------------------------------------------
pure module function csr_size(x, dim) result(rst)
    ! Arguments
    class(csr_matrix), intent(in) :: x
    integer(int32), intent(in) :: dim
    integer(int32) :: rst

    ! Process
    select case (dim)
    case (1)
        if (allocated(x%ia)) then
            rst = size(x%ia) - 1
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

    ! Local Variables
    integer(int32) :: m1

    ! Process
    if (allocated(x%ia)) then
        m1 = size(x, 1) + 1
        rst = x%ia(m1) - 1
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

    ! Parameters
    integer(int32), parameter :: zero = 0

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
    allocate(rst%ia(m1), rst%ja(nnz), source = zero, stat = flag)
    if (flag == 0) allocate(rst%v(nnz), source = 0.0d0, stat = flag)
    if (flag /= 0) then
        call errmgr%report_error("create_empty_csr_matrix", &
            "Memory allocation error.", LA_OUT_OF_MEMORY_ERROR)
        return
    end if
    rst%ia(m1) = nnz
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
    rst%ia(1) = 1
    do i = 1, m
        inner_loop : do j = 1, n
            if (abs(a(i,j)) < t) cycle inner_loop
            rst%ja(k) = j
            rst%v(k) = a(i,j)
            k = k + 1
        end do inner_loop
        rst%ia(i+1) = k
    end do
end function

! ------------------------------------------------------------------------------
module function csr_to_dense(a, err) result(rst)
    ! Arguments
    class(csr_matrix), intent(in) :: a
    class(errors), intent(inout), optional, target :: err
    real(real64), allocatable, dimension(:,:) :: rst

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
    allocate(rst(m, n), source = 0.0d0, stat = flag)
    if (flag /= 0) then
        call errmgr%report_error("csr_to_dense", "Memory allocation error.", &
            LA_OUT_OF_MEMORY_ERROR)
        return
    end if

    ! Process
    do i = 1, m
        do k = a%ia(i), a%ia(i+1) - 1
            j = a%ja(k)
            rst(i,j) = a%v(k)
        end do
    end do
end function

! ------------------------------------------------------------------------------
module function csr_mtx_mtx_mult(a, b, err) result(rst)
    ! Arguments
    class(csr_matrix), intent(in) :: a, b
    class(errors), intent(inout), optional, target :: err
    type(csr_matrix) :: rst

    ! Local Variables
    integer(int32), parameter :: sym_mult = 0
    integer(int32), parameter :: full_mult = 1
    integer(int32) :: flag, m, n, k, nnza, nnzb, nnzc, ierr
    integer(int32), allocatable, dimension(:) :: ic, jc, iw
    real(real64) :: dummy(1)
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
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
    call amub(m, n, sym_mult, a%v, a%ja, a%ia, b%v, b%ja, b%ia, dummy, jc, ic, &
        nnzc, iw, ierr)
    if (ierr /= 0) then
        ! NNZC was too small - try increasing it
        do while (ierr /= 0)
            deallocate(jc)
            nnzc = nnzc + nnza + nnzb
            allocate(jc(nnzc), stat = flag)
            if (flag /= 0) go to 10
            call amub(m, n, sym_mult, a%v, a%ja, a%ia, b%v, b%ja, b%ia, dummy, &
                jc, ic, nnzc, iw, ierr)
        end do
    end if

    ! Determine the actual NNZ for C & allocate space for the output
    nnzc = ic(m + 1) - 1
    deallocate(ic)
    deallocate(jc)
    rst = create_empty_csr_matrix(m, n, nnzc, errmgr)
    if (errmgr%has_error_occurred()) return

    ! Compute the actual product
    call amub(m, n, full_mult, a%v, a%ja, a%ia, b%v, b%ja, b%ia, rst%v, &
        rst%ja, rst%ia, nnzc, iw, ierr)

    ! End
    return

    ! Memory Error
10  continue
    call errmgr%report_error("csr_mtx_mtx_mult", &
        "Memory allocation error.", LA_OUT_OF_MEMORY_ERROR)
    return
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
    call aplb(m, n, sym_add, a%v, a%ja, a%ia, b%v, b%ja, b%ia, dummy, jc, ic, &
        nnzc, iw, ierr)
    if (ierr /= 0) then
        ! NNZC was too small - try increasing it
        do while (ierr /= 0)
            deallocate(jc)
            nnzc = nnzc + nnza + nnzb
            allocate(jc(nnzc), stat = flag)
            if (flag /= 0) go to 10
            call aplb(m, n, sym_add, a%v, a%ja, a%ia, b%v, b%ja, b%ia, dummy, &
                jc, ic, nnzc, iw, ierr)
        end do
    end if

    ! Determine the actuall NNZ for C & allocate space for the output
    nnzc = ic(m + 1) - 1
    deallocate(ic)
    deallocate(jc)
    rst = create_empty_csr_matrix(m, n, nnzc, errmgr)
    if (errmgr%has_error_occurred()) return

    ! Compute the actual sum
    call aplb(m, n, full_add, a%v, a%ja, a%ia, b%v, b%ja, b%ia, rst%v, rst%ja, &
        rst%ia, nnzc, iw, ierr)

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
    call aplb(m, n, sym_add, a%v, a%ja, a%ia, b%v, b%ja, b%ia, dummy, jc, ic, &
        nnzc, iw, ierr)
    if (ierr /= 0) then
        ! NNZC was too small - try increasing it
        do while (ierr /= 0)
            deallocate(jc)
            nnzc = nnzc + nnza + nnzb
            allocate(jc(nnzc), stat = flag)
            if (flag /= 0) go to 10
            call aplb(m, n, sym_add, a%v, a%ja, a%ia, b%v, b%ja, b%ia, dummy, &
                jc, ic, nnzc, iw, ierr)
        end do
    end if

    ! Determine the actuall NNZ for C & allocate space for the output
    nnzc = ic(m + 1) - 1
    deallocate(ic)
    deallocate(jc)
    rst = create_empty_csr_matrix(m, n, nnzc, errmgr)
    if (errmgr%has_error_occurred()) return

    ! Compute the actual sum
    call aplsb(m, n, a%v, a%ja, a%ia, -1.0d0, b%v, b%ja, b%ia, rst%v, rst%ja, &
        rst%ia, nnzc, iw, ierr)

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
    rst%ia = a%ia
    rst%ja = a%ja
    rst%v = b * a%v
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
    rst%ia = b%ia
    rst%ja = b%ja
    rst%v = a * b%v
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
    rst%ia = a%ia
    rst%ja = a%ja
    rst%v = a%v / b
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
    call csrcsc2(m, n, 1, 1, a%v, a%ja, a%ia, rst%v, rst%ja, rst%ia)
end function

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end submodule