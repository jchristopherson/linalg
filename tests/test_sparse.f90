module test_sparse
    use iso_fortran_env
    use linalg
    use fortran_test_helper
    implicit none

contains
! ------------------------------------------------------------------------------
function test_csr_1() result(rst)
    ! Arguments
    logical :: rst

    ! Local Variables
    integer(int32), parameter :: m = 4
    integer(int32), parameter :: n = 6
    integer(int32), parameter :: nnz = 8
    real(real64) :: dense(m, n), v(nnz)
    integer(int32) :: ja(nnz), ia(m + 1)
    real(real64), allocatable, dimension(:,:) :: check
    type(csr_matrix) :: sparse

    ! Initialization
    rst = .true.

    ! Construct a small sparse matrix, but in dense form
    dense = reshape([ &
        1.0d0, 0.0d0, 0.0d0, 0.0d0, &
        2.0d0, 3.0d0, 0.0d0, 0.0d0, &
        0.0d0, 0.0d0, 5.0d0, 0.0d0, &
        0.0d0, 4.0d0, 6.0d0, 0.0d0, &
        0.0d0, 0.0d0, 7.0d0, 0.0d0, &
        0.0d0, 0.0d0, 0.0d0, 8.0d0], [m, n])

    ! Define the sparse format - CSR
    v = [1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0, 6.0d0, 7.0d0, 8.0d0]
    ja = [1, 2, 2, 4, 3, 4, 5, 6]
    ia = [1, 3, 5, 8, 9]

    ! Form the sparse matrix
    sparse = dense_to_csr(dense)

    ! Test
    if (.not.assert(v, sparse%v)) then
        rst = .false.
        print "(A)", "Test Failed: test_csr_1 -1"
    end if
    if (.not.assert(ja, sparse%ja)) then
        rst = .false.
        print "(A)", "Test Failed: test_csr_1 -2"
    end if
    if (.not.assert(ia, sparse%ia)) then
        rst = .false.
        print "(A)", "Test Failed: test_csr_1 -3"
    end if
    if (sparse%n /= n) then
        rst = .false.
        print "(A)", "Test Failed: test_csr_1 -4"
    end if
    if (nonzero_count(sparse) /= nnz) then
        rst = .false.
        print "(A)", "Test Failed: test_csr_1 -5"
    end if

    ! Convert back to dense
    check = csr_to_dense(sparse)

    ! Test
    if (.not.assert(dense, check)) then
        rst = .false.
        print "(A)", "Test Failed: test_csr_1 -6"
    end if
end function

! ------------------------------------------------------------------------------
function test_csr_mult_1() result(rst)
    ! Arguments
    logical :: rst

    ! Local Variables
    real(real64) :: a(4, 4), ans(4, 4), dense(4, 4)
    type(csr_matrix) :: sa, sc

    ! Initialization
    rst = .true.

    ! Initialize the dense matrices
    a = reshape([ &
        5.0d0, 0.0d0, 0.0d0, 0.0d0, &
        0.0d0, 8.0d0, 0.0d0, 6.0d0, &
        0.0d0, 0.0d0, 3.0d0, 0.0d0, &
        0.0d0, 0.0d0, 0.0d0, 0.0d0], [4, 4])

    ! Construct a sparse matrix from the dense matrix
    sa = dense_to_csr(a)

    ! Multiply
    sc = matmul(sa, sa)
    ans = matmul(a, a)

    ! Convert the sparse matrix to a dense matrix for testing purposes
    dense = csr_to_dense(sc)

    ! Test
    if (.not.assert(dense, ans)) then
        rst = .false.
        print "(A)", "Test Failed: test_csr_mult_1 -1"
    end if
end function

! ------------------------------------------------------------------------------
function test_csr_mult_2() result(rst)
    ! Arguments
    logical :: rst

    ! Local Variables
    real(real64) :: a(4, 4), b(4, 6), ans(4, 6), dense(4, 6)
    type(csr_matrix) :: sa, sb, sc

    ! Initialization
    rst = .true.

    ! Initialize the dense matrices
    a = reshape([ &
        5.0d0, 0.0d0, 0.0d0, 0.0d0, &
        0.0d0, 8.0d0, 0.0d0, 6.0d0, &
        0.0d0, 0.0d0, 3.0d0, 0.0d0, &
        0.0d0, 0.0d0, 0.0d0, 0.0d0], [4, 4])
    b = reshape([ &
        1.0d0, 0.0d0, 0.0d0, 0.0d0, &
        2.0d0, 3.0d0, 0.0d0, 0.0d0, &
        0.0d0, 0.0d0, 5.0d0, 0.0d0, &
        0.0d0, 4.0d0, 6.0d0, 0.0d0, &
        0.0d0, 0.0d0, 7.0d0, 0.0d0, &
        0.0d0, 0.0d0, 0.0d0, 8.0d0], [4, 6])

    ! Construct a sparse matrix from the dense matrix
    sa = dense_to_csr(a)
    sb = dense_to_csr(b)

    ! Multiply
    sc = matmul(sa, sb)
    ans = matmul(a, b)

    ! Convert the sparse matrix to a dense matrix for testing purposes
    dense = csr_to_dense(sc)

    ! Test
    if (.not.assert(dense, ans)) then
        rst = .false.
        print "(A)", "Test Failed: test_csr_mult_2 -1"
    end if
end function

! ------------------------------------------------------------------------------
function test_csr_add_1() result(rst)
    ! Arguments
    logical :: rst

    ! Local Variables
    real(real64) :: a(4, 4), ans(4, 4), dense(4, 4)
    type(csr_matrix) :: sa, sc

    ! Initialization
    rst = .true.

    ! Initialize the dense matrices
    a = reshape([ &
        5.0d0, 0.0d0, 0.0d0, 0.0d0, &
        0.0d0, 8.0d0, 0.0d0, 6.0d0, &
        0.0d0, 0.0d0, 3.0d0, 0.0d0, &
        0.0d0, 0.0d0, 0.0d0, 0.0d0], [4, 4])

    ! Construct a sparse matrix from the dense matrix
    sa = dense_to_csr(a)

    ! Add
    sc = sa + sa
    ans = a + a

    ! Convert the sparse matrix to a dense matrix for testing purposes
    dense = csr_to_dense(sc)

    ! Test
    if (.not.assert(dense, ans)) then
        rst = .false.
        print "(A)", "Test Failed: test_csr_add_1 -1"
    end if
end function

! ------------------------------------------------------------------------------
function test_csr_subtract_1() result(rst)
    ! Arguments
    logical :: rst

    ! Local Variables
    real(real64) :: a(4, 4), dense(4, 4)
    type(csr_matrix) :: sa, sc

    ! Initialization
    rst = .true.

    ! Initialize the dense matrices
    a = reshape([ &
        5.0d0, 0.0d0, 0.0d0, 0.0d0, &
        0.0d0, 8.0d0, 0.0d0, 6.0d0, &
        0.0d0, 0.0d0, 3.0d0, 0.0d0, &
        0.0d0, 0.0d0, 0.0d0, 0.0d0], [4, 4])

    ! Construct a sparse matrix from the dense matrix
    sa = dense_to_csr(a)

    ! Add
    sc = sa + sa - sa

    ! Convert the sparse matrix to a dense matrix for testing purposes
    dense = csr_to_dense(sc)

    ! Test
    if (.not.assert(dense, a)) then
        rst = .false.
        print "(A)", "Test Failed: test_csr_subtract_1 -1"
    end if
end function

! ------------------------------------------------------------------------------
function test_csr_scalar_mult_1() result(rst)
    ! Arguments
    logical :: rst

    ! Local Variables
    real(real64) :: a(4, 4), ans(4, 4), dense(4, 4), s
    type(csr_matrix) :: sa, sc

    ! Initialization
    rst = .true.
    call random_number(s)

    ! Initialize the dense matrices
    a = reshape([ &
        5.0d0, 0.0d0, 0.0d0, 0.0d0, &
        0.0d0, 8.0d0, 0.0d0, 6.0d0, &
        0.0d0, 0.0d0, 3.0d0, 0.0d0, &
        0.0d0, 0.0d0, 0.0d0, 0.0d0], [4, 4])

    ! Construct a sparse matrix from the dense matrix
    sa = dense_to_csr(a)

    ! Multiply
    sc = s * sa
    ans = s * a

    ! Convert the sparse matrix to a dense matrix for testing purposes
    dense = csr_to_dense(sc)

    ! Test 1
    if (.not.assert(dense, ans)) then
        rst = .false.
        print "(A)", "Test Failed: test_csr_scalar_mult_1 -1"
    end if

    ! Test 2
    sc = sa * s
    dense = csr_to_dense(sc)
    if (.not.assert(dense, ans)) then
        rst = .false.
        print "(A)", "Test Failed: test_csr_scalar_mult_1 -2"
    end if
end function

! ------------------------------------------------------------------------------
function test_csr_scalar_divide_1() result(rst)
    ! Arguments
    logical :: rst

    ! Local Variables
    real(real64) :: a(4, 4), ans(4, 4), dense(4, 4), s
    type(csr_matrix) :: sa, sc

    ! Initialization
    rst = .true.
    call random_number(s)

    ! Initialize the dense matrices
    a = reshape([ &
        5.0d0, 0.0d0, 0.0d0, 0.0d0, &
        0.0d0, 8.0d0, 0.0d0, 6.0d0, &
        0.0d0, 0.0d0, 3.0d0, 0.0d0, &
        0.0d0, 0.0d0, 0.0d0, 0.0d0], [4, 4])

    ! Construct a sparse matrix from the dense matrix
    sa = dense_to_csr(a)

    ! Multiply
    sc = sa / s
    ans = a / s

    ! Convert the sparse matrix to a dense matrix for testing purposes
    dense = csr_to_dense(sc)

    ! Test 1
    if (.not.assert(dense, ans)) then
        rst = .false.
        print "(A)", "Test Failed: test_csr_scalar_divide_1 -1"
    end if
end function

! ------------------------------------------------------------------------------
function test_csr_transpose_1() result(rst)
    ! Arguments
    logical :: rst

    ! Local Variables
    real(real64) :: a(4, 6), ans(6, 4), dense(6, 4)
    type(csr_matrix) :: sa, sat

    ! Initialization
    rst = .true.
    a = reshape([ &
        1.0d0, 0.0d0, 0.0d0, 0.0d0, &
        2.0d0, 3.0d0, 0.0d0, 0.0d0, &
        0.0d0, 0.0d0, 5.0d0, 0.0d0, &
        0.0d0, 4.0d0, 6.0d0, 0.0d0, &
        0.0d0, 0.0d0, 7.0d0, 0.0d0, &
        0.0d0, 0.0d0, 0.0d0, 8.0d0], [4, 6])
    ans = transpose(a)
    sa = dense_to_csr(a)

    ! Operation
    sat = transpose(sa)
    dense = csr_to_dense(sat)

    ! Test 1
    if (.not.assert(dense, ans)) then
        rst = .false.
        print "(A)", "Test Failed: test_csr_transpose_1 -1"
    end if
end function

! ------------------------------------------------------------------------------
function test_csr_diag_mult_1() result(rst)
    ! Arguments
    logical :: rst

    ! Local Variables
    real(real64) :: alpha, a(4, 6), diag(4, 4), d(4), ans(4, 6), dense(4, 6)
    type(csr_matrix) :: sa

    ! Initialization
    rst = .true.
    call random_number(alpha)
    a = reshape([ &
        1.0d0, 0.0d0, 0.0d0, 0.0d0, &
        2.0d0, 3.0d0, 0.0d0, 0.0d0, &
        0.0d0, 0.0d0, 5.0d0, 0.0d0, &
        0.0d0, 4.0d0, 6.0d0, 0.0d0, &
        0.0d0, 0.0d0, 7.0d0, 0.0d0, &
        0.0d0, 0.0d0, 0.0d0, 8.0d0], [4, 6])
    diag = reshape([ &
        5.0d0, 0.0d0, 0.0d0, 0.0d0, &
        0.0d0, 8.0d0, 0.0d0, 0.0d0, &
        0.0d0, 0.0d0, 3.0d0, 0.0d0, &
        0.0d0, 0.0d0, 0.0d0, 6.0d0], [4, 4])
    d = [diag(1,1), diag(2,2), diag(3,3), diag(4,4)]
    sa = dense_to_csr(a)

    ! Test
    ans = alpha * matmul(diag, a)
    call diag_mtx_mult(.true., alpha, d, sa)
    dense = csr_to_dense(sa)
    if (.not.assert(dense, ans)) then
        rst = .false.
        print "(A)", "Test Failed: test_csr_diag_mult_1 -1"
    end if
end function

! ------------------------------------------------------------------------------
function test_csr_diag_mult_2() result(rst)
    ! Arguments
    logical :: rst

    ! Local Variables
    real(real64) :: alpha, a(4, 6), diag(6), dense(4, 6)
    type(csr_matrix) :: sa

    ! Initialization
    rst = .true.
    call random_number(alpha)
    a = reshape([ &
        1.0d0, 0.0d0, 0.0d0, 0.0d0, &
        2.0d0, 3.0d0, 0.0d0, 0.0d0, &
        0.0d0, 0.0d0, 5.0d0, 0.0d0, &
        0.0d0, 4.0d0, 6.0d0, 0.0d0, &
        0.0d0, 0.0d0, 7.0d0, 0.0d0, &
        0.0d0, 0.0d0, 0.0d0, 8.0d0], [4, 6])
    call random_number(diag)
    sa = dense_to_csr(a)

    ! Test
    call diag_mtx_mult(.false., alpha, diag, a) ! solution
    call diag_mtx_mult(.false., alpha, diag, sa)
    dense = csr_to_dense(sa)
    if (.not.assert(dense, a)) then
        rst = .false.
        print "(A)", "Test Failed: test_csr_diag_mult_2 -1"
    end if
end function

! ------------------------------------------------------------------------------
end module