! linalg_basic.f90

module linalg_basic
    use iso_fortran_env, only: int32, real64
    use blas
    use lapack
    use linalg_sparse
    use linalg_errors
    implicit none
    private
    public :: LA_NO_OPERATION
    public :: LA_TRANSPOSE
    public :: LA_HERMITIAN_TRANSPOSE
    public :: mtx_mult
    public :: rank1_update
    public :: diag_mtx_mult
    public :: trace
    public :: mtx_rank
    public :: det
    public :: swap
    public :: recip_mult_array
    public :: tri_mtx_mult
    public :: band_mtx_mult
    public :: band_mtx_to_full_mtx
    public :: band_diag_mtx_mult
    public :: banded_to_dense
    public :: dense_to_banded
    public :: extract_diagonal
    public :: extract_upper_triangular
    public :: extract_lower_triangular
    public :: identity

    integer(int32), parameter :: LA_NO_OPERATION = 0
        !! Defines no operation should be performed on the matrix.
    integer(int32), parameter :: LA_TRANSPOSE = 1
        !! Defines a transpose operation.
    integer(int32), parameter :: LA_HERMITIAN_TRANSPOSE = 2
        !! Defines a Hermitian transpose operation for a complex-valued matrix.

    interface mtx_mult
        !! An interface to the matrix multiplication routines.
        module procedure :: mtx_mult_mtx
        module procedure :: mtx_mult_vec
        module procedure :: cmtx_mult_mtx
        module procedure :: cmtx_mult_vec
    end interface

    interface rank1_update
        !! An interface to the rank-1 update routines.
        module procedure :: rank1_update_dbl
        module procedure :: rank1_update_cmplx
    end interface

    interface diag_mtx_mult
        !! An interface to the diagonal matrix multiplication routines.
        module procedure :: diag_mtx_mult_mtx
        module procedure :: diag_mtx_mult_mtx2
        module procedure :: diag_mtx_mult_mtx3
        module procedure :: diag_mtx_mult_mtx4
        module procedure :: diag_mtx_mult_mtx_cmplx
        module procedure :: diag_mtx_mult_mtx2_cmplx
        module procedure :: diag_mtx_mult_mtx_mix
        module procedure :: diag_mtx_mult_mtx2_mix
        module procedure :: diag_mtx_sparse_mult
    end interface

    interface trace
        !! An interface to the trace routines.
        module procedure :: trace_dbl
        module procedure :: trace_cmplx
    end interface

    interface mtx_rank
        !! An interface to the matrix rank routines.
        module procedure :: mtx_rank_dbl
        module procedure :: mtx_rank_cmplx
    end interface

    interface det
        !! An interface to the determinant routines.
        module procedure :: det_dbl
        module procedure :: det_cmplx
    end interface

    interface swap
        !! An interface to the swap routines.
        module procedure :: swap_dbl
        module procedure :: swap_cmplx
    end interface

    interface recip_mult_array
        !! An interface to the reciprocal multiplication routines.
        module procedure :: recip_mult_array_dbl
    end interface

    interface tri_mtx_mult
        !! An interface to the triangular matrix multiplication routines.
        module procedure :: tri_mtx_mult_dbl
        module procedure :: tri_mtx_mult_cmplx
    end interface

    interface band_mtx_mult
        !! An interface to the banded matrix multiplication routines.
        module procedure :: band_mtx_vec_mult_dbl
        module procedure :: band_mtx_vec_mult_cmplx
    end interface

    interface band_mtx_to_full_mtx
        !! An interface to the banded matrix to full matrix conversion routines.
        module procedure :: band_to_full_mtx_dbl
        module procedure :: band_to_full_mtx_cmplx
    end interface

    interface band_diag_mtx_mult
        !! An interface to the banded diagonal matrix multiplication routines.
        module procedure :: band_diag_mtx_mult_dbl
        module procedure :: band_diag_mtx_mult_cmplx
    end interface

    interface banded_to_dense
        !! An interface to the banded to dense matrix conversion routines.
        module procedure :: banded_to_dense_dbl
        module procedure :: banded_to_dense_cmplx
    end interface

    interface dense_to_banded
        !! An interface to the dense to banded matrix conversion routines.
        module procedure :: dense_to_banded_dbl
        module procedure :: dense_to_banded_cmplx
    end interface

    interface extract_diagonal
        !! An interface to the diagonal extraction routines.
        module procedure :: extract_diagonal_dbl
        module procedure :: extract_diagonal_cmplx
        module procedure :: extract_diagonal_csr
    end interface

    interface extract_upper_triangular
        !! An interface to the upper triangular matrix extraction routines.
        module procedure :: extract_upper_triangular_dbl
        module procedure :: extract_upper_triangular_cmplx
    end interface

    interface extract_lower_triangular
        !! An interface to the lower triangular matrix extraction routines.
        module procedure :: extract_lower_triangular_dbl
        module procedure :: extract_lower_triangular_cmplx
    end interface

contains
! ******************************************************************************
! MATRIX MULTIPLICATION ROUTINES
! ------------------------------------------------------------------------------
pure subroutine mtx_mult_mtx(transa, transb, alpha, a, b, beta, c)
    !! Performs the matrix operation \(C = \alpha A B + \beta C \).
    logical, intent(in) :: transa
        !! A logical flag indicating if the matrix \(A\) should be transposed.
    logical, intent(in) :: transb
        !! A logical flag indicating if the matrix \(B\) should be transposed.
    real(real64), intent(in) :: alpha
        !! The scalar \(\alpha\) to multiply the product of \(A\) and \(B\).
    real(real64), intent(in) :: beta
        !! The scalar \(\beta\) to multiply the matrix \(C\).
    real(real64), intent(in), dimension(:,:) :: a
        !! The matrix \(A\) in the operation.
    real(real64), intent(in), dimension(:,:) :: b
        !! The matrix \(B\) in the operation.
    real(real64), intent(inout), dimension(:,:) :: c
        !! The matrix \(C\) in the operation.

    ! Parameters
    real(real64), parameter :: zero = 0.0d0
    real(real64), parameter :: one = 1.0d0

    ! Local Variables
    character :: ta, tb
    integer(int32) :: m, n, k, lda, ldb, flag

    ! Initialization
    m = size(c, 1)
    n = size(c, 2)
    if (transa) then ! K = # of columns in op(A) (# of rows in op(B))
        k = size(a, 1)
        ta = 'T'
        lda = k
    else
        k = size(a, 2)
        ta = 'N'
        lda = m
    end if
    if (transb) then
        tb = 'T'
        ldb = n
    else
        tb = 'N'
        ldb = k
    end if

    ! Input Check
    flag = 0
    if (transa) then
        if (size(a, 2) /= m) flag = 4
    else
        if (size(a, 1) /= m) flag = 4
    end if
    if (transb) then
        if (size(b, 2) /= k .or. size(b, 1) /= n) flag = 5
    else
        if (size(b, 1) /= k .or. size(b, 2) /= n) flag = 5
    end if
    if (flag /= 0) then
        error stop flag
    end if

    ! Call DGEMM
    call DGEMM(ta, tb, m, n, k, alpha, a, lda, b, ldb, beta, c, m)
end subroutine

! ------------------------------------------------------------------------------
pure subroutine mtx_mult_vec(trans, alpha, a, b, beta, c)
    !! Performs the matrix-vector operation \(C = \alpha A B + \beta C \).
    logical, intent(in) :: trans
        !! A logical flag indicating if the matrix \(A\) should be transposed.
    real(real64), intent(in) :: alpha
        !! The scalar \(\alpha\) to multiply the product of \(A\) and \(B\).
    real(real64), intent(in) :: beta
        !! The scalar \(\beta\) to multiply the vector \(C\).
    real(real64), intent(in), dimension(:,:) :: a
        !! The matrix \(A\) in the operation.
    real(real64), intent(in), dimension(:) :: b
        !! The vector \(B\) in the operation.
    real(real64), intent(inout), dimension(:) :: c
        !! The vector \(C\) in the operation.

    ! Local Variables
    character :: t
    integer(int32) :: m, n, flag

    ! Initialization
    m = size(a, 1)
    n = size(a, 2)
    t = 'N'
    if (trans) t = 'T'

    ! Input Check
    flag = 0
    if (trans) then
        if (size(b) /= m) then
            flag = 4
        else if (size(c) /= n) then
            flag = 6
        end if
    else
        if (size(b) /= n) then
            flag = 4
        else if (size(c) /= m) then
            flag = 6
        end if
    end if
    if (flag /= 0) then
        ! ERROR: Matrix dimensions mismatch
        error stop flag
    end if

    ! Call DGEMV
    call DGEMV(t, m, n, alpha, a, m, b, 1, beta, c, 1)
end subroutine

! xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx !
!                           COMPLEX VALUED VERSIONS                            !
! xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx !
pure subroutine cmtx_mult_mtx(opa, opb, alpha, a, b, beta, c)
    !! Performs the matrix operation \(C = \alpha A B + \beta C \).
    integer(int32), intent(in) :: opa
        !! An integer flag indicating the operation to perform on matrix \(A\).
        !! Possible options are:
        !!
        !! - LA_NO_OPERATION: No operation is performed on matrix.
        !!
        !! - LA_TRANSPOSE: The transpose of matrix is used.
        !!
        !! - LA_HERMITIAN_TRANSPOSE: The Hermitian transpose of matrix is used.
    integer(int32), intent(in) :: opb
        !! An integer flag indicating the operation to perform on matrix \(B\).
        !! Possible options are:
        !!
        !! - LA_NO_OPERATION: No operation is performed on matrix.
        !!
        !! - LA_TRANSPOSE: The transpose of matrix is used.
        !!
        !! - LA_HERMITIAN_TRANSPOSE: The Hermitian transpose of matrix is used.
    complex(real64), intent(in) :: alpha
        !! The scalar \(\alpha\) to multiply the product of \(A\) and \(B\).
    complex(real64), intent(in) :: beta
        !! The scalar \(\beta\) to multiply the matrix \(C\).
    complex(real64), intent(in), dimension(:,:) :: a
        !! The matrix \(A\) in the operation.
    complex(real64), intent(in), dimension(:,:) :: b
        !! The matrix \(B\) in the operation.
    complex(real64), intent(inout), dimension(:,:) :: c
        !! The matrix \(C\) in the operation.

    ! Parameters
    real(real64), parameter :: zero = 0.0d0
    real(real64), parameter :: one = 1.0d0

    ! Local Variables
    character :: ta, tb
    integer(int32) :: m, n, k, lda, ldb, flag

    ! Initialization
    m = size(c, 1)
    n = size(c, 2)
    if (opa == LA_TRANSPOSE) then ! K = # of columns in op(A) (# of rows in op(B))
        k = size(a, 1)
        ta = 'T'
        lda = k
    else if (opa == LA_HERMITIAN_TRANSPOSE) then
        k = size(a, 1)
        ta = 'C'
        lda = k
    else
        k = size(a, 2)
        ta = 'N'
        lda = m
    end if
    if (opb == LA_TRANSPOSE) then
        tb = 'T'
        ldb = n
    else if (opb ==  LA_HERMITIAN_TRANSPOSE) then
        tb = 'C'
        ldb = n
    else
        tb = 'N'
        ldb = k
    end if

    ! Input Check
    flag = 0
    if (opa == LA_TRANSPOSE .or. opa ==  LA_HERMITIAN_TRANSPOSE) then
        if (size(a, 2) /= m) flag = 4
    else
        if (size(a, 1) /= m) flag = 4
    end if
    if (opb == LA_TRANSPOSE .or. opb ==  LA_HERMITIAN_TRANSPOSE) then
        if (size(b, 2) /= k .or. size(b, 1) /= n) flag = 5
    else
        if (size(b, 1) /= k .or. size(b, 2) /= n) flag = 5
    end if
    if (flag /= 0) then
        ! ERROR: Matrix dimensions mismatch
        error stop flag
    end if

    ! Call ZGEMM
    call ZGEMM(ta, tb, m, n, k, alpha, a, lda, b, ldb, beta, c, m)
end subroutine

! ------------------------------------------------------------------------------
pure subroutine cmtx_mult_vec(opa, alpha, a, b, beta, c)
    !! Performs the matrix-vector operation \(C = \alpha A B + \beta C \).
    integer(int32), intent(in) :: opa
        !! An integer flag indicating the operation to perform on matrix \(A\).
        !! Possible options are:
        !!
        !! - LA_NO_OPERATION: No operation is performed on matrix.
        !!
        !! - LA_TRANSPOSE: The transpose of matrix is used.
        !!
        !! - LA_HERMITIAN_TRANSPOSE: The Hermitian transpose of matrix is used.
    complex(real64), intent(in) :: alpha
        !! The scalar \(\alpha\) to multiply the product of \(A\) and \(B\).
    complex(real64), intent(in) :: beta
        !! The scalar \(\beta\) to multiply the vector \(C\).
    complex(real64), intent(in), dimension(:,:) :: a
        !! The matrix \(A\) in the operation.
    complex(real64), intent(in), dimension(:) :: b
        !! The vector \(B\) in the operation.
    complex(real64), intent(inout), dimension(:) :: c
        !! The vector \(C\) in the operation.

    ! Local Variables
    character :: t
    integer(int32) :: m, n, flag

    ! Initialization
    m = size(a, 1)
    n = size(a, 2)
    if (opa == LA_TRANSPOSE) then
        t = 'T'
    else if (opa ==  LA_HERMITIAN_TRANSPOSE) then
        t = 'C'
    else
        t = 'N'
    end if

    ! Input Check
    flag = 0
    if (opa == LA_TRANSPOSE .or. opa ==  LA_HERMITIAN_TRANSPOSE) then
        if (size(b) /= m) then
            flag = 4
        else if (size(c) /= n) then
            flag = 6
        end if
    else
        if (size(b) /= n) then
            flag = 4
        else if (size(c) /= m) then
            flag = 6
        end if
    end if
    if (flag /= 0) then
        ! ERROR: Matrix dimensions mismatch
        error stop flag
    end if

    ! Call ZGEMV
    call ZGEMV(t, m, n, alpha, a, m, b, 1, beta, c, 1)
end subroutine

! ******************************************************************************
! RANK 1 UPDATE
! ------------------------------------------------------------------------------
pure subroutine rank1_update_dbl(alpha, x, y, a)
    !! Performs a rank-1 update of a matrix of the form \(A = \alpha x y^T + A\).
    real(real64), intent(in) :: alpha
        !! The scalar \(\alpha\) to multiply the outer product of \(x\) and \(y\).
    real(real64), intent(in), dimension(:) :: x
        !! The vector \(x\) in the outer product.
    real(real64), intent(in), dimension(:) :: y
        !! The vector \(y\) in the outer product.
    real(real64), intent(inout), dimension(:,:) :: a
        !! The matrix \(A\) to update.

    ! Parameters
    real(real64), parameter :: zero = 0.0d0

    ! Local Variables
    integer(int32) :: j, m, n
    real(real64) :: temp

    ! Initialization
    m = size(x)
    n = size(y)

    ! Input Check
    if (size(a, 1) /= m .or. size(a, 2) /= n) then
        ! ERROR: Matrix dimension array
        error stop 4
    end if

    ! Process
    do j = 1, n
        if (y(j) /= zero) then
            temp = alpha * y(j)
            a(:,j) = a(:,j) + temp * x
        end if
    end do
end subroutine

! xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx !
!                           COMPLEX VALUED VERSIONS                            !
! xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx !
pure subroutine rank1_update_cmplx(alpha, x, y, a)
    !! Performs a rank-1 update of a matrix of the form \(A = \alpha x y^H + A\).
    complex(real64), intent(in) :: alpha
        !! The scalar \(\alpha\) to multiply the outer product of \(x\) and \(y\).
    complex(real64), intent(in), dimension(:) :: x
        !! The vector \(x\) in the outer product.
    complex(real64), intent(in), dimension(:) :: y
        !! The vector \(y\) in the outer product.
    complex(real64), intent(inout), dimension(:,:) :: a
        !! The matrix \(A\) to update.

    ! Parameters
    complex(real64), parameter :: zero = (0.0d0, 0.0d0)

    ! Local Variables
    integer(int32) :: j, m, n
    complex(real64) :: temp

    ! Initialization
    m = size(x)
    n = size(y)

    ! Input Check
    if (size(a, 1) /= m .or. size(a, 2) /= n) then
        ! ERROR: Matrix dimension array
        error stop 4
    end if

    ! Process
    do j = 1, n
        if (y(j) /= zero) then
            temp = alpha * conjg(y(j))
            a(:,j) = a(:,j) + temp * x
        end if
    end do
end subroutine

! ******************************************************************************
! DIAGONAL MATRIX MULTIPLICATION ROUTINES
! ------------------------------------------------------------------------------
pure subroutine diag_mtx_mult_mtx(lside, trans, alpha, a, b, beta, c)
    !! Performs the matrix operation \(C = \alpha A op(B) + \beta C \) or
    !! \(C = \alpha op(B) A + \beta C \) where \(A\) is a diagonal matrix.
    logical, intent(in) :: lside
        !! A logical flag indicating if the diagonal matrix is on the left.
    logical, intent(in) :: trans
        !! A logical flag indicating if the matrix \(B\) should be transposed.
    real(real64), intent(in) :: alpha
        !! The scalar \(\alpha\) to multiply the product of \(A\) and \(B\).
    real(real64), intent(in) :: beta
        !! The scalar \(\beta\) to multiply the matrix \(C\).
    real(real64), intent(in), dimension(:) :: a
        !! The diagonal matrix \(A\) in the operation.
    real(real64), intent(in), dimension(:,:) :: b
        !! The matrix \(B\) in the operation.
    real(real64), intent(inout), dimension(:,:) :: c
        !! The matrix \(C\) in the operation.

    ! Parameters
    real(real64), parameter :: zero = 0.0d0
    real(real64), parameter :: one = 1.0d0

    ! Local Variables
    integer(int32) :: i, m, n, k, nrowb, ncolb, flag
    real(real64) :: temp

    ! Initialization
    m = size(c, 1)
    n = size(c, 2)
    k = size(a)
    nrowb = size(b, 1)
    ncolb = size(b, 2)

    ! Input Check
    flag = 0
    if (lside) then
        if (k > m) then
            flag = 4
        else
            if (trans) then
                ! Compute C = alpha * A * B**T + beta * C
                if (nrowb /= n .or. ncolb < k) flag = 5
            else
                ! Compute C = alpha * A * B + beta * C
                if (nrowb < k .or. ncolb /= n) flag = 5
            end if
        end if
    else
        if (k > n) then
            flag = 4
        else
            if (trans) then
                ! Compute C = alpha * B**T * A + beta * C
                if (ncolb /= m .or. nrowb < k) flag = 5
            else
                ! Compute C = alpha * B * A + beta * C
                if (nrowb /= m .or. ncolb < k) flag = 5
            end if
        end if
    end if
    if (flag /= 0) then
        ! ERROR: One of the input arrays is not sized correctly
        error stop flag
    end if

    ! Deal with ALPHA == 0
    if (alpha == 0) then
        if (beta == zero) then
            c = zero
        else if (beta /= one) then
            c = beta * c
        end if
        return
    end if

    ! Process
    if (lside) then
        if (trans) then
            ! Compute C = alpha * A * B**T + beta * C
            do i = 1, k
                if (beta == zero) then
                    c(i,:) = zero
                else if (beta /= one) then
                    c(i,:) = beta * c(i,:)
                end if
                temp = alpha * a(i)
                c(i,:) = c(i,:) + temp * b(:,i)
            end do
        else
            ! Compute C = alpha * A * B + beta * C
            do i = 1, k
                if (beta == zero) then
                    c(i,:) = zero
                else if (beta /= one) then
                    c(i,:) = beta * c(i,:)
                end if
                temp = alpha * a(i)
                c(i,:) = c(i,:) + temp * b(i,:)
            end do
        end if

        ! Handle extra rows
        if (m > k) then
            if (beta == zero) then
                c(k+1:m,:) = zero
            else
                c(k+1:m,:) = beta * c(k+1:m,:)
            end if
        end if
    else
        if (trans) then
            ! Compute C = alpha * B**T * A + beta * C
            do i = 1, k
                if (beta == zero) then
                    c(:,i) = zero
                else if (beta /= one) then
                    c(:,i) = beta * c(:,i)
                end if
                temp = alpha * a(i)
                c(:,i) = c(:,i) + temp * b(i,:)
            end do
        else
            ! Compute C = alpha * B * A + beta * C
            do i = 1, k
                if (beta == zero) then
                    c(:,i) = zero
                else if (beta /= one) then
                    c(:,i) = beta * c(:,i)
                end if
                temp = alpha * a(i)
                c(:,i) = c(:,i) + temp * b(:,i)
            end do
        end if

        ! Handle extra columns
        if (n > k) then
            if (beta == zero) then
                c(:,k+1:m) = zero
            else if (beta /= one) then
                c(:,k+1:m) = beta * c(:,k+1:m)
            end if
        end if
    end if
end subroutine

! ------------------------------------------------------------------------------
pure subroutine diag_mtx_mult_mtx2(lside, alpha, a, b)
    !! Performs the matrix operation \(B = \alpha A B \) or \(B = \alpha B A \)
    !! where \(A\) is a diagonal matrix.
    logical, intent(in) :: lside
        !! A logical flag indicating if the diagonal matrix is on the left.
    real(real64), intent(in) :: alpha
        !! The scalar \(\alpha\) to multiply the product of \(A\) and \(B\).
    real(real64), intent(in), dimension(:) :: a
        !! The diagonal matrix \(A\) in the operation.
    real(real64), intent(inout), dimension(:,:) :: b
        !! The matrix \(B\) in the operation.

    ! Parameters
    real(real64), parameter :: zero = 0.0d0
    real(real64), parameter :: one = 1.0d0

    ! Local Variables
    integer(int32) :: i, m, n, k
    real(real64) :: temp

    ! Initialization
    m = size(b, 1)
    n = size(b, 2)
    k = size(a)

    ! Input Check
    if ((lside .and. k > m) .or. (.not.lside .and. k > n)) then
        ! ERROR: One of the input arrays is not sized correctly
        error stop 3
    end if

    ! Process
    if (lside) then
        ! Compute B = alpha * A * B
        do i = 1, k
            temp = alpha * a(i)
            b(i,:) = temp * b(i,:)
        end do
        if (m > k) b(k+1:m,:) = zero
    else
        ! Compute B = alpha * B * A
        do i = 1, k
            temp = alpha * a(i)
            b(:,i) = temp * b(:,i)
        end do
        if (n > k) b(:,k+1:n) = zero
    end if
end subroutine

! ------------------------------------------------------------------------------
pure subroutine diag_mtx_mult_mtx3(lside, trans, alpha, a, b, beta, c)
    !! Performs the matrix operation \(C = \alpha A op(B) + \beta C \) or
    !! \(C = \alpha B A + \beta C \) where \(A\) is a diagonal matrix.
    logical, intent(in) :: lside
        !! A logical flag indicating if the diagonal matrix is on the left.
    logical, intent(in) :: trans
        !! A logical flag indicating if the matrix \(B\) should be transposed.
    real(real64), intent(in) :: alpha
        !! The scalar \(\alpha\) to multiply the product of \(A\) and \(B\).
    real(real64), intent(in) :: beta
        !! The scalar \(\beta\) to multiply the matrix \(C\).
    complex(real64), intent(in), dimension(:) :: a
        !! The diagonal matrix \(A\) in the operation.
    real(real64), intent(in), dimension(:,:) :: b
        !! The matrix \(B\) in the operation.
    complex(real64), intent(inout), dimension(:,:) :: c
        !! The matrix \(C\) in the operation.

    ! Parameters
    complex(real64), parameter :: zero = (0.0d0, 0.0d0)
    complex(real64), parameter :: one = (1.0d0, 0.0d0)

    ! Local Variables
    integer(int32) :: i, m, n, k, nrowb, ncolb, flag
    complex(real64) :: temp

    ! Initialization
    m = size(c, 1)
    n = size(c, 2)
    k = size(a)
    nrowb = size(b, 1)
    ncolb = size(b, 2)

    ! Input Check
    flag = 0
    if (lside) then
        if (k > m) then
            flag = 4
        else
            if (trans) then
                ! Compute C = alpha * A * B**T + beta * C
                if (nrowb /= n .or. ncolb < k) flag = 5
            else
                ! Compute C = alpha * A * B + beta * C
                if (nrowb < k .or. ncolb /= n) flag = 5
            end if
        end if
    else
        if (k > n) then
            flag = 4
        else
            if (trans) then
                ! Compute C = alpha * B**T * A + beta * C
                if (ncolb /= m .or. nrowb < k) flag = 5
            else
                ! Compute C = alpha * B * A + beta * C
                if (nrowb /= m .or. ncolb < k) flag = 5
            end if
        end if
    end if
    if (flag /= 0) then
        ! ERROR: One of the input arrays is not sized correctly
        error stop flag
    end if

    ! Deal with ALPHA == 0
    if (alpha == 0) then
        if (beta == zero) then
            c = zero
        else if (beta /= one) then
            c = beta * c
        end if
        return
    end if

    ! Process
    if (lside) then
        if (trans) then
            ! Compute C = alpha * A * B**T + beta * C
            do i = 1, k
                if (beta == zero) then
                    c(i,:) = zero
                else if (beta /= one) then
                    c(i,:) = beta * c(i,:)
                end if
                temp = alpha * a(i)
                c(i,:) = c(i,:) + temp * b(:,i)
            end do
        else
            ! Compute C = alpha * A * B + beta * C
            do i = 1, k
                if (beta == zero) then
                    c(i,:) = zero
                else if (beta /= one) then
                    c(i,:) = beta * c(i,:)
                end if
                temp = alpha * a(i)
                c(i,:) = c(i,:) + temp * b(i,:)
            end do
        end if

        ! Handle extra rows
        if (m > k) then
            if (beta == zero) then
                c(k+1:m,:) = zero
            else
                c(k+1:m,:) = beta * c(k+1:m,:)
            end if
        end if
    else
        if (trans) then
            ! Compute C = alpha * B**T * A + beta * C
            do i = 1, k
                if (beta == zero) then
                    c(:,i) = zero
                else if (beta /= one) then
                    c(:,i) = beta * c(:,i)
                end if
                temp = alpha * a(i)
                c(:,i) = c(:,i) + temp * b(i,:)
            end do
        else
            ! Compute C = alpha * B * A + beta * C
            do i = 1, k
                if (beta == zero) then
                    c(:,i) = zero
                else if (beta /= one) then
                    c(:,i) = beta * c(:,i)
                end if
                temp = alpha * a(i)
                c(:,i) = c(:,i) + temp * b(:,i)
            end do
        end if

        ! Handle extra columns
        if (n > k) then
            if (beta == zero) then
                c(:,k+1:m) = zero
            else if (beta /= one) then
                c(:,k+1:m) = beta * c(:,k+1:m)
            end if
        end if
    end if
end subroutine

! ------------------------------------------------------------------------------
pure subroutine diag_mtx_mult_mtx4(lside, opb, alpha, a, b, beta, c)
    !! Performs the matrix operation \(C = \alpha A op(B) + \beta C \) or
    !! \(C = \alpha op(B) A + \beta C \) where \(A\) is a diagonal matrix.
    logical, intent(in) :: lside
        !! A logical flag indicating if the diagonal matrix is on the left.
    integer(int32), intent(in) :: opb
        !! An integer flag indicating the operation to perform on matrix \(B\).
        !! Possible options are:
        !!
        !! - LA_NO_OPERATION: No operation is performed on matrix.
        !!
        !! - LA_TRANSPOSE: The transpose of matrix is used.
        !!
        !! - LA_HERMITIAN_TRANSPOSE: The Hermitian transpose of matrix is used.
    real(real64), intent(in) :: alpha
        !! The scalar \(\alpha\) to multiply the product of \(A\) and \(B\).
    real(real64), intent(in) :: beta
        !! The scalar \(\beta\) to multiply the matrix \(C\).
    complex(real64), intent(in), dimension(:) :: a
        !! The diagonal matrix \(A\) in the operation.
    complex(real64), intent(in), dimension(:,:) :: b
        !! The matrix \(B\) in the operation.
    complex(real64), intent(inout), dimension(:,:) :: c
        !! The matrix \(C\) in the operation.

    ! Parameters
    complex(real64), parameter :: zero = (0.0d0, 0.0d0)
    complex(real64), parameter :: one = (1.0d0, 0.0d0)

    ! Local Variables
    integer(int32) :: i, m, n, k, nrowb, ncolb, flag
    complex(real64) :: temp

    ! Initialization
    m = size(c, 1)
    n = size(c, 2)
    k = size(a)
    nrowb = size(b, 1)
    ncolb = size(b, 2)

    ! Input Check
    flag = 0
    if (lside) then
        if (k > m) then
            flag = 4
        else
            if (opb == LA_TRANSPOSE .or. opb ==  LA_HERMITIAN_TRANSPOSE) then
                ! Compute C = alpha * A * B**T + beta * C
                if (nrowb /= n .or. ncolb < k) flag = 5
            else
                ! Compute C = alpha * A * B + beta * C
                if (nrowb < k .or. ncolb /= n) flag = 5
            end if
        end if
    else
        if (k > n) then
            flag = 4
        else
            if (opb == LA_TRANSPOSE .or. opb ==  LA_HERMITIAN_TRANSPOSE) then
                ! Compute C = alpha * B**T * A + beta * C
                if (ncolb /= m .or. nrowb < k) flag = 5
            else
                ! Compute C = alpha * B * A + beta * C
                if (nrowb /= m .or. ncolb < k) flag = 5
            end if
        end if
    end if
    if (flag /= 0) then
        ! ERROR: One of the input arrays is not sized correctly
        error stop flag
    end if

    ! Deal with ALPHA == 0
    if (alpha == 0) then
        if (beta == zero) then
            c = zero
        else if (beta /= one) then
            c = beta * c
        end if
        return
    end if

    ! Process
    if (lside) then
        if (opb == LA_TRANSPOSE) then
            ! Compute C = alpha * A * B**T + beta * C
            do i = 1, k
                if (beta == zero) then
                    c(i,:) = zero
                else if (beta /= one) then
                    c(i,:) = beta * c(i,:)
                end if
                temp = alpha * a(i)
                c(i,:) = c(i,:) + temp * b(:,i)
            end do
        else if (opb ==  LA_HERMITIAN_TRANSPOSE) then
            ! Compute C = alpha * A * B**H + beta * C
            do i = 1, k
                if (beta == zero) then
                    c(i,:) = zero
                else if (beta /= one) then
                    c(i,:) = beta * c(i,:)
                end if
                temp = alpha * a(i)
                c(i,:) = c(i,:) + temp * conjg(b(:,i))
            end do
        else
            ! Compute C = alpha * A * B + beta * C
            do i = 1, k
                if (beta == zero) then
                    c(i,:) = zero
                else if (beta /= one) then
                    c(i,:) = beta * c(i,:)
                end if
                temp = alpha * a(i)
                c(i,:) = c(i,:) + temp * b(i,:)
            end do
        end if

        ! Handle extra rows
        if (m > k) then
            if (beta == zero) then
                c(k+1:m,:) = zero
            else
                c(k+1:m,:) = beta * c(k+1:m,:)
            end if
        end if
    else
        if (opb == LA_TRANSPOSE) then
            ! Compute C = alpha * B**T * A + beta * C
            do i = 1, k
                if (beta == zero) then
                    c(:,i) = zero
                else if (beta /= one) then
                    c(:,i) = beta * c(:,i)
                end if
                temp = alpha * a(i)
                c(:,i) = c(:,i) + temp * b(i,:)
            end do
        else if (opb ==  LA_HERMITIAN_TRANSPOSE) then
            ! Compute C = alpha * B**H * A + beta * C
            do i = 1, k
                if (beta == zero) then
                    c(:,i) = zero
                else if (beta /= one) then
                    c(:,i) = beta * c(:,i)
                end if
                temp = alpha * a(i)
                c(:,i) = c(:,i) + temp * conjg(b(i,:))
            end do
        else
            ! Compute C = alpha * B * A + beta * C
            do i = 1, k
                if (beta == zero) then
                    c(:,i) = zero
                else if (beta /= one) then
                    c(:,i) = beta * c(:,i)
                end if
                temp = alpha * a(i)
                c(:,i) = c(:,i) + temp * b(:,i)
            end do
        end if

        ! Handle extra columns
        if (n > k) then
            if (beta == zero) then
                c(:,k+1:m) = zero
            else if (beta /= one) then
                c(:,k+1:m) = beta * c(:,k+1:m)
            end if
        end if
    end if
end subroutine

! ------------------------------------------------------------------------------
pure subroutine diag_mtx_mult_mtx_cmplx(lside, opb, alpha, a, b, beta, c)
    !! Performs the matrix operation \(C = \alpha A op(B) + \beta C \) or
    !! \(C = \alpha op(B) A + \beta C \) where \(A\) is a diagonal matrix.
    logical, intent(in) :: lside
        !! A logical flag indicating if the diagonal matrix is on the left.
    integer(int32), intent(in) :: opb
        !! An integer flag indicating the operation to perform on matrix \(B\).
        !! Possible options are:
        !!
        !! - LA_NO_OPERATION: No operation is performed on matrix.
        !!
        !! - LA_TRANSPOSE: The transpose of matrix is used.
        !!
        !! - LA_HERMITIAN_TRANSPOSE: The Hermitian transpose of matrix is used.
    complex(real64), intent(in) :: alpha
        !! The scalar \(\alpha\) to multiply the product of \(A\) and \(B\).
    complex(real64), intent(in) :: beta
        !! The scalar \(\beta\) to multiply the matrix \(C\).
    complex(real64), intent(in), dimension(:) :: a
        !! The diagonal matrix \(A\) in the operation.
    complex(real64), intent(in), dimension(:,:) :: b
        !! The matrix \(B\) in the operation.
    complex(real64), intent(inout), dimension(:,:) :: c
        !! The matrix \(C\) in the operation.

    ! Parameters
    complex(real64), parameter :: zero = (0.0d0, 0.0d0)
    complex(real64), parameter :: one = (1.0d0, 0.0d0)

    ! Local Variables
    integer(int32) :: i, m, n, k, nrowb, ncolb, flag
    complex(real64) :: temp

    ! Initialization
    m = size(c, 1)
    n = size(c, 2)
    k = size(a)
    nrowb = size(b, 1)
    ncolb = size(b, 2)

    ! Input Check
    flag = 0
    if (lside) then
        if (k > m) then
            flag = 4
        else
            if (opb == LA_TRANSPOSE .or. opb ==  LA_HERMITIAN_TRANSPOSE) then
                ! Compute C = alpha * A * B**T + beta * C
                if (nrowb /= n .or. ncolb < k) flag = 5
            else
                ! Compute C = alpha * A * B + beta * C
                if (nrowb < k .or. ncolb /= n) flag = 5
            end if
        end if
    else
        if (k > n) then
            flag = 4
        else
            if (opb == LA_TRANSPOSE .or. opb ==  LA_HERMITIAN_TRANSPOSE) then
                ! Compute C = alpha * B**T * A + beta * C
                if (ncolb /= m .or. nrowb < k) flag = 5
            else
                ! Compute C = alpha * B * A + beta * C
                if (nrowb /= m .or. ncolb < k) flag = 5
            end if
        end if
    end if
    if (flag /= 0) then
        ! ERROR: One of the input arrays is not sized correctly
        error stop flag
    end if

    ! Deal with ALPHA == 0
    if (alpha == 0) then
        if (beta == zero) then
            c = zero
        else if (beta /= one) then
            c = beta * c
        end if
        return
    end if

    ! Process
    if (lside) then
        if (opb == LA_TRANSPOSE) then
            ! Compute C = alpha * A * B**T + beta * C
            do i = 1, k
                if (beta == zero) then
                    c(i,:) = zero
                else if (beta /= one) then
                    c(i,:) = beta * c(i,:)
                end if
                temp = alpha * a(i)
                c(i,:) = c(i,:) + temp * b(:,i)
            end do
        else if (opb ==  LA_HERMITIAN_TRANSPOSE) then
            ! Compute C = alpha * A * B**H + beta * C
            do i = 1, k
                if (beta == zero) then
                    c(i,:) = zero
                else if (beta /= one) then
                    c(i,:) = beta * c(i,:)
                end if
                temp = alpha * a(i)
                c(i,:) = c(i,:) + temp * conjg(b(:,i))
            end do
        else
            ! Compute C = alpha * A * B + beta * C
            do i = 1, k
                if (beta == zero) then
                    c(i,:) = zero
                else if (beta /= one) then
                    c(i,:) = beta * c(i,:)
                end if
                temp = alpha * a(i)
                c(i,:) = c(i,:) + temp * b(i,:)
            end do
        end if

        ! Handle extra rows
        if (m > k) then
            if (beta == zero) then
                c(k+1:m,:) = zero
            else
                c(k+1:m,:) = beta * c(k+1:m,:)
            end if
        end if
    else
        if (opb == LA_TRANSPOSE) then
            ! Compute C = alpha * B**T * A + beta * C
            do i = 1, k
                if (beta == zero) then
                    c(:,i) = zero
                else if (beta /= one) then
                    c(:,i) = beta * c(:,i)
                end if
                temp = alpha * a(i)
                c(:,i) = c(:,i) + temp * b(i,:)
            end do
        else if (opb ==  LA_HERMITIAN_TRANSPOSE) then
            ! Compute C = alpha * B**H * A + beta * C
            do i = 1, k
                if (beta == zero) then
                    c(:,i) = zero
                else if (beta /= one) then
                    c(:,i) = beta * c(:,i)
                end if
                temp = alpha * a(i)
                c(:,i) = c(:,i) + temp * conjg(b(i,:))
            end do
        else
            ! Compute C = alpha * B * A + beta * C
            do i = 1, k
                if (beta == zero) then
                    c(:,i) = zero
                else if (beta /= one) then
                    c(:,i) = beta * c(:,i)
                end if
                temp = alpha * a(i)
                c(:,i) = c(:,i) + temp * b(:,i)
            end do
        end if

        ! Handle extra columns
        if (n > k) then
            if (beta == zero) then
                c(:,k+1:m) = zero
            else if (beta /= one) then
                c(:,k+1:m) = beta * c(:,k+1:m)
            end if
        end if
    end if
end subroutine

! ------------------------------------------------------------------------------
pure subroutine diag_mtx_mult_mtx2_cmplx(lside, alpha, a, b)
    !! Performs the matrix operation \(B = \alpha A B \) or \(B = \alpha B A \)
    !! where \(A\) is a diagonal matrix.
    logical, intent(in) :: lside
        !! A logical flag indicating if the diagonal matrix is on the left.
    complex(real64), intent(in) :: alpha
        !! The scalar \(\alpha\) to multiply the product of \(A\) and \(B\).
    complex(real64), intent(in), dimension(:) :: a
        !! The diagonal matrix \(A\) in the operation.
    complex(real64), intent(inout), dimension(:,:) :: b
        !! The matrix \(B\) in the operation.

    ! Parameters
    complex(real64), parameter :: zero = (0.0d0, 0.0d0)
    complex(real64), parameter :: one = (1.0d0, 0.0d0)

    ! Local Variables
    integer(int32) :: i, m, n, k
    complex(real64) :: temp

    ! Initialization
    m = size(b, 1)
    n = size(b, 2)
    k = size(a)

    ! Input Check
    if ((lside .and. k > m) .or. (.not.lside .and. k > n)) then
        ! ERROR: One of the input arrays is not sized correctly
        error stop 3
    end if

    ! Process
    if (lside) then
        ! Compute B = alpha * A * B
        do i = 1, k
            temp = alpha * a(i)
            b(i,:) = temp * b(i,:)
        end do
        if (m > k) b(k+1:m,:) = zero
    else
        ! Compute B = alpha * B * A
        do i = 1, k
            temp = alpha * a(i)
            b(:,i) = temp * b(:,i)
        end do
        if (n > k) b(:,k+1:n) = zero
    end if
end subroutine

! ------------------------------------------------------------------------------
pure subroutine diag_mtx_mult_mtx_mix(lside, opb, alpha, a, b, beta, c)
    !! Performs the matrix operation \(C = \alpha A op(B) + \beta C \) or
    !! \(C = \alpha op(B) A + \beta C \) where \(A\) is a diagonal matrix.
    logical, intent(in) :: lside
        !! A logical flag indicating if the diagonal matrix is on the left.
    integer(int32), intent(in) :: opb
        !! An integer flag indicating the operation to perform on matrix \(B\).
        !! Possible options are:
        !!
        !! - LA_NO_OPERATION: No operation is performed on matrix.
        !!
        !! - LA_TRANSPOSE: The transpose of matrix is used.
        !!
        !! - LA_HERMITIAN_TRANSPOSE: The Hermitian transpose of matrix is used.
    complex(real64), intent(in) :: alpha
        !! The scalar \(\alpha\) to multiply the product of \(A\) and \(B\).
    complex(real64), intent(in) :: beta
        !! The scalar \(\beta\) to multiply the matrix \(C\).
    real(real64), intent(in), dimension(:) :: a
        !! The diagonal matrix \(A\) in the operation.
    complex(real64), intent(in), dimension(:,:) :: b
        !! The matrix \(B\) in the operation.
    complex(real64), intent(inout), dimension(:,:) :: c
        !! The matrix \(C\) in the operation.

    ! Parameters
    complex(real64), parameter :: zero = (0.0d0, 0.0d0)
    complex(real64), parameter :: one = (1.0d0, 0.0d0)

    ! Local Variables
    integer(int32) :: i, m, n, k, nrowb, ncolb, flag
    complex(real64) :: temp

    ! Initialization
    m = size(c, 1)
    n = size(c, 2)
    k = size(a)
    nrowb = size(b, 1)
    ncolb = size(b, 2)

    ! Input Check
    flag = 0
    if (lside) then
        if (k > m) then
            flag = 4
        else
            if (opb == LA_TRANSPOSE .or. opb ==  LA_HERMITIAN_TRANSPOSE) then
                ! Compute C = alpha * A * B**T + beta * C
                if (nrowb /= n .or. ncolb < k) flag = 5
            else
                ! Compute C = alpha * A * B + beta * C
                if (nrowb < k .or. ncolb /= n) flag = 5
            end if
        end if
    else
        if (k > n) then
            flag = 4
        else
            if (opb == LA_TRANSPOSE .or. opb ==  LA_HERMITIAN_TRANSPOSE) then
                ! Compute C = alpha * B**T * A + beta * C
                if (ncolb /= m .or. nrowb < k) flag = 5
            else
                ! Compute C = alpha * B * A + beta * C
                if (nrowb /= m .or. ncolb < k) flag = 5
            end if
        end if
    end if
    if (flag /= 0) then
        ! ERROR: One of the input arrays is not sized correctly
        error stop flag
    end if

    ! Deal with ALPHA == 0
    if (alpha == 0) then
        if (beta == zero) then
            c = zero
        else if (beta /= one) then
            c = beta * c
        end if
        return
    end if

    ! Process
    if (lside) then
        if (opb == LA_TRANSPOSE) then
            ! Compute C = alpha * A * B**T + beta * C
            do i = 1, k
                if (beta == zero) then
                    c(i,:) = zero
                else if (beta /= one) then
                    c(i,:) = beta * c(i,:)
                end if
                temp = alpha * a(i)
                c(i,:) = c(i,:) + temp * b(:,i)
            end do
        else if (opb ==  LA_HERMITIAN_TRANSPOSE) then
            ! Compute C = alpha * A * B**H + beta * C
            do i = 1, k
                if (beta == zero) then
                    c(i,:) = zero
                else if (beta /= one) then
                    c(i,:) = beta * c(i,:)
                end if
                temp = alpha * a(i)
                c(i,:) = c(i,:) + temp * conjg(b(:,i))
            end do
        else
            ! Compute C = alpha * A * B + beta * C
            do i = 1, k
                if (beta == zero) then
                    c(i,:) = zero
                else if (beta /= one) then
                    c(i,:) = beta * c(i,:)
                end if
                temp = alpha * a(i)
                c(i,:) = c(i,:) + temp * b(i,:)
            end do
        end if

        ! Handle extra rows
        if (m > k) then
            if (beta == zero) then
                c(k+1:m,:) = zero
            else
                c(k+1:m,:) = beta * c(k+1:m,:)
            end if
        end if
    else
        if (opb == LA_TRANSPOSE) then
            ! Compute C = alpha * B**T * A + beta * C
            do i = 1, k
                if (beta == zero) then
                    c(:,i) = zero
                else if (beta /= one) then
                    c(:,i) = beta * c(:,i)
                end if
                temp = alpha * a(i)
                c(:,i) = c(:,i) + temp * b(i,:)
            end do
        else if (opb ==  LA_HERMITIAN_TRANSPOSE) then
            ! Compute C = alpha * B**H * A + beta * C
            do i = 1, k
                if (beta == zero) then
                    c(:,i) = zero
                else if (beta /= one) then
                    c(:,i) = beta * c(:,i)
                end if
                temp = alpha * a(i)
                c(:,i) = c(:,i) + temp * conjg(b(i,:))
            end do
        else
            ! Compute C = alpha * B * A + beta * C
            do i = 1, k
                if (beta == zero) then
                    c(:,i) = zero
                else if (beta /= one) then
                    c(:,i) = beta * c(:,i)
                end if
                temp = alpha * a(i)
                c(:,i) = c(:,i) + temp * b(:,i)
            end do
        end if

        ! Handle extra columns
        if (n > k) then
            if (beta == zero) then
                c(:,k+1:m) = zero
            else if (beta /= one) then
                c(:,k+1:m) = beta * c(:,k+1:m)
            end if
        end if
    end if
end subroutine

! ------------------------------------------------------------------------------
pure subroutine diag_mtx_mult_mtx2_mix(lside, alpha, a, b)
    !! Performs the matrix operation \(B = \alpha A B \) or \(B = \alpha B A \)
    !! where \(A\) is a diagonal matrix.
    logical, intent(in) :: lside
        !! A logical flag indicating if the diagonal matrix is on the left.
    complex(real64), intent(in) :: alpha
        !! The scalar \(\alpha\) to multiply the product of \(A\) and \(B\).
    real(real64), intent(in), dimension(:) :: a
        !! The diagonal matrix \(A\) in the operation.
    complex(real64), intent(inout), dimension(:,:) :: b
        !! The matrix \(B\) in the operation.

    ! Parameters
    complex(real64), parameter :: zero = (0.0d0, 0.0d0)
    complex(real64), parameter :: one = (1.0d0, 0.0d0)

    ! Local Variables
    integer(int32) :: i, m, n, k
    complex(real64) :: temp

    ! Initialization
    m = size(b, 1)
    n = size(b, 2)
    k = size(a)

    ! Input Check
    if ((lside .and. k > m) .or. (.not.lside .and. k > n)) then
        ! ERROR: One of the input arrays is not sized correctly
        error stop 3
    end if

    ! Process
    if (lside) then
        ! Compute B = alpha * A * B
        do i = 1, k
            temp = alpha * a(i)
            b(i,:) = temp * b(i,:)
        end do
        if (m > k) b(k+1:m,:) = zero
    else
        ! Compute B = alpha * B * A
        do i = 1, k
            temp = alpha * a(i)
            b(:,i) = temp * b(:,i)
        end do
        if (n > k) b(:,k+1:n) = zero
    end if
end subroutine

! ------------------------------------------------------------------------------
pure subroutine diag_mtx_sparse_mult(lside, alpha, a, b)
    !! Performs the matrix operation \(B = \alpha A B \) or \(B = \alpha B A \)
    !! where \(A\) is a diagonal matrix and \(B\) is a sparse matrix.
    logical, intent(in) :: lside
    real(real64), intent(in) :: alpha
    real(real64), intent(in), dimension(:) :: a
    class(csr_matrix), intent(inout) :: b

    ! Local Variables
    integer(int32) :: ii, k, k1, k2, nrow
    real(real64) :: scal
    
    ! Initialization
    nrow = size(b, 1)

    ! Input Check
    if (lside) then
        if (size(a) /= nrow) then
            error stop 3
        end if
    else
        if (size(a) /= size(b, 2)) then
            error stop 4
        end if
    end if

    ! Process
    if (lside) then
        ! Compute B = DIAG * B
        do ii = 1, nrow
            k1 = b%row_indices(ii)
            k2 = b%row_indices(ii+1) - 1
            if (alpha == 1.0d0) then
                scal = a(ii)
            else
                scal = alpha * a(ii)
            end if
            do k = k1, k2
                b%values(k) = b%values(k) * scal
            end do
        end do
    else
        ! Compute B = B * DIAG
        do ii = 1, nrow
            k1 = b%row_indices(ii)
            k2 = b%row_indices(ii+1) - 1
            if (alpha == 1.0d0) then
                do k = k1, k2
                    b%values(k) = b%values(k) * a(b%column_indices(k))
                end do
            else
                do k = k1, k2
                    b%values(k) = alpha * b%values(k) * a(b%column_indices(k))
                end do
            end if
        end do
    end if
end subroutine

! ******************************************************************************
! BASIC OPERATION ROUTINES
! ------------------------------------------------------------------------------
pure function trace_dbl(x) result(y)
    !! Computes the trace of a matrix.
    real(real64), intent(in), dimension(:,:) :: x
        !! The matrix.
    real(real64) :: y
        !! The trace of the matrix.

    ! Parameters
    real(real64), parameter :: zero = 0.0d0

    ! Local Variables
    integer(int32) :: i, m, n, mn

    ! Initialization
    y = zero
    m = size(x, 1)
    n = size(x, 2)
    mn = min(m, n)

    ! Process
    do i = 1, mn
        y = y + x(i,i)
    end do
end function

! ------------------------------------------------------------------------------
pure function trace_cmplx(x) result(y)
    !! Computes the trace of a matrix.
    complex(real64), intent(in), dimension(:,:) :: x
        !! The matrix.
    complex(real64) :: y
        !! The trace of the matrix.

    ! Parameters
    complex(real64), parameter :: zero = (0.0d0, 0.0d0)

    ! Local Variables
    integer(int32) :: i, m, n, mn

    ! Initialization
    y = zero
    m = size(x, 1)
    n = size(x, 2)
    mn = min(m, n)

    ! Process
    do i = 1, mn
        y = y + x(i,i)
    end do
end function

! ------------------------------------------------------------------------------
pure function mtx_rank_dbl(a, tol) result(rnk)
    !! Computes the rank of a matrix.
    real(real64), intent(in), dimension(:,:) :: a
        !! The matrix.
    real(real64), intent(in), optional :: tol
        !! An optional input, that if supplied, overrides the default
        !! tolerance on singular values such that singular values less than 
        !! this tolerance are treated as zero.  The default tolerance is:
        !! MAX(M, N) * EPS * MAX(S).  If the supplied value is less than the
        !! smallest value that causes an overflow if inverted, the tolerance
        !! reverts back to its default value, and the operation continues; 
        !! however, a warning message is issued.
    integer(int32) :: rnk
        !! The rank of the matrix.

    ! Local Variables
    integer(int32) :: i, m, n, mn, istat, lwork, flag
    real(real64), allocatable, dimension(:) :: s, w
    real(real64), allocatable, dimension(:,:) :: ac
    real(real64) :: t, tref, smlnum
    real(real64), dimension(1) :: dummy, temp

    ! Initialization
    m = size(a, 1)
    n = size(a, 2)
    mn = min(m, n)
    smlnum = DLAMCH('s')
    allocate(ac(m, n), source = a)
    rnk = 0

    ! Workspace Query
    call DGESVD('N', 'N', m, n, ac, m, dummy, dummy, m, dummy, n, temp, &
        -1, flag)
    lwork = int(temp(1), int32)

    ! Local Memory Allocation
    allocate(s(mn), w(lwork))

    ! Compute the singular values of A
    call DGESVD('N', 'N', m, n, ac, m, s, dummy, m, dummy, n, w, &
        lwork - mn, flag)
    if (flag > 0) then
        error stop LA_CONVERGENCE_ERROR
    end if

    ! Determine the threshold tolerance for the singular values such that
    ! singular values less than the threshold result in zero when inverted.
    tref = max(m, n) * epsilon(t) * s(1)
    if (present(tol)) then
        t = tol
    else
        t = tref
    end if
    if (t < smlnum) then
        ! ! The supplied tolerance is too small, simply fall back to the
        ! ! default, but issue a warning to the user
        ! t = tref
        ! call report_warning("mtx_rank", "The supplied tolerance was " // &
        !     "smaller than a value that would result in an overflow " // &
        !     "condition, or is negative; therefore, the tolerance has " // &
        !     "been reset to its default value.")
    end if

    ! Count the singular values that are larger than the tolerance value
    do i = 1, mn
        if (s(i) < t) exit
        rnk = rnk + 1
    end do
end function

! ------------------------------------------------------------------------------
pure function mtx_rank_cmplx(a, tol) result(rnk)
    !! Computes the rank of a matrix.
    complex(real64), intent(in), dimension(:,:) :: a
        !! The matrix.
    real(real64), intent(in), optional :: tol
        !! An optional input, that if supplied, overrides the default
        !! tolerance on singular values such that singular values less than 
        !! this tolerance are treated as zero.  The default tolerance is:
        !! MAX(M, N) * EPS * MAX(S).  If the supplied value is less than the
        !! smallest value that causes an overflow if inverted, the tolerance
        !! reverts back to its default value, and the operation continues; 
        !! however, a warning message is issued.
    integer(int32) :: rnk
        !! The rank of the matrix.

    ! External Function Interfaces
    interface
        pure function DLAMCH(cmach) result(x)
            use, intrinsic :: iso_fortran_env, only : real64
            character, intent(in) :: cmach
            real(real64) :: x
        end function
    end interface

    ! Local Variables
    integer(int32) :: i, m, n, mn, istat, lwork, flag, lrwork
    real(real64), allocatable, dimension(:) :: s, rw
    complex(real64), allocatable, dimension(:) :: w
    complex(real64), allocatable, dimension(:,:) :: ac
    real(real64) :: t, tref, smlnum
    real(real64), dimension(1) :: dummy
    complex(real64), dimension(1) :: cdummy, temp

    ! Initialization
    m = size(a, 1)
    n = size(a, 2)
    mn = min(m, n)
    lrwork = 6 * mn
    smlnum = DLAMCH('s')
    allocate(ac(m, n), source = a)
    rnk = 0

    ! Workspace Query
    call ZGESVD('N', 'N', m, n, ac, m, dummy, cdummy, m, cdummy, n, temp, &
        -1, dummy, flag)
    lwork = int(temp(1), int32)

    ! Local Memory Allocation
    allocate(s(mn), rw(lrwork), w(lwork))

    ! Compute the singular values of A
    call ZGESVD('N', 'N', m, n, ac, m, s, cdummy, m, cdummy, n, w, &
        lwork - mn, rw, flag)
    if (flag > 0) then
        error stop LA_CONVERGENCE_ERROR
    end if

    ! Determine the threshold tolerance for the singular values such that
    ! singular values less than the threshold result in zero when inverted.
    tref = max(m, n) * epsilon(t) * s(1)
    if (present(tol)) then
        t = tol
    else
        t = tref
    end if
    if (t < smlnum) then
        ! ! The supplied tolerance is too small, simply fall back to the
        ! ! default, but issue a warning to the user
        ! t = tref
        ! call report_warning("mtx_rank", "The supplied tolerance was " // &
        !     "smaller than a value that would result in an overflow " // &
        !     "condition, or is negative; therefore, the tolerance has " // &
        !     "been reset to its default value.")
    end if

    ! Count the singular values that are larger than the tolerance value
    do i = 1, mn
        if (s(i) < t) exit
        rnk = rnk + 1
    end do
end function

! ------------------------------------------------------------------------------
pure function det_dbl(a) result(x)
    !! Computes the determinant of a matrix.
    real(real64), intent(in), dimension(:,:) :: a
        !! The N-by-N matrix on which to operate.
    real(real64) :: x
        !! The determinant of the matrix.

    ! Parameters
    real(real64), parameter :: zero = 0.0d0
    real(real64), parameter :: one = 1.0d0
    real(real64), parameter :: ten = 1.0d1
    real(real64), parameter :: p1 = 1.0d-1

    ! Local Variables
    integer(int32) :: i, ep, n, flag
    real(real64) :: temp
    real(real64), allocatable, dimension(:,:) :: ac
    integer(int32), allocatable, dimension(:) :: ipvt

    ! Input Check
    if (size(a, 1) /= size(a, 2)) then
        error stop 1
    end if

    ! Initialization
    n = size(a, 1)
    allocate(ac(n, n), source = a)
    allocate(ipvt(n))
    x = zero

    ! Compute the LU factorization of A
    call DGETRF(n, n, ac, n, ipvt, flag)
    if (flag > 0) then
        ! A singular matrix has a determinant of zero
        x = zero
        return
    end if

    ! Compute the product of the diagonal of A
    temp = one
    ep = 0
    do i = 1, n
        if (ipvt(i) /= i) temp = -temp

        temp = ac(i,i) * temp
        if (temp == zero) then
            x = zero
            exit
        end if

        do while (abs(temp) < one)
            temp = ten * temp
            ep = ep - 1
        end do

        do while (abs(temp) > ten)
            temp = p1 * temp
            ep = ep + 1
        end do
    end do
    x = temp * ten**ep
end function

! ------------------------------------------------------------------------------
pure function det_cmplx(a) result(x)
    !! Computes the determinant of a matrix.
    complex(real64), intent(in), dimension(:,:) :: a
        !! The matrix on which to operate.
    complex(real64) :: x
        !! The determinant of the matrix.

    ! Parameters
    complex(real64), parameter :: zero = (0.0d0, 0.0d0)
    complex(real64), parameter :: one = (1.0d0, 0.0d0)
    complex(real64), parameter :: ten = (1.0d1, 0.0d0)
    complex(real64), parameter :: p1 = (1.0d-1, 0.0d0)
    real(real64), parameter :: real_one = 1.0d0
    real(real64), parameter :: real_ten = 1.0d1

    ! Local Variables
    integer(int32) :: i, ep, n, flag
    complex(real64) :: temp
    complex(real64), allocatable, dimension(:,:) :: ac
    integer(int32), allocatable, dimension(:) :: ipvt

    ! Input Check
    if (size(a, 1) /= size(a, 2)) then
        error stop 1
    end if

    ! Initialization
    n = size(a, 1)
    allocate(ac(n, n), source = a)
    allocate(ipvt(n))
    x = zero

    ! Compute the LU factorization of A
    call ZGETRF(n, n, ac, n, ipvt, flag)
    if (flag > 0) then
        ! A singular matrix has a determinant of zero
        x = zero
        return
    end if

    ! Compute the product of the diagonal of A
    temp = one
    ep = 0
    do i = 1, n
        if (ipvt(i) /= i) temp = -temp

        temp = ac(i,i) * temp
        if (temp == zero) then
            x = zero
            exit
        end if

        do while (abs(temp) < real_one)
            temp = ten * temp
            ep = ep - 1
        end do

        do while (abs(temp) > real_ten)
            temp = p1 * temp
            ep = ep + 1
        end do
    end do
    x = temp * ten**ep
end function

! ******************************************************************************
! ARRAY SWAPPING ROUTINE
! ------------------------------------------------------------------------------
pure subroutine swap_dbl(x, y)
    !! Swaps the contents of two arrays.
    real(real64), intent(inout), dimension(:) :: x
        !! On input, the first array to swap.  On output, the contents of the 
        !! first array are copied to the second array.
    real(real64), intent(inout), dimension(:) :: y
        !! On input, the second array to swap.  On output, the contents of the 
        !! second array are copied to the first array.

    ! Local Variables
    integer(int32) :: i, n
    real(real64) :: temp

    ! Initialization
    n = size(x)

    ! Input Check
    if (size(y) /= n) then
        error stop 2
    end if

    ! Process
    do i = 1, n
        temp = x(i)
        x(i) = y(i)
        y(i) = temp
    end do
end subroutine

! ------------------------------------------------------------------------------
pure subroutine swap_cmplx(x, y)
    !! Swaps the contents of two arrays.
    complex(real64), intent(inout), dimension(:) :: x
        !! On input, the first array to swap.  On output, the contents of the
        !! first array are copied to the second array.
    complex(real64), intent(inout), dimension(:) :: y
        !! On input, the second array to swap.  On output, the contents of the
        !! second array are copied to the first array.

    ! Local Variables
    integer(int32) :: i, n
    complex(real64) :: temp

    ! Initialization
    n = size(x)

    ! Input Check
    if (size(y) /= n) then
        error stop 2
    end if

    ! Process
    do i = 1, n
        temp = x(i)
        x(i) = y(i)
        y(i) = temp
    end do
end subroutine

! ******************************************************************************
! ARRAY MULTIPLICIATION ROUTINES
! ------------------------------------------------------------------------------
pure subroutine recip_mult_array_dbl(a, x)
    !! Computes the product of a scalar and a vector, where the scalar is 
    !! the reciprocal of the scalar A.
    real(real64), intent(in) :: a
        !! The scalar A, which is the reciprocal of the scalar to multiply by.
    real(real64), intent(inout), dimension(:) :: x
        !! On input, the vector to multiply.  On output, the product of the
        !! vector and the scalar reciprocal.

    ! External Function Interfaces
    interface
        pure function DLAMCH(cmach) result(x)
            use, intrinsic :: iso_fortran_env, only : real64
            character, intent(in) :: cmach
            real(real64) :: x
        end function
    end interface

    ! Parameters
    real(real64), parameter :: zero = 0.0d0
    real(real64), parameter :: one = 1.0d0
    real(real64), parameter :: twotho = 2.0d3

    ! Local Variables
    logical :: done
    real(real64) :: bignum, cden, cden1, cnum, cnum1, mul, smlnum

    ! Initialization
    smlnum = DLAMCH('s')
    bignum = one / smlnum
    if (log10(bignum) > twotho) then
        smlnum = sqrt(smlnum)
        bignum = sqrt(bignum)
    end if

    ! Initialize the denominator to A, and the numerator to ONE
    cden = a
    cnum = one

    ! Process
    do
        cden1 = cden * smlnum
        cnum1 = cnum / bignum
        if (abs(cden1) > abs(cnum) .and. cnum /= zero) then
            mul = smlnum
            done = .false.
            cden = cden1
        else if (abs(cnum1) > abs(cden)) then
            mul = bignum
            done = .false.
            cnum = cnum1
        else
            mul = cnum / cden
            done = .true.
        end if

        ! Scale the vector X by MUL
        x = mul * x

        ! Exit if done
        if (done) exit
    end do
end subroutine

! ******************************************************************************
! TRIANGULAR MATRIX MULTIPLICATION ROUTINES
! ------------------------------------------------------------------------------
pure subroutine tri_mtx_mult_dbl(upper, alpha, a, beta, b)
    !! Performs the matrix operation \(B = \alpha A^T A + \beta B\) or 
    !! \(B = \alpha A A^T + \beta B\) where \(A\) is a triangular matrix.
    logical, intent(in) :: upper
        !! A logical flag indicating whether the matrix A is upper triangular 
        !! (TRUE) or lower triangular (FALSE).
    real(real64), intent(in) :: alpha
        !! The scalar \(\alpha\) to multiply by.
    real(real64), intent(in) :: beta
        !! The scalar \(\beta\) to multiply by.
    real(real64), intent(in), dimension(:,:) :: a
        !! The triangular matrix \(A\) to multiply by.
    real(real64), intent(inout), dimension(:,:) :: b
        !! On input, the matrix \(B\) to multiply.  On output, the result of the
        !! operation.

    ! Parameters
    real(real64), parameter :: zero = 0.0d0

    ! Local Variables
    integer(int32) :: i, j, k, n
    real(real64) :: temp

    ! Initialization
    n = size(a, 1)

    ! Input Check
    if (size(a, 2) /= n) then
        error stop 2
    else if (size(b, 1) /= n .or. size(b, 2) /= n) then
        error stop 4
    end if

    ! Process
    if (upper) then
        ! Form: B = alpha * A**T * A + beta * B
        if (beta == zero) then
            do j = 1, n
                do i = 1, j
                    temp = zero
                    do k = 1, j
                        temp = temp + a(k,i) * a(k,j)
                    end do
                    temp = alpha * temp
                    b(i,j) = temp
                    if (i /= j) b(j,i) = temp
                end do
            end do
        else
            do j = 1, n
                do i = 1, j
                    temp = zero
                    do k = 1, j
                        temp = temp + a(k,i) * a(k,j)
                    end do
                    temp = alpha * temp
                    b(i,j) = temp + beta * b(i,j)
                    if (i /= j) b(j,i) = temp + beta * b(j,i)
                end do
            end do
        end if
    else
        ! Form: B = alpha * A * A**T + beta * B
        if (beta == zero) then
            do j = 1, n
                do i = j, n
                    temp = zero
                    do k = 1, j
                        temp = temp + a(i,k) * a(j,k)
                    end do
                    temp = alpha * temp
                    b(i,j) = temp
                    if (i /= j) b(j,i) = temp
                end do
            end do
        else
            do j = 1, n
                do i = j, n
                    temp = zero
                    do k = 1, j
                        temp = temp + a(i,k) * a(j,k)
                    end do
                    temp = alpha * temp
                    b(i,j) = temp + beta * b(i,j)
                    if (i /= j) b(j,i) = temp + beta * b(j,i)
                end do
            end do
        end if
    end if
end subroutine

! ------------------------------------------------------------------------------
pure subroutine tri_mtx_mult_cmplx(upper, alpha, a, beta, b)
    !! Performs the matrix operation \(B = \alpha A^T A + \beta B\) or
    !! \(B = \alpha A A^T + \beta B\) where \(A\) is a triangular matrix.
    logical, intent(in) :: upper
        !! A logical flag indicating whether the matrix A is upper triangular
        !! (TRUE) or lower triangular (FALSE).
    complex(real64), intent(in) :: alpha
    !! The scalar \(\alpha\) to multiply by.
    complex(real64), intent(in) :: beta
        !! The scalar \(\beta\) to multiply by.
    complex(real64), intent(in), dimension(:,:) :: a
        !! The triangular matrix \(A\) to multiply by.
    complex(real64), intent(inout), dimension(:,:) :: b
        !! On input, the matrix \(B\) to multiply.  On output, the result of the
        !! operation.

    ! Parameters
    complex(real64), parameter :: zero = (0.0d0, 0.0d0)

    ! Local Variables
    integer(int32) :: i, j, k, n
    complex(real64) :: temp

    ! Initialization
    n = size(a, 1)

    ! Input Check
    if (size(a, 2) /= n) then
        error stop 2
    else if (size(b, 1) /= n .or. size(b, 2) /= n) then
        error stop 4
    end if

    ! Process
    if (upper) then
        ! Form: B = alpha * A**T * A + beta * B
        if (beta == zero) then
            do j = 1, n
                do i = 1, j
                    temp = zero
                    do k = 1, j
                        temp = temp + a(k,i) * a(k,j)
                    end do
                    temp = alpha * temp
                    b(i,j) = temp
                    if (i /= j) b(j,i) = temp
                end do
            end do
        else
            do j = 1, n
                do i = 1, j
                    temp = zero
                    do k = 1, j
                        temp = temp + a(k,i) * a(k,j)
                    end do
                    temp = alpha * temp
                    b(i,j) = temp + beta * b(i,j)
                    if (i /= j) b(j,i) = temp + beta * b(j,i)
                end do
            end do
        end if
    else
        ! Form: B = alpha * A * A**T + beta * B
        if (beta == zero) then
            do j = 1, n
                do i = j, n
                    temp = zero
                    do k = 1, j
                        temp = temp + a(i,k) * a(j,k)
                    end do
                    temp = alpha * temp
                    b(i,j) = temp
                    if (i /= j) b(j,i) = temp
                end do
            end do
        else
            do j = 1, n
                do i = j, n
                    temp = zero
                    do k = 1, j
                        temp = temp + a(i,k) * a(j,k)
                    end do
                    temp = alpha * temp
                    b(i,j) = temp + beta * b(i,j)
                    if (i /= j) b(j,i) = temp + beta * b(j,i)
                end do
            end do
        end if
    end if
end subroutine

! ******************************************************************************
! BANDED MATRIX MULTIPLICATION ROUTINES
! ------------------------------------------------------------------------------
pure subroutine band_mtx_vec_mult_dbl(trans, kl, ku, alpha, a, x, beta, &
    y)
    !! Performs the matrix operation \(y = \alpha A x + \beta y\) or
    !! \(y = \alpha A^T A + \beta y\) where \(A\) is a banded matrix.
    !!
    !! The banded matrix is stored in a compressed form supplied column by 
    !! column.  The following code segment transfers between a full matrix
    !! to the bonded matrix storage scheme.
    !! \code{fortran}
    !! do j = 1, n
    !!    k = ku + 1 - j
    !!    do i = max(1, j - ku), min(n, j + kl)
    !!       a(k + i, j) = matrix(i, j)
    !!    end do
    !! end do
    !! \endcode
    logical, intent(in) :: trans
    !! A logical flag indicating whether to perform the operation
    !! \(y = \alpha A x + \beta y\) (FALSE) or \(y = \alpha A^T x + \beta y\)
    !! (TRUE).
    integer(int32), intent(in) :: kl
        !! The number of subdiagonals in the banded matrix \(A\).
    integer(int32), intent(in) :: ku
        !! The number of superdiagonals in the banded matrix \(A\).
    real(real64), intent(in) :: alpha
        !! The scalar \(\alpha\) to multiply by.
    real(real64), intent(in) :: beta
        !! The scalar \(\beta\) to multiply by.
    real(real64), intent(in), dimension(:,:) :: a
        !! The banded matrix \(A\) to multiply by.
    real(real64), intent(in), dimension(:) :: x
        !! The vector \(x\) to multiply by.
    real(real64), intent(inout), dimension(:) :: y
        !! On input, the vector \(y\) to multiply.  On output, the result of the
        !! operation.

    ! Local Variables
    integer(int32) :: m, n
    
    ! Initialization
    if (trans) then
        m = size(x)
        n = size(y)
    else
        m = size(y)
        n = size(x)
    end if

    ! Input Checking
    if (kl < 0) error stop 2
    if (ku < 0) error stop 3
    if (size(a, 1) /= kl + ku + 1) error stop 5
    if (size(a, 2) /= n) error stop 5

    ! Process
    if (trans) then
        call DGBMV("T", m, n, kl, ku, alpha, a, size(a, 1), x, 1, beta, y, 1)
    else
        call DGBMV("N", m, n, kl, ku, alpha, a, size(a, 1), x, 1, beta, y, 1)
    end if
end subroutine

! ------------------------------------------------------------------------------
pure subroutine band_mtx_vec_mult_cmplx(trans, kl, ku, alpha, a, x, &
    beta, y)
    !! Performs the matrix operation \(y = \alpha op(A) x + \beta y\)  where 
    !! \(A\) is a banded matrix.
    !!
    !! The banded matrix is stored in a compressed form supplied column by 
    !! column.  The following code segment transfers between a full matrix
    !! to the bonded matrix storage scheme.
    !! \code{fortran}
    !! do j = 1, n
    !!    k = ku + 1 - j
    !!    do i = max(1, j - ku), min(n, j + kl)
    !!       a(k + i, j) = matrix(i, j)
    !!    end do
    !! end do
    !! \endcode
    integer(int32), intent(in) :: trans
        !! An integer flag indicating the operation to perform on matrix \(A\).
        !! Possible options are:
        !!
        !! - LA_NO_OPERATION: No operation is performed on matrix.
        !!
        !! - LA_TRANSPOSE: The transpose of matrix is used.
        !!
        !! - LA_HERMITIAN_TRANSPOSE: The Hermitian transpose of matrix is used.
    integer(int32), intent(in) :: kl
        !! The number of subdiagonals in the banded matrix \(A\).
    integer(int32), intent(in) :: ku
        !! The number of superdiagonals in the banded matrix \(A\).
    complex(real64), intent(in) :: alpha
        !! The scalar \(\alpha\) to multiply by.
    complex(real64), intent(in) :: beta
        !! The scalar \(\beta\) to multiply by.
    complex(real64), intent(in), dimension(:,:) :: a
        !! The banded matrix \(A\) to multiply by.
    complex(real64), intent(in), dimension(:) :: x
        !! The vector \(x\) to multiply by.
    complex(real64), intent(inout), dimension(:) :: y
        !! On input, the vector \(y\) to multiply.  On output, the result of the
        !! operation.

    ! Local Variables
    character :: op
    logical :: trns
    integer(int32) :: m, n
    
    ! Initialization
    if (trans == LA_TRANSPOSE) then
        op = "T"
        trns = .true.
    else if (trans == LA_HERMITIAN_TRANSPOSE) then
        op = "C"
        trns = .true.
    else
        op = "N"
        trns = .false.
    end if
    if (trns) then
        m = size(x)
        n = size(y)
    else
        m = size(y)
        n = size(x)
    end if

    ! Input Checking
    if (kl < 0) error stop 2
    if (ku < 0) error stop 3
    if (size(a, 1) /= kl + ku + 1) error stop 5
    if (size(a, 2) /= n) error stop 5

    ! Process
    call ZGBMV(op, m, n, kl, ku, alpha, a, size(a, 1), x, 1, beta, y, 1)
end subroutine

! ------------------------------------------------------------------------------
pure subroutine band_to_full_mtx_dbl(kl, ku, b, f)
    !! Converts a banded matrix to a full matrix.
    !!
    !! The banded matrix is stored in a compressed form supplied column by 
    !! column.  The following code segment transfers between a full matrix
    !! to the bonded matrix storage scheme.
    !! \code{fortran}
    !! do j = 1, n
    !!    k = ku + 1 - j
    !!    do i = max(1, j - ku), min(n, j + kl)
    !!       a(k + i, j) = matrix(i, j)
    !!    end do
    !! end do
    !! \endcode
    integer(int32), intent(in) :: kl
        !! The number of subdiagonals in the banded matrix.
    integer(int32), intent(in) :: ku
        !! The number of superdiagonals in the banded matrix.
    real(real64), intent(in), dimension(:,:) :: b
        !! The banded matrix to convert.
    real(real64), intent(out), dimension(:,:) :: f
        !! The full matrix to store the result in.

    ! Parameters
    real(real64), parameter :: zero = 0.0d0

    ! Local Variables
    integer(int32) :: i, j, k, m, n, i1, i2
    
    ! Initialization
    m = size(f, 1)
    n = size(f, 2)

    ! Input Check
    if (kl < 0) error stop 1
    if (ku < 0) error stop 2
    if (size(b, 2) /= n) error stop 4
    if (size(b, 1) /= kl + ku + 1) error stop 4

    ! Process
    do j = 1, n
        k = ku + 1 - j
        i1 = max(1, j - ku)
        i2 = min(m, j + kl)
        do i = 1, i1 - 1
            f(i,j) = zero
        end do
        do i = i1, i2
            f(i,j) = b(k+i,j)
        end do
        do i = i2 + 1, m
            f(i,j) = zero
        end do
    end do
end subroutine

! ------------------------------------------------------------------------------
pure subroutine band_to_full_mtx_cmplx(kl, ku, b, f)
    !! Converts a banded matrix to a full matrix.
    !!
    !! The banded matrix is stored in a compressed form supplied column by 
    !! column.  The following code segment transfers between a full matrix
    !! to the bonded matrix storage scheme.
    !! \code{fortran}
    !! do j = 1, n
    !!    k = ku + 1 - j
    !!    do i = max(1, j - ku), min(n, j + kl)
    !!       a(k + i, j) = matrix(i, j)
    !!    end do
    !! end do
    !! \endcode
    integer(int32), intent(in) :: kl
        !! The number of subdiagonals in the banded matrix.
    integer(int32), intent(in) :: ku
        !! The number of superdiagonals in the banded matrix.
    complex(real64), intent(in), dimension(:,:) :: b
        !! The banded matrix to convert.
    complex(real64), intent(out), dimension(:,:) :: f
        !! The full matrix to store the result in.

    ! Parameters
    complex(real64), parameter :: zero = (0.0d0, 0.0d0)

    ! Local Variables
    integer(int32) :: i, j, k, m, n, i1, i2
    
    ! Initialization
    m = size(f, 1)
    n = size(f, 2)

    ! Input Check
    if (kl < 0) error stop 1
    if (ku < 0) error stop 2
    if (size(b, 2) /= n) error stop 4
    if (size(b, 1) /= kl + ku + 1) error stop 4

    ! Process
    do j = 1, n
        k = ku + 1 - j
        i1 = max(1, j - ku)
        i2 = min(m, j + kl)
        do i = 1, i1 - 1
            f(i,j) = zero
        end do
        do i = i1, i2
            f(i,j) = b(k+i,j)
        end do
        do i = i2 + 1, m
            f(i,j) = zero
        end do
    end do
end subroutine

! ------------------------------------------------------------------------------
pure subroutine band_diag_mtx_mult_dbl(left, m, kl, ku, alpha, a, b)
    !! Performs the matrix operation \(A = \alpha A B\) or \(A = \alpha B A\) 
    !! where \(A\) is a banded matrix and \(B\) is a diagonal matrix.
    !!
    !! The banded matrix is stored in a compressed form supplied column by 
    !! column.  The following code segment transfers between a full matrix
    !! to the bonded matrix storage scheme.
    !! \code{fortran}
    !! do j = 1, n
    !!    k = ku + 1 - j
    !!    do i = max(1, j - ku), min(n, j + kl)
    !!       a(k + i, j) = matrix(i, j)
    !!    end do
    !! end do
    !! \endcode
    logical, intent(in) :: left
        !! A logical flag indicating whether to perform the operation
        !! \(A = \alpha A B\) (TRUE) or \(A = \alpha B A\) (FALSE).
    integer(int32), intent(in) :: m
        !! The number of rows in the banded matrix \(A\).
    integer(int32), intent(in) :: kl
        !! The number of subdiagonals in the banded matrix.
    integer(int32), intent(in) :: ku
        !! The number of superdiagonals in the banded matrix.
    real(real64), intent(in) :: alpha
        !! The scalar \(\alpha\) to multiply by.
    real(real64), intent(inout), dimension(:,:) :: a
        !! The banded matrix to multiply.
    real(real64), intent(in), dimension(:) :: b
        !! The diagonal matrix to multiply by.

    ! Parameters
    real(real64), parameter :: one = 1.0d0

    ! Local Variables
    integer(int32) :: i, i1, i2, j, k, n
    real(real64) :: temp
    
    ! Initialization
    n = size(a, 2)

    ! Input Checking
    if (kl < 0) error stop 3
    if (ku < 0) error stop 4
    if (left) then
        if (size(b) /= n) error stop 7
    else
        if (size(b) < m) error stop 7
    end if

    ! Process
    if (left) then
        ! Compute A = A * B
        do j = 1, n
            k = ku + 1 - j
            i1 = max(1, j - ku) + k
            i2 = min(m, j + kl) + k
            if (alpha == one) then
                temp = b(j)
            else
                temp = alpha * b(j)
            end if
            do i = i1, i2
                a(i,j) = a(i,j) * temp
            end do
        end do
    else
        ! Compute A = B * A
        do j = 1, n
            k = ku + 1 - j
            i1 = max(1, j - ku)
            i2 = min(m, j + kl)
            if (alpha == 1.0d0) then
                do i = i1, i2
                    a(i+k,j) = a(i+k,j) * b(i)
                end do
            else
                do i = i1, i2
                    a(i+k,j) = alpha * a(i+k,j) * b(i)
                end do
            end if
        end do
    end if
end subroutine

! ------------------------------------------------------------------------------
pure subroutine band_diag_mtx_mult_cmplx(left, m, kl, ku, alpha, a, b)
    !! Performs the matrix operation \(A = \alpha A B\) or \(A = \alpha B A\) 
    !! where \(A\) is a banded matrix and \(B\) is a diagonal matrix.
    !!
    !! The banded matrix is stored in a compressed form supplied column by 
    !! column.  The following code segment transfers between a full matrix
    !! to the bonded matrix storage scheme.
    !! \code{fortran}
    !! do j = 1, n
    !!    k = ku + 1 - j
    !!    do i = max(1, j - ku), min(n, j + kl)
    !!       a(k + i, j) = matrix(i, j)
    !!    end do
    !! end do
    !! \endcode
    logical, intent(in) :: left
        !! A logical flag indicating whether to perform the operation
        !! \(A = \alpha A B\) (TRUE) or \(A = \alpha B A\) (FALSE).
    integer(int32), intent(in) :: m
        !! The number of rows in the banded matrix \(A\).
    integer(int32), intent(in) :: kl
        !! The number of subdiagonals in the banded matrix.
    integer(int32), intent(in) :: ku
        !! The number of superdiagonals in the banded matrix.
    complex(real64), intent(in) :: alpha
        !! The scalar \(\alpha\) to multiply by.
    complex(real64), intent(inout), dimension(:,:) :: a
        !! The banded matrix to multiply.
    complex(real64), intent(in), dimension(:) :: b
        !! The diagonal matrix to multiply by.

    ! Parameters
    complex(real64), parameter :: one = (1.0d0, 0.0d0)

    ! Local Variables
    integer(int32) :: i, i1, i2, j, k, n
    complex(real64) :: temp
    
    ! Initialization
    n = size(a, 2)

    ! Input Checking
    if (kl < 0) error stop 3
    if (ku < 0) error stop 4
    if (left) then
        if (size(b) /= n) error stop 7
    else
        if (size(b) < m) error stop 7
    end if

    ! Process
    if (left) then
        ! Compute A = A * B
        do j = 1, n
            k = ku + 1 - j
            i1 = max(1, j - ku) + k
            i2 = min(m, j + kl) + k
            if (alpha == one) then
                temp = b(j)
            else
                temp = alpha * b(j)
            end if
            do i = i1, i2
                a(i,j) = a(i,j) * temp
            end do
        end do
    else
        ! Compute A = B * A
        do j = 1, n
            k = ku + 1 - j
            i1 = max(1, j - ku)
            i2 = min(m, j + kl)
            if (alpha == 1.0d0) then
                do i = i1, i2
                    a(i+k,j) = a(i+k,j) * b(i)
                end do
            else
                do i = i1, i2
                    a(i+k,j) = alpha * a(i+k,j) * b(i)
                end do
            end if
        end do
    end if
end subroutine

! ******************************************************************************
! MATRIX CONVERSION ROUTINES
! ------------------------------------------------------------------------------
pure function banded_to_dense_dbl(m, kl, ku, a) result(x)
    !! Converts a banded matrix to a dense matrix.
    !!
    !! The banded matrix is stored in a compressed form supplied column by 
    !! column.  The following code segment transfers between a full matrix
    !! to the bonded matrix storage scheme.
    !! \code{fortran}
    !! do j = 1, n
    !!    k = ku + 1 - j
    !!    do i = max(1, j - ku), min(n, j + kl)
    !!       a(k + i, j) = matrix(i, j)
    !!    end do
    !! end do
    !! \endcode
    integer(int32), intent(in) :: m
        !! The number of rows in the M-by-N dense matrix.
    integer(int32), intent(in) :: kl
        !! The number of subdiagonals.  Must be at least 0.
    integer(int32), intent(in) :: ku
        !! The number of superdiagonals.  Must be at least 0.
    real(real64), intent(in), dimension(:,:) :: a
        !! The (KL+KU+1)-by-N banded matrix.
    real(real64), allocatable, dimension(:,:) :: x
        !! The M-by-N dense matrix.

    ! Parameters
    real(real64), parameter :: zero = 0.0d0

    ! Local Variables
    integer(int32) :: i, j, k, n, i1, i2
    
    ! Initialization
    n = size(a, 2)

    ! Input Checking
    if (kl < 0) then
        error stop 2
    end if
    if (ku < 0) then
        error stop 3
    end if
    if (size(a, 1) /= kl + ku + 1) then
        error stop 4
    end if

    ! Process
    allocate(x(m, n), source = zero)
    do j = 1, n
        k = ku + 1 - j
        i1 = max(1, j - ku)
        i2 = min(m, j + kl)
        do i = i1, i2
            x(i, j) = a(k + i, j)
        end do
    end do
end function

! ------------------------------------------------------------------------------
pure function banded_to_dense_cmplx(m, kl, ku, a) result(x)
    !! Converts a banded matrix to a dense matrix.
    !!
    !! The banded matrix is stored in a compressed form supplied column by 
    !! column.  The following code segment transfers between a full matrix
    !! to the bonded matrix storage scheme.
    !! \code{fortran}
    !! do j = 1, n
    !!    k = ku + 1 - j
    !!    do i = max(1, j - ku), min(n, j + kl)
    !!       a(k + i, j) = matrix(i, j)
    !!    end do
    !! end do
    !! \endcode
    integer(int32), intent(in) :: m
        !! The M-by-N dense matrix.
    integer(int32), intent(in) :: kl
        !! The number of subdiagonals.  Must be at least 0.
    integer(int32), intent(in) :: ku
        !! The number of superdiagonals.  Must be at least 0.
    complex(real64), intent(in), dimension(:,:) :: a
        !! The (KL+KU+1)-by-N banded matrix.
    complex(real64), allocatable, dimension(:,:) :: x
        !! The M-by-N dense matrix.

    ! Parameters
    complex(real64), parameter :: zero = (0.0d0, 0.0d0)

    ! Local Variables
    integer(int32) :: i, j, k, n, i1, i2
    
    ! Initialization
    n = size(a, 2)

    ! Input Checking
    if (kl < 0) then
        error stop 2
    end if
    if (ku < 0) then
        error stop 3
    end if
    if (size(a, 1) /= kl + ku + 1) then
        error stop 4
    end if

    ! Process
    allocate(x(m, n), source = zero)
    do j = 1, n
        k = ku + 1 - j
        i1 = max(1, j - ku)
        i2 = min(m, j + kl)
        do i = i1, i2
            x(i, j) = a(k + i, j)
        end do
    end do
end function

! ------------------------------------------------------------------------------
pure function dense_to_banded_dbl(a, kl, ku) result(x)
    !! Converts a banded matrix stored in dense format to a compressed form.
    !!
    !! The banded matrix is stored in a compressed form supplied column by 
    !! column.  The following code segment transfers between a full matrix
    !! to the bonded matrix storage scheme.
    !! \code{fortran}
    !! do j = 1, n
    !!    k = ku + 1 - j
    !!    do i = max(1, j - ku), min(n, j + kl)
    !!       a(k + i, j) = matrix(i, j)
    !!    end do
    !! end do
    !! \endcode
    real(real64), intent(in), dimension(:,:) :: a
        !! The matrix to convert.
    integer(int32), intent(in) :: kl
        !! The number of subdiagonals.  Must be at least 0.
    integer(int32), intent(in) :: ku
        !! The number of superdiagonals.  Must be at least 0.
    real(real64), allocatable, dimension(:,:) :: x
        !! The (KL+KU+1)-by-N banded matrix.

    ! Local Variables
    integer(int32) :: i, j, k, m, n, mm, flag
    
    ! Initialization
    m = size(a, 1)
    n = size(a, 2)
    mm = kl + ku + 1

    ! Input Check
    if (kl < 0) then
        error stop 2
    end if
    if (ku < 0) then
        error stop 3
    end if

    ! Process
    allocate(x(mm, n), source = 0.0d0)
    do j = 1, n
        k = ku + 1 - j
        do i = max(1, j - ku), min(m, j + kl)
            x(k + i, j) = a(i,j)
        end do
    end do
end function

! ------------------------------------------------------------------------------
pure function dense_to_banded_cmplx(a, kl, ku) result(x)
    !! Converts a banded matrix stored in dense format to a compressed form.
    !!
    !! The banded matrix is stored in a compressed form supplied column by 
    !! column.  The following code segment transfers between a full matrix
    !! to the bonded matrix storage scheme.
    !! \code{fortran}
    !! do j = 1, n
    !!    k = ku + 1 - j
    !!    do i = max(1, j - ku), min(n, j + kl)
    !!       a(k + i, j) = matrix(i, j)
    !!    end do
    !! end do
    !! \endcode
    complex(real64), intent(in), dimension(:,:) :: a
        !! The matrix to convert.
    integer(int32), intent(in) :: kl
        !! The number of subdiagonals.  Must be at least 0.
    integer(int32), intent(in) :: ku
        !! The number of superdiagonals.  Must be at least 0.
    complex(real64), allocatable, dimension(:,:) :: x
        !! The (KL+KU+1)-by-N banded matrix.

    ! Local Variables
    integer(int32) :: i, j, k, m, n, mm, flag
    
    ! Initialization
    m = size(a, 1)
    n = size(a, 2)
    mm = kl + ku + 1

    ! Input Check
    if (kl < 0) then
        error stop 2
    end if
    if (ku < 0) then
        error stop 3
    end if

    ! Process
    allocate(x(mm, n), source = (0.0d0, 0.0d0))
    do j = 1, n
        k = ku + 1 - j
        do i = max(1, j - ku), min(m, j + kl)
            x(k + i, j) = a(i,j)
        end do
    end do
end function

! ------------------------------------------------------------------------------
pure function extract_diagonal_dbl(a) result(diag)
    !! Extracts the diagonal of a matrix.
    real(real64), intent(in), dimension(:,:) :: a
        !! The M-by-N matrix.
    real(real64), allocatable, dimension(:) :: diag
        !! The MIN(M, N) element array for the diagonal elements.

    ! Local Variables
    integer(int32) :: i, m, n, mn
    
    ! Initialization
    m = size(a, 1)
    n = size(a, 2)
    mn = min(m, n)
    allocate(diag(mn))

    ! Process
    do i = 1, mn
        diag(i) = a(i,i)
    end do
end function

! ------------------------------------------------------------------------------
pure function extract_diagonal_cmplx(a) result(diag)
    !! Extracts the diagonal of a matrix.
    complex(real64), intent(in), dimension(:,:) :: a
        !! The M-by-N matrix.
    complex(real64), allocatable, dimension(:) :: diag
        !! The MIN(M, N) element array for the diagonal elements.

    ! Local Variables
    integer(int32) :: i, m, n, mn
    
    ! Initialization
    m = size(a, 1)
    n = size(a, 2)
    mn = min(m, n)
    allocate(diag(mn))

    ! Process
    do i = 1, mn
        diag(i) = a(i,i)
    end do
end function

! ------------------------------------------------------------------------------
pure function extract_diagonal_csr(a) result(diag)
    !! Extracts the diagonal of a matrix.
    class(csr_matrix), intent(in) :: a
        !! The M-by-N matrix.
    real(real64), allocatable, dimension(:) :: diag
        !! The MIN(M, N) element array for the diagonal elements.

    ! Local Variables
    integer(int32) :: i, m, n, mn
    
    ! Initialization
    m = size(a, 1)
    n = size(a, 2)
    mn = min(m, n)
    allocate(diag(mn))

    ! Process
    diag = a%extract_diagonal()
end function

! ------------------------------------------------------------------------------
pure function extract_upper_triangular_dbl(x) result(rst)
    !! Extracts the upper triangular portion of a matrix.
    real(real64), intent(in), dimension(:,:) :: x
        !! The M-by-N matrix.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The K-by-N upper triangular matrix where K = MIN(M, N).

    ! Local Variables
    integer(int32) :: j, m, n, k, kj

    ! Process
    m = size(x, 1)
    n = size(x, 2)
    k = min(m, n)
    allocate(rst(k, n))
    do j = 1, n
        kj = min(k, j)
        rst(1:kj,j) = x(1:kj,j)
        if (j < k) rst(j+1:k,j) = 0.0d0
    end do
end function

! ------------------------------------------------------------------------------
pure function extract_upper_triangular_cmplx(x) result(rst)
    !! Extracts the upper triangular portion of a matrix.
    complex(real64), intent(in), dimension(:,:) :: x
        !! The M-by-N matrix.
    complex(real64), allocatable, dimension(:,:) :: rst
        !! The K-by-N upper triangular matrix where K = MIN(M, N).

    ! Local Variables
    integer(int32) :: j, m, n, k, kj

    ! Process
    m = size(x, 1)
    n = size(x, 2)
    k = min(m, n)
    allocate(rst(k, n))
    do j = 1, n
        kj = min(k, j)
        rst(1:kj,j) = x(1:kj,j)
        if (j < k) rst(j+1:,j) = (0.0d0, 0.0d0)
    end do
end function

! ------------------------------------------------------------------------------
pure function extract_lower_triangular_dbl(x) result(rst)
    !! Extracts the lower triangular portion of a matrix.
    real(real64), intent(in), dimension(:,:) :: x
        !! The M-by-N matrix.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The M-by-K lower triangular matrix where K = MIN(M, N).

    ! Local Variables
    integer(int32) :: j, m, n, k

    ! Process
    m = size(x, 1)
    n = size(x, 2)
    k = min(m, n)
    allocate(rst(m, k))
    do j = 1, k
        if (j > 1) rst(1:j-1,j) = 0.0d0
        rst(j:,j) = x(j:,j)
    end do
end function

! ------------------------------------------------------------------------------
pure function extract_lower_triangular_cmplx(x) result(rst)
    !! Extracts the lower triangular portion of a matrix.
    complex(real64), intent(in), dimension(:,:) :: x
        !! The M-by-N matrix.
    complex(real64), allocatable, dimension(:,:) :: rst
        !! The M-by-K lower triangular matrix where K = MIN(M, N).

    ! Local Variables
    integer(int32) :: j, m, n, k

    ! Process
    m = size(x, 1)
    n = size(x, 2)
    k = min(m, n)
    allocate(rst(m, k))
    do j = 1, k
        if (j > 1) rst(1:j-1,j) = (0.0d0, 0.0d0)
        rst(j:,j) = x(j:,j)
    end do
end function

! ******************************************************************************
! MISC. OPERATIONS
! ------------------------------------------------------------------------------
pure function identity(n) result(rst)
    !! Constructs an N-by-N identity matrix.
    integer(int32), intent(in) :: n
        !! The size of the matrix.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The resulting matrix.

    ! Local Variables
    integer(int32) :: i

    ! Process
    if (n < 1) return
    allocate(rst(n, n), source = 0.0d0)
    do i = 1, n
        rst(i,i) = 1.0d0
    end do
end function

! ------------------------------------------------------------------------------
end module
