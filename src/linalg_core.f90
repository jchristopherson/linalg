! linalg_core.f90


!> @mainpage
!!
!! @section intro_sec Introduction
!! LINALG is a linear algebra library that provides a user-friendly interface 
!! to several BLAS and LAPACK routines.
!!
!! @author Jason Christopherson
!! @version 1.0


!> @brief \b linalg_core
!!
!! @par Purpose
!! Provides common "core" linear algebra routines.
module linalg_core
    use ferror, only : errors
    use lapack
    use linalg_constants
    implicit none
    private
    public :: solve_triangular_system
    public :: mtx_mult
    public :: rank1_update
    public :: diag_mtx_mult
    public :: trace
    public :: mtx_rank
    public :: det
    public :: swap
    public :: recip_mult_array

! ******************************************************************************
! INTERFACES
! ------------------------------------------------------------------------------
    !> @brief Solves a triangular system of equations.
    interface solve_triangular_system
        module procedure :: solve_tri_mtx
        module procedure :: solve_tri_vec
    end interface

! ------------------------------------------------------------------------------
    !> @brief Performs the matrix operation: 
    !!  C = alpha * op(A) * op(B) + beta * C.
    interface mtx_mult
        module procedure :: mtx_mult_mtx
        module procedure :: mtx_mult_vec
    end interface

! ------------------------------------------------------------------------------
!> @brief Multiplies a diagonal matrix with another matrix or array.
interface diag_mtx_mult
    module procedure :: diag_mtx_mult_mtx
    module procedure :: diag_mtx_mult_mtx2
    module procedure :: diag_mtx_mult_mtx3
    module procedure :: diag_mtx_mult_mtx4
end interface

contains
! ******************************************************************************
! TRIANGULAR MATRIX SOLUTION ROUTINES
! ------------------------------------------------------------------------------
    !> @brief Solves one of the matrix equations: op(A) * X = alpha * B, or
    !! X * op(A) = alpha * B, where A is a triangular matrix.
    !!
    !! @param[in] lside Set to true to solve op(A) * X = alpha * B; else, set to
    !!  false to solve X * op(A) = alpha * B.
    !! @param[in] upper Set to true if A is an upper triangular matrix; else,
    !!  set to false if A is a lower triangular matrix.
    !! @param[in] trans Set to true if op(A) = A**T; else, set to false if
    !!  op(A) = A.
    !! @param[in] nounit Set to true if A is not a unit-diagonal matrix (ones on
    !!  every diagonal element); else, set to false if A is a unit-diagonal
    !!  matrix.
    !! @param[in] alpha The scalar multiplier to B.
    !! @param[in] a If @p lside is true, the M-by-M triangular matrix on which
    !!  to operate; else, if @p lside is false, the N-by-N triangular matrix on
    !!  which to operate.
    !! @param[in,out] b On input, the M-by-N right-hand-side.  On output, the
    !!  M-by-N solution.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - LA_ARRAY_SIZE_ERROR: Occurs if @p a is not square, or if the sizes of
    !!      @p a and @p b are not compatible.
    !!
    !! @par Usage
    !! To solve a triangular system of N equations of N unknowns A*X = B, where
    !! A is an N-by-N upper triangular matrix, and B and X are N-by-NRHS
    !! matrices, the following code will suffice.
    !!
    !! @code{.f90}
    !! ! Solve the system: A*X = B, where A is an upper triangular N-by-N
    !! ! matrix, and B and X are N-by-NRHS in size.
    !!
    !! ! Variables
    !! integer(i32) :: info
    !! real(dp), dimension(n, n) :: a
    !! real(dp), dimension(n, nrhs) :: b
    !!
    !! ! Initialize A and B...
    !!
    !! ! Solve A*X = B for X - Note: X overwrites B.
    !! call solve_triangular_system(.true., .true., .false., .true., &
    !!      1.0d0, a, b)
    !! @endcode
    !!
    !! @par Notes
    !! This routine is based upon the BLAS routine DTRSM.
    subroutine solve_tri_mtx(lside, upper, trans, nounit, alpha, a, b, err)
        ! Arguments
        logical, intent(in) :: lside, upper, trans, nounit
        real(dp), intent(in) :: alpha
        real(dp), intent(in), dimension(:,:) :: a
        real(dp), intent(inout), dimension(:,:) :: b
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        character :: side, uplo, transa, diag
        real(dp), parameter :: zero = 0.0d0
        real(dp), parameter :: one = 1.0d0

        ! Local Variables
        integer(i32) :: m, n, nrowa
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
            ! ERROR: A must be square
            call errmgr%report_error("solve_tri_mtx", &
                "The input matrix must be square.", LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Call DTRSM
        call DTRSM(side, uplo, transa, diag, m, n, alpha, a, nrowa, b, m)
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Solves the system of equations: op(A) * X = B, where A is a
    !!  triangular matrix.
    !!
    !! @param[in] upper Set to true if A is an upper triangular matrix; else,
    !!  set to false if A is a lower triangular matrix.
    !! @param[in] trans Set to true if op(A) = A**T; else, set to false if
    !!  op(A) = A.
    !! @param[in] nounit Set to true if A is not a unit-diagonal matrix (ones on
    !!  every diagonal element); else, set to false if A is a unit-diagonal
    !!  matrix.
    !! @param[in] a The N-by-N triangular matrix.
    !! @param[in,out] x On input, the N-element right-hand-side array.  On
    !!  output, the N-element solution array.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - LA_ARRAY_SIZE_ERROR: Occurs if @p a is not square, or if the sizes of
    !!      @p a and @p b are not compatible.
    !!
    !!
    !! @par Usage
    !! To solve a triangular system of N equations of N unknowns A*X = B, where
    !! A is an N-by-N upper triangular matrix, and B and X are N-element
    !! arrays, the following code will suffice.
    !!
    !! @code{.f90}
    !! ! Solve the system: A*X = B, where A is an upper triangular N-by-N
    !! ! matrix, and B and X are N-elements in size.
    !!
    !! ! Variables
    !! integer(i32) :: info
    !! real(dp), dimension(n, n) :: a
    !! real(dp), dimension(n) :: b
    !!
    !! ! Initialize A and B...
    !!
    !! ! Solve A*X = B for X - Note: X overwrites B.
    !! call solve_triangular_system(.true., .false., a, b)
    !! @endcode
    !!
    !! @par Notes
    !! This routine is based upon the BLAS routine DTRSV.
    subroutine solve_tri_vec(upper, trans, nounit, a, x, err)
        ! Arguments
        logical, intent(in) :: upper, trans, nounit
        real(dp), intent(in), dimension(:,:) :: a
        real(dp), intent(inout), dimension(:) :: x
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        real(dp), parameter :: zero = 0.0d0

        ! Local Variables
        character :: uplo, t, diag
        integer(i32) :: n
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
            ! ERROR: A must be square
            call errmgr%report_error("solve_tri_vec", &
                "The input matrix must be square.", LA_ARRAY_SIZE_ERROR)
            return
        else if (size(x) /= n) then
            ! ERROR: Inner matrix dimensions must agree
            call errmgr%report_error("solve_tri_vec", &
                "The inner matrix dimensions must be equal.", &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Call DTRSV
        call DTRSV(uplo, t, diag, n, a, n, x, 1)
    end subroutine

! ******************************************************************************
! MATRIX MULTIPLICATION ROUTINES
! ------------------------------------------------------------------------------
    !> @brief Performs the matrix operation: C = alpha * op(A) * op(B) +
    !! beta * C.
    !!
    !! @param[in] transa Set to true if op(A) = A**T; else, set to false for
    !!  op(A) = A.
    !! @param[in] transb Set to true if op(B) = B**T; else, set to false for
    !!  op(B) = B.
    !! @param[in] alpha A scalar multiplier.
    !! @param[in] a If @p transa is set to true, an K-by-M matrix; else, if
    !!  @p transa is set to false, an M-by-K matrix.
    !! @param[in] b If @p transb is set to true, an N-by-K matrix; else, if
    !!  @p transb is set to false, a K-by-N matrix.
    !! @param[in] beta A scalar multiplier.
    !! @param[in,out] c On input, the M-by-N matrix C.  On output, the M-by-N
    !!  result.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - LA_ARRAY_SIZE_ERROR: Occurs if any of the input array sizes are 
    !!      incorrect.
    !!
    !! @par Notes
    !! This routine utilizes the BLAS routine DGEMM.
    subroutine mtx_mult_mtx(transa, transb, alpha, a, b, beta, c, err)
        ! Arguments
        logical, intent(in) :: transa, transb
        real(dp), intent(in) :: alpha, beta
        real(dp), intent(in), dimension(:,:) :: a, b
        real(dp), intent(inout), dimension(:,:) :: c
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        real(dp), parameter :: zero = 0.0d0
        real(dp), parameter :: one = 1.0d0

        ! Local Variables
        character :: ta, tb
        integer(i32) :: m, n, k, lda, ldb, flag
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg

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
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
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
            ! ERROR: Matrix dimensions mismatch
            write(errmsg, '(AI0A)') &
                "Matrix dimension mismatch.  Input number ", flag, &
                " was not sized correctly."
            call errmgr%report_error("mtx_mult_mtx", errmsg, &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Call DGEMM
        call DGEMM(ta, tb, m, n, k, alpha, a, lda, b, ldb, beta, c, m)
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Performs the matrix-vector operation: c = alpha * op(A) * b +
    !! beta * c.
    !!
    !! @param[in] trans Set to true if op(A) = A**T; else, set to false for
    !!  op(A) = A.
    !! @param[in] alpha A scalar multiplier.
    !! @param[in] a The M-by-N matrix A.
    !! @param[in] b If @p trans is set to true, an M-element array; else, if
    !!  @p trans is set to false, an N-element array.
    !! @param[in] beta A scalar multiplier.
    !! @param[in,out] c On input, if @p trans is set to true, an N-element
    !!  array; else, if @p trans is set to false, an M-element array.  On
    !!  output, the results of the operation.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - LA_ARRAY_SIZE_ERROR: Occurs if any of the input array sizes are 
    !!      incorrect.
    !!
    !! @par Notes
    !! This routine utilizes the BLAS routine DGEMV.
    subroutine mtx_mult_vec(trans, alpha, a, b, beta, c, err)
        ! Arguments
        logical, intent(in) :: trans
        real(dp), intent(in) :: alpha, beta
        real(dp), intent(in), dimension(:,:) :: a
        real(dp), intent(in), dimension(:) :: b
        real(dp), intent(inout), dimension(:) :: c
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        character :: t
        integer(i32) :: m, n, flag
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg

        ! Initialization
        m = size(a, 1)
        n = size(a, 2)
        t = 'N'
        if (trans) t = 'T'
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

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
            write(errmsg, '(AI0A)') &
                "Matrix dimension mismatch.  Input number ", flag, &
                " was not sized correctly."
            call errmgr%report_error("mtx_mult_vec", errmsg, &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Call DGEMV
        call DGEMV(t, m, n, alpha, a, m, b, 1, beta, c, 1)
    end subroutine

! ******************************************************************************
! RANK 1 UPDATE
! ------------------------------------------------------------------------------
    !> @brief Performs the rank-1 update to matrix A such that:
    !! A = alpha * X * Y**T + A, where A is an M-by-N matrix, alpha is a scalar,
    !! X is an M-element array, and N is an N-element array.
    !!
    !! @param[in] alpha The scalar multiplier.
    !! @param[in] x An M-element array.
    !! @param[in] y An N-element array.
    !! @param[in,out] a On input, the M-by-N matrix to update.  On output, the
    !!  updated M-by-N matrix.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - LA_ARRAY_SIZE_ERROR: Occurs if the size of @p a does not match with 
    !!      @p x and @p y.
    !!
    !! @par Notes
    !! This routine is based upon the BLAS routine DGER.
    subroutine rank1_update(alpha, x, y, a, err)
        ! Arguments
        real(dp), intent(in) :: alpha
        real(dp), intent(in), dimension(:) :: x, y
        real(dp), intent(inout), dimension(:,:) :: a
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        real(dp), parameter :: zero = 0.0d0

        ! Local Variables
        integer(i32) :: j, m, n
        real(dp) :: temp
        class(errors), pointer :: errmgr
        type(errors), target :: deferr

        ! Initialization
        m = size(x)
        n = size(y)
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        if (size(a, 1) /= m .or. size(a, 2) /= n) then
            ! ERROR: Matrix dimension array
            call errmgr%report_error("rank1_update", &
                "Matrix dimension mismatch.", LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Process
        do j = 1, n
            if (y(j) /= zero) then
                temp = alpha * y(j)
                a(:,j) = a(:,j) + temp * x
            end if
        end do
    end subroutine

! ******************************************************************************
! DIAGONAL MATRIX MULTIPLICATION ROUTINES
! ------------------------------------------------------------------------------
    !> @brief Computes the matrix operation: C = alpha * A * op(B) + beta * C,
    !! or C = alpha * op(B) * A + beta * C.
    !!
    !! @param[in] lside Set to true to apply matrix A from the left; else, set
    !!  to false to apply matrix A from the left.
    !! @param[in] trans Set to true if op(B) == B**T; else, set to false if
    !!  op(B) == B.
    !! @param[in] alpha A scalar multiplier.
    !! @param[in] a A K-element array containing the diagonal elements of A
    !!  where MIN(M,P) >= K >= 0 if @p lside is true; else, if @p lside is
    !!  false, MIN(N,P) >= K >= 0.
    !! @param[in] b The LDB-by-TDB matrix B where:
    !!  - @p lside == true & @p trans == true: LDA = N, TDB = P
    !!  - @p lside == true & @p trans == false: LDA = P, TDB = N
    !!  - @p lside == false & @p trans == true: LDA = P, TDB = M
    !!  - @p lside == false & @p trans == false: LDA = M, TDB = P
    !! @param[in] beta A scalar multiplier.
    !! @param[in,out] c On input, the M-by-N matrix C.  On output, the resulting
    !!  M-by-N matrix.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - LA_ARRAY_SIZE_ERROR: Occurs if any of the input array sizes are 
    !!      incorrect.
    subroutine diag_mtx_mult_mtx(lside, trans, alpha, a, b, beta, c, err)
        ! Arguments
        logical, intent(in) :: lside, trans
        real(dp) :: alpha, beta
        real(dp), intent(in), dimension(:) :: a
        real(dp), intent(in), dimension(:,:) :: b
        real(dp), intent(inout), dimension(:,:) :: c
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        real(dp), parameter :: zero = 0.0d0
        real(dp), parameter :: one = 1.0d0

        ! Local Variables
        integer(i32) :: i, m, n, k, nrowb, ncolb, flag
        real(dp) :: temp
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg

        ! Initialization
        m = size(c, 1)
        n = size(c, 2)
        k = size(a)
        nrowb = size(b, 1)
        ncolb = size(b, 2)
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

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
            write(errmsg, '(AI0A)') "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("diag_mtx_mult_mtx", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
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
                    if (temp /= one) c(i,:) = c(i,:) + temp * b(:,i)
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
                    if (temp /= one) c(i,:) = c(i,:) + temp * b(i,:)
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
                    if (temp /= one) c(:,i) = c(:,i) + temp * b(i,:)
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
                    if (temp /= one) c(:,i) = c(:,i) + temp * b(:,i)
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
    !> @brief Computes the matrix operation: B = alpha * A * op(B), or
    !! B = alpha * op(B) * A.
    !!
    !! @param[in] lside Set to true to apply matrix A from the left; else, set
    !!  to false to apply matrix A from the left.
    !! @param[in] alpha A scalar multiplier.
    !! @param[in] a A K-element array containing the diagonal elements of A
    !!  where MIN(M,P) >= K >= 0 if @p lside is true; else, if @p lside is
    !!  false, MIN(N,P) >= K >= 0.
    !! @param[in] b On input, the M-by-N matrix B.  On output, the resulting
    !!  M-by-N matrix.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - LA_ARRAY_SIZE_ERROR: Occurs if any of the input array sizes are 
    !!      incorrect.
    subroutine diag_mtx_mult_mtx2(lside, alpha, a, b, err)
        ! Arguments
        logical, intent(in) :: lside
        real(dp), intent(in) :: alpha
        real(dp), intent(in), dimension(:) :: a
        real(dp), intent(inout), dimension(:,:) :: b
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        real(dp), parameter :: zero = 0.0d0
        real(dp), parameter :: one = 1.0d0

        ! Local Variables
        integer(i32) :: i, m, n, k
        real(dp) :: temp
        class(errors), pointer :: errmgr
        type(errors), target :: deferr

        ! Initialization
        m = size(b, 1)
        n = size(b, 2)
        k = size(a)
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        if ((lside .and. k > m) .or. (.not.lside .and. k > n)) then
            ! ERROR: One of the input arrays is not sized correctly
            call errmgr%report_error("diag_mtx_mult_mtx", &
                "Input number 3 is not sized correctly.", &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Process
        if (lside) then
            ! Compute B = alpha * A * B
            do i = 1, k
                temp = alpha * a(i)
                if (temp /= one) b(i,:) = temp * b(i,:)
            end do
            if (m > k) b(k+1:m,:) = zero
        else
            ! Compute B = alpha * B * A
            do i = 1, k
                temp = alpha * a(i)
                if (temp /= one) b(:,i) = temp * b(:,i)
            end do
            if (n > k) b(:,k+1:n) = zero
        end if
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Computes the matrix operation: C = alpha * A * op(B) + beta * C,
    !! or C = alpha * op(B) * A + beta * C, where A and C are complex-valued.
    !!
    !! @param[in] lside Set to true to apply matrix A from the left; else, set
    !!  to false to apply matrix A from the left.
    !! @param[in] trans Set to true if op(B) == B**T; else, set to false if
    !!  op(B) == B.
    !! @param[in] alpha A scalar multiplier.
    !! @param[in] a A K-element array containing the diagonal elements of A
    !!  where MIN(M,P) >= K >= 0 if @p lside is true; else, if @p lside is
    !!  false, MIN(N,P) >= K >= 0.
    !! @param[in] b The LDB-by-TDB matrix B where:
    !!  - @p lside == true & @p trans == true: LDA = N, TDB = P
    !!  - @p lside == true & @p trans == false: LDA = P, TDB = N
    !!  - @p lside == false & @p trans == true: LDA = P, TDB = M
    !!  - @p lside == false & @p trans == false: LDA = M, TDB = P
    !! @param[in] beta A scalar multiplier.
    !! @param[in,out] c On input, the M-by-N matrix C.  On output, the resulting
    !!  M-by-N matrix.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - LA_ARRAY_SIZE_ERROR: Occurs if any of the input array sizes are 
    !!      incorrect.
    subroutine diag_mtx_mult_mtx3(lside, trans, alpha, a, b, beta, c, err)
        ! Arguments
        logical, intent(in) :: lside, trans
        real(dp) :: alpha, beta
        complex(dp), intent(in), dimension(:) :: a
        real(dp), intent(in), dimension(:,:) :: b
        complex(dp), intent(inout), dimension(:,:) :: c
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        complex(dp), parameter :: zero = (0.0d0, 0.0d0)
        complex(dp), parameter :: one = (1.0d0, 0.0d0)

        ! Local Variables
        integer(i32) :: i, m, n, k, nrowb, ncolb, flag
        complex(dp) :: temp
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg

        ! Initialization
        m = size(c, 1)
        n = size(c, 2)
        k = size(a)
        nrowb = size(b, 1)
        ncolb = size(b, 2)
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

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
            write(errmsg, '(AI0A)') "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("diag_mtx_mult_mtx3", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
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
                    if (temp /= one) c(i,:) = c(i,:) + temp * b(:,i)
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
                    if (temp /= one) c(i,:) = c(i,:) + temp * b(i,:)
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
                    if (temp /= one) c(:,i) = c(:,i) + temp * b(i,:)
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
                    if (temp /= one) c(:,i) = c(:,i) + temp * b(:,i)
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
    !> @brief Computes the matrix operation: C = alpha * A * op(B) + beta * C,
    !! or C = alpha * op(B) * A + beta * C, where A, B,  and C are
    !! complex-valued.
    !!
    !! @param[in] lside Set to true to apply matrix A from the left; else, set
    !!  to false to apply matrix A from the left.
    !! @param[in] trans Set to true if op(B) == B**T; else, set to false if
    !!  op(B) == B.
    !! @param[in] alpha A scalar multiplier.
    !! @param[in] a A K-element array containing the diagonal elements of A
    !!  where MIN(M,P) >= K >= 0 if @p lside is true; else, if @p lside is
    !!  false, MIN(N,P) >= K >= 0.
    !! @param[in] b The LDB-by-TDB matrix B where:
    !!  - @p lside == true & @p trans == true: LDA = N, TDB = P
    !!  - @p lside == true & @p trans == false: LDA = P, TDB = N
    !!  - @p lside == false & @p trans == true: LDA = P, TDB = M
    !!  - @p lside == false & @p trans == false: LDA = M, TDB = P
    !! @param[in] beta A scalar multiplier.
    !! @param[in,out] c On input, the M-by-N matrix C.  On output, the resulting
    !!  M-by-N matrix.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - LA_ARRAY_SIZE_ERROR: Occurs if any of the input array sizes are 
    !!      incorrect.
    subroutine diag_mtx_mult_mtx4(lside, trans, alpha, a, b, beta, c, err)
        ! Arguments
        logical, intent(in) :: lside, trans
        real(dp) :: alpha, beta
        complex(dp), intent(in), dimension(:) :: a
        complex(dp), intent(in), dimension(:,:) :: b
        complex(dp), intent(inout), dimension(:,:) :: c
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        complex(dp), parameter :: zero = (0.0d0, 0.0d0)
        complex(dp), parameter :: one = (1.0d0, 0.0d0)

        ! Local Variables
        integer(i32) :: i, m, n, k, nrowb, ncolb, flag
        complex(dp) :: temp
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg

        ! Initialization
        m = size(c, 1)
        n = size(c, 2)
        k = size(a)
        nrowb = size(b, 1)
        ncolb = size(b, 2)
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

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
            write(errmsg, '(AI0A)') "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("diag_mtx_mult_mtx4", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
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
                    if (temp /= one) c(i,:) = c(i,:) + temp * b(:,i)
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
                    if (temp /= one) c(i,:) = c(i,:) + temp * b(i,:)
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
                    if (temp /= one) c(:,i) = c(:,i) + temp * b(i,:)
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
                    if (temp /= one) c(:,i) = c(:,i) + temp * b(:,i)
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

! ******************************************************************************
! BASIC OPERATION ROUTINES
! ------------------------------------------------------------------------------
    !> @brief Computes the trace of a matrix (the sum of the main diagonal
    !! elements).
    !!
    !! @param[in] x The matrix on which to operate.
    !!
    !! @return The trace of @p x.
    pure function trace(x) result(y)
        ! Arguments
        real(dp), intent(in), dimension(:,:) :: x
        real(dp) :: y

        ! Parameters
        real(dp), parameter :: zero = 0.0d0

        ! Local Variables
        integer(i32) :: i, m, n, mn

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
    !> @brief Computes the rank of a matrix.
    !!
    !! @param[in,out] a On input, the M-by-N matrix of interest.  On output, the
    !!  contents of the matrix are overwritten.
    !! @param[in] tol An optional input, that if supplied, overrides the default
    !!  tolerance on singular values such that singular values less than this
    !!  tolerance are treated as zero.  The default tolerance is:
    !!  MAX(M, N) * EPS * MAX(S).  If the supplied value is less than the
    !!  smallest value that causes an overflow if inverted, the tolerance
    !!  reverts back to its default value, and the operation continues; however,
    !!  a warning message is issued.
    !! @param[out] work An optional input, that if provided, prevents any local
    !!  memory allocation.  If not provided, the memory required is allocated
    !!  within.  If provided, the length of the array must be at least
    !!  @p olwork.
    !! @param[out] olwork An optional output used to determine workspace size.
    !!  If supplied, the routine determines the optimal size for @p work, and
    !!  returns without performing any actual calculations.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - LA_ARRAY_SIZE_ERROR: Occurs if any of the input arrays are not sized
    !!      appropriately.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
    !!      there is insufficient memory available.
    !!  - LA_CONVERGENCE_ERROR: Occurs as a warning if the QR iteration process
    !!      could not converge to a zero value.
    !!
    !! @par See Also
    !! - [Wolfram MathWorld](http://mathworld.wolfram.com/MatrixRank.html)
    function mtx_rank(a, tol, work, olwork, err) result(rnk)
        ! Arguments
        real(dp), intent(inout), dimension(:,:) :: a
        real(dp), intent(in), optional :: tol
        real(dp), intent(out), pointer, optional, dimension(:) :: work
        integer(i32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err
        integer(i32) :: rnk

        ! Local Variables
        integer(i32) :: i, m, n, mn, istat, lwork, flag
        real(dp), pointer, dimension(:) :: wptr, s, w
        real(dp), allocatable, target, dimension(:) :: wrk
        real(dp) :: t, tref, smlnum
        real(dp), dimension(1) :: dummy, temp
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg

        ! Initialization
        m = size(a, 1)
        n = size(a, 2)
        mn = min(m, n)
        smlnum = DLAMCH('s')
        rnk = 0
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Workspace Query
        !call svd(a, a(1:mn,1), olwork = lwork)
        call DGESVD('N', 'N', m, n, a, m, dummy, dummy, m, dummy, n, temp, &
            -1, flag)
        lwork = int(temp(1), i32) + mn
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (size(work) < lwork) then
                ! ERROR: WORK not sized correctly
                call errmgr%report_error("mtx_rank", &
                    "Incorrectly sized input array WORK, argument 5.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lwork)
        else
            allocate(wrk(lwork), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("mtx_rank", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if
        s => wptr(1:mn)
        w => wptr(mn+1:lwork)

        ! Compute the singular values of A
        call DGESVD('N', 'N', m, n, a, m, s, dummy, m, dummy, n, w, &
            lwork - mn, flag)
        if (flag > 0) then
            write(errmsg, '(I0A)') flag, " superdiagonals could not " // &
                "converge to zero as part of the QR iteration process."
            call errmgr%report_warning("mtx_rank", errmsg, LA_CONVERGENCE_ERROR)
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
    !> @brief Computes the determinant of a square matrix.
    !!
    !! @param[in,out] a On input, the N-by-N matrix on which to operate.  On
    !! output the contents are overwritten by the LU factorization of the
    !! original matrix.
    !! @param[out] iwork An optional input, that if provided, prevents any local
    !!  memory allocation.  If not provided, the memory required is allocated
    !!  within.  If provided, the length of the array must be at least
    !!  N-elements.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - LA_ARRAY_SIZE_ERROR: Occurs if any of the input arrays are not sized
    !!      appropriately.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
    !!      there is insufficient memory available.
    !!
    !! @return The determinant of @p a.
    function det(a, iwork, err) result(x)
        ! Arguments
        real(dp), intent(inout), dimension(:,:) :: a
        integer(i32), intent(out), pointer, optional, dimension(:) :: iwork
        class(errors), intent(inout), optional, target :: err
        real(dp) :: x

        ! Parameters
        real(dp), parameter :: zero = 0.0d0
        real(dp), parameter :: one = 1.0d0
        real(dp), parameter :: ten = 1.0d1
        real(dp), parameter :: p1 = 1.0d-1

        ! Local Variables
        integer(i32) :: i, ep, n, istat, flag
        integer(i32), pointer, dimension(:) :: ipvt
        integer(i32), allocatable, target, dimension(:) :: iwrk
        real(dp) :: temp
        class(errors), pointer :: errmgr
        type(errors), target :: deferr

        ! Initialization
        n = size(a, 1)
        x = zero
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        if (size(a, 2) /= n) then
            call errmgr%report_error("det", &
                "The supplied matrix must be square.", LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Local Memory Allocation
        if (present(iwork)) then
            if (size(iwork) < n) then
                ! ERROR: WORK not sized correctly
                call errmgr%report_error("det", &
                    "Incorrectly sized input array IWORK, argument 2.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            ipvt => iwork(1:n)
        else
            allocate(iwrk(n), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("det", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            ipvt => iwrk
        end if

        ! Compute the LU factorization of A
        call DGETRF(n, n, a, n, ipvt, flag)
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

            temp = a(i,i) * temp
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

! ******************************************************************************
! ARRAY SWAPPING ROUTINE
! ------------------------------------------------------------------------------
    !> @brief Swaps the contents of two arrays.
    !!
    !! @param[in,out] x One of the N-element arrays.
    !! @param[in,out] y The other N-element array.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - LA_ARRAY_SIZE_ERROR: Occurs if @p x and @p y are not the same size.
    subroutine swap(x, y, err)
        ! Arguments
        real(dp), intent(inout), dimension(:) :: x, y
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        integer(i32) :: i, n
        real(dp) :: temp
        class(errors), pointer :: errmgr
        type(errors), target :: deferr

        ! Initialization
        n = size(x)
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        if (size(y) /= n) then
            call errmgr%report_error("swap", &
                "The input arrays are not the same size.", &
                LA_ARRAY_SIZE_ERROR)
            return
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
    !> @brief Multiplies a vector by the reciprocal of a real scalar.
    !!
    !! @param[in] a The scalar which is used to divide each component of @p X.
    !!  The value must be >= 0, or the subroutine will divide by zero.
    !! @param[in,out] x The vector.
    !!
    !! @par Notes
    !! This routine is based upon the LAPACK routine DRSCL.
    subroutine recip_mult_array(a, x)
        ! Arguments
        real(dp), intent(in) :: a
        real(dp), intent(inout), dimension(:) :: x

        ! Parameters
        real(dp), parameter :: zero = 0.0d0
        real(dp), parameter :: one = 1.0d0
        real(dp), parameter :: twotho = 2.0d3

        ! Local Variables
        logical :: done
        real(dp) :: bignum, cden, cden1, cnum, cnum1, mul, smlnum

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


end module
