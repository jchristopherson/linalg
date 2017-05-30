! linalg_core.f90

!> @brief \b linalg_core
!!
!! @par Purpose
!! Provides common "core" linear algebra routines.
module linalg_core
    use ferror, only : errors
    use linalg_constants
    implicit none
    private
    public :: solve_triangular_system
    public :: mtx_mult
    public :: rank1_update

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

! ------------------------------------------------------------------------------
end module
