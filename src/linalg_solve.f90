! linalg_solve.f90

!> @brief \b linalg_solve
!!
!! @par Purpose
!! Provides a set of routines for solving systems of linear equations.
module linalg_solve
    use, intrinsic :: iso_fortran_env, only : int32, real64
    use ferror, only : errors
    use linalg_constants
    use linalg_factor, only : rz_factor, mult_rz, mult_qr
    use linalg_core, only : mtx_mult, recip_mult_array
    implicit none
    private
    public :: solve_triangular_system
    public :: solve_lu
    public :: solve_qr
    public :: solve_cholesky
    public :: mtx_inverse
    public :: mtx_pinverse
    public :: solve_least_squares
    public :: solve_least_squares_full
    public :: solve_least_squares_svd

! ******************************************************************************
! INTERFACES
! ------------------------------------------------------------------------------
    !> @brief Solves a triangular system of equations.
    !!
    !! @par Usage
    !! The following example illustrates the solution of two triangular systems
    !! to solve a system of LU factored equations.
    !! @code{.f90}
    !! program example
    !!     use iso_fortran_env, only : real64, int32
    !!     use linalg_factor, only : lu_factor, form_lu
    !!     use linalg_solve, only : solve_triangular_system
    !!     implicit none
    !!
    !!     ! Variables
    !!     real(real64) :: a(3,3), b(3), u(3,3), p(3,3)
    !!     integer(int32) :: i, pvt(3)
    !!
    !!     ! Build the 3-by-3 matrix A.
    !!     !     | 1   2   3 |
    !!     ! A = | 4   5   6 |
    !!     !     | 7   8   0 |
    !!     a = reshape( &
    !!         [1.0d0, 4.0d0, 7.0d0, 2.0d0, 5.0d0, 8.0d0, 3.0d0, 6.0d0, 0.0d0], &
    !!         [3, 3])
    !!
    !!     ! Build the right-hand-side vector B.
    !!     !     | -1 |
    !!     ! b = | -2 |
    !!     !     | -3 |
    !!     b = [-1.0d0, -2.0d0, -3.0d0]
    !!
    !!     ! The solution is:
    !!     !     |  1/3 |
    !!     ! x = | -2/3 |
    !!     !     |   0  |
    !!
    !!     ! Compute the LU factorization
    !!     call lu_factor(a, pvt)
    !!
    !!     ! Extract the L and U matrices. A is overwritten with L.
    !!     call form_lu(a, pvt, u, p)
    !!
    !!     ! Solve the lower triangular system L * Y = P * B for Y, but first compute
    !!     ! P * B, and store the results in B
    !!     b = matmul(p, b)
    !!
    !!     ! Now, compute the solution to the lower triangular system.  Store the
    !!     ! result in B.  Remember, L is unit diagonal (ones on its diagonal)
    !!     call solve_triangular_system(.false., .false., .false., a, b)
    !! 
    !!     ! Solve the upper triangular system U * X = Y for X.
    !!     call solve_triangular_system(.true., .false., .true., u, b)
    !!
    !!     ! Display the results.
    !!     print '(A)', "LU Solution: X = "
    !!     print '(F8.4)', (b(i), i = 1, size(b))
    !! end program
    !! @endcode
    !! The above program produces the following output.
    !! @code{.txt}
    !! LU Solution: X =
    !! 0.3333
    !! -0.6667
    !! 0.0000
    !! @endcode
    interface solve_triangular_system
        module procedure :: solve_tri_mtx
        module procedure :: solve_tri_vec
    end interface

! ------------------------------------------------------------------------------
    !> @brief Solves a system of LU-factored equations.
    !!
    !! @par Usage
    !! To solve a system of 3 equations of 3 unknowns using LU factorization,
    !! the following code will suffice.
    !! @code{.f90}
    !! program example
    !!     use iso_fortran_env
    !!     use linalg_factor, only : lu_factor
    !!     use linalg_solve, only : solve_lu
    !!     implicit none
    !!
    !!     ! Local Variables
    !!     real(real64) :: a(3,3), b(3)
    !!     integer(int32) :: i, pvt(3)
    !!
    !!     ! Build the 3-by-3 matrix A.
    !!     !     | 1   2   3 |
    !!     ! A = | 4   5   6 |
    !!     !     | 7   8   0 |
    !!     a = reshape( &
    !!         [1.0d0, 4.0d0, 7.0d0, 2.0d0, 5.0d0, 8.0d0, 3.0d0, 6.0d0, 0.0d0], &
    !!         [3, 3])
    !!
    !!     ! Build the right-hand-side vector B.
    !!     !     | -1 |
    !!     ! b = | -2 |
    !!     !     | -3 |
    !!     b = [-1.0d0, -2.0d0, -3.0d0]
    !!
    !!     ! The solution is:
    !!     !     |  1/3 |
    !!     ! x = | -2/3 |
    !!     !     |   0  |
    !!
    !!     ! Compute the LU factorization
    !!     call lu_factor(a, pvt)
    !!
    !!     ! Compute the solution.  The results overwrite b.
    !!     call solve_lu(a, pvt, b)
    !!
    !!     ! Display the results.
    !!     print '(A)', "LU Solution: X = "
    !!     print '(F8.4)', (b(i), i = 1, size(b))
    !! end program
    !! @endcode
    !! The program generates the following output.
    !! @code{.txt}
    !!  LU Solution: X =
    !!   0.3333
    !!  -0.6667
    !!   0.0000
    !! @endcode
    !!
    !! @par See Also
    !! - [Wikipedia](https://en.wikipedia.org/wiki/LU_decomposition)
    !! - [Wolfram MathWorld](http://mathworld.wolfram.com/LUDecomposition.html)
    interface solve_lu
        module procedure :: solve_lu_mtx
        module procedure :: solve_lu_vec
    end interface

! ------------------------------------------------------------------------------
    !> @brief Solves a system of M QR-factored equations of N unknowns.
    !!
    !! @par Usage
    !! The following example illustrates the solution of a system of equations
    !! using QR factorization.
    !! @code{.f90}
    !! program example
    !!     use iso_fortran_env, only : real64, int32
    !!     use linalg_factor, only : qr_factor
    !!     use linalg_solve, only : solve_qr
    !!
    !!     ! Local Variables
    !!     real(real64) :: a(3,3), tau(3), b(3)
    !!     integer(int32) :: i, pvt(3)
    !!
    !!     ! Build the 3-by-3 matrix A.
    !!     !     | 1   2   3 |
    !!     ! A = | 4   5   6 |
    !!     !     | 7   8   0 |
    !!     a = reshape( &
    !!         [1.0d0, 4.0d0, 7.0d0, 2.0d0, 5.0d0, 8.0d0, 3.0d0, 6.0d0, 0.0d0], &
    !!         [3, 3])
    !!
    !!     ! Build the right-hand-side vector B.
    !!     !     | -1 |
    !!     ! b = | -2 |
    !!     !     | -3 |
    !!     b = [-1.0d0, -2.0d0, -3.0d0]
    !!
    !!     ! The solution is:
    !!     !     |  1/3 |
    !!     ! x = | -2/3 |
    !!     !     |   0  |
    !!
    !!     ! Compute the QR factorization, using pivoting
    !!     pvt = 0     ! Zero every entry in order not to lock any column in place
    !!     call qr_factor(a, tau, pvt)
    !!
    !!     ! Compute the solution.  The results overwrite b.
    !!     call solve_qr(a, tau, pvt, b)
    !!
    !!     ! Display the results.
    !!     print '(A)', "QR Solution: X = "
    !!     print '(F8.4)', (b(i), i = 1, size(b))
    !!
    !!     ! Notice, QR factorization without pivoting could be accomplished in the
    !!     ! same manner.  The only difference is to omit the PVT array (column pivot
    !!     ! tracking array).
    !! end program
    !! @endcode
    !! The above program produces the following output.
    !! @code{.txt}
    !! QR Solution: X =
    !! 0.3333
    !! -0.6667
    !! 0.0000
    !! @endcode
    !!
    !! @par See Also
    !! - [Wikipedia](https://en.wikipedia.org/wiki/QR_decomposition)
    !! - [LAPACK Users Manual](http://netlib.org/lapack/lug/node39.html)
    interface solve_qr
        module procedure :: solve_qr_no_pivot_mtx
        module procedure :: solve_qr_no_pivot_vec
        module procedure :: solve_qr_pivot_mtx
        module procedure :: solve_qr_pivot_vec
    end interface

! ------------------------------------------------------------------------------
    !> @brief Solves a system of Cholesky factored equations.
    !!
    !! @par Usage
    !! The following example illustrates the solution of a positive-definite 
    !! system of equations via Cholesky factorization.
    !! @code{.f90}
    !! program example
    !!     use iso_fortran_env, only : real64, int32
    !!     use linalg_factor, only : cholesky_factor
    !!     use linalg_solve, only : solve_cholesky, solve_triangular_system
    !!     implicit none
    !!
    !!     ! Variables
    !!     real(real64) :: a(3, 3), b(3), bu(3)
    !!     integer(int32) :: i
    !!
    !!     ! Build the 3-by-3 positive-definite matrix A.
    !!     !     | 4   12   -16 |
    !!     ! A = | 12  37   -43 |
    !!     !     |-16 -43    98 |
    !!     a = reshape([4.0d0, 12.0d0, -16.0d0, 12.0d0, 37.0d0, -43.0d0, -16.0d0, &
    !!         -43.0d0, 98.0d0], [3, 3])
    !!
    !!     ! Build the 3-element array B
    !!     !     | 5 |
    !!     ! b = | 1 |
    !!     !     | 3 |
    !!     b = [5.0d0, 1.0d0, 3.0d0]
    !!
    !!     ! Make a copy of B for later use - not necessary, but just for example to
    !!     ! illustrate the long or manual method of solving a Cholesky factored system
    !!     bu = b
    !!
    !!     ! Compute the Cholesky factorization of A considering only the upper 
    !!     ! triangular portion of A (the default configuration).
    !!     call cholesky_factor(a)
    !!
    !!     ! Compute the solution
    !!     call solve_cholesky(.true., a, b)
    !!
    !!     ! Display the results
    !!     print '(A)', "Cholesky Solution: X = "
    !!     print '(F8.4)', (b(i), i = 1, size(b))
    !!
    !!     ! The solution could also be computed manually noting the Cholesky 
    !!     ! factorization causes A = U**T * U.  Then U**T * U * X = B.  
    !!
    !!     ! Step 1 would then be to solve the problem U**T * Y = B, for Y.
    !!     call solve_triangular_system(.true., .true., .true., a, bu)
    !!
    !!     ! Now, solve the problem U * X = Y, for X
    !!     call solve_triangular_system(.true., .false., .true., a, bu)
    !!
    !!     ! Display the results
    !!     print '(A)', "Cholesky Solution (Manual Approach): X = "
    !!     print '(F8.4)', (bu(i), i = 1, size(bu))
    !! end program
    !! @endcode
    !! The above program produces the following output.
    !! @code{.txt}
    !! Cholesky Solution: X =
    !!  239.5833
    !!  -65.6667
    !!  10.3333
    !! Cholesky Solution (Manual Approach): X =
    !!  239.5833
    !!  -65.6667
    !!  10.3333
    !! @endcode
    interface solve_cholesky
        module procedure :: solve_cholesky_mtx
        module procedure :: solve_cholesky_vec
    end interface

! ------------------------------------------------------------------------------
    !> @brief Solves the overdetermined or underdetermined system (A*X = B) of
    !! M equations of N unknowns.
    interface solve_least_squares
        module procedure :: solve_least_squares_mtx
        module procedure :: solve_least_squares_vec
    end interface

! ------------------------------------------------------------------------------
    !> @brief Solves the overdetermined or underdetermined system (A*X = B) of
    !! M equations of N unknowns, but uses a full orthogonal factorization of
    !! the system.
    interface solve_least_squares_full
        module procedure :: solve_least_squares_mtx_pvt
        module procedure :: solve_least_squares_vec_pvt
    end interface

! ------------------------------------------------------------------------------
    !> @brief Solves the overdetermined or underdetermined system (A*X = B) of
    !! M equations of N unknowns using a singular value decomposition of
    !! matrix A.
    interface solve_least_squares_svd
        module procedure :: solve_least_squares_mtx_svd
        module procedure :: solve_least_squares_vec_svd
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
    !! @par Notes
    !! This routine is based upon the BLAS routine DTRSM.
    subroutine solve_tri_mtx(lside, upper, trans, nounit, alpha, a, b, err)
        ! Arguments
        logical, intent(in) :: lside, upper, trans, nounit
        real(real64), intent(in) :: alpha
        real(real64), intent(in), dimension(:,:) :: a
        real(real64), intent(inout), dimension(:,:) :: b
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        character :: side, uplo, transa, diag
        real(real64), parameter :: zero = 0.0d0
        real(real64), parameter :: one = 1.0d0

        ! Local Variables
        integer(int32) :: m, n, nrowa
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
    !! integer(int32) :: info
    !! real(real64), dimension(n, n) :: a
    !! real(real64), dimension(n) :: b
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
        real(real64), intent(in), dimension(:,:) :: a
        real(real64), intent(inout), dimension(:) :: x
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        real(real64), parameter :: zero = 0.0d0

        ! Local Variables
        character :: uplo, t, diag
        integer(int32) :: n
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
! LU SOLUTION
! ------------------------------------------------------------------------------
    !> @brief Solves a system of LU-factored equations.
    !!
    !! @param[in] a The N-by-N LU factored matrix as output by lu_factor.
    !! @param[in] ipvt The N-element pivot array as output by lu_factor.
    !! @param[in,out] b On input, the N-by-NRHS right-hand-side matrix.  On
    !!  output, the N-by-NRHS solution matrix.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - LA_ARRAY_SIZE_ERROR: Occurs if any of the input array sizes are 
    !!      incorrect.
    !!
    !! @par Notes
    !! The routine is based upon the LAPACK routine DGETRS.
    subroutine solve_lu_mtx(a, ipvt, b, err)
        ! Arguments
        real(real64), intent(in), dimension(:,:) :: a
        integer(int32), intent(in), dimension(:) :: ipvt
        real(real64), intent(inout), dimension(:,:) :: b
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        integer(int32) :: n, nrhs, flag
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg

        ! Initialization
        n = size(a, 1)
        nrhs = size(b, 2)
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        flag = 0
        if (size(a, 2) /= n) then
            flag = 1
        else if (size(ipvt) /= n) then
            flag = 2
        else if (size(b, 1) /= n) then
            flag = 3
        end if
        if (flag /= 0) then
            ! One of the input arrays is not sized correctly
            write(errmsg, '(AI0A)') "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("solve_lu_mtx", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Call DGETRS
        call DGETRS("N", n, nrhs, a, n, ipvt, b, n, flag)
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Solves a system of LU-factored equations.
    !!
    !! @param[in] a The N-by-N LU factored matrix as output by lu_factor.
    !! @param[in] ipvt The N-element pivot array as output by lu_factor.
    !! @param[in,out] b On input, the N-element right-hand-side array.  On
    !!  output, the N-element solution array.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - LA_ARRAY_SIZE_ERROR: Occurs if any of the input array sizes are 
    !!      incorrect.
    !!
    !! @par Notes
    !! The routine is based upon the LAPACK routine DGETRS.
    subroutine solve_lu_vec(a, ipvt, b, err)
        ! Arguments
        real(real64), intent(in), dimension(:,:) :: a
        integer(int32), intent(in), dimension(:) :: ipvt
        real(real64), intent(inout), dimension(:) :: b
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        integer(int32) :: n, flag
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg

        ! Initialization
        n = size(a, 1)
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        flag = 0
        if (size(a, 2) /= n) then
            flag = 1
        else if (size(ipvt) /= n) then
            flag = 2
        else if (size(b) /= n) then
            flag = 3
        end if
        if (flag /= 0) then
            ! One of the input arrays is not sized correctly
            write(errmsg, '(AI0A)') "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("solve_lu_vec", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Call DGETRS
        call DGETRS("N", n, 1, a, n, ipvt, b, n, flag)
    end subroutine

! ******************************************************************************
! QR SOLUTION
! ------------------------------------------------------------------------------
    !> @brief Solves a system of M QR-factored equations of N unknowns where
    !! M >= N.
    !!
    !! @param[in] a On input, the M-by-N QR factored matrix as returned by
    !!  qr_factor.  On output, the contents of this matrix are restored.
    !!  Notice, M must be greater than or equal to N.
    !! @param[in] tau A MIN(M, N)-element array containing the scalar factors of
    !!  the elementary reflectors as returned by qr_factor.
    !! @param[in] b On input, the M-by-NRHS right-hand-side matrix.  On output,
    !!  the first N columns are overwritten by the solution matrix X.
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
    !!
    !! @par Notes
    !! This routine is based upon a subset of the LAPACK routine DGELS.
    subroutine solve_qr_no_pivot_mtx(a, tau, b, work, olwork, err)
        ! Arguments
        real(real64), intent(inout), dimension(:,:) :: a, b
        real(real64), intent(in), dimension(:) :: tau
        real(real64), intent(out), target, optional, dimension(:) :: work
        integer(int32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        real(real64), parameter :: one = 1.0d0

        ! Local Variables
        integer(int32) :: m, n, nrhs, k, lwork, flag, istat
        real(real64), pointer, dimension(:) :: wptr
        real(real64), allocatable, target, dimension(:) :: wrk
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg

        ! Initialization
        m = size(a, 1)
        n = size(a, 2)
        nrhs = size(b, 2)
        k = min(m, n)
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        flag = 0
        if (m < n) then
            flag = 1
        else if (size(tau) /= k) then
            flag = 2
        else if (size(b, 1) /= m) then
            flag = 3
        end if
        if (flag /= 0) then
            ! ERROR: One of the input arrays is not sized correctly
            write(errmsg, '(AI0A)') "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("solve_qr_no_pivot_mtx", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Workspace Query
        call mult_qr(.true., .true., a, tau, b, olwork = lwork)
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (size(work) < lwork) then
                ! ERROR: WORK not sized correctly
                call errmgr%report_error("solve_qr_no_pivot_mtx", &
                    "Incorrectly sized input array WORK, argument 4.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lwork)
        else
            allocate(wrk(lwork), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("solve_qr_no_pivot_mtx", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if

        ! Compute Q**T * B, and store in B
        call mult_qr(.true., .true., a, tau, b, wptr)

        ! Solve the triangular system: A(1:N,1:N)*X = B(1:N,:)
        call solve_triangular_system(.true., .true., .false., .true., one, &
            a(1:n,1:n), b(1:n,:))
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Solves a system of M QR-factored equations of N unknowns where
    !! M >= N.
    !!
    !! @param[in] a On input, the M-by-N QR factored matrix as returned by
    !!  qr_factor.  On output, the contents of this matrix are restored.
    !!  Notice, M must be greater than or equal to N.
    !! @param[in] tau A MIN(M, N)-element array containing the scalar factors of
    !!  the elementary reflectors as returned by qr_factor.
    !! @param[in] b On input, the M-element right-hand-side vector.  On output,
    !!  the first N elements are overwritten by the solution vector X.
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
    !!
    !! @par Notes
    !! This routine is based upon a subset of the LAPACK routine DGELS.
    subroutine solve_qr_no_pivot_vec(a, tau, b, work, olwork, err)
        ! Arguments
        real(real64), intent(inout), dimension(:,:) :: a
        real(real64), intent(in), dimension(:) :: tau
        real(real64), intent(inout), dimension(:) :: b
        real(real64), intent(out), target, optional, dimension(:) :: work
        integer(int32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        integer(int32) :: m, n, k, flag, lwork, istat
        real(real64), pointer, dimension(:) :: wptr
        real(real64), allocatable, target, dimension(:) :: wrk
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg

        ! Initialization
        m = size(a, 1)
        n = size(a, 2)
        k = min(m, n)
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        flag = 0
        if (m < n) then
            flag = 1
        else if (size(tau) /= k) then
            flag = 2
        else if (size(b) /= m) then
            flag = 3
        end if
        if (flag /= 0) then
            ! ERROR: One of the input arrays is not sized correctly
            write(errmsg, '(AI0A)') "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("solve_qr_no_pivot_vec", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Workspace Query
        call mult_qr(.true., a, tau, b, olwork = lwork)
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (size(work) < lwork) then
                ! ERROR: WORK not sized correctly
                call errmgr%report_error("solve_qr_no_pivot_vec", &
                    "Incorrectly sized input array WORK, argument 4.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lwork)
        else
            allocate(wrk(lwork), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("solve_qr_no_pivot_vec", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if

        ! Compute Q**T * B, and store in B
        call mult_qr(.true., a, tau, b, work = wptr)

        ! Solve the triangular system: A(1:N,1:N)*X = B(1:N)
        call solve_triangular_system(.true., .false., .true., a(1:n,1:n), b)
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Solves a system of M QR-factored equations of N unknowns where the
    !! QR factorization made use of column pivoting.
    !!
    !! @param[in] a On input, the M-by-N QR factored matrix as returned by
    !!  qr_factor.  On output, the contents of this matrix are altered.
    !! @param[in] tau A MIN(M, N)-element array containing the scalar factors of
    !!  the elementary reflectors as returned by qr_factor.
    !! @param[in] jpvt An N-element array, as output by qr_factor, used to
    !!  track the column pivots.
    !! @param[in] b On input, the MAX(M, N)-by-NRHS matrix where the first M
    !!  rows contain the right-hand-side matrix B.  On output, the first N rows
    !!  are overwritten by the solution matrix X.
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
    !!
    !! @par Notes
    !! This routine is based upon a subset of the LAPACK routine DGELSY.
    subroutine solve_qr_pivot_mtx(a, tau, jpvt, b, work, olwork, err)
        ! Arguments
        real(real64), intent(inout), dimension(:,:) :: a
        real(real64), intent(in), dimension(:) :: tau
        integer(int32), intent(in), dimension(:) :: jpvt
        real(real64), intent(inout), dimension(:,:) :: b
        real(real64), intent(out), target, optional, dimension(:) :: work
        integer(int32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        integer(int32), parameter :: imin = 2
        integer(int32), parameter :: imax = 1
        real(real64), parameter :: zero = 0.0d0
        real(real64), parameter :: one = 1.0d0

        ! Local Variables
        integer(int32) :: i, j, m, n, mn, nrhs, lwork, ismin, ismax, &
            rnk, maxmn, flag, istat, lwork1, lwork2, lwork3
        real(real64) :: rcond, smax, smin, smaxpr, sminpr, s1, c1, s2, c2
        real(real64), pointer, dimension(:) :: wptr, w, tau2
        real(real64), allocatable, target, dimension(:) :: wrk
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg

        ! Initialization
        m = size(a, 1)
        n = size(a, 2)
        mn = min(m, n)
        maxmn = max(m, n)
        nrhs = size(b, 2)
        ismin = mn + 1
        ismax = 2 * mn + 1
        rcond = epsilon(rcond)
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        flag = 0
        if (size(tau) /= mn) then
            flag = 2
        else if (size(jpvt) /= n) then
            flag = 3
        else if (size(b, 1) /= maxmn) then
            flag = 4
        end if
        if (flag /= 0) then
            ! ERROR: One of the input arrays is not sized correctly
            write(errmsg, '(AI0A)') "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("solve_qr_pivot_mtx", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Workspace Query
        call rz_factor(a(1:mn,:), a(1:mn,1), olwork = lwork1)
        call mult_qr(.true., .true., a, tau, b(1:m,:), olwork = lwork2)
        call mult_rz(.true., .true., n, a(1:mn,:), a(1:mn,1), b(1:n,:), &
            olwork = lwork3)
        lwork = max(lwork1, lwork2, lwork3, 2 * mn + 1) + mn
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (size(work) < lwork) then
                ! ERROR: WORK not sized correctly
                call errmgr%report_error("solve_qr_no_pivot_mtx", &
                    "Incorrectly sized input array WORK, argument 5.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lwork)
        else
            allocate(wrk(lwork), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("solve_qr_pivot_mtx", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if

        ! Determine the rank of R11 using an incremental condition estimation
        wptr(ismin) = one
        wptr(ismax) = one
        smax = abs(a(1,1))
        smin = smax
        if (abs(a(1,1)) == zero) then
            rnk = 0
            b(1:maxmn,:) = zero
            return
        else
            rnk = 1
        end if
        do
            if (rnk < mn) then
                i = rnk + 1
                call DLAIC1(imin, rnk, wptr(ismin:ismin+rnk-1), smin, &
                    a(1:rnk-1,i), a(i,i), sminpr, s1, c1)
                call DLAIC1(imax, rnk, wptr(ismax:ismax+rnk-1), smax, &
                    a(1:rnk-1,i), a(i,i), smaxpr, s2, c2)
                if (smaxpr * rcond <= sminpr) then
                    do i = 1, rnk
                        wptr(ismin+i-1) = s1 * wptr(ismin+i-1)
                        wptr(ismax+i-1) = s2 * wptr(ismax+i-1)
                    end do
                    wptr(ismin+rnk) = c1
                    wptr(ismax+rnk) = c2
                    smin = sminpr
                    smax = smaxpr
                    rnk = rnk + 1
                    cycle
                end if
            end if
            exit
        end do

        ! Partition R = [R11 R12]
        !               [ 0  R22]
        tau2 => wptr(1:rnk)
        w => wptr(rnk+1:lwork)
        if (rnk < n) call rz_factor(a(1:rnk,:), tau2, w)

        ! Compute B(1:m,1:NRHS) = Q**T * B(1:M,1:NRHS)
        call mult_qr(.true., .true., a, tau, b(1:m,:), w)

        ! Solve the triangular system T11 * B(1:rnk,1:nrhs) = B(1:rnk,1:nrhs)
        call solve_triangular_system(.true., .true., .false., .true., one, &
            a(1:rnk,1:rnk), b(1:rnk,:))
        if (n > rnk) b(rnk+1:n,:) = zero

        ! Compute B(1:n,1:nrhs) = Y**T * B(1:n,1:nrhs)
        if (rnk < n) then
            call mult_rz(.true., .true., n - rnk, a(1:rnk,:), tau2, b(1:n,:), w)
        end if

        ! Apply the pivoting: B(1:N,1:NRHS) = P * B(1:N,1:NRHS)
        do j = 1, nrhs
            do i = 1, n
                wptr(jpvt(i)) = b(i,j)
            end do
            b(:,j) = wptr(1:n)
        end do
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Solves a system of M QR-factored equations of N unknowns where the
    !! QR factorization made use of column pivoting.
    !!
    !! @param[in] a On input, the M-by-N QR factored matrix as returned by
    !!  qr_factor.  On output, the contents of this matrix are altered.
    !! @param[in] tau A MIN(M, N)-element array containing the scalar factors of
    !!  the elementary reflectors as returned by qr_factor.
    !! @param[in] jpvt An N-element array, as output by qr_factor, used to
    !!  track the column pivots.
    !! @param[in] b On input, the MAX(M, N)-element array where the first M
    !!  elements contain the right-hand-side vector B.  On output, the first N
    !!  elements are overwritten by the solution vector X.
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
    !!
    !! @par Notes
    !! This routine is based upon a subset of the LAPACK routine DGELSY.
    subroutine solve_qr_pivot_vec(a, tau, jpvt, b, work, olwork, err)
        ! Arguments
        real(real64), intent(inout), dimension(:,:) :: a
        real(real64), intent(in), dimension(:) :: tau
        integer(int32), intent(in), dimension(:) :: jpvt
        real(real64), intent(inout), dimension(:) :: b
        real(real64), intent(out), target, optional, dimension(:) :: work
        integer(int32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        integer(int32), parameter :: imin = 2
        integer(int32), parameter :: imax = 1
        real(real64), parameter :: zero = 0.0d0
        real(real64), parameter :: one = 1.0d0

        ! Local Variables
        integer(int32) :: i, m, n, mn, lwork, ismin, ismax, rnk, maxmn, flag, &
            istat, lwork1, lwork2
        real(real64) :: rcond, smax, smin, smaxpr, sminpr, s1, c1, s2, c2
        real(real64), pointer, dimension(:) :: wptr, w, tau2
        real(real64), allocatable, target, dimension(:) :: wrk
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg

        ! Initialization
        m = size(a, 1)
        n = size(a, 2)
        mn = min(m, n)
        maxmn = max(m, n)
        ismin = mn + 1
        ismax = 2 * mn + 1
        rcond = epsilon(rcond)
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        flag = 0
        if (size(tau) /= mn) then
            flag = 2
        else if (size(jpvt) /= n) then
            flag = 3
        else if (size(b) /= maxmn) then
            flag = 4
        end if
        if (flag /= 0) then
            ! ERROR: One of the input arrays is not sized correctly
            write(errmsg, '(AI0A)') "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("solve_qr_pivot_vec", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Workspace Query
        call rz_factor(a(1:mn,:), a(1:mn,1), olwork = lwork1)
        call mult_rz(.true., n, a(1:mn,:), a(1:mn,1), b(1:n), olwork = lwork2)
        lwork = max(lwork1, lwork2, 2 * mn + 1) + mn
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (size(work) < lwork) then
                ! ERROR: WORK not sized correctly
                call errmgr%report_error("solve_qr_no_pivot_mtx", &
                    "Incorrectly sized input array WORK, argument 5.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lwork)
        else
            allocate(wrk(lwork), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("solve_qr_pivot_vec", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if

        ! Determine the rank of R11 using an incremental condition estimation
        wptr(ismin) = one
        wptr(ismax) = one
        smax = abs(a(1,1))
        smin = smax
        if (abs(a(1,1)) == zero) then
            rnk = 0
            b(maxmn) = zero
            return
        else
            rnk = 1
        end if
        do
            if (rnk < mn) then
                i = rnk + 1
                call DLAIC1(imin, rnk, wptr(ismin:ismin+rnk-1), smin, &
                    a(1:rnk-1,i), a(i,i), sminpr, s1, c1)
                call DLAIC1(imax, rnk, wptr(ismax:ismax+rnk-1), smax, &
                    a(1:rnk-1,i), a(i,i), smaxpr, s2, c2)
                if (smaxpr * rcond <= sminpr) then
                    do i = 1, rnk
                        wptr(ismin+i-1) = s1 * wptr(ismin+i-1)
                        wptr(ismax+i-1) = s2 * wptr(ismax+i-1)
                    end do
                    wptr(ismin+rnk) = c1
                    wptr(ismax+rnk) = c2
                    smin = sminpr
                    smax = smaxpr
                    rnk = rnk + 1
                    cycle
                end if
            end if
            exit
        end do

        ! Partition R = [R11 R12]
        !               [ 0  R22]
        tau2 => wptr(1:rnk)
        w => wptr(rnk+1:lwork)
        if (rnk < n) call rz_factor(a(1:rnk,:), tau2, w)

        ! Compute B(1:m,1:NRHS) = Q**T * B(1:M,1:NRHS)
        call mult_qr(.true., a, tau, b(1:m))

        ! Solve the triangular system T11 * B(1:rnk) = B(1:rnk)
        call solve_triangular_system(.true., .false., .true., a(1:rnk,1:rnk), &
            b(1:rnk))
        if (n > rnk) b(rnk+1:n) = zero

        ! Compute B(1:n) = Y**T * B(1:n)
        if (rnk < n) then
            call mult_rz(.true., n - rnk, a(1:rnk,:), tau2, b(1:n), w)
        end if

        ! Apply the pivoting: B(1:N) = P * B(1:N)
        do i = 1, n
            wptr(jpvt(i)) = b(i)
        end do
        b = wptr(1:n)
    end subroutine

! ******************************************************************************
! CHOLESKY SOLVE
! ------------------------------------------------------------------------------
    !> @brief Solves a system of Cholesky factored equations.
    !!
    !! @param[in] upper Set to true if the original matrix A was factored such
    !!  that A = U**T * U; else, set to false if the factorization of A was
    !!  A = L**T * L.
    !! @param[in] a The N-by-N Cholesky factored matrix.
    !! @param[in,out] b On input, the N-by-NRHS right-hand-side matrix B.  On
    !!  output, the solution matrix X.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - LA_ARRAY_SIZE_ERROR: Occurs if any of the input array sizes are 
    !!      incorrect.
    !!
    !! @par Notes
    !! This routine utilizes the LAPACK routine DPOTRS.
    subroutine solve_cholesky_mtx(upper, a, b, err)
        ! Arguments
        logical, intent(in) :: upper
        real(real64), intent(in), dimension(:,:) :: a
        real(real64), intent(inout), dimension(:,:) :: b
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        character :: uplo
        integer(int32) :: n, nrhs, flag
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg

        ! Initialization
        n = size(a, 1)
        nrhs = size(b, 2)
        if (upper) then
            uplo = 'U'
        else
            uplo = 'L'
        end if
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        flag = 0
        if (size(a, 2) /= n) then
            flag = 2
        else if (size(b, 1) /= n) then
            flag = 3
        end if
        if (flag /= 0) then
            write(errmsg, '(AI0A)') "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("solve_cholesky_mtx", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Process
        call DPOTRS(uplo, n, nrhs, a, n, b, n, flag)
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Solves a system of Cholesky factored equations.
    !!
    !! @param[in] upper Set to true if the original matrix A was factored such
    !!  that A = U**T * U; else, set to false if the factorization of A was
    !!  A = L**T * L.
    !! @param[in] a The N-by-N Cholesky factored matrix.
    !! @param[in,out] b On input, the N-element right-hand-side vector B.  On
    !!  output, the solution vector X.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - LA_ARRAY_SIZE_ERROR: Occurs if any of the input array sizes are 
    !!      incorrect.
    !!
    !! @par Notes
    !! This routine utilizes the LAPACK routine DPOTRS.
    subroutine solve_cholesky_vec(upper, a, b, err)
        ! Arguments
        logical, intent(in) :: upper
        real(real64), intent(in), dimension(:,:) :: a
        real(real64), intent(inout), dimension(:) :: b
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        character :: uplo
        integer(int32) :: n, flag
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg

        ! Initialization
        n = size(a, 1)
        if (upper) then
            uplo = 'U'
        else
            uplo = 'L'
        end if
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        flag = 0
        if (size(a, 2) /= n) then
            flag = 2
        else if (size(b) /= n) then
            flag = 3
        end if
        if (flag /= 0) then
            write(errmsg, '(AI0A)') "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("solve_cholesky_vec", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Process
        call DPOTRS(uplo, n, 1, a, n, b, n, flag)
    end subroutine

! ******************************************************************************
! MATRIX INVERSION ROUTINES
! ------------------------------------------------------------------------------
    !> @brief Computes the inverse of a square matrix.
    !!
    !! @param[in,out] a On input, the N-by-N matrix to invert.  On output, the
    !!  inverted matrix.
    !! @param[out] iwork An optional N-element integer workspace array.
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
    !!  - LA_ARRAY_SIZE_ERROR: Occurs if @p a is not square.  Will also occur if
    !!      incorrectly sized workspace arrays are provided.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
    !!      there is insufficient memory available.
    !!  - LA_SINGULAR_MATRIX_ERROR: Occurs if the input matrix is singular.
    !!
    !! @par Usage
    !! @code {.f90}
    !! ! The following example illustrates how to solve a system of linear
    !! ! equations by matrix inversion.  Notice, this is not a preferred
    !! ! solution technique (use LU factorization instead), but is merely a
    !! ! means of illustrating how to compute the inverse of a square matrix.
    !!
    !! ! Variables
    !! real(real64), dimension(n, n) :: a
    !! real(real64), dimension(n, nrhs) :: b, x
    !!
    !! ! Initialize A and B...
    !!
    !! ! Compute the inverse of A.  The inverse will overwrite the original
    !! ! matrix.
    !! call mtx_inverse(a)
    !!
    !! ! Solve A*X = B as X = inv(A) * B.
    !! x = matmul(a, b)
    !! @endcode
    !!
    !! @par Notes
    !! This routine utilizes the LAPACK routines DGETRF to perform an LU
    !! factorization of the matrix, and DGETRI to invert the LU factored
    !! matrix.
    !!
    !! @par See Also
    !! - [Wikipedia](https://en.wikipedia.org/wiki/Invertible_matrix)
    !! - [Wolfram MathWorld](http://mathworld.wolfram.com/MatrixInverse.html)
    subroutine mtx_inverse(a, iwork, work, olwork, err)
        ! Arguments
        real(real64), intent(inout), dimension(:,:) :: a
        integer(int32), intent(out), target, optional, dimension(:) :: iwork
        real(real64), intent(out), target, optional, dimension(:) :: work
        integer(int32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        integer(int32) :: n, liwork, lwork, istat, flag
        integer(int32), pointer, dimension(:) :: iptr
        integer(int32), allocatable, target, dimension(:) :: iwrk
        real(real64), pointer, dimension(:) :: wptr
        real(real64), allocatable, target, dimension(:) :: wrk
        real(real64), dimension(1) :: temp
        class(errors), pointer :: errmgr
        type(errors), target :: deferr

        ! Initialization
        n = size(a, 1)
        liwork = n
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        if (size(a, 2) /= n) then
            call errmgr%report_error("mtx_inverse", &
                "The matrix must be squre to invert.", LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Workspace Query
        call DGETRI(n, a, n, istat, temp, -1, flag)
        lwork = int(temp(1), int32)
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Workspace Allocation
        if (present(work)) then
            if (size(work) < lwork) then
                ! ERROR: WORK not sized correctly
                call errmgr%report_error("svd", &
                    "Incorrectly sized input array WORK, argument 3.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lwork)
        else
            allocate(wrk(lwork), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("mtx_inverse", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if

        ! Integer Workspace Allocation
        if (present(iwork)) then
            if (size(iwork) < liwork) then
                ! ERROR: IWORK not sized correctly
                call errmgr%report_error("svd", &
                    "Incorrectly sized input array IWORK, argument 2.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            iptr => iwork(1:liwork)
        else
            allocate(iwrk(liwork), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("mtx_inverse", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            iptr => iwrk
        end if

        ! Compute the LU factorization of A
        call DGETRF(n, n, a, n, iptr, flag)

        ! Compute the inverse of the LU factored matrix
        call DGETRI(n, a, n, iptr, wptr, lwork, flag)

        ! Check for a singular matrix
        if (flag > 0) then
            call errmgr%report_error("mtx_inverse", &
                "The matrix is singular; therefore, the inverse could " // &
                "not be computed.", LA_SINGULAR_MATRIX_ERROR)
        end if
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Computes the Moore-Penrose pseudo-inverse of a M-by-N matrix
    !! using the singular value decomposition of the matrix.
    !!
    !! @param[in,out] a On input, the M-by-N matrix to invert.  The matrix is
    !!  overwritten on output.
    !! @param[out] ainv The N-by-M matrix where the pseudo-inverse of @p a
    !!  will be written.
    !! @param[in] tol An optional input, that if supplied, overrides the default
    !!  tolerance on singular values such that singular values less than this
    !!  tolerance are forced to have a reciprocal of zero, as opposed to 1/S(I).
    !!  The default tolerance is: MAX(M, N) * EPS * MAX(S).  If the supplied
    !!  value is less than a value that causes an overflow, the tolerance
    !!  reverts back to its default value, and the operation continues;
    !!  however, a warning message is issued.
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
    !! @par Usage
    !! @code{.f90}
    !! ! Use the pseudo-inverse to obtain a least-squares solution to the
    !! ! overdetermined problem A*X = B, where A is an M-by-N matrix (M >= N),
    !! ! B is an M-by-NRHS matrix, and X is an N-by-NRHS matrix.
    !!
    !! ! Variables
    !! real(real64), dimension(m, n) :: a
    !! real(real64), dimension(n, m) :: ainv
    !! real(real64), dimension(m, nrhs) :: b
    !! real(real64), dimension(n, nrhs) :: x
    !!
    !! ! Initialize A, and B...
    !!
    !! ! Compute the pseudo-inverse of A.  Let the subroutine allocate its
    !! ! own workspace array.
    !! call mtx_pinverse(a, ainv)
    !!
    !! ! Compute X = AINV * B to obtain the solution.
    !! x = matmul(ainv, b)
    !! @endcode
    !!
    !! @par See Also
    !! - [Wikipedia](https://en.wikipedia.org/wiki/Moore%E2%80%93Penrose_pseudoinverse)
    !! - [Wolfram MathWorld](http://mathworld.wolfram.com/Moore-PenroseMatrixInverse.html)
    !! - [MathWorks](http://www.mathworks.com/help/matlab/ref/pinv.html?s_tid=srchtitle)
    subroutine mtx_pinverse(a, ainv, tol, work, olwork, err)
        ! Arguments
        real(real64), intent(inout), dimension(:,:) :: a
        real(real64), intent(out), dimension(:,:) :: ainv
        real(real64), intent(in), optional :: tol
        real(real64), intent(out), target, dimension(:), optional :: work
        integer(int32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! External Function Interfaces
        interface
            function DLAMCH(cmach) result(x)
                use, intrinsic :: iso_fortran_env, only : real64
                character, intent(in) :: cmach
                real(real64) :: x
            end function
        end interface

        ! Parameters
        real(real64), parameter :: zero = 0.0d0
        real(real64), parameter :: one = 1.0d0

        ! Local Variables
        integer(int32) :: i, m, n, mn, lwork, istat, flag, i1, i2a, i2b, i3a, &
            i3b, i4
        real(real64), pointer, dimension(:) :: s, wptr, w
        real(real64), pointer, dimension(:,:) :: u, vt
        real(real64), allocatable, target, dimension(:) :: wrk
        real(real64), dimension(1) :: temp
        real(real64) :: t, tref, tolcheck
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg

        ! Initialization
        m = size(a, 1)
        n = size(a, 2)
        mn = min(m, n)
        i1 = m * mn
        i2a = i1 + 1
        i2b = i2a + n * n - 1
        i3a = i2b + 1
        i3b = i3a + mn - 1
        i4 = i3b + 1
        tolcheck = dlamch('s')
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        if (size(ainv, 1) /= n .or. size(ainv, 2) /= m) then
            write(errmsg, '(AI0AI0A)') &
                "The output matrix AINV is not sized appropriately.  " // &
                "It is expected to be ", n, "-by-", m, "."
            call errmgr%report_error("mtx_pinverse", errmsg, &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Workspace Query
        call DGESVD('S', 'A', m, n, a, m, a(1:mn,:), a, m, a, n, temp, -1, flag)
        lwork = int(temp(1), int32)
        lwork = lwork + m * mn + n * n + mn
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (size(work) < lwork) then
                ! ERROR: WORK not sized correctly
                call errmgr%report_error("mtx_pinverse", &
                    "Incorrectly sized input array WORK, argument 4.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lwork)
        else
            allocate(wrk(lwork), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("mtx_pinverse", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if
        u(1:m,1:mn) => wptr(1:i1)
        vt(1:n,1:n) => wptr(i2a:i2b)
        s => wptr(i3a:i3b)
        w => wptr(i4:lwork)

        ! Compute the SVD of A
        call DGESVD('S', 'A', m, n, a, m, s, u, m, vt, n, w, size(w), flag)

        ! Check for convergence
        if (flag > 0) then
            write(errmsg, '(I0A)') flag, " superdiagonals could not " // &
                "converge to zero as part of the QR iteration process."
            call errmgr%report_warning("mtx_pinverse", errmsg, &
                LA_CONVERGENCE_ERROR)
            return
        end if

        ! Determine the threshold tolerance for the singular values such that
        ! singular values less than the threshold result in zero when inverted.
        tref = max(m, n) * epsilon(t) * s(1)
        if (present(tol)) then
            t = tol
        else
            t = tref
        end if
        !if (t < safe_denom(t)) then
        if (t < tolcheck) then
            ! The supplied tolerance is too small, simply fall back to the
            ! default, but issue a warning to the user
            t = tref
            ! call errmgr%report_warning("pinverse_1", "The supplied tolerance was " // &
            !     "smaller than a value that would result in an overflow " // &
            !     "condition, or is negative; therefore, the tolerance has " // &
            !     "been reset to its default value.")
        end if

        ! Compute the pseudoinverse such that pinv(A) = V * inv(S) * U**T by
        ! first computing V * inv(S) (result is N-by-M), and store in the first
        ! MN rows of VT in a transposed manner.
        do i = 1, mn
            ! Apply 1 / S(I) to VT(I,:)
            if (s(i) < t) then
                vt(i,:) = zero
            else
                call recip_mult_array(s(i), vt(i,1:n))
            end if
        end do

        ! Compute (VT**T * inv(S)) * U**T
        call mtx_mult(.true., .true., one, vt(1:mn,:), u, zero, ainv)
    end subroutine

! ******************************************************************************
! LEAST SQUARES SOLUTION ROUTINES
! ------------------------------------------------------------------------------
    !> @brief Solves the overdetermined or underdetermined system (A*X = B) of
    !! M equations of N unknowns using a QR or LQ factorization of the matrix A.
    !! Notice, it is assumed that matrix A has full rank.
    !!
    !! @param[in,out] a On input, the M-by-N matrix A.  On output, if M >= N,
    !!  the QR factorization of A in the form as output by qr_factor; else,
    !!  if M < N, the LQ factorization of A.
    !! @param[in,out] b If M >= N, the M-by-NRHS matrix B.  On output, the first
    !!  N rows contain the N-by-NRHS solution matrix X.  If M < N, an
    !!  N-by-NRHS matrix with the first M rows containing the matrix B.  On
    !!  output, the N-by-NRHS solution matrix X.
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
    !!  - LA_INVALID_OPERATION_ERROR: Occurs if @p a is not of full rank.
    !!
    !! @par Notes
    !! This routine utilizes the LAPACK routine DGELS.
    subroutine solve_least_squares_mtx(a, b, work, olwork, err)
        ! Arguments
        real(real64), intent(inout), dimension(:,:) :: a, b
        real(real64), intent(out), target, optional, dimension(:) :: work
        integer(int32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        integer(int32) :: m, n, maxmn, nrhs, lwork, istat, flag
        real(real64), pointer, dimension(:) :: wptr
        real(real64), allocatable, target, dimension(:) :: wrk
        real(real64), dimension(1) :: temp
        class(errors), pointer :: errmgr
        type(errors), target :: deferr

        ! Initialization
        m = size(a, 1)
        n = size(a, 2)
        maxmn = max(m, n)
        nrhs = size(b, 2)
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        if (size(b, 1) /= maxmn) then
            call errmgr%report_error("solve_least_squares_mtx", &
                "Input 2 is not sized correctly.", LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Workspace Query
        call DGELS('N', m, n, nrhs, a, m, b, maxmn, temp, -1, flag)
        lwork = int(temp(1), int32)
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (size(work) < lwork) then
                ! ERROR: WORK not sized correctly
                call errmgr%report_error("solve_least_squares_mtx", &
                    "Incorrectly sized input array WORK, argument 3.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lwork)
        else
            allocate(wrk(lwork), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("solve_least_squares_mtx", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if

        ! Process
        call DGELS('N', m, n, nrhs, a, m, b, maxmn, wptr, lwork, flag)
        if (flag > 0) then
            call errmgr%report_error("solve_least_squares_mtx", &
                "The supplied matrix is not of full rank; therefore, " // &
                "the solution could not be computed via this routine.  " // &
                "Try a routine that utilizes column pivoting.", &
                LA_INVALID_OPERATION_ERROR)
        end if
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Solves the overdetermined or underdetermined system (A*X = B) of
    !! M equations of N unknowns using a QR or LQ factorization of the matrix A.
    !! Notice, it is assumed that matrix A has full rank.
    !!
    !! @param[in,out] a On input, the M-by-N matrix A.  On output, if M >= N,
    !!  the QR factorization of A in the form as output by qr_factor; else,
    !!  if M < N, the LQ factorization of A.
    !! @param[in,out] b If M >= N, the M-element array B.  On output, the first
    !!  N elements contain the N-element solution array X.  If M < N, an
    !!  N-element array with the first M elements containing the array B.  On
    !!  output, the N-element solution array X.
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
    !!  - LA_INVALID_OPERATION_ERROR: Occurs if @p a is not of full rank.
    !!
    !! @par Notes
    !! This routine utilizes the LAPACK routine DGELS.
    subroutine solve_least_squares_vec(a, b, work, olwork, err)
        ! Arguments
        real(real64), intent(inout), dimension(:,:) :: a
        real(real64), intent(inout), dimension(:) :: b
        real(real64), intent(out), target, optional, dimension(:) :: work
        integer(int32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        integer(int32) :: m, n, maxmn, lwork, istat, flag
        real(real64), pointer, dimension(:) :: wptr
        real(real64), allocatable, target, dimension(:) :: wrk
        real(real64), dimension(1) :: temp
        class(errors), pointer :: errmgr
        type(errors), target :: deferr

        ! Initialization
        m = size(a, 1)
        n = size(a, 2)
        maxmn = max(m, n)
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        if (size(b) /= maxmn) then
            call errmgr%report_error("solve_least_squares_vec", &
                "Input 2 is not sized correctly.", LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Workspace Query
        call DGELS('N', m, n, 1, a, m, b, maxmn, temp, -1, flag)
        lwork = int(temp(1), int32)
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (size(work) < lwork) then
                ! ERROR: WORK not sized correctly
                call errmgr%report_error("solve_least_squares_vec", &
                    "Incorrectly sized input array WORK, argument 3.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lwork)
        else
            allocate(wrk(lwork), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("solve_least_squares_vec", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if

        ! Process
        call DGELS('N', m, n, 1, a, m, b, maxmn, wptr, lwork, flag)
        if (flag > 0) then
            call errmgr%report_error("solve_least_squares_mtx", &
                "The supplied matrix is not of full rank; therefore, " // &
                "the solution could not be computed via this routine.  " // &
                "Try a routine that utilizes column pivoting.", &
                LA_INVALID_OPERATION_ERROR)
        end if
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Solves the overdetermined or underdetermined system (A*X = B) of
    !! M equations of N unknowns using a complete orthogonal factorization of
    !! matrix A.
    !!
    !! @param[in,out] a On input, the M-by-N matrix A.  On output, the matrix
    !!  is overwritten by the details of its complete orthogonal factorization.
    !! @param[in,out] b If M >= N, the M-by-NRHS matrix B.  On output, the first
    !!  N rows contain the N-by-NRHS solution matrix X.  If M < N, an
    !!  N-by-NRHS matrix with the first M rows containing the matrix B.  On
    !!  output, the N-by-NRHS solution matrix X.
    !! @param[out] ipvt An optional input that on input, an N-element array 
    !!  that if IPVT(I) .ne. 0, the I-th column of A is permuted to the front 
    !!  of A * P; if IPVT(I) = 0, the I-th column of A is a free column.  On 
    !!  output, if IPVT(I) = K, then the I-th column of A * P was the K-th 
    !!  column of A.  If not supplied, memory is allocated internally, and IPVT
    !!  is set to all zeros such that all columns are treated as free.
    !! @param[out] arnk An optional output, that if provided, will return the
    !!  rank of @p a.
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
    !!
    !! @par Notes
    !! This routine utilizes the LAPACK routine DGELSY.
    subroutine solve_least_squares_mtx_pvt(a, b, ipvt, arnk, work, olwork, err)
        ! Arguments
        real(real64), intent(inout), dimension(:,:) :: a, b
        integer(int32), intent(inout), target, optional, dimension(:) :: ipvt
        integer(int32), intent(out), optional :: arnk
        real(real64), intent(out), target, optional, dimension(:) :: work
        integer(int32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        integer(int32) :: m, n, maxmn, nrhs, lwork, istat, flag, rnk
        real(real64), pointer, dimension(:) :: wptr
        real(real64), allocatable, target, dimension(:) :: wrk
        integer(int32), allocatable, target, dimension(:) :: iwrk
        integer(int32), pointer, dimension(:) :: iptr
        real(real64), dimension(1) :: temp
        integer(int32), dimension(1) :: itemp
        real(real64) :: rc
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg

        ! Initialization
        m = size(a, 1)
        n = size(a, 2)
        maxmn = max(m, n)
        nrhs = size(b, 2)
        rc = epsilon(rc)
        if (present(arnk)) arnk = 0
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        flag = 0
        if (size(b, 1) /= maxmn) then
            flag = 2
        else if (size(ipvt) /= n) then
            flag = 3
        end if
        if (flag /= 0) then
            write(errmsg, '(AI0A)') "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("solve_least_squares_mtx_pvt", &
                trim(errmsg), LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Workspace Query
        call DGELSY(m, n, nrhs, a, m, b, maxmn, itemp, rc, rnk, temp, -1, flag)
        lwork = int(temp(1), int32)
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Local Memory Allocation
        if (present(ipvt)) then
            if (size(ipvt) < n) then
                ! ERROR: IPVT is not big enough
                call errmgr%report_error("solve_least_squares_mtx_pvt", &
                    "Incorrectly sized pivot array, argument 3.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            iptr => ipvt(1:n)
        else
            allocate(iwrk(n), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("solve_least_squares_mtx_pvt", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            iptr => iwrk
            iptr = 0
        end if

        if (present(work)) then
            if (size(work) < lwork) then
                ! ERROR: WORK not sized correctly
                call errmgr%report_error("solve_least_squares_mtx_pvt", &
                    "Incorrectly sized input array WORK, argument 5.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lwork)
        else
            allocate(wrk(lwork), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("solve_least_squares_mtx_pvt", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if

        ! Process
        call DGELSY(m, n, nrhs, a, m, b, maxmn, iptr, rc, rnk, wptr, lwork, &
            flag)
        if (present(arnk)) arnk = rnk
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Solves the overdetermined or underdetermined system (A*X = B) of
    !! M equations of N unknowns using a complete orthogonal factorization of
    !! matrix A.
    !!
    !! @param[in,out] a On input, the M-by-N matrix A.  On output, the matrix
    !!  is overwritten by the details of its complete orthogonal factorization.
    !! @param[in,out] b If M >= N, the M-element array B.  On output, the first
    !!  N elements contain the N-element solution array X.  If M < N, an
    !!  N-element array with the first M elements containing the array B.  On
    !!  output, the N-element solution array X.
    !! @param[out] ipvt An optional input that on input, an N-element array 
    !!  that if IPVT(I) .ne. 0, the I-th column of A is permuted to the front 
    !!  of A * P; if IPVT(I) = 0, the I-th column of A is a free column.  On 
    !!  output, if IPVT(I) = K, then the I-th column of A * P was the K-th 
    !!  column of A.  If not supplied, memory is allocated internally, and IPVT
    !!  is set to all zeros such that all columns are treated as free.
    !! @param[out] arnk An optional output, that if provided, will return the
    !!  rank of @p a.
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
    !!
    !! @par Notes
    !! This routine utilizes the LAPACK routine DGELSY.
    subroutine solve_least_squares_vec_pvt(a, b, ipvt, arnk, work, olwork, err)
        ! Arguments
        real(real64), intent(inout), dimension(:,:) :: a
        real(real64), intent(inout), dimension(:) :: b
        integer(int32), intent(inout), target, optional, dimension(:) :: ipvt
        integer(int32), intent(out), optional :: arnk
        real(real64), intent(out), target, optional, dimension(:) :: work
        integer(int32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        integer(int32) :: m, n, maxmn, lwork, istat, flag, rnk
        real(real64), pointer, dimension(:) :: wptr
        real(real64), allocatable, target, dimension(:) :: wrk
        integer(int32), allocatable, target, dimension(:) :: iwrk
        integer(int32), pointer, dimension(:) :: iptr
        real(real64), dimension(1) :: temp
        integer(int32), dimension(1) :: itemp
        real(real64) :: rc
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg

        ! Initialization
        m = size(a, 1)
        n = size(a, 2)
        maxmn = max(m, n)
        rc = epsilon(rc)
        if (present(arnk)) arnk = 0
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        flag = 0
        if (size(b, 1) /= maxmn) then
            flag = 2
        else if (size(ipvt) /= n) then
            flag = 3
        end if
        if (flag /= 0) then
            write(errmsg, '(AI0A)') "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("solve_least_squares_vec_pvt", &
                trim(errmsg), LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Workspace Query
        call DGELSY(m, n, 1, a, m, b, maxmn, itemp, rc, rnk, temp, -1, flag)
        lwork = int(temp(1), int32)
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Local Memory Allocation
        if (present(ipvt)) then
            if (size(ipvt) < n) then
                ! ERROR: IPVT is not big enough
                call errmgr%report_error("solve_least_squares_mtx_pvt", &
                    "Incorrectly sized pivot array, argument 3.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            iptr => ipvt(1:n)
        else
            allocate(iwrk(n), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("solve_least_squares_mtx_pvt", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            iptr => iwrk
            iptr = 0
        end if

        if (present(work)) then
            if (size(work) < lwork) then
                ! ERROR: WORK not sized correctly
                call errmgr%report_error("solve_least_squares_vec_pvt", &
                    "Incorrectly sized input array WORK, argument 5.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lwork)
        else
            allocate(wrk(lwork), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("solve_least_squares_vec_pvt", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if

        ! Process
        call DGELSY(m, n, 1, a, m, b, maxmn, ipvt, rc, rnk, wptr, lwork, flag)
        if (present(arnk)) arnk = rnk
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Solves the overdetermined or underdetermined system (A*X = B) of
    !! M equations of N unknowns using a singular value decomposition of
    !! matrix A.
    !!
    !! @param[in,out] a On input, the M-by-N matrix A.  On output, the matrix
    !!  is overwritten by the details of its complete orthogonal factorization.
    !! @param[in,out] b If M >= N, the M-by-NRHS matrix B.  On output, the first
    !!  N rows contain the N-by-NRHS solution matrix X.  If M < N, an
    !!  N-by-NRHS matrix with the first M rows containing the matrix B.  On
    !!  output, the N-by-NRHS solution matrix X.
    !! @param[out] arnk An optional output, that if provided, will return the
    !!  rank of @p a.
    !! @param[out] s A MIN(M, N)-element array that on output contains the
    !!  singular values of @p a in descending order.  Notice, the condition
    !!  number of @p a can be determined by S(1) / S(MIN(M, N)).
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
    !! @par Notes
    !! This routine utilizes the LAPACK routine DGELSS.
    subroutine solve_least_squares_mtx_svd(a, b, arnk, s, work, olwork, err)
        ! Arguments
        real(real64), intent(inout), dimension(:,:) :: a, b
        real(real64), intent(out), dimension(:) :: s
        integer(int32), intent(out), optional :: arnk
        real(real64), intent(out), target, optional, dimension(:) :: work
        integer(int32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        integer(int32) :: m, n, nrhs, mn, maxmn, istat, flag, lwork, rnk
        real(real64), pointer, dimension(:) :: wptr
        real(real64), allocatable, target, dimension(:) :: wrk
        real(real64), dimension(1) :: temp
        real(real64) :: rcond
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg

        ! Initialization
        m = size(a, 1)
        n = size(a, 2)
        nrhs = size(b, 2)
        mn = min(m, n)
        maxmn = max(m, n)
        rcond = epsilon(rcond)
        if (present(arnk)) arnk = 0
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        flag = 0
        if (size(b, 1) /= maxmn) then
            flag = 2
        else if (size(s) /= mn) then
            flag = 4
        end if
        if (flag /= 0) then
            ! ERROR: One of the input arrays is not sized correctly
            write(errmsg, '(AI0A)') "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("solve_least_squares_mtx_svd", &
                trim(errmsg), LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Workspace Query
        call DGELSS(m, n, nrhs, a, m, b, maxmn, s, rcond, rnk, temp, -1, flag)
        lwork = int(temp(1), int32)
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (size(work) < lwork) then
                ! ERROR: WORK not sized correctly
                call errmgr%report_error("solve_least_squares_mtx_svd", &
                    "Incorrectly sized input array WORK, argument 5.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lwork)
        else
            allocate(wrk(lwork), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("solve_least_squares_mtx_svd", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if

        ! Process
        call DGELSS(m, n, nrhs, a, m, b, maxmn, s, rcond, rnk, wptr, lwork, &
            flag)
        if (present(arnk)) arnk = rnk
        if (flag > 0) then
            write(errmsg, '(I0A)') flag, " superdiagonals could not " // &
                "converge to zero as part of the QR iteration process."
            call errmgr%report_warning("solve_least_squares_mtx_svd", errmsg, &
                LA_CONVERGENCE_ERROR)
        end if
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Solves the overdetermined or underdetermined system (A*X = B) of
    !! M equations of N unknowns using a singular value decomposition of
    !! matrix A.
    !!
    !! @param[in,out] a On input, the M-by-N matrix A.  On output, the matrix
    !!  is overwritten by the details of its complete orthogonal factorization.
    !! @param[in,out] b If M >= N, the M-by-NRHS matrix B.  On output, the first
    !!  N rows contain the N-by-NRHS solution matrix X.  If M < N, an
    !!  N-by-NRHS matrix with the first M rows containing the matrix B.  On
    !!  output, the N-by-NRHS solution matrix X.
    !! @param[out] arnk An optional output, that if provided, will return the
    !!  rank of @p a.
    !! @param[out] s A MIN(M, N)-element array that on output contains the
    !!  singular values of @p a in descending order.  Notice, the condition
    !!  number of @p a can be determined by S(1) / S(MIN(M, N)).
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
    !! @par Notes
    !! This routine utilizes the LAPACK routine DGELSS.
    subroutine solve_least_squares_vec_svd(a, b, arnk, s, work, olwork, err)
        ! Arguments
        real(real64), intent(inout), dimension(:,:) :: a
        real(real64), intent(inout), dimension(:) :: b
        real(real64), intent(out), dimension(:) :: s
        integer(int32), intent(out), optional :: arnk
        real(real64), intent(out), target, optional, dimension(:) :: work
        integer(int32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        integer(int32) :: m, n, mn, maxmn, istat, flag, lwork, rnk
        real(real64), pointer, dimension(:) :: wptr
        real(real64), allocatable, target, dimension(:) :: wrk
        real(real64), dimension(1) :: temp
        real(real64) :: rcond
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg

        ! Initialization
        m = size(a, 1)
        n = size(a, 2)
        mn = min(m, n)
        maxmn = max(m, n)
        rcond = epsilon(rcond)
        if (present(arnk)) arnk = 0
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        flag = 0
        if (size(b) /= maxmn) then
            flag = 2
        else if (size(s) /= mn) then
            flag = 4
        end if
        if (flag /= 0) then
            ! ERROR: One of the input arrays is not sized correctly
            write(errmsg, '(AI0A)') "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("solve_least_squares_vec_svd", &
                trim(errmsg), LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Workspace Query
        call DGELSS(m, n, 1, a, m, b, maxmn, s, rcond, rnk, temp, -1, flag)
        lwork = int(temp(1), int32)
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (size(work) < lwork) then
                ! ERROR: WORK not sized correctly
                call errmgr%report_error("solve_least_squares_vec_svd", &
                    "Incorrectly sized input array WORK, argument 5.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lwork)
        else
            allocate(wrk(lwork), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("solve_least_squares_vec_svd", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if

        ! Process
        call DGELSS(m, n, 1, a, m, b, maxmn, s, rcond, rnk, wptr, lwork, flag)
        if (present(arnk)) arnk = rnk
        if (flag > 0) then
            write(errmsg, '(I0A)') flag, " superdiagonals could not " // &
                "converge to zero as part of the QR iteration process."
            call errmgr%report_warning("solve_least_squares_vec_svd", errmsg, &
                LA_CONVERGENCE_ERROR)
        end if
    end subroutine


end module
