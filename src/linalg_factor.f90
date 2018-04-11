! linalg_factor.f90

!> @brief \b linalg_factor
!!
!! @par Purpose
!! Provides a set of matrix factorization routines.
module linalg_factor
    use ferror, only : errors
    use linalg_constants
    use linalg_core, only : swap
    implicit none
    private
    public :: lu_factor
    public :: form_lu
    public :: qr_factor
    public :: form_qr
    public :: mult_qr
    public :: qr_rank1_update
    public :: cholesky_factor
    public :: cholesky_rank1_update
    public :: cholesky_rank1_downdate
    public :: rz_factor
    public :: mult_rz
    public :: svd

! ******************************************************************************
! INTERFACES
! ------------------------------------------------------------------------------
    !> @brief Extracts the L and U matrices from the condensed [L\\U] storage 
    !! format used by the @ref lu_factor.
    !!
    !! @par Usage
    !! The following example illustrates how to extract the L, U, and P matrices
    !! in order to solve a system of LU factored equations.
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
    interface form_lu
        module procedure :: form_lu_all
        module procedure :: form_lu_only
    end interface

! ------------------------------------------------------------------------------
    !> @brief Computes the QR factorization of an M-by-N matrix.
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
    !! - [Wolfram MathWorld](http://mathworld.wolfram.com/QRDecomposition.html)
    !! - [LAPACK Users Manual](http://netlib.org/lapack/lug/node39.html)
    interface qr_factor
        module procedure :: qr_factor_no_pivot
        module procedure :: qr_factor_pivot
    end interface

! ------------------------------------------------------------------------------
    !> @brief Forms the full M-by-M orthogonal matrix Q from the elementary
    !! reflectors returned by the base QR factorization algorithm.
    !!
    !! @par Usage
    !! The following example illustrates how to explicitly form the Q and R
    !! matrices from the output of qr_factor, and then use the resulting 
    !! matrices to solve a system of linear equations.
    !! @code{.f90}
    !! program example
    !!     use iso_fortran_env, only : real64, int32
    !!     use linalg_factor, only : qr_factor, form_qr
    !!     use linalg_solve, only : solve_triangular_system
    !!     implicit none
    !!
    !!     ! Variables
    !!     real(real64) :: a(3,3), b(3), q(3,3), tau(3)
    !!     integer(int32) :: i
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
    !!     ! Compute the QR factorization without column pivoting
    !!     call qr_factor(a, tau)
    !!
    !!     ! Build Q and R.  A is overwritten with R
    !!     call form_qr(a, tau, q)
    !!
    !!     ! As this system is square, matrix R is upper triangular.  Also, Q is
    !!     ! always orthogonal such that it's inverse and transpose are equal.  As the
    !!     ! system is now factored, its form is: Q * R * X = B.  Solving this system
    !!     ! is then as simple as solving the upper triangular system: 
    !!     ! R * X = Q**T * B.
    !!
    !!     ! Compute Q**T * B, and store the results in B
    !!     b = matmul(transpose(q), b)
    !!
    !!     ! Solve the upper triangular system R * X = Q**T * B for X
    !!     call solve_triangular_system(.true., .false., .true., a, b)
    !!
    !!     ! Display the results
    !!     print '(A)', "QR Solution: X = "
    !!     print '(F8.4)', (b(i), i = 1, size(b))
    !!
    !!     ! Notice, QR factorization with column pivoting could be accomplished via
    !!     ! a similar approach, but the column pivoting would need to be accounted
    !!     ! for by noting that Q * R = A * P, where P is an N-by-N matrix describing
    !!     ! the column pivoting operations.
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
    interface form_qr
        module procedure :: form_qr_no_pivot
        module procedure :: form_qr_pivot
    end interface

! ------------------------------------------------------------------------------
    !> @brief Multiplies a general matrix by the orthogonal matrix Q from a QR
    !! factorization.
    !!
    !! @par Usage
    !! The following example illustrates how to perform the multiplication 
    !! Q**T * B when solving a system of QR factored equations without 
    !! explicitly forming the matrix Q.
    !! @code{.f90}
    !! program example
    !!     use iso_fortran_env, only : real64, int32
    !!     use linalg_factor, only : qr_factor, mult_qr
    !!     use linalg_solve, only : solve_triangular_system
    !!     implicit none
    !!
    !!     ! Variables
    !!     real(real64) :: a(3,3), b(3), tau(3)
    !!     integer(int32) :: i
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
    !!     ! Compute the QR factorization without column pivoting
    !!     call qr_factor(a, tau)
    !!
    !!     ! As this system is square, matrix R is upper triangular.  Also, Q is
    !!     ! always orthogonal such that it's inverse and transpose are equal.  As the
    !!     ! system is now factored, its form is: Q * R * X = B.  Solving this system
    !!     ! is then as simple as solving the upper triangular system: 
    !!     ! R * X = Q**T * B.
    !!
    !!     ! Compute Q**T * B, and store the results in B.  Notice, using mult_qr
    !!     ! avoids direct construction of the full Q and R matrices.
    !!     call mult_qr(.true., a, tau, b)
    !!
    !!     ! Solve the upper triangular system R * X = Q**T * B for X
    !!     call solve_triangular_system(.true., .false., .true., a, b)
    !!
    !!     ! Display the results
    !!     print '(A)', "QR Solution: X = "
    !!     print '(F8.4)', (b(i), i = 1, size(b))
    !!
    !!     ! Notice, QR factorization with column pivoting could be accomplished via
    !!     ! a similar approach, but the column pivoting would need to be accounted
    !!     ! for by noting that Q * R = A * P, where P is an N-by-N matrix describing
    !!     ! the column pivoting operations.
    !! end program
    !! @endcode
    !! The above program produces the following output.
    !! @code{.txt}
    !! QR Solution: X =
    !! 0.3333
    !! -0.6667
    !! 0.0000
    !! @endcode
    interface mult_qr
        module procedure :: mult_qr_mtx
        module procedure :: mult_qr_vec
    end interface

! ------------------------------------------------------------------------------
    !> @brief Multiplies a general matrix by the orthogonal matrix Z from an
    !! RZ factorization.
    interface mult_rz
        module procedure :: mult_rz_mtx
        module procedure :: mult_rz_vec
    end interface

contains
! ******************************************************************************
! LU FACTORIZATION
! ------------------------------------------------------------------------------
    !> @brief Computes the LU factorization of an M-by-N matrix.
    !!
    !! @param[in,out] a On input, the M-by-N matrix on which to operate.  On
    !! output, the LU factored matrix in the form [L\\U] where the unit diagonal
    !! elements of L are not stored.
    !! @param[out] ipvt An MIN(M, N)-element array used to track row-pivot
    !!  operations.  The array stored pivot information such that row I is
    !!  interchanged with row IPVT(I).
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - LA_ARRAY_SIZE_ERROR: Occurs if @p ipvt is not sized appropriately.
    !!  - LA_SINGULAR_MATRIX_ERROR: Occurs as a warning if @p a is found to be
    !!      singular.
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
    !! @par Notes
    !! This routine utilizes the LAPACK routine DGETRF.
    !!
    !! @par See Also
    !! - [Wikipedia](https://en.wikipedia.org/wiki/LU_decomposition)
    !! - [Wolfram MathWorld](http://mathworld.wolfram.com/LUDecomposition.html)
    subroutine lu_factor(a, ipvt, err)
        ! Arguments
        real(dp), intent(inout), dimension(:,:) :: a
        integer(i32), intent(out), dimension(:) :: ipvt
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        real(dp), parameter :: zero = 0.0d0
        real(dp), parameter :: one = 1.0d0

        ! Local Variables
        integer(i32) :: m, n, mn, flag
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg

        ! Initialization
        m = size(a, 1)
        n = size(a, 2)
        mn = min(m, n)
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        flag = 0
        if (size(ipvt) /= mn) then
            ! ERROR: IPVT not sized correctly
            call errmgr%report_error("lu_factor", &
                "Incorrectly sized input array IPVT, argument 2.", &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Compute the LU factorization by calling the LAPACK routine DGETRF
        call DGETRF(m, n, a, m, ipvt, flag)

        ! If flag > 0, the matrix is singular.  Notice, flag should not be
        ! able to be < 0 as we've already verrified inputs prior to making the
        ! call to LAPACK
        if (flag > 0) then
            ! WARNING: Singular matrix
            write(errmsg, '(AI0A)') &
                "Singular matrix encountered (row ", flag, ")"
            call errmgr%report_warning("lu_factor", trim(errmsg), &
                LA_SINGULAR_MATRIX_ERROR)
        end if
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Extracts the L, U, and P matrices from the output of the
    !! @ref lu_factor routine.
    !!
    !! @param[in,out] lu On input, the N-by-N matrix as output by
    !!  @ref lu_factor.  On output, the N-by-N lower triangular matrix L.
    !! @param[in] ipvt The N-element pivot array as output by
    !!  @ref lu_factor.
    !! @param[out] u An N-by-N matrix where the U matrix will be written.
    !! @param[out] p An N-by-N matrix where the row permutation matrix will be
    !!  written.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - LA_ARRAY_SIZE_ERROR: Occurs if any of the input array sizes are 
    !!      incorrect.
    !!
    !! @par Remarks
    !! This routine allows extraction of the actual "L", "U", and "P" matrices
    !! of the decomposition.  To use these matrices to solve the system A*X = B,
    !! the following approach is used.
    !!
    !! 1. First, solve the linear system: L*Y = P*B for Y.
    !! 2. Second, solve the linear system: U*X = Y for X.
    !!
    !! Notice, as both L and U are triangular in structure, the above equations
    !! can be solved by forward and backward substitution.
    !!
    !! @par See Also
    !! - [Wikipedia](https://en.wikipedia.org/wiki/LU_decomposition)
    !! - [Wolfram MathWorld](http://mathworld.wolfram.com/LUDecomposition.html)
    subroutine form_lu_all(lu, ipvt, u, p, err)
        ! Arguments
        real(dp), intent(inout), dimension(:,:) :: lu
        integer(i32), intent(in), dimension(:) :: ipvt
        real(dp), intent(out), dimension(:,:) :: u, p
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        integer(i32) :: j, jp, n, flag
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg

        ! Parameters
        real(dp), parameter :: zero = 0.0d0
        real(dp), parameter :: one = 1.0d0

        ! Initialization
        n = size(lu, 1)
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        flag = 0
        if (size(lu, 2) /= n) then
            flag = 1
        else if (size(ipvt) /= n) then
            flag = 2
        else if (size(u, 1) /= n .or. size(u, 2) /= n) then
            flag = 3
        else if (size(p, 1) /= n .or. size(p, 2) /= n) then
            flag = 4
        end if
        if (flag /= 0) then
            ! One of the input arrays is not sized correctly
            write(errmsg, '(AI0A)') "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("form_lu_all", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
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
    !> @brief Extracts the L, and U matrices from the output of the
    !! @ref lu_factor routine.
    !!
    !! @param[in,out] lu On input, the N-by-N matrix as output by
    !!  @ref lu_factor.  On output, the N-by-N lower triangular matrix L.
    !! @param[out] u An N-by-N matrix where the U matrix will be written.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - LA_ARRAY_SIZE_ERROR: Occurs if any of the input array sizes are 
    !!      incorrect.
    subroutine form_lu_only(lu, u, err)
        ! Arguments
        real(dp), intent(inout), dimension(:,:) :: lu
        real(dp), intent(out), dimension(:,:) :: u
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        integer(i32) :: j, n, flag
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg

        ! Parameters
        real(dp), parameter :: zero = 0.0d0
        real(dp), parameter :: one = 1.0d0

        ! Initialization
        n = size(lu, 1)
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        flag = 0
        if (size(lu, 2) /= n) then
            flag = 2
        else if (size(u, 1) /= n .or. size(u, 2) /= n) then
            flag = 3
        end if
        if (flag /= 0) then
            ! One of the input arrays is not sized correctly
            write(errmsg, '(AI0A)') "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("form_lu_only", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
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

! ******************************************************************************
! QR FACTORIZATION
! ------------------------------------------------------------------------------
    !> @brief Computes the QR factorization of an M-by-N matrix without
    !! pivoting.
    !!
    !! @param[in,out] a On input, the M-by-N matrix to factor.  On output, the
    !!  elements on and above the diagonal contain the MIN(M, N)-by-N upper
    !!  trapezoidal matrix R (R is upper triangular if M >= N).  The elements
    !!  below the diagonal, along with the array @p tau, represent the
    !!  orthogonal matrix Q as a product of elementary reflectors.
    !! @param[out] tau A MIN(M, N)-element array used to store the scalar
    !!  factors of the elementary reflectors.
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
    !!  - LA_ARRAY_SIZE_ERROR: Occurs if @p tau or @p work are not sized 
    !!      appropriately.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
    !!      there is insufficient memory available.
    !!
    !! @remarks
    !! QR factorization without pivoting is best suited to solving an
    !! overdetermined system in least-squares terms, or to solve a normally
    !! defined system.  To solve an underdetermined system, it is recommended to
    !! use either LQ factorization, or a column-pivoting based QR factorization.
    !!
    !! @par Notes
    !! This routine utilizes the LAPACK routine DGEQRF.
    subroutine qr_factor_no_pivot(a, tau, work, olwork, err)
        ! Arguments
        real(dp), intent(inout), dimension(:,:) :: a
        real(dp), intent(out), dimension(:) :: tau
        real(dp), intent(out), target, dimension(:), optional :: work
        integer(i32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        integer(i32) :: m, n, mn, istat, lwork, flag
        real(dp), dimension(1) :: temp
        real(dp), pointer, dimension(:) :: wptr
        real(dp), allocatable, target, dimension(:) :: wrk
        class(errors), pointer :: errmgr
        type(errors), target :: deferr

        ! Initialization
        m = size(a, 1)
        n = size(a, 2)
        mn = min(m, n)
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        if (size(tau) /= mn) then
            ! ERROR: TAU not sized correctly
            call errmgr%report_error("qr_factor_no_pivot", &
                "Incorrectly sized input array TAU, argument 2.", &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Workspace Query
        call DGEQRF(m, n, a, m, tau, temp, -1, flag)
        lwork = int(temp(1), i32)
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (size(work) < lwork) then
                ! ERROR: WORK not sized correctly
                call errmgr%report_error("qr_factor_no_pivot", &
                    "Incorrectly sized input array WORK, argument 3.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lwork)
        else
            allocate(wrk(lwork), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("qr_factor_no_pivot", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if

        ! Call DGEQRF
        call DGEQRF(m, n, a, m, tau, wptr, lwork, flag)
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Computes the QR factorization of an M-by-N matrix with column
    !! pivoting such that A * P = Q * R.
    !!
    !! @param[in,out] a On input, the M-by-N matrix to factor.  On output, the
    !!  elements on and above the diagonal contain the MIN(M, N)-by-N upper
    !!  trapezoidal matrix R (R is upper triangular if M >= N).  The elements
    !!  below the diagonal, along with the array @p tau, represent the
    !!  orthogonal matrix Q as a product of elementary reflectors.
    !! @param[out] tau A MIN(M, N)-element array used to store the scalar
    !!  factors of the elementary reflectors.
    !! @param[in,out] jpvt On input, an N-element array that if JPVT(I) .ne. 0,
    !!  the I-th column of A is permuted to the front of A * P; if JPVT(I) = 0,
    !!  the I-th column of A is a free column.  On output, if JPVT(I) = K, then
    !!  the I-th column of A * P was the K-th column of A.
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
    !! This routine utilizes the LAPACK routine DGEQP3.
    subroutine qr_factor_pivot(a, tau, jpvt, work, olwork, err)
        ! Arguments
        real(dp), intent(inout), dimension(:,:) :: a
        real(dp), intent(out), dimension(:) :: tau
        integer(i32), intent(inout), dimension(:) :: jpvt
        real(dp), intent(out), target, dimension(:), optional :: work
        integer(i32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        integer(i32) :: m, n, mn, istat, lwork, flag
        real(dp), dimension(1) :: temp
        real(dp), pointer, dimension(:) :: wptr
        real(dp), allocatable, target, dimension(:) :: wrk
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg

        ! Initialization
        m = size(a, 1)
        n = size(a, 2)
        mn = min(m, n)
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
        end if
        if (flag /= 0) then
            ! ERROR: One of the input arrays is not sized correctly
            write(errmsg, '(AI0A)') "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("qr_factor_pivot", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Workspace Query
        call DGEQP3(m, n, a, m, jpvt, tau, temp, -1, flag)
        lwork = int(temp(1), i32)
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (size(work) < lwork) then
                ! ERROR: WORK not sized correctly
                call errmgr%report_error("qr_factor_pivot", &
                    "Incorrectly sized input array WORK, argument 4.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lwork)
        else
            allocate(wrk(lwork), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("qr_factor_pivot", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if

        ! Call DGEQP3
        call DGEQP3(m, n, a, m, jpvt, tau, wptr, lwork, flag)

        ! End
        if (allocated(wrk)) deallocate(wrk)
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Forms the full M-by-M orthogonal matrix Q from the elementary
    !! reflectors returned by the base QR factorization algorithm.
    !!
    !! @param[in,out] r On input, an M-by-N matrix where the elements below the
    !!  diagonal contain the elementary reflectors generated from the QR
    !!  factorization.  On and above the diagonal, the matrix contains the
    !!  matrix R.  On output, the elements below the diagonal are zeroed such
    !!  that the remaining matrix is simply the M-by-N matrix R.
    !! @param[in] tau A MIN(M, N)-element array containing the scalar factors of
    !!  each elementary reflector defined in @p r.
    !! @param[out] q An M-by-M matrix where the full orthogonal matrix Q will be
    !!  written.  In the event that M > N, Q may be supplied as M-by-N, and
    !!  therefore only return the useful submatrix Q1 (Q = [Q1, Q2]) as the
    !!  factorization can be written as Q * R = [Q1, Q2] * [R1; 0].
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
    !! This routine utilizes the LAPACK routine DORGQR.
    subroutine form_qr_no_pivot(r, tau, q, work, olwork, err)
        ! Arguments
        real(dp), intent(inout), dimension(:,:) :: r
        real(dp), intent(in), dimension(:) :: tau
        real(dp), intent(out), dimension(:,:) :: q
        real(dp), intent(out), target, dimension(:), optional :: work
        integer(i32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        real(dp), parameter :: zero = 0.0d0

        ! Local Variables
        integer(i32) :: j, m, n, mn, qcol, istat, flag, lwork
        real(dp), pointer, dimension(:) :: wptr
        real(dp), allocatable, target, dimension(:) :: wrk
        real(dp), dimension(1) :: temp
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg

        ! Initialization
        m = size(r, 1)
        n = size(r, 2)
        mn = min(m, n)
        qcol = size(q, 2)
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        flag = 0
        if (size(tau) /= mn) then
            flag = 2
        else if (size(q, 1) /= m .or. (qcol /= m .and. qcol /= n)) then
            flag = 3
        else if (qcol == n .and. m < n) then
            flag = 3
        end if
        if (flag /= 0) then
            ! ERROR: One of the input arrays is not sized correctly
            write(errmsg, '(AI0A)') "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("form_qr_no_pivot", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Workspace Query
        call DORGQR(m, qcol, mn, q, m, tau, temp, -1, flag)
        lwork = int(temp(1), i32)
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (size(work) < lwork) then
                ! ERROR: WORK not sized correctly
                call errmgr%report_error("form_qr_no_pivot", &
                    "Incorrectly sized input array WORK, argument 4.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lwork)
        else
            allocate(wrk(lwork), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("form_qr_no_pivot", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if

        ! Copy the sub-diagonal portion of R to Q, and then zero out the
        ! sub-diagonal portion of R
        do j = 1, mn
            q(j+1:m,j) = r(j+1:m,j)
            r(j+1:m,j) = zero
        end do

        ! Build Q - Build M-by-M or M-by-N, but M-by-N only for M >= N
        call DORGQR(m, qcol, mn, q, m, tau, wptr, lwork, flag)

        ! End
        if (allocated(wrk)) deallocate(wrk)
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Forms the full M-by-M orthogonal matrix Q from the elementary
    !! reflectors returned by the base QR factorization algorithm.
    !!
    !! @param[in,out] r On input, an M-by-N matrix where the elements below the
    !!  diagonal contain the elementary reflectors generated from the QR
    !!  factorization.  On and above the diagonal, the matrix contains the
    !!  matrix R.  On output, the elements below the diagonal are zeroed such
    !!  that the remaining matrix is simply the M-by-N matrix R.
    !! @param[in] tau A MIN(M, N)-element array containing the scalar factors of
    !!  each elementary reflector defined in @p r.
    !! @param[in] pvt An N-element column pivot array as returned by the QR
    !!  factorization.
    !! @param[out] q An M-by-M matrix where the full orthogonal matrix Q will be
    !!  written.  In the event that M > N, Q may be supplied as M-by-N, and
    !!  therefore only return the useful submatrix Q1 (Q = [Q1, Q2]) as the
    !!  factorization can be written as Q * R = [Q1, Q2] * [R1; 0].
    !! @param[out] p An N-by-N matrix where the pivot matrix will be written.
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
    !! This routine utilizes the LAPACK routine DORGQR.
    subroutine form_qr_pivot(r, tau, pvt, q, p, work, olwork, err)
        ! Arguments
        real(dp), intent(inout), dimension(:,:) :: r
        real(dp), intent(in), dimension(:) :: tau
        integer(i32), intent(in), dimension(:) :: pvt
        real(dp), intent(out), dimension(:,:) :: q, p
        real(dp), intent(out), target, dimension(:), optional :: work
        integer(i32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        real(dp), parameter :: zero = 0.0d0
        real(dp), parameter :: one = 1.0d0

        ! Local Variables
        integer(i32) :: j, jp, m, n, mn, flag
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg

        ! Initialization
        m = size(r, 1)
        n = size(r, 2)
        mn = min(m, n)
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        flag = 0
        if (size(tau) /= mn) then
            flag = 2
        else if (size(pvt) /= n) then
            flag = 3
        else if (size(q, 1) /= m .or. &
            (size(q, 2) /= m .and. size(q, 2) /= n)) then
            flag = 4
        else if (size(q, 2) == n .and. m < n) then
            flag = 4
        else if (size(p, 1) /= n .or. size(p, 2) /= n) then
            flag = 5
        end if
        if (flag /= 0) then
            ! ERROR: One of the input arrays is not sized correctly
            write(errmsg, '(AI0A)') "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("form_qr_pivot", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Generate Q and R
        call form_qr_no_pivot(r, tau, q, work = work, olwork = olwork, &
            err = errmgr)
        if (present(olwork)) return ! Just a workspace query
        if (errmgr%has_error_occurred()) return

        ! Form P
        do j = 1, n
            jp = pvt(j)
            p(:,j) = zero
            p(jp,j) = one
        end do
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Multiplies a general matrix by the orthogonal matrix Q from a QR
    !! factorization such that: C = op(Q) * C, or C = C * op(Q).
    !!
    !! @param[in] lside Set to true to apply Q or Q**T from the left; else, set
    !!  to false to apply Q or Q**T from the right.
    !! @param[in] trans Set to true to apply Q**T; else, set to false.
    !! @param[in] a On input, an LDA-by-K matrix containing the elementary
    !!  reflectors output from the QR factorization.  If @p lside is set to
    !!  true, LDA = M, and M >= K >= 0; else, if @p lside is set to false,
    !!  LDA = N, and N >= K >= 0.  Notice, the contents of this matrix are
    !!  restored on exit.
    !! @param[in] tau A K-element array containing the scalar factors of each
    !!  elementary reflector defined in @p a.
    !! @param[in,out] c On input, the M-by-N matrix C.  On output, the product
    !!  of the orthogonal matrix Q and the original matrix C.
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
    !! This routine utilizes the LAPACK routine DORMQR.
    subroutine mult_qr_mtx(lside, trans, a, tau, c, work, olwork, err)
        ! Arguments
        logical, intent(in) :: lside, trans
        real(dp), intent(in), dimension(:) :: tau
        real(dp), intent(inout), dimension(:,:) :: a, c
        real(dp), intent(out), target, dimension(:), optional :: work
        integer(i32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        real(dp), parameter :: one = 1.0d0

        ! Local Variables
        character :: side, t
        integer(i32) :: m, n, k, nrowa, istat, flag, lwork
        real(dp), pointer, dimension(:) :: wptr
        real(dp), allocatable, target, dimension(:) :: wrk
        real(dp), dimension(1) :: temp
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg

        ! Initialization
        m = size(c, 1)
        n = size(c, 2)
        k = size(tau)
        if (lside) then
            side = 'L'
            nrowa = m
        else
            side = 'R'
            nrowa = n
        end if
        if (trans) then
            t = 'T'
        else
            t = 'N'
        end if
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        flag = 0
        if (lside) then
            ! A is M-by-K, M >= K >= 0
            if (size(a, 1) /= m .or. size(a, 2) < k) flag = 3
        else
            ! A is N-by-K, N >= K >= 0
            if (size(a, 1) /= n .or. size(a, 2) < k) flag = 3
        end if
        if (flag /= 0) then
            ! ERROR: One of the input arrays is not sized correctly
            write(errmsg, '(AI0A)') "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("mult_qr_mtx", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Workspace Query
        call DORMQR(side, t, m, n, k, a, nrowa, tau, c, m, temp, -1, flag)
        lwork = int(temp(1), i32)
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (size(work) < lwork) then
                ! ERROR: WORK not sized correctly
                call errmgr%report_error("mult_qr_mtx", &
                    "Incorrectly sized input array WORK, argument 6.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lwork)
        else
            allocate(wrk(lwork), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("mult_qr_mtx", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if

        ! Call DORMQR
        call DORMQR(side, t, m, n, k, a, nrowa, tau, c, m, wptr, lwork, flag)
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Multiplies a vector by the orthogonal matrix Q from a QR
    !! factorization such that: C = op(Q) * C.
    !!
    !! @param[in] trans Set to true to apply Q**T; else, set to false.
    !! @param[in] a On input, an M-by-K matrix containing the elementary
    !!  reflectors output from the QR factorization.  Notice, the contents of
    !!  this matrix are restored on exit.
    !! @param[in] tau A K-element array containing the scalar factors of each
    !!  elementary reflector defined in @p a.
    !! @param[in,out] c On input, the M-element vector C.  On output, the
    !!  product of the orthogonal matrix Q and the original vector C.
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
    !! This routine is based upon the LAPACK routine DORM2R.
    subroutine mult_qr_vec(trans, a, tau, c, work, olwork, err)
        ! Arguments
        logical, intent(in) :: trans
        real(dp), intent(inout), dimension(:,:) :: a
        real(dp), intent(in), dimension(:) :: tau
        real(dp), intent(inout), dimension(:) :: c
        real(dp), intent(out), target, dimension(:), optional :: work
        integer(i32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        real(dp), parameter :: one = 1.0d0

        ! Local Variables
        character :: side, t
        integer(i32) :: m, k, nrowa, istat, flag, lwork
        real(dp), pointer, dimension(:) :: wptr
        real(dp), allocatable, target, dimension(:) :: wrk
        real(dp), dimension(1) :: temp
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg

        ! Initialization
        m = size(c)
        k = size(tau)
        side = 'L'
        nrowa = m
        if (trans) then
            t = 'T'
        else
            t = 'N'
        end if
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        flag = 0
        if (size(a, 1) /= m .or. size(a, 2) < k) flag = 3
        if (flag /= 0) then
            ! ERROR: One of the input arrays is not sized correctly
            write(errmsg, '(AI0A)') "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("mult_qr_vec", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Workspace Query
        call DORMQR(side, t, m, 1, k, a, nrowa, tau, c, m, temp, -1, flag)
        lwork = int(temp(1), i32)
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (size(work) < lwork) then
                ! ERROR: WORK not sized correctly
                call errmgr%report_error("mult_qr_vec", &
                    "Incorrectly sized input array WORK, argument 6.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lwork)
        else
            allocate(wrk(lwork), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("mult_qr_vec", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if

        ! Call DORMQR
        call DORMQR(side, t, m, 1, k, a, nrowa, tau, c, m, wptr, lwork, flag)
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Computes the rank 1 update to an M-by-N QR factored matrix A
    !! (M >= N) where A = Q * R, and A1 = A + U * V**T such that A1 = Q1 * R1.
    !!
    !! @param[in,out] q On input, the original M-by-K orthogonal matrix Q.  On
    !!  output, the updated matrix Q1.
    !! @param[in,out] r On input, the M-by-N matrix R.  On output, the updated
    !!  matrix R1.
    !! @param[in,out] u On input, the M-element U update vector.  On output,
    !!  the original content of the array is overwritten.
    !! @param[in,out] v On input, the N-element V update vector.  On output,
    !!  the original content of the array is overwritten.
    !! @param[out] work An optional argument that if supplied prevents local
    !!  memory allocation.  If provided, the array must have at least 2*K
    !!  elements.
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
    !! @par Remarks
    !! @verbatim
    !! Notice, K must either be equal to M, or to N.  In the event that K = N,
    !! only the submatrix Qa is updated.  This is appropriate as the QR
    !! factorization for an overdetermined system can be written as follows:
    !!  A = Q * R = [Qa, Qb] * [Ra]
    !!                         [0 ]
    !!
    !! Note: Ra is upper triangular of dimension N-by-N.
    !! @endverbatim
    !!
    !! @par Usage
    !! The following example illustrates a rank 1 update to a QR factored 
    !! system.  The results are compared to updating the original matrix, and
    !! then performing the factorization.
    !! @code{.f90}
    !! program example
    !!     use iso_fortran_env
    !!     use linalg_factor
    !!     use linalg_solve
    !!     use linalg_core
    !!     implicit none
    !!
    !!     ! Variables
    !!     real(real64) :: a(3,3), u(3), v(3), r(3,3), tau(3), q(3,3), qu(3,3)
    !!     integer(int32) :: i
    !!
    !!     ! Build the 3-by-3 matrix A.
    !!     !     | 1   2   3 |
    !!     ! A = | 4   5   6 |
    !!     !     | 7   8   0 |
    !!     a = reshape( &
    !!         [1.0d0, 4.0d0, 7.0d0, 2.0d0, 5.0d0, 8.0d0, 3.0d0, 6.0d0, 0.0d0], &
    !!         [3, 3])
    !!
    !!     ! Build the update vectors
    !!     !     | 1/2 |      | 1 |
    !!     ! u = | 3/2 |, v = | 5 |
    !!     !     |  3  |      | 2 |
    !!     u = [0.5d0, 1.5d0, 3.0d0]
    !!     v = [1.0d0, 5.0d0, 2.0d0]
    !!
    !!     ! Compute the QR factorization of the original matrix
    !!     r = a   ! Making a copy as the matrix will be overwritten by qr_factor
    !!     call qr_factor(r, tau)
    !!
    !!     ! Form Q & R
    !!     call form_qr(r, tau, q)
    !!
    !!     ! Compute the rank 1 update to the original matrix such that: 
    !!     ! A = A + u * v**T
    !!     call rank1_update(1.0d0, u, v, a)
    !!
    !!     ! Compute the rank 1 update to the factorization.  Notice, the contents 
    !!     ! of U & V are destroyed as part of this process.
    !!     call qr_rank1_update(q, r, u, v)
    !!
    !!     ! As comparison, compute the QR factorization of the rank 1 updated matrix
    !!     call qr_factor(a, tau)
    !!     call form_qr(a, tau, qu)
    !!
    !!     ! Display the matrices
    !!     print '(A)', "Updating the Factored Form:"
    !!     print '(A)', "Q = "
    !!     do i = 1, size(q, 1)
    !!         print *, q(i,:)
    !!     end do
    !!     print '(A)', "R = "
    !!     do i = 1, size(r, 1)
    !!         print *, r(i,:)
    !!     end do
    !!
    !!     print '(A)', "Updating A Directly:"
    !!     print '(A)', "Q = "
    !!     do i = 1, size(qu, 1)
    !!         print *, qu(i,:)
    !!     end do
    !!     print '(A)', "R = "
    !!     do i = 1, size(a, 1)
    !!         print *, a(i,:)
    !!     end do
    !! end program 
    !! @endcode
    !! The above program produces the following output.
    !! @code{.txt}
    !! Updating the Factored Form:
    !! Q =
    !!  -0.13031167282892092       0.98380249683206911      -0.12309149097933236
    !!  -0.47780946703937632      -0.17109608640557677      -0.86164043685532932
    !!  -0.86874448552613881       -5.3467527001743037E-002  0.49236596391733078
    !! R =
    !!  -11.510864433221338       -26.540144032823541       -10.033998807826904
    !!  0.0000000000000000        1.0586570346345126        2.0745400476676279
    !!  0.0000000000000000        0.0000000000000000       -5.2929341121113067
    !! Updating A Directly:
    !! Q =
    !!  -0.13031167282892087       0.98380249683206955      -0.12309149097933178
    !!  -0.47780946703937643      -0.17109608640557616      -0.86164043685532943
    !!  -0.86874448552613903       -5.3467527001742954E-002  0.49236596391733084
    !! R =
    !!  -11.510864433221336       -26.540144032823545       -10.033998807826906
    !!  0.0000000000000000        1.0586570346345205        2.0745400476676350
    !!  0.0000000000000000        0.0000000000000000       -5.2929341121113058
    !! @endcode
    !!
    !! @par Notes
    !! This routine utilizes the QRUPDATE routine DQR1UP.
    !!
    !! @par See Also
    !! [Source](https://sourceforge.net/projects/qrupdate/)
    subroutine qr_rank1_update(q, r, u, v, work, err)
        ! Arguments
        real(dp), intent(inout), dimension(:,:) :: q, r
        real(dp), intent(inout), dimension(:) :: u, v
        real(dp), intent(out), target, optional, dimension(:) :: work
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        logical :: full
        integer(i32) :: m, n, k, lwork, istat, flag
        real(dp), pointer, dimension(:) :: wptr
        real(dp), allocatable, target, dimension(:) :: wrk
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg

        ! Initialization
        m = size(u, 1)
        n = size(r, 2)
        k = min(m, n)
        full = size(q, 2) == m
        lwork = 2 * k
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        flag = 0
        if (m < n) then
            flag = 1
        else if (.not.full .and. size(q, 2) /= k) then
            flag = 1
        else if (size(r, 1) /= m) then
            flag = 2
        else if (size(u) /= m) then
            flag = 3
        else if (size(v) /= n) then
            flag = 4
        end if
        if (flag /= 0) then
            ! ERROR: One of the input arrays is not sized correctly
            write(errmsg, '(AI0A)') "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("qr_rank1_update", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (size(work) < lwork) then
                ! ERROR: WORK not sized correctly
                call errmgr%report_error("qr_rank1_update", &
                    "Incorrectly sized input array WORK, argument 5.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lwork)
        else
            allocate(wrk(lwork), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("qr_rank1_update", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if

        ! Process
        call DQR1UP(m, n, k, q, m, r, m, u, v, wptr)

        ! End
        if (allocated(wrk)) deallocate(wrk)
    end subroutine

! ******************************************************************************
! CHOLESKY FACTORIZATION
! ------------------------------------------------------------------------------
    !> @brief Computes the Cholesky factorization of a symmetric, positive
    !! definite matrix.
    !!
    !! @param[in,out] a On input, the N-by-N matrix to factor.  On output, the
    !!  factored matrix is returned in either the upper or lower triangular
    !!  portion of the matrix, dependent upon the value of @p upper.
    !! @param[in] upper An optional input that, if specified, provides control
    !!  over whether the factorization is computed as A = U**T * U (set to
    !!  true), or as A = L * L**T (set to false).  The default value is true
    !!  such that A = U**T * U.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - LA_ARRAY_SIZE_ERROR: Occurs if @p a is not square.
    !!  - LA_MATRIX_FORMAT_ERROR: Occurs if @p a is not positive definite.
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
    !!
    !! @par Notes
    !! This routine utilizes the LAPACK routine DPOTRF.
    subroutine cholesky_factor(a, upper, err)
        ! Arguments
        real(dp), intent(inout), dimension(:,:) :: a
        logical, intent(in), optional :: upper
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        real(dp), parameter :: zero = 0.0d0

        ! Local Variables
        character :: uplo
        integer(i32) :: i, n, flag
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg

        ! Initialization
        n = size(a, 1)
        if (present(upper)) then
            if (upper) then
                uplo = 'U'
            else
                uplo = 'L'
            end if
        else
            uplo = 'U'
        end if
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        if (size(a, 2) /= n) then
            ! ERROR: A must be square
            call errmgr%report_error("cholesky_factor", &
                "The input matrix must be square.", LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Process
        call DPOTRF(uplo, n, a, n, flag)
        if (flag > 0) then
            ! ERROR: Matrix is not positive definite
            write(errmsg, '(AI0A)') "The leading minor of order ", flag, &
                " is not positive definite."
            call errmgr%report_error("cholesky_factor", trim(errmsg), &
                LA_MATRIX_FORMAT_ERROR)
        end if

        ! Zero out the non-used upper or lower diagonal
        if (uplo == 'U') then
            ! Zero out the lower
            do i = 1, n - 1
                a(i+1:n,i) = zero
            end do
        else
            ! Zero out the upper
            do i = 2, n
                a(1:i-1,i) = zero
            end do
        end if
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Computes the rank 1 update to a Cholesky factored matrix (upper
    !! triangular).
    !!
    !! @param[in,out] r On input, the N-by-N upper triangular matrix R.  On
    !!  output, the updated matrix R1.
    !! @param[in,out] u On input, the N-element update vector U.  On output,
    !!  the rotation sines used to transform R to R1.
    !! @param[out] work An optional argument that if supplied prevents local
    !!  memory allocation.  If provided, the array must have at least N
    !!  elements.  Additionally, this workspace array is used to contain the
    !!  rotation cosines used to transform R to R1.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - LA_ARRAY_SIZE_ERROR: Occurs if any of the input array sizes are 
    !!      incorrect.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
    !!      there is insufficient memory available.
    !!
    !! @par Usage
    !! The following example illustrates the use of the rank 1 Cholesky update,
    !! and compares the results to factoring the original rank 1 updated matrix.
    !! @code{.f90}
    !! program example
    !!     use iso_fortran_env, only : real64, int32
    !!     use linalg_factor, only : cholesky_factor, cholesky_rank1_update
    !!     use linalg_core, only : rank1_update
    !!     implicit none
    !!
    !!     ! Variables
    !!     real(real64) :: a(3,3), u(3), au(3,3)
    !!     integer(int32) :: i
    !!
    !!     ! Build the 3-by-3 positive-definite matrix A.
    !!     !     | 4   12   -16 |
    !!     ! A = | 12  37   -43 |
    !!     !     |-16 -43    98 |
    !!     a = reshape([4.0d0, 12.0d0, -16.0d0, 12.0d0, 37.0d0, -43.0d0, -16.0d0, &
    !!         -43.0d0, 98.0d0], [3, 3])
    !!
    !!     ! Build the update vector U
    !!     u = [0.5d0, -1.5d0, 2.0d0]
    !!
    !!     ! Compute the rank 1 update of A
    !!     au = a
    !!     call rank1_update(1.0d0, u, u, au)
    !!
    !!     ! Compute the Cholesky factorization of the original matrix
    !!     call cholesky_factor(a)
    !!
    !!     ! Apply the rank 1 update to the factored matrix
    !!     call cholesky_rank1_update(a, u)
    !!
    !!     ! Compute the Cholesky factorization of the update of the original matrix
    !!     call cholesky_factor(au)
    !!
    !!     ! Display the matrices
    !!     print '(A)', "Updating the Factored Form:"
    !!     do i = 1, size(a, 1)
    !!         print *, a(i,:)
    !!     end do
    !!
    !!     print '(A)', "Updating A Directly:"
    !!     do i = 1, size(au, 1)
    !!         print *, au(i,:)
    !!     end do
    !! end program
    !! @endcode
    !! The above program produces the following output.
    !! @code{.txt}
    !! Updating the Factored Form:
    !!  2.0615528128088303        5.4570515633174921       -7.2760687510899889
    !!  0.0000000000000000        3.0774320845949008       -2.0452498947307731
    !!  0.0000000000000000        0.0000000000000000        6.6989384530323566
    !! Updating A Directly:
    !!  2.0615528128088303        5.4570515633174921       -7.2760687510899889
    !!  0.0000000000000000        3.0774320845949008       -2.0452498947307736
    !!  0.0000000000000000        0.0000000000000000        6.6989384530323557
    !! @endcode
    !!
    !! @par Notes
    !! This routine utilizes the QRUPDATE routine DCH1UP.
    !!
    !! @par See Also
    !! [Source](https://sourceforge.net/projects/qrupdate/)
    subroutine cholesky_rank1_update(r, u, work, err)
        ! Arguments
        real(dp), intent(inout), dimension(:,:) :: r
        real(dp), intent(inout), dimension(:) :: u
        real(dp), intent(out), target, optional, dimension(:) :: work
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        integer(i32) :: n, lwork, istat, flag
        real(dp), pointer, dimension(:) :: wptr
        real(dp), allocatable, target, dimension(:) :: wrk
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg

        ! Initialization
        n = size(r, 1)
        lwork = n
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        flag = 0
        if (size(r, 2) /= n) then
            flag = 1
        else if (size(u) /= n) then
            flag = 2
        end if
        if (flag /= 0) then
            write(errmsg, '(AI0A)') "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("cholesky_rank1_update", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (size(work) < lwork) then
                ! ERROR: Workspace array is not sized correctly
                call errmgr%report_error("cholesky_rank1_update", &
                    "The workspace array is too short.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lwork)
        else
            allocate(wrk(lwork), stat = istat)
            if (istat /= 0) then
                call errmgr%report_error("cholesky_rank1_update", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if

        ! Process
        call DCH1UP(n, r, n, u, wptr)
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Computes the rank 1 downdate to a Cholesky factored matrix (upper
    !! triangular).
    !!
    !! @param[in,out] r On input, the N-by-N upper triangular matrix R.  On
    !!  output, the updated matrix R1.
    !! @param[in,out] u On input, the N-element update vector U.  On output,
    !!  the rotation sines used to transform R to R1.
    !! @param[out] work An optional argument that if supplied prevents local
    !!  memory allocation.  If provided, the array must have at least N
    !!  elements.  Additionally, this workspace array is used to contain the
    !!  rotation cosines used to transform R to R1.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - LA_ARRAY_SIZE_ERROR: Occurs if any of the input array sizes are 
    !!      incorrect.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
    !!      there is insufficient memory available.
    !!  - LA_MATRIX_FORMAT_ERROR: Occurs if the downdated matrix is not 
    !!      positive definite.
    !!  - LA_SINGULAR_MATRIX_ERROR: Occurs if @p r is singular.
    !!
    !! @par Notes
    !! This routine utilizes the QRUPDATE routine DCH1DN.
    !!
    !! @par See Also
    !! [Source](https://sourceforge.net/projects/qrupdate/)
    subroutine cholesky_rank1_downdate(r, u, work, err)
        ! Arguments
        real(dp), intent(inout), dimension(:,:) :: r
        real(dp), intent(inout), dimension(:) :: u
        real(dp), intent(out), target, optional, dimension(:) :: work
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        integer(i32) :: n, lwork, istat, flag
        real(dp), pointer, dimension(:) :: wptr
        real(dp), allocatable, target, dimension(:) :: wrk
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg

        ! Initialization
        n = size(r, 1)
        lwork = n
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        flag = 0
        if (size(r, 2) /= n) then
            flag = 1
        else if (size(u) /= n) then
            flag = 2
        end if
        if (flag /= 0) then
            write(errmsg, '(AI0A)') "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("cholesky_rank1_downdate", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (size(work) < lwork) then
                ! ERROR: Workspace array is not sized correctly
                call errmgr%report_error("cholesky_rank1_downdate", &
                    "The workspace array is too short.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lwork)
        else
            allocate(wrk(lwork), stat = istat)
            if (istat /= 0) then
                call errmgr%report_error("cholesky_rank1_downdate", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if

        ! Process
        call DCH1DN(n, r, n, u, wptr, flag)
        if (flag == 1) then
            ! ERROR: The matrix is not positive definite
            call errmgr%report_error("cholesky_rank1_downdate", &
                "The downdated matrix is not positive definite.", &
                LA_MATRIX_FORMAT_ERROR)
        else if (flag == 2) then
            ! ERROR: The matrix is singular
            call errmgr%report_error("cholesky_rank1_downdate", &
                "The input matrix is singular.", LA_SINGULAR_MATRIX_ERROR)
        end if
    end subroutine

! ******************************************************************************
! RZ FACTORIZATION ROUTINES
! ------------------------------------------------------------------------------
    !> @brief Factors an upper trapezoidal matrix by means of orthogonal
    !! transformations such that A = R * Z = (R 0) * Z.  Z is an orthogonal
    !! matrix of dimension N-by-N, and R is an M-by-M upper triangular
    !! matrix.
    !!
    !! @param[in,out] a On input, the M-by-N upper trapezoidal matrix to factor.
    !!  On output, the leading M-by-M upper triangular part of the matrix
    !!  contains the upper triangular matrix R, and elements N-L+1 to N of the
    !!  first M rows of A, with the array @p tau, represent the orthogonal
    !!  matrix Z as a product of M elementary reflectors.
    !! @param[out] tau An M-element array used to store the scalar
    !!  factors of the elementary reflectors.
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
    !! @par Further Details
    !! @verbatim
    !!  The factorization is obtained by Householder's method.  The kth
    !!  transformation matrix, Z( k ), which is used to introduce zeros into
    !!  the ( m - k + 1 )th row of A, is given in the form
    !!
    !!     Z( k ) = ( I     0   ),
    !!              ( 0  T( k ) )
    !!
    !!  where
    !!
    !!     T( k ) = I - tau*u( k )*u( k )**T,   u( k ) = (   1    ),
    !!                                                   (   0    )
    !!                                                   ( z( k ) )
    !!
    !!  tau is a scalar and z( k ) is an l element vector. tau and z( k )
    !!  are chosen to annihilate the elements of the kth row of A2.
    !!
    !!  The scalar tau is returned in the kth element of TAU and the vector
    !!  u( k ) in the kth row of A2, such that the elements of z( k ) are
    !!  in  a( k, l + 1 ), ..., a( k, n ). The elements of R are returned in
    !!  the upper triangular part of A1.
    !!
    !!  Z is given by
    !!
    !!     Z =  Z( 1 ) * Z( 2 ) * ... * Z( m ).
    !! @endverbatim
    !!
    !! @par Notes
    !! This routine is based upon the LAPACK routine DTZRZF.
    !!
    !! @par See Also
    !! - [LAPACK Users Manual](http://netlib.org/lapack/lug/node44.html)
    subroutine rz_factor(a, tau, work, olwork, err)
        ! Arguments
        real(dp), intent(inout), dimension(:,:) :: a
        real(dp), intent(out), dimension(:) :: tau
        real(dp), intent(out), target, optional, dimension(:) :: work
        integer(i32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        integer(i32) :: m, n, lwork, flag, istat
        real(dp), pointer, dimension(:) :: wptr
        real(dp), allocatable, target, dimension(:) :: wrk
        real(dp), dimension(1) :: temp
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg

        ! Initialization
        m = size(a, 1)
        n = size(a, 2)
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        flag = 0
        if (size(tau) /= m) then
            flag = 3
        end if
        if (flag /= 0) then
            ! ERROR: One of the input arrays is not sized correctly
            write(errmsg, '(AI0A)') "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("rz_factor", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Workspace Query
        call DTZRZF(m, n, a, m, tau, temp, -1, flag)
        lwork = int(temp(1), i32)
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (size(work) < lwork) then
                ! ERROR: WORK not sized correctly
                call errmgr%report_error("rz_factor", &
                    "Incorrectly sized input array WORK, argument 3.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lwork)
        else
            allocate(wrk(lwork), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("rz_factor", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if

        ! Call DTZRZF
        call DTZRZF(m, n, a, m, tau, wptr, lwork, flag)
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Multiplies a general matrix by the orthogonal matrix Z from an
    !! RZ factorization such that: C = op(Z) * C, or C = C * op(Z).
    !!
    !! @param[in] lside Set to true to apply Z or Z**T from the left; else, set
    !!  to false to apply Z or Z**T from the right.
    !! @param[in] trans Set to true to apply Z**T; else, set to false.
    !! @param[in] l The number of columns in matrix @p a containing the
    !!  meaningful part of the Householder vectors.  If @p lside is true,
    !!  M >= L >= 0; else, if @p lside is false, N >= L >= 0.
    !! @param[in,out] a On input the K-by-LTA matrix Z, where LTA = M if
    !!  @p lside is true; else, LTA = N if @p lside is false.  The I-th row must
    !!  contain the Householder vector in the last k rows. Notice, the contents
    !!  of this matrix are restored on exit.
    !! @param[in] tau A K-element array containing the scalar factors of the
    !!  elementary reflectors, where M >= K >= 0 if @p lside is true; else,
    !!  N >= K >= 0 if @p lside is false.
    !! @param[in,out] c On input, the M-by-N matrix C.  On output, the product
    !!  of the orthogonal matrix Z and the original matrix C.
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
    !! This routine utilizes the LAPACK routine DORMRZ.
    subroutine mult_rz_mtx(lside, trans, l, a, tau, c, work, olwork, err)
        ! Arguments
        logical, intent(in) :: lside, trans
        integer(i32), intent(in) :: l
        real(dp), intent(inout), dimension(:,:) :: a, c
        real(dp), intent(in), dimension(:) :: tau
        real(dp), intent(out), target, optional, dimension(:) :: work
        integer(i32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        character :: side, t
        integer(i32) :: m, n, k, lwork, flag, istat, lda
        real(dp), pointer, dimension(:) :: wptr
        real(dp), allocatable, target, dimension(:) :: wrk
        real(dp), dimension(1) :: temp
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg

        ! Initialization
        m = size(c, 1)
        n = size(c, 2)
        k = size(tau)
        lda = size(a, 1)
        if (lside) then
            side = 'L'
        else
            side = 'R'
        end if
        if (trans) then
            t = 'T'
        else
            t = 'N'
        end if
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        flag = 0
        if (lside) then
            if (l > m .or. l < 0) then
               flag = 3
            else if (k > m) then
                flag = 5
            else if (size(a, 1) < k .or. size(a, 2) /= m) then
                flag = 4
            end if
        else
            if (l > n .or. l < 0) then
                flag = 3
            else if (k > n) then
                flag = 5
            else if (size(a, 1) < k .or. size(a, 2) /= n) then
                flag = 4
            end if
        end if
        if (flag /= 0) then
            ! ERROR: One of the input arrays is not sized correctly
            write(errmsg, '(AI0A)') "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("mult_rz_mtx", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Workspace Query
        call DORMRZ(side, t, m, n, k, l, a, lda, tau, c, m, temp, -1, flag)
        lwork = int(temp(1), i32)
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (size(work) < lwork) then
                ! ERROR: WORK not sized correctly
                call errmgr%report_error("mult_rz_mtx", &
                    "Incorrectly sized input array WORK, argument 7.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lwork)
        else
            allocate(wrk(lwork), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("mult_rz_mtx", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if

        ! Call DORMRZ
        call DORMRZ(side, t, m, n, k, l, a, lda, tau, c, m, wptr, lwork, flag)

        ! End
        if (allocated(wrk)) deallocate(wrk)
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Multiplies a vector by the orthogonal matrix Z from an
    !! RZ factorization such that: C = op(Z) * C.
    !!
    !! @param[in] trans Set to true to apply Z**T; else, set to false.
    !! @param[in] l The number of columns in matrix @p a containing the
    !!  meaningful part of the Householder vectors.  If @p lside is true,
    !!  M >= L >= 0; else, if @p lside is false, N >= L >= 0.
    !! @param[in,out] a On input the K-by-LTA matrix Z, where LTA = M if
    !!  @p lside is true; else, LTA = N if @p lside is false.  The I-th row must
    !!  contain the Householder vector in the last k rows. Notice, the contents
    !!  of this matrix are restored on exit.
    !! @param[in] tau A K-element array containing the scalar factors of the
    !!  elementary reflectors, where M >= K >= 0 if @p lside is true; else,
    !!  N >= K >= 0 if @p lside is false.
    !! @param[in,out] c On input, the M-element array C.  On output, the product
    !!  of the orthogonal matrix Z and the original array C.
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
    !! This routine utilizes the LAPACK routine DORMRZ.
    subroutine mult_rz_vec(trans, l, a, tau, c, work, olwork, err)
        ! Arguments
        logical, intent(in) :: trans
        integer(i32), intent(in) :: l
        real(dp), intent(inout), dimension(:,:) :: a
        real(dp), intent(in), dimension(:) :: tau
        real(dp), intent(inout), dimension(:) :: c
        real(dp), intent(out), target, optional, dimension(:) :: work
        integer(i32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        character :: side, t
        integer(i32) :: m, k, lwork, flag, istat, lda
        real(dp), pointer, dimension(:) :: wptr
        real(dp), allocatable, target, dimension(:) :: wrk
        real(dp), dimension(1) :: temp
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg

        ! Initialization
        m = size(c)
        k = size(tau)
        lda = size(a, 1)
        side = 'L'
        if (trans) then
            t = 'T'
        else
            t = 'N'
        end if
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        flag = 0
        if (l > m .or. l < 0) then
            flag = 2
        else if (k > m) then
            flag = 4
        else if (size(a, 1) < k .or. size(a, 2) /= m) then
            flag = 3
        end if
        if (flag /= 0) then
            ! ERROR: One of the input arrays is not sized correctly
            write(errmsg, '(AI0A)') "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("mult_rz_vec", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Workspace Query
        call DORMRZ(side, t, m, 1, k, l, a, lda, tau, c, m, temp, -1, flag)
        lwork = int(temp(1), i32)
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (size(work) < lwork) then
                ! ERROR: WORK not sized correctly
                call errmgr%report_error("mult_rz_vec", &
                    "Incorrectly sized input array WORK, argument 6.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lwork)
        else
            allocate(wrk(lwork), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("mult_rz_vec", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if

        ! Call DORMRZ
        call DORMRZ(side, t, m, 1, k, l, a, lda, tau, c, m, wptr, lwork, flag)
    end subroutine

! ******************************************************************************
! SVD ROUTINES
! ------------------------------------------------------------------------------
    !> @brief Computes the singular value decomposition of a matrix A.  The
    !!  SVD is defined as: A = U * S * V**T, where U is an M-by-M orthogonal
    !!  matrix, S is an M-by-N diagonal matrix, and V is an N-by-N orthogonal
    !!  matrix.
    !!
    !! @param[in,out] a On input, the M-by-N matrix to factor.  The matrix is
    !!  overwritten on output.
    !! @param[out] s A MIN(M, N)-element array containing the singular values
    !!  of @p a sorted in descending order.
    !! @param[out] u An optional argument, that if supplied, is used to contain
    !!  the orthogonal matrix U from the decomposition.  The matrix U contains
    !!  the left singular vectors, and can be either M-by-M (all left singular
    !!  vectors are computed), or M-by-MIN(M,N) (only the first MIN(M, N) left
    !!  singular vectors are computed).
    !! @param[out] vt An optional argument, that if supplied, is used to contain
    !!  the transpose of the N-by-N orthogonal matrix V.  The matrix V contains
    !!  the right singular vectors.
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
    !! ! Decompose matrix the M-by-N matrix A such that A = U * S * V**T with
    !! ! M >= N.
    !!
    !! ! Variables
    !! real(dp), dimension(m, n) :: a
    !! real(dp), dimension(m, m) :: u
    !! real(dp), dimension(n, n) :: vt
    !! real(dp), dimension(n) :: s
    !!
    !! ! Initialize A...
    !!
    !! ! Compute the SVD of A. On output, S contains the MIN(M,N) singular
    !! ! values of A in descending order, U contains the left singular vectors
    !! ! (one per column), and VT contains the right singular vectors (one per
    !! ! row).
    !! call svd(a, s, u, vt)
    !!
    !! ! Note: If M > N, then we can make U M-by-N, and compute the N
    !! ! left singular vectors of A, as there are at most N singular values
    !! ! of A.  Also, if M < N, then there are at most M singular values of A,
    !! ! and as such, the length of the array S should be M.
    !! @endcode
    !!
    !! @par Notes
    !! This routine utilizes the LAPACK routine DGESVD.
    !!
    !! @par See Also
    !! - [Wikipedia](https://en.wikipedia.org/wiki/Singular_value_decomposition)
    !! - [Wolfram MathWorld](http://mathworld.wolfram.com/SingularValueDecomposition.html)
    subroutine svd(a, s, u, vt, work, olwork, err)
        ! Arguments
        real(dp), intent(inout), dimension(:,:) :: a
        real(dp), intent(out), dimension(:) :: s
        real(dp), intent(out), optional, dimension(:,:) :: u, vt
        real(dp), intent(out), target, optional, dimension(:) :: work
        integer(i32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        character :: jobu, jobvt
        integer(i32) :: m, n, mn, istat, lwork, flag
        real(dp), pointer, dimension(:) :: wptr
        real(dp), allocatable, target, dimension(:) :: wrk
        real(dp), dimension(1) :: temp
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg

        ! Initialization
        m = size(a, 1)
        n = size(a, 2)
        mn = min(m, n)
        if (present(u)) then
            if (size(u, 2) == m) then
                jobu = 'A'
            else if (size(u, 2) == mn) then
                jobu = 'S'
            end if
        else
            jobu = 'N'
        end if
        if (present(vt)) then
            jobvt = 'A'
        else
            jobvt = 'N'
        end if
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        flag = 0
        if (size(s) /= mn) then
            flag = 2
        else if (present(u)) then
            if (size(u, 1) /= m) flag = 3
            if (size(u, 2) /= m .and. size(u, 2) /= mn) flag = 3
        else if (present(vt)) then
            if (size(vt, 1) /= n .or. size(vt, 2) /= n) flag = 4
        end if
        if (flag /= 0) then
            ! ERROR: One of the input arrays is not sized correctly
            write(errmsg, '(AI0A)') "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("svd", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Workspace Query
        call DGESVD(jobu, jobvt, m, n, a, m, s, temp, m, temp, n, temp, -1, &
            flag)
        lwork = int(temp(1), i32)
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (size(work) < lwork) then
                ! ERROR: WORK not sized correctly
                call errmgr%report_error("svd", &
                    "Incorrectly sized input array WORK, argument 5.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lwork)
        else
            allocate(wrk(lwork), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("svd", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if

        ! Call DGESVD
        if (present(u) .and. present(vt)) then
            call DGESVD(jobu, jobvt, m, n, a, m, s, u, m, vt, n, wptr, lwork, &
                flag)
        else if (present(u) .and. .not.present(vt)) then
            call DGESVD(jobu, jobvt, m, n, a, m, s, u, m, temp, n, wptr, &
                lwork, flag)
        else if (.not.present(u) .and. present(vt)) then
            call DGESVD(jobu, jobvt, m, n, a, m, s, temp, m, vt, n, wptr, &
                lwork, flag)
        else
            call DGESVD(jobu, jobvt, m, n, a, m, s, temp, m, temp, n, wptr, &
                lwork, flag)
        end if

        ! Check for convergence
        if (flag > 0) then
            write(errmsg, '(I0A)') flag, " superdiagonals could not " // &
                "converge to zero as part of the QR iteration process."
            call errmgr%report_warning("svd", errmsg, LA_CONVERGENCE_ERROR)
        end if
    end subroutine


end module
