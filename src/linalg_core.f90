! linalg_core.f90


!> @mainpage
!!
!! @section intro_sec Introduction
!! LINALG is a linear algebra library that provides a user-friendly interface
!! to several BLAS and LAPACK routines.
!!
!! @author Jason Christopherson
!! @version 1.5.0


!> @brief \b linalg_core
!!
!! @par Purpose
!! Provides common "core" linear algebra routines.
module linalg_core
    use, intrinsic :: iso_fortran_env, only : int32, real64
    use ferror, only : errors
    use linalg_constants
    implicit none

    private
    public :: mtx_mult
    public :: rank1_update
    public :: diag_mtx_mult
    public :: trace
    public :: mtx_rank
    public :: det
    public :: swap
    public :: recip_mult_array
    public :: tri_mtx_mult
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
    public :: solve_triangular_system
    public :: solve_lu
    public :: solve_qr
    public :: solve_cholesky
    public :: mtx_inverse
    public :: mtx_pinverse
    public :: solve_least_squares
    public :: solve_least_squares_full
    public :: solve_least_squares_svd
    public :: eigen
    public :: sort

! ******************************************************************************
! INTERFACES
! ------------------------------------------------------------------------------
!> @brief Performs the matrix operation:
!!  C = alpha * op(A) * op(B) + beta * C.
interface mtx_mult
    module procedure :: mtx_mult_mtx
    module procedure :: mtx_mult_vec
end interface

! ------------------------------------------------------------------------------
!> @brief Performs the rank-1 update to matrix A such that:
!! A = alpha * X * Y**T + A, where A is an M-by-N matrix, alpha is a scalar,
!! X is an M-element array, and N is an N-element array.
interface rank1_update
    module procedure :: rank1_update_dbl
end interface

! ------------------------------------------------------------------------------
!> @brief Multiplies a diagonal matrix with another matrix or array.
!!
!! @par Usage
!! The following example illustrates the use of the diagonal matrix
!! multiplication routine to compute the S * V**T component of a singular
!! value decomposition.
!! @code{.f90}
!! program example
!!     use iso_fortran_env, only : int32, real64
!!     use linalg_core
!!     implicit none
!!
!!     ! Variables
!!     real(real64) :: a(3,2), s(2), u(3,3), vt(2,2), ac(3,2)
!!     integer(int32) :: i
!!
!!     ! Initialize the 3-by-2 matrix A
!!     !     | 2   1 |
!!     ! A = |-3   1 |
!!     !     |-1   1 |
!!     a = reshape([2.0d0, -3.0d0, -1.0d0, 1.0d0, 1.0d0, 1.0d0], [3, 2])
!!
!!     ! Compute the singular value decomposition of A.  Notice, V**T is returned
!!     ! instead of V.  Also note, A is overwritten.
!!     call svd(a, s, u, vt)
!!
!!     ! Display the results
!!     print '(A)', "U ="
!!     do i = 1, size(u, 1)
!!         print *, u(i,:)
!!     end do
!!
!!     print '(A)', "S ="
!!     print '(F9.5)', (s(i), i = 1, size(a, 2))
!!
!!     print '(A)', "V**T ="
!!     do i = 1, size(vt, 1)
!!         print *, vt(i,:)
!!     end do
!!
!!     ! Compute U * S * V**T, but first establish S in full form
!!     call diag_mtx_mult(.true., 1.0d0, s, vt) ! Compute: VT = S * V**T
!!     ac = matmul(u(:,1:2), vt)
!!     print '(A)', "U * S * V**T ="
!!     do i = 1, size(ac, 1)
!!         print *, ac(i,:)
!!     end do
!! end program
!! @endcode
!! The above program produces the following output.
!! @code{.txt}
!! U =
!!  -0.47411577501825380      -0.81850539032073777      -0.32444284226152509
!!   0.82566838523833064      -0.28535874325972488      -0.48666426339228758
!!   0.30575472113569685      -0.49861740208412991       0.81110710565381272
!! S =
!!   3.78845
!!   1.62716
!! V**T =
!!  -0.98483334211643059       0.17350299206578967
!!  -0.17350299206578967      -0.98483334211643059
!! U * S * V**T =
!!    1.9999999999999993       0.99999999999999956
!!   -3.0000000000000000        1.0000000000000000
!!   -1.0000000000000000       0.99999999999999967
!! @endcode
interface diag_mtx_mult
    module procedure :: diag_mtx_mult_mtx
    module procedure :: diag_mtx_mult_mtx2
    module procedure :: diag_mtx_mult_mtx3
    module procedure :: diag_mtx_mult_mtx4
    module procedure :: diag_mtx_mult_mtx_cmplx
    module procedure :: diag_mtx_mult_mtx2_cmplx
end interface

! ------------------------------------------------------------------------------
!> @brief Computes the trace of a matrix (the sum of the main diagonal
!! elements).
interface trace
    module procedure :: trace_dbl
end interface

! ------------------------------------------------------------------------------
!> @brief Computes the rank of a matrix.
interface mtx_rank
    module procedure :: mtx_rank_dbl
end interface

! ------------------------------------------------------------------------------
!> @brief Computes the determinant of a square matrix.
interface det
    module procedure :: det_dbl
end interface

! ------------------------------------------------------------------------------
!> @brief Swaps the contents of two arrays.
interface swap
    module procedure :: swap_dbl
end interface

! ------------------------------------------------------------------------------
!> @brief Multiplies a vector by the reciprocal of a real scalar.
interface recip_mult_array
    module procedure :: recip_mult_array_dbl
end interface

! ------------------------------------------------------------------------------
!> @brief Computes the triangular matrix operation:
!! B = alpha * A**T * A + beta * B, or B = alpha * A * A**T + beta * B,
!! where A is a triangular matrix.
interface tri_mtx_mult
    module procedure :: tri_mtx_mult_dbl
    module procedure :: tri_mtx_mult_cmplx
end interface

! ------------------------------------------------------------------------------
!> @brief Computes the LU factorization of an M-by-N matrix.
!!
!! @par Usage
!! To solve a system of 3 equations of 3 unknowns using LU factorization,
!! the following code will suffice.
!! @code{.f90}
!! program example
!!     use iso_fortran_env
!!     use linalg_core
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
interface lu_factor
    module procedure :: lu_factor_dbl
    module procedure :: lu_factor_cmplx
end interface

!> @brief Extracts the L and U matrices from the condensed [L\\U] storage
!! format used by the @ref lu_factor.
!!
!! @par Usage
!! The following example illustrates how to extract the L, U, and P matrices
!! in order to solve a system of LU factored equations.
!! @code{.f90}
!! program example
!!     use iso_fortran_env, only : real64, int32
!!     use linalg_core
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
    module procedure :: form_lu_all_cmplx
    module procedure :: form_lu_only
    module procedure :: form_lu_only_cmplx
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
!!     use linalg_core
!!     implicit none
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
!!     use linalg_core
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
!!     use linalg_core
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
!> @brief Computes the rank 1 update to an M-by-N QR factored matrix A
!! (M >= N) where A = Q * R, and A1 = A + U * V**T such that A1 = Q1 * R1.
!!
!! @par Usage
!! The following example illustrates a rank 1 update to a QR factored
!! system.  The results are compared to updating the original matrix, and
!! then performing the factorization.
!! @code{.f90}
!! program example
!!     use iso_fortran_env
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
interface qr_rank1_update
    module procedure :: qr_rank1_update_dbl
end interface

! ------------------------------------------------------------------------------
!> @brief Computes the Cholesky factorization of a symmetric, positive
!! definite matrix.
!!
!! @par Usage
!! The following example illustrates the solution of a positive-definite
!! system of equations via Cholesky factorization.
!! @code{.f90}
!! program example
!!     use iso_fortran_env, only : real64, int32
!!     use linalg_core
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
interface cholesky_factor
    module procedure :: cholesky_factor_dbl
end interface

! ------------------------------------------------------------------------------
!> @brief Computes the rank 1 update to a Cholesky factored matrix (upper
!! triangular).
!!
!! @par Usage
!! The following example illustrates the use of the rank 1 Cholesky update,
!! and compares the results to factoring the original rank 1 updated matrix.
!! @code{.f90}
!! program example
!!     use iso_fortran_env, only : real64, int32
!!     use linalg_core
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
interface cholesky_rank1_update
    module procedure :: cholesky_rank1_update_dbl
end interface

! ------------------------------------------------------------------------------
!> @brief Computes the rank 1 downdate to a Cholesky factored matrix (upper
!! triangular).
!!
!! @par Usage
!! The following example illustrates the use of the rank 1 Cholesky
!! downdate, and compares the results to factoring the original rank 1
!! downdated matrix.
!! @code{.f90}
!! program example
!!     use iso_fortran_env, only : real64, int32
!!     use linalg_factor, only : cholesky_factor, cholesky_rank1_downdate
!!     use linalg_core, only : rank1_update
!!     implicit none
!!
!!     ! Variables
!!     real(real64) :: a(3,3), u(3), ad(3,3)
!!     integer(int32) :: i
!!
!!     ! Build the 3-by-3 matrix A.
!!     !     |  4.25   11.25   -15 |
!!     ! A = | 11.25   39.25   -46 |
!!     !     |  -15     -46    102 |
!!     a = reshape([4.25d0, 11.25d0, -15.0d0, 11.25d0, 39.25d0, -46.0d0, &
!!         -15.0d0, -46.0d0, 102.0d0], [3, 3])
!!
!!     ! The downdate vector
!!     !     |  0.5 |
!!     ! u = | -1.5 |
!!     !     |   2  |
!!     u = [0.5d0, -1.5d0, 2.0d0]
!!
!!     ! Compute the rank 1 downdate of A
!!     ad = a
!!     call rank1_update(-1.0d0, u, u, ad)
!!
!!     ! Compute the Cholesky factorization of the original matrix
!!     call cholesky_factor(a)
!!
!!     ! Apply the rank 1 downdate to the factored matrix
!!     call cholesky_rank1_downdate(a, u)
!!
!!     ! Compute the Cholesky factorization of the downdate to the original matrix
!!     call cholesky_factor(ad)
!!
!!     ! Display the matrices
!!     print '(A)', "Downdating the Factored Form:"
!!     do i = 1, size(a, 1)
!!         print *, a(i,:)
!!     end do
!!
!!     print '(A)', "Downdating A Directly:"
!!     do i = 1, size(ad, 1)
!!         print *, ad(i,:)
!!     end do
!! end program
!! @endcode
!! The above program produces the following output.
!! @code{.txt}
!! Downdating the Factored Form:
!!  2.0000000000000000        6.0000000000000000       -8.0000000000000000
!!  0.0000000000000000        1.0000000000000000        4.9999999999999973
!!  0.0000000000000000        0.0000000000000000        3.0000000000000049
!! Downdating A Directly:
!!  2.0000000000000000        6.0000000000000000       -8.0000000000000000
!!  0.0000000000000000        1.0000000000000000        5.0000000000000000
!!  0.0000000000000000        0.0000000000000000        3.0000000000000000
!! @endcode
interface cholesky_rank1_downdate
    module procedure :: cholesky_rank1_downdate_dbl
end interface

! ------------------------------------------------------------------------------
!> @brief Factors an upper trapezoidal matrix by means of orthogonal
!! transformations such that A = R * Z = (R 0) * Z.  Z is an orthogonal
!! matrix of dimension N-by-N, and R is an M-by-M upper triangular
!! matrix.
interface rz_factor
    module procedure :: rz_factor_dbl
end interface

! ------------------------------------------------------------------------------
!> @brief Multiplies a general matrix by the orthogonal matrix Z from an
!! RZ factorization.
interface mult_rz
    module procedure :: mult_rz_mtx
    module procedure :: mult_rz_vec
end interface

! ------------------------------------------------------------------------------
!> @brief Computes the singular value decomposition of a matrix A.  The
!!  SVD is defined as: A = U * S * V**T, where U is an M-by-M orthogonal
!!  matrix, S is an M-by-N diagonal matrix, and V is an N-by-N orthogonal
!!  matrix.
!!
!! @par Usage
!! The following example illustrates the calculation of the singular value
!! decomposition of an overdetermined system.
!! @code{.f90}
!! program example
!!     use iso_fortran_env, only : int32, real64
!!     use linalg_core
!!     implicit none
!!
!!     ! Variables
!!     real(real64) :: a(3,2), s(2), u(3,3), vt(2,2), ac(3,2)
!!     integer(int32) :: i
!!
!!     ! Initialize the 3-by-2 matrix A
!!     !     | 2   1 |
!!     ! A = |-3   1 |
!!     !     |-1   1 |
!!     a = reshape([2.0d0, -3.0d0, -1.0d0, 1.0d0, 1.0d0, 1.0d0], [3, 2])
!!
!!     ! Compute the singular value decomposition of A.  Notice, V**T is returned
!!     ! instead of V.  Also note, A is overwritten.
!!     call svd(a, s, u, vt)
!!
!!     ! Display the results
!!     print '(A)', "U ="
!!     do i = 1, size(u, 1)
!!         print *, u(i,:)
!!     end do
!!
!!     print '(A)', "S ="
!!     print '(F9.5)', (s(i), i = 1, size(a, 2))
!!
!!     print '(A)', "V**T ="
!!     do i = 1, size(vt, 1)
!!         print *, vt(i,:)
!!     end do
!!
!!     ! Compute U * S * V**T, but first establish S in full form
!!     call diag_mtx_mult(.true., 1.0d0, s, vt) ! Compute: VT = S * V**T
!!     ac = matmul(u(:,1:2), vt)
!!     print '(A)', "U * S * V**T ="
!!     do i = 1, size(ac, 1)
!!         print *, ac(i,:)
!!     end do
!! end program
!! @endcode
!! The above program produces the following output.
!! @code{.txt}
!! U =
!!  -0.47411577501825380      -0.81850539032073777      -0.32444284226152509
!!   0.82566838523833064      -0.28535874325972488      -0.48666426339228758
!!   0.30575472113569685      -0.49861740208412991       0.81110710565381272
!! S =
!!   3.78845
!!   1.62716
!! V**T =
!!  -0.98483334211643059       0.17350299206578967
!!  -0.17350299206578967      -0.98483334211643059
!! U * S * V**T =
!!    1.9999999999999993       0.99999999999999956
!!   -3.0000000000000000        1.0000000000000000
!!   -1.0000000000000000       0.99999999999999967
!! @endcode
interface svd
    module procedure :: svd_dbl
end interface

! ------------------------------------------------------------------------------
!> @brief Solves a triangular system of equations.
!!
!! @par Usage
!! The following example illustrates the solution of two triangular systems
!! to solve a system of LU factored equations.
!! @code{.f90}
!! program example
!!     use iso_fortran_env, only : real64, int32
!!     use linalg_core
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
    module procedure :: solve_tri_mtx_cmplx
    module procedure :: solve_tri_vec
    module procedure :: solve_tri_vec_cmplx
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
!!     use linalg_core
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
    module procedure :: solve_lu_mtx_cmplx
    module procedure :: solve_lu_vec
    module procedure :: solve_lu_vec_cmplx
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
!!     use linalg_core
!!     implicit none
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
!!     use linalg_core
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
!!
!! @par Usage
!! The following example illustrates the least squares solution of an
!! overdetermined system of linear equations.
!! @code{.f90}
!! program example
!!     use iso_fortran_env, only : real64, int32
!!     use linalg_core
!!     implicit none
!!
!!     ! Local Variables
!!     real(real64) :: a(3,2), b(3)
!!     integer(int32) :: i
!!
!!     ! Build the 3-by-2 matrix A
!!     !     | 2   1 |
!!     ! A = |-3   1 |
!!     !     |-1   1 |
!!     a = reshape([2.0d0, -3.0d0, -1.0d0, 1.0d0, 1.0d0, 1.0d0], [3, 2])
!!
!!     ! Build the right-hand-side vector B.
!!     !     |-1 |
!!     ! b = |-2 |
!!     !     | 1 |
!!     b = [-1.0d0, -2.0d0, 1.0d0]
!!
!!     ! The solution is:
!!     ! x = [0.13158, -0.57895]**T
!!
!!     ! Compute the solution via a least-squares approach.  The results overwrite
!!     ! the first 2 elements in b.
!!     call solve_least_squares(a, b)
!!
!!     ! Display the results
!!     print '(A)', "Least Squares Solution: X = "
!!     print '(F9.5)', (b(i), i = 1, size(a, 2))
!! end program
!! @endcode
!! The above program produces the following output.
!! @code{.txt}
!! Least Squares Solution: X =
!!  0.13158
!! -0.57895
!! @endcode
interface solve_least_squares
    module procedure :: solve_least_squares_mtx
    module procedure :: solve_least_squares_vec
end interface

! ------------------------------------------------------------------------------
!> @brief Solves the overdetermined or underdetermined system (A*X = B) of
!! M equations of N unknowns, but uses a full orthogonal factorization of
!! the system.
!!
!! @par Usage
!! The following example illustrates the least squares solution of an
!! overdetermined system of linear equations.
!! @code{.f90}
!! program example
!!     use iso_fortran_env, only : real64, int32
!!     use linalg_core
!!     implicit none
!!
!!     ! Local Variables
!!     real(real64) :: a(3,2), b(3)
!!     integer(int32) :: i
!!
!!     ! Build the 3-by-2 matrix A
!!     !     | 2   1 |
!!     ! A = |-3   1 |
!!     !     |-1   1 |
!!     a = reshape([2.0d0, -3.0d0, -1.0d0, 1.0d0, 1.0d0, 1.0d0], [3, 2])
!!
!!     ! Build the right-hand-side vector B.
!!     !     |-1 |
!!     ! b = |-2 |
!!     !     | 1 |
!!     b = [-1.0d0, -2.0d0, 1.0d0]
!!
!!     ! The solution is:
!!     ! x = [0.13158, -0.57895]**T
!!
!!     ! Compute the solution via a least-squares approach.  The results overwrite
!!     ! the first 2 elements in b.
!!     call solve_least_squares_full(a, b)
!!
!!     ! Display the results
!!     print '(A)', "Least Squares Solution: X = "
!!     print '(F9.5)', (b(i), i = 1, size(a, 2))
!! end program
!! @endcode
!! The above program produces the following output.
!! @code{.txt}
!! Least Squares Solution: X =
!!  0.13158
!! -0.57895
!! @endcode
interface solve_least_squares_full
    module procedure :: solve_least_squares_mtx_pvt
    module procedure :: solve_least_squares_vec_pvt
end interface

! ------------------------------------------------------------------------------
!> @brief Solves the overdetermined or underdetermined system (A*X = B) of
!! M equations of N unknowns using a singular value decomposition of
!! matrix A.
!!
!! @par Usage
!! The following example illustrates the least squares solution of an
!! overdetermined system of linear equations.
!! @code{.f90}
!! program example
!!     use iso_fortran_env, only : real64, int32
!!     use linalg_core
!!     implicit none
!!
!!     ! Local Variables
!!     real(real64) :: a(3,2), b(3)
!!     integer(int32) :: i
!!
!!     ! Build the 3-by-2 matrix A
!!     !     | 2   1 |
!!     ! A = |-3   1 |
!!     !     |-1   1 |
!!     a = reshape([2.0d0, -3.0d0, -1.0d0, 1.0d0, 1.0d0, 1.0d0], [3, 2])
!!
!!     ! Build the right-hand-side vector B.
!!     !     |-1 |
!!     ! b = |-2 |
!!     !     | 1 |
!!     b = [-1.0d0, -2.0d0, 1.0d0]
!!
!!     ! The solution is:
!!     ! x = [0.13158, -0.57895]**T
!!
!!     ! Compute the solution via a least-squares approach.  The results overwrite
!!     ! the first 2 elements in b.
!!     call solve_least_squares_svd(a, b)
!!
!!     ! Display the results
!!     print '(A)', "Least Squares Solution: X = "
!!     print '(F9.5)', (b(i), i = 1, size(a, 2))
!! end program
!! @endcode
!! The above program produces the following output.
!! @code{.txt}
!! Least Squares Solution: X =
!!  0.13158
!! -0.57895
!! @endcode
interface solve_least_squares_svd
    module procedure :: solve_least_squares_mtx_svd
    module procedure :: solve_least_squares_vec_svd
end interface

! ------------------------------------------------------------------------------
!> @brief Computes the inverse of a square matrix.
!!
!! @par Usage
!! The following example illustrates the inversion of a 3-by-3 matrix.
!! @code{.f90}
!! program example
!!     use iso_fortran_env, only : real64, int32
!!     use linalg_core
!!     implicit none
!!
!!     ! Variables
!!     real(real64) :: a(3,3), ai(3,3), c(3,3)
!!     integer(int32) :: i
!!
!!     ! Construct the 3-by-3 matrix A to invert
!!     !     | 1   2   3 |
!!     ! A = | 4   5   6 |
!!     !     | 7   8   0 |
!!     a = reshape([1.0d0, 4.0d0, 7.0d0, 2.0d0, 5.0d0, 8.0d0, 3.0d0, 6.0d0, &
!!         0.0d0], [3, 3])
!!
!!     ! Compute the inverse of A.  Notice, the original matrix is overwritten
!!     ! with it's inverse.
!!     ai = a
!!     call mtx_inverse(ai)
!!
!!     ! Show that A * inv(A) = I
!!     c = matmul(a, ai)
!!
!!     ! Display the inverse
!!     print '(A)', "Inverse:"
!!     do i = 1, size(ai, 1)
!!         print *, ai(i,:)
!!     end do
!!
!!     ! Display the result of A * inv(A)
!!     print '(A)', "A * A**-1:"
!!     do i = 1, size(c, 1)
!!         print *, c(i,:)
!!     end do
!! end program
!! @endcode
!! The above program produces the following output.
!! @code{.txt}
!! Inverse:
!!  -1.7777777777777777       0.88888888888888884      -0.11111111111111110
!!   1.5555555555555556      -0.77777777777777779       0.22222222222222221
!!  -0.11111111111111119      0.22222222222222227      -0.11111111111111112
!! A * A**-1:
!!   0.99999999999999989       5.5511151231257827E-017  -4.1633363423443370E-017
!!   5.5511151231257827E-017   1.0000000000000000       -8.3266726846886741E-017
!!   1.7763568394002505E-015  -8.8817841970012523E-016   1.0000000000000000
!! @endcode
interface mtx_inverse
    module procedure :: mtx_inverse_dbl
end interface

! ------------------------------------------------------------------------------
!> @brief Computes the Moore-Penrose pseudo-inverse of a M-by-N matrix
!! using the singular value decomposition of the matrix.
!!
!! @par Usage
!! The following example illustrates how to compute the Moore-Penrose
!! pseudo-inverse of a matrix.
!! @code{.f90}
!! program example
!!     use iso_fortran_env, only : int32, real64
!!     use linalg_core
!!     implicit none
!!
!!     ! Variables
!!     real(real64) :: a(3,2), ai(2,3), ao(3,2), c(2,2)
!!     integer(int32) :: i
!!
!!     ! Create the 3-by-2 matrix A
!!     !     | 1   0 |
!!     ! A = | 0   1 |
!!     !     | 0   1 |
!!     a = reshape([1.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0, 1.0d0], [3, 2])
!!     ao = a  ! Just making a copy for later as mtx_pinverse will destroy the
!!             ! contents of the original matrix
!!
!!     ! The Moore-Penrose pseudo-inverse of this matrix is:
!!     !         | 1   0    0  |
!!     ! A**-1 = |             |
!!     !         | 0  1/2  1/2 |
!!     call mtx_pinverse(a, ai)
!!
!!     ! Notice, A**-1 * A is an identity matrix.
!!     c = matmul(ai, ao)
!!
!!     ! Display the inverse
!!     print '(A)', "Inverse:"
!!     do i = 1, size(ai, 1)
!!         print *, ai(i,:)
!!     end do
!!
!!     ! Display the result of inv(A) * A
!!     print '(A)', "A**-1 * A:"
!!     do i = 1, size(c, 1)
!!         print *, c(i,:)
!!     end do
!! end program
!! @endcode
!! The above program produces the following output.
!! @code{.txt}
!! Inverse:
!!  1.0000000000000000        0.0000000000000000        0.0000000000000000
!!  0.0000000000000000       0.49999999999999978       0.49999999999999989
!! A**-1 * A:
!!  1.0000000000000000        0.0000000000000000
!!  0.0000000000000000       0.99999999999999967
!! @endcode
interface mtx_pinverse
    module procedure :: mtx_pinverse_dbl
end interface

! ------------------------------------------------------------------------------
!> @brief Computes the eigenvalues, and optionally the eigenvectors, of a
!! matrix.
!!
!! @par Usage
!! As an example, consider the eigenvalue problem arising from a mechanical
!! system of masses and springs such that the masses are described by
!! a mass matrix M, and the arrangement of springs are described by a
!! stiffness matrix K.
!! @code {.f90}
!! ! This is an example illustrating the use of the eigenvalue and eigenvector
!! ! routines to solve a free vibration problem of 3 masses connected by springs.
!! !
!! !     k1           k2           k3           k4
!! ! |-\/\/\-| m1 |-\/\/\-| m2 |-\/\/\-| m3 |-\/\/\-|
!! !
!! ! As illustrated above, the system consists of 3 masses connected by springs.
!! ! Spring k1 and spring k4 connect the end masses to ground.  The equations of
!! ! motion for this system are as follows.
!! !
!! ! | m1  0   0 | |x1"|   | k1+k2  -k2      0  | |x1|   |0|
!! ! | 0   m2  0 | |x2"| + |  -k2  k2+k3    -k3 | |x2| = |0|
!! ! | 0   0   m3| |x3"|   |   0    -k3    k3+k4| |x3|   |0|
!! !
!! ! Notice: x1" = the second time derivative of x1.
!! program example
!!     use iso_fortran_env, only : int32, real64
!!     use linalg_core
!!     implicit none
!!
!!     ! Define the model parameters
!!     real(real64), parameter :: pi = 3.14159265359d0
!!     real(real64), parameter :: m1 = 0.5d0
!!     real(real64), parameter :: m2 = 2.5d0
!!     real(real64), parameter :: m3 = 0.75d0
!!     real(real64), parameter :: k1 = 5.0d6
!!     real(real64), parameter :: k2 = 10.0d6
!!     real(real64), parameter :: k3 = 10.0d6
!!     real(real64), parameter :: k4 = 5.0d6
!!
!!     ! Local Variables
!!     integer(int32) :: i, j
!!     real(real64) :: m(3,3), k(3,3), natFreq(3)
!!    complex(real64) :: vals(3), modeShapes(3,3)
!!
!!     ! Define the mass matrix
!!     m = reshape([m1, 0.0d0, 0.0d0, 0.0d0, m2, 0.0d0, 0.0d0, 0.0d0, m3], [3, 3])
!!
!!     ! Define the stiffness matrix
!!     k = reshape([k1 + k2, -k2, 0.0d0, -k2, k2 + k3, -k3, 0.0d0, -k3, k3 + k4], &
!!         [3, 3])
!!
!!     ! Compute the eigenvalues and eigenvectors.
!!     call eigen(k, m, vals, vecs = modeShapes)
!!
!!     ! Compute the natural frequency values, and return them with units of Hz.
!!     ! Notice, all eigenvalues and eigenvectors are real for this example.
!!     natFreq = sqrt(real(vals)) / (2.0d0 * pi)
!!
!!     ! Display the natural frequency and mode shape values.  Notice, the eigen
!!     ! routine does not necessarily sort the values.
!!     print '(A)', "Modal Information (Not Sorted):"
!!     do i = 1, size(natFreq)
!!         print '(AI0AF8.4A)', "Mode ", i, ": (", natFreq(i), " Hz)"
!!         print '(F10.3)', (real(modeShapes(j,i)), j = 1, size(natFreq))
!!     end do
!! end program
!! @endcode
!! The above program produces the following output.
!! @code{.txt}
!! Modal Information:
!! Mode 1: (232.9225 Hz)
!!     -0.718
!!     -1.000
!!     -0.747
!! Mode 2: (749.6189 Hz)
!!     -0.419
!!     -0.164
!!      1.000
!! Mode 3: (923.5669 Hz)
!!      1.000
!!     -0.184
!!      0.179
!! @endcode
!!
!! @par See Also
!! - [Wikipedia](https://en.wikipedia.org/wiki/Eigenvalues_and_eigenvectors)
!! - [Wolfram MathWorld](http://mathworld.wolfram.com/Eigenvalue.html)
!! - [LAPACK Users Manual](http://netlib.org/lapack/lug/node56.html)
interface eigen
    module procedure :: eigen_symm
    module procedure :: eigen_asymm
    module procedure :: eigen_gen
end interface

! ------------------------------------------------------------------------------
!> @brief Sorts an array.
interface sort
    module procedure :: sort_dbl_array
    module procedure :: sort_dbl_array_ind
    module procedure :: sort_cmplx_array
    module procedure :: sort_cmplx_array_ind
    module procedure :: sort_eigen_cmplx
    module procedure :: sort_eigen_dbl
end interface


! ******************************************************************************
! LINALG_BASIC.F90
! ------------------------------------------------------------------------------
interface
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
    module subroutine mtx_mult_mtx(transa, transb, alpha, a, b, beta, c, err)
        logical, intent(in) :: transa, transb
        real(real64), intent(in) :: alpha, beta
        real(real64), intent(in), dimension(:,:) :: a, b
        real(real64), intent(inout), dimension(:,:) :: c
        class(errors), intent(inout), optional, target :: err
    end subroutine

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
    module subroutine mtx_mult_vec(trans, alpha, a, b, beta, c, err)
        logical, intent(in) :: trans
        real(real64), intent(in) :: alpha, beta
        real(real64), intent(in), dimension(:,:) :: a
        real(real64), intent(in), dimension(:) :: b
        real(real64), intent(inout), dimension(:) :: c
        class(errors), intent(inout), optional, target :: err
    end subroutine

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
    module subroutine rank1_update_dbl(alpha, x, y, a, err)
        real(real64), intent(in) :: alpha
        real(real64), intent(in), dimension(:) :: x, y
        real(real64), intent(inout), dimension(:,:) :: a
        class(errors), intent(inout), optional, target :: err
    end subroutine

    !> @brief Computes the matrix operation: C = alpha * A * op(B) + beta * C,
    !! or C = alpha * op(B) * A + beta * C.
    !!
    !! @param[in] lside Set to true to apply matrix A from the left; else, set
    !!  to false to apply matrix A from the left.
    !! @param[in] trans Set to true if op(B) == B**T; else, set to false if
    !!  op(B) == B.
    !! @param[in] alpha A scalar multiplier.
    !! @param[in] a A K-element array containing the diagonal elements of A
    !!  where K = MIN(M,P) if @p lside is true; else, if @p lside is
    !!  false, K = MIN(N,P).
    !! @param[in] b The LDB-by-TDB matrix B where (LDB = leading dimension of B,
    !!  and TDB = trailing dimension of B):
    !!  - @p lside == true & @p trans == true: LDB = N, TDB = P
    !!  - @p lside == true & @p trans == false: LDB = P, TDB = N
    !!  - @p lside == false & @p trans == true: LDB = P, TDB = M
    !!  - @p lside == false & @p trans == false: LDB = M, TDB = P
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
    module subroutine diag_mtx_mult_mtx(lside, trans, alpha, a, b, beta, c, err)
        logical, intent(in) :: lside, trans
        real(real64) :: alpha, beta
        real(real64), intent(in), dimension(:) :: a
        real(real64), intent(in), dimension(:,:) :: b
        real(real64), intent(inout), dimension(:,:) :: c
        class(errors), intent(inout), optional, target :: err
    end subroutine

    !> @brief Computes the matrix operation: B = alpha * A * op(B), or
    !! B = alpha * op(B) * A.
    !!
    !! @param[in] lside Set to true to apply matrix A from the left; else, set
    !!  to false to apply matrix A from the left.
    !! @param[in] alpha A scalar multiplier.
    !! @param[in] a A K-element array containing the diagonal elements of A
    !!  where K = MIN(M,P) if @p lside is true; else, if @p lside is
    !!  false, K = MIN(N,P).
    !! @param[in] b On input, the M-by-N matrix B.  On output, the resulting
    !!  M-by-N matrix.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - LA_ARRAY_SIZE_ERROR: Occurs if any of the input array sizes are
    !!      incorrect.
    module subroutine diag_mtx_mult_mtx2(lside, alpha, a, b, err)
        logical, intent(in) :: lside
        real(real64), intent(in) :: alpha
        real(real64), intent(in), dimension(:) :: a
        real(real64), intent(inout), dimension(:,:) :: b
        class(errors), intent(inout), optional, target :: err
    end subroutine

    !> @brief Computes the matrix operation: C = alpha * A * op(B) + beta * C,
    !! or C = alpha * op(B) * A + beta * C, where A and C are complex-valued.
    !!
    !! @param[in] lside Set to true to apply matrix A from the left; else, set
    !!  to false to apply matrix A from the left.
    !! @param[in] trans Set to true if op(B) == B**T; else, set to false if
    !!  op(B) == B.
    !! @param[in] alpha A scalar multiplier.
    !! @param[in] a A K-element array containing the diagonal elements of A
    !!  where K = MIN(M,P) if @p lside is true; else, if @p lside is
    !!  false, K = MIN(N,P).
    !! @param[in] b The LDB-by-TDB matrix B where (LDB = leading dimension of B,
    !!  and TDB = trailing dimension of B):
    !!  - @p lside == true & @p trans == true: LDB = N, TDB = P
    !!  - @p lside == true & @p trans == false: LDB = P, TDB = N
    !!  - @p lside == false & @p trans == true: LDB = P, TDB = M
    !!  - @p lside == false & @p trans == false: LDB = M, TDB = P
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
    module subroutine diag_mtx_mult_mtx3(lside, trans, alpha, a, b, beta, c, err)
        logical, intent(in) :: lside, trans
        real(real64) :: alpha, beta
        complex(real64), intent(in), dimension(:) :: a
        real(real64), intent(in), dimension(:,:) :: b
        complex(real64), intent(inout), dimension(:,:) :: c
        class(errors), intent(inout), optional, target :: err
    end subroutine

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
    !!  where K = MIN(M,P) if @p lside is true; else, if @p lside is
    !!  false, K = MIN(N,P).
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
    module subroutine diag_mtx_mult_mtx4(lside, trans, alpha, a, b, beta, c, err)
        logical, intent(in) :: lside, trans
        real(real64) :: alpha, beta
        complex(real64), intent(in), dimension(:) :: a
        complex(real64), intent(in), dimension(:,:) :: b
        complex(real64), intent(inout), dimension(:,:) :: c
        class(errors), intent(inout), optional, target :: err
    end subroutine

    !> @brief Computes the matrix operation: C = alpha * A * op(B) + beta * C,
    !! or C = alpha * op(B) * A + beta * C.
    !!
    !! @param[in] lside Set to true to apply matrix A from the left; else, set
    !!  to false to apply matrix A from the left.
    !! @param[in] trans Set to true if op(B) == B**T; else, set to false if
    !!  op(B) == B.
    !! @param[in] alpha A scalar multiplier.
    !! @param[in] a A K-element array containing the diagonal elements of A
    !!  where K = MIN(M,P) if @p lside is true; else, if @p lside is
    !!  false, K = MIN(N,P).
    !! @param[in] b The LDB-by-TDB matrix B where (LDB = leading dimension of B,
    !!  and TDB = trailing dimension of B):
    !!  - @p lside == true & @p trans == true: LDB = N, TDB = P
    !!  - @p lside == true & @p trans == false: LDB = P, TDB = N
    !!  - @p lside == false & @p trans == true: LDB = P, TDB = M
    !!  - @p lside == false & @p trans == false: LDB = M, TDB = P
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
    module subroutine diag_mtx_mult_mtx_cmplx(lside, trans, alpha, a, b, beta, c, err)
        logical, intent(in) :: lside, trans
        complex(real64) :: alpha, beta
        complex(real64), intent(in), dimension(:) :: a
        complex(real64), intent(in), dimension(:,:) :: b
        complex(real64), intent(inout), dimension(:,:) :: c
        class(errors), intent(inout), optional, target :: err
    end subroutine

    !> @brief Computes the matrix operation: B = alpha * A * op(B), or
    !! B = alpha * op(B) * A.
    !!
    !! @param[in] lside Set to true to apply matrix A from the left; else, set
    !!  to false to apply matrix A from the left.
    !! @param[in] alpha A scalar multiplier.
    !! @param[in] a A K-element array containing the diagonal elements of A
    !!  where K = MIN(M,P) if @p lside is true; else, if @p lside is
    !!  false, K = MIN(N,P).
    !! @param[in] b On input, the M-by-N matrix B.  On output, the resulting
    !!  M-by-N matrix.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - LA_ARRAY_SIZE_ERROR: Occurs if any of the input array sizes are
    !!      incorrect.
    module subroutine diag_mtx_mult_mtx2_cmplx(lside, alpha, a, b, err)
        logical, intent(in) :: lside
        complex(real64), intent(in) :: alpha
        complex(real64), intent(in), dimension(:) :: a
        complex(real64), intent(inout), dimension(:,:) :: b
        class(errors), intent(inout), optional, target :: err
    end subroutine

    !> @brief Computes the trace of a matrix (the sum of the main diagonal
    !! elements).
    !!
    !! @param[in] x The matrix on which to operate.
    !!
    !! @return The trace of @p x.
    pure module function trace_dbl(x) result(y)
        real(real64), intent(in), dimension(:,:) :: x
        real(real64) :: y
    end function

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
    module function mtx_rank_dbl(a, tol, work, olwork, err) result(rnk)
        real(real64), intent(inout), dimension(:,:) :: a
        real(real64), intent(in), optional :: tol
        real(real64), intent(out), pointer, optional, dimension(:) :: work
        integer(int32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err
        integer(int32) :: rnk
    end function

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
    module function det_dbl(a, iwork, err) result(x)
        real(real64), intent(inout), dimension(:,:) :: a
        integer(int32), intent(out), pointer, optional, dimension(:) :: iwork
        class(errors), intent(inout), optional, target :: err
        real(real64) :: x
    end function

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
    module subroutine swap_dbl(x, y, err)
        real(real64), intent(inout), dimension(:) :: x, y
        class(errors), intent(inout), optional, target :: err
    end subroutine

    !> @brief Multiplies a vector by the reciprocal of a real scalar.
    !!
    !! @param[in] a The scalar which is used to divide each component of @p X.
    !!  The value must be >= 0, or the subroutine will divide by zero.
    !! @param[in,out] x The vector.
    !!
    !! @par Notes
    !! This routine is based upon the LAPACK routine DRSCL.
    module subroutine recip_mult_array_dbl(a, x)
        real(real64), intent(in) :: a
        real(real64), intent(inout), dimension(:) :: x
    end subroutine

    !> @brief Computes the triangular matrix operation:
    !! B = alpha * A**T * A + beta * B, or B = alpha * A * A**T + beta * B,
    !! where A is a triangular matrix.
    !!
    !! @param[in] upper Set to true if matrix A is upper triangular, and
    !!  B = alpha * A**T * A + beta * B is to be calculated; else, set to false
    !!  if A is lower triangular, and B = alpha * A * A**T + beta * B is to
    !!  be computed.
    !! @param[in] alpha A scalar multiplier.
    !! @param[in] a The N-by-N triangular matrix.  Notice, if @p upper is true
    !!  only the upper triangular portion of this matrix is referenced; else,
    !!  if @p upper is false, only the lower triangular portion of this matrix
    !!  is referenced.
    !! @param[in] beta A scalar multiplier.
    !! @param[in,out] b On input, the N-by-N matrix B.  On output, the N-by-N
    !!  solution matrix.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - LA_ARRAY_SIZE_ERROR: Occurs if any of the input arrays are not sized
    !!      appropriately.
    module subroutine tri_mtx_mult_dbl(upper, alpha, a, beta, b, err)
        logical, intent(in) :: upper
        real(real64), intent(in) :: alpha, beta
        real(real64), intent(in), dimension(:,:) :: a
        real(real64), intent(inout), dimension(:,:) :: b
        class(errors), intent(inout), optional, target :: err
    end subroutine

    !> @brief Computes the triangular matrix operation:
    !! B = alpha * A**T * A + beta * B, or B = alpha * A * A**T + beta * B,
    !! where A is a triangular matrix.
    !!
    !! @param[in] upper Set to true if matrix A is upper triangular, and
    !!  B = alpha * A**T * A + beta * B is to be calculated; else, set to false
    !!  if A is lower triangular, and B = alpha * A * A**T + beta * B is to
    !!  be computed.
    !! @param[in] alpha A scalar multiplier.
    !! @param[in] a The N-by-N triangular matrix.  Notice, if @p upper is true
    !!  only the upper triangular portion of this matrix is referenced; else,
    !!  if @p upper is false, only the lower triangular portion of this matrix
    !!  is referenced.
    !! @param[in] beta A scalar multiplier.
    !! @param[in,out] b On input, the N-by-N matrix B.  On output, the N-by-N
    !!  solution matrix.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - LA_ARRAY_SIZE_ERROR: Occurs if any of the input arrays are not sized
    !!      appropriately.
    module subroutine tri_mtx_mult_cmplx(upper, alpha, a, beta, b, err)
        logical, intent(in) :: upper
        complex(real64), intent(in) :: alpha, beta
        complex(real64), intent(in), dimension(:,:) :: a
        complex(real64), intent(inout), dimension(:,:) :: b
        class(errors), intent(inout), optional, target :: err
    end subroutine

end interface

! ******************************************************************************
! LINALG_FACTOR.F90
! ------------------------------------------------------------------------------
interface
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
    !! @par Notes
    !! This routine utilizes the LAPACK routine DGETRF.
    !!
    !! @par See Also
    !! - [Wikipedia](https://en.wikipedia.org/wiki/LU_decomposition)
    !! - [Wolfram MathWorld](http://mathworld.wolfram.com/LUDecomposition.html)
    module subroutine lu_factor_dbl(a, ipvt, err)
        real(real64), intent(inout), dimension(:,:) :: a
        integer(int32), intent(out), dimension(:) :: ipvt
        class(errors), intent(inout), optional, target :: err
    end subroutine

    !> @brief Computes the LU factorization of a complex-valued M-by-N matrix.
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
    !! @par Notes
    !! This routine utilizes the LAPACK routine ZGETRF.
    !!
    !! @par See Also
    !! - [Wikipedia](https://en.wikipedia.org/wiki/LU_decomposition)
    !! - [Wolfram MathWorld](http://mathworld.wolfram.com/LUDecomposition.html)
    module subroutine lu_factor_cmplx(a, ipvt, err)
        complex(real64), intent(inout), dimension(:,:) :: a
        integer(int32), intent(out), dimension(:) :: ipvt
        class(errors), intent(inout), optional, target :: err
    end subroutine

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
    module subroutine form_lu_all(lu, ipvt, u, p, err)
        real(real64), intent(inout), dimension(:,:) :: lu
        integer(int32), intent(in), dimension(:) :: ipvt
        real(real64), intent(out), dimension(:,:) :: u, p
        class(errors), intent(inout), optional, target :: err
    end subroutine

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
    module subroutine form_lu_all_cmplx(lu, ipvt, u, p, err)
        complex(real64), intent(inout), dimension(:,:) :: lu
        integer(int32), intent(in), dimension(:) :: ipvt
        complex(real64), intent(out), dimension(:,:) :: u
        real(real64), intent(out), dimension(:,:) :: p
        class(errors), intent(inout), optional, target :: err
    end subroutine

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
    module subroutine form_lu_only(lu, u, err)
        real(real64), intent(inout), dimension(:,:) :: lu
        real(real64), intent(out), dimension(:,:) :: u
        class(errors), intent(inout), optional, target :: err
    end subroutine

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
    module subroutine form_lu_only_cmplx(lu, u, err)
        complex(real64), intent(inout), dimension(:,:) :: lu
        complex(real64), intent(out), dimension(:,:) :: u
        class(errors), intent(inout), optional, target :: err
    end subroutine

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
    module subroutine qr_factor_no_pivot(a, tau, work, olwork, err)
        real(real64), intent(inout), dimension(:,:) :: a
        real(real64), intent(out), dimension(:) :: tau
        real(real64), intent(out), target, dimension(:), optional :: work
        integer(int32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err
    end subroutine

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
    module subroutine qr_factor_pivot(a, tau, jpvt, work, olwork, err)
        real(real64), intent(inout), dimension(:,:) :: a
        real(real64), intent(out), dimension(:) :: tau
        integer(int32), intent(inout), dimension(:) :: jpvt
        real(real64), intent(out), target, dimension(:), optional :: work
        integer(int32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err
    end subroutine

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
    module subroutine form_qr_no_pivot(r, tau, q, work, olwork, err)
        real(real64), intent(inout), dimension(:,:) :: r
        real(real64), intent(in), dimension(:) :: tau
        real(real64), intent(out), dimension(:,:) :: q
        real(real64), intent(out), target, dimension(:), optional :: work
        integer(int32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err
    end subroutine

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
    module subroutine form_qr_pivot(r, tau, pvt, q, p, work, olwork, err)
        real(real64), intent(inout), dimension(:,:) :: r
        real(real64), intent(in), dimension(:) :: tau
        integer(int32), intent(in), dimension(:) :: pvt
        real(real64), intent(out), dimension(:,:) :: q, p
        real(real64), intent(out), target, dimension(:), optional :: work
        integer(int32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err
    end subroutine

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
    module subroutine mult_qr_mtx(lside, trans, a, tau, c, work, olwork, err)
        logical, intent(in) :: lside, trans
        real(real64), intent(in), dimension(:) :: tau
        real(real64), intent(inout), dimension(:,:) :: a, c
        real(real64), intent(out), target, dimension(:), optional :: work
        integer(int32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err
    end subroutine

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
    module subroutine mult_qr_vec(trans, a, tau, c, work, olwork, err)
        logical, intent(in) :: trans
        real(real64), intent(inout), dimension(:,:) :: a
        real(real64), intent(in), dimension(:) :: tau
        real(real64), intent(inout), dimension(:) :: c
        real(real64), intent(out), target, dimension(:), optional :: work
        integer(int32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err
    end subroutine

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
    !! @par Notes
    !! This routine utilizes the QRUPDATE routine DQR1UP.
    !!
    !! @par See Also
    !! [Source](https://sourceforge.net/projects/qrupdate/)
    module subroutine qr_rank1_update_dbl(q, r, u, v, work, err)
        real(real64), intent(inout), dimension(:,:) :: q, r
        real(real64), intent(inout), dimension(:) :: u, v
        real(real64), intent(out), target, optional, dimension(:) :: work
        class(errors), intent(inout), optional, target :: err
    end subroutine

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
    !! @par Notes
    !! This routine utilizes the LAPACK routine DPOTRF.
    module subroutine cholesky_factor_dbl(a, upper, err)
        real(real64), intent(inout), dimension(:,:) :: a
        logical, intent(in), optional :: upper
        class(errors), intent(inout), optional, target :: err
    end subroutine

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
    !! @par Notes
    !! This routine utilizes the QRUPDATE routine DCH1UP.
    !!
    !! @par See Also
    !! [Source](https://sourceforge.net/projects/qrupdate/)
    module subroutine cholesky_rank1_update_dbl(r, u, work, err)
        real(real64), intent(inout), dimension(:,:) :: r
        real(real64), intent(inout), dimension(:) :: u
        real(real64), intent(out), target, optional, dimension(:) :: work
        class(errors), intent(inout), optional, target :: err
    end subroutine

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
    module subroutine cholesky_rank1_downdate_dbl(r, u, work, err)
        real(real64), intent(inout), dimension(:,:) :: r
        real(real64), intent(inout), dimension(:) :: u
        real(real64), intent(out), target, optional, dimension(:) :: work
        class(errors), intent(inout), optional, target :: err
    end subroutine

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
    module subroutine rz_factor_dbl(a, tau, work, olwork, err)
        real(real64), intent(inout), dimension(:,:) :: a
        real(real64), intent(out), dimension(:) :: tau
        real(real64), intent(out), target, optional, dimension(:) :: work
        integer(int32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err
    end subroutine

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
    module subroutine mult_rz_mtx(lside, trans, l, a, tau, c, work, olwork, err)
        logical, intent(in) :: lside, trans
        integer(int32), intent(in) :: l
        real(real64), intent(inout), dimension(:,:) :: a, c
        real(real64), intent(in), dimension(:) :: tau
        real(real64), intent(out), target, optional, dimension(:) :: work
        integer(int32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err
    end subroutine

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
    module subroutine mult_rz_vec(trans, l, a, tau, c, work, olwork, err)
        logical, intent(in) :: trans
        integer(int32), intent(in) :: l
        real(real64), intent(inout), dimension(:,:) :: a
        real(real64), intent(in), dimension(:) :: tau
        real(real64), intent(inout), dimension(:) :: c
        real(real64), intent(out), target, optional, dimension(:) :: work
        integer(int32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err
    end subroutine

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
    !! @par Notes
    !! This routine utilizes the LAPACK routine DGESVD.
    !!
    !! @par See Also
    !! - [Wikipedia](https://en.wikipedia.org/wiki/Singular_value_decomposition)
    !! - [Wolfram MathWorld](http://mathworld.wolfram.com/SingularValueDecomposition.html)
    module subroutine svd_dbl(a, s, u, vt, work, olwork, err)
        real(real64), intent(inout), dimension(:,:) :: a
        real(real64), intent(out), dimension(:) :: s
        real(real64), intent(out), optional, dimension(:,:) :: u, vt
        real(real64), intent(out), target, optional, dimension(:) :: work
        integer(int32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err
    end subroutine

end interface

! ******************************************************************************
! LINALG_SOLVE.F90
! ------------------------------------------------------------------------------
interface
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
    module subroutine solve_tri_mtx(lside, upper, trans, nounit, alpha, a, b, err)
        logical, intent(in) :: lside, upper, trans, nounit
        real(real64), intent(in) :: alpha
        real(real64), intent(in), dimension(:,:) :: a
        real(real64), intent(inout), dimension(:,:) :: b
        class(errors), intent(inout), optional, target :: err
    end subroutine

    !> @brief Solves one of the matrix equations: op(A) * X = alpha * B, or
    !! X * op(A) = alpha * B, where A is a triangular matrix.
    !!
    !! @param[in] lside Set to true to solve op(A) * X = alpha * B; else, set to
    !!  false to solve X * op(A) = alpha * B.
    !! @param[in] upper Set to true if A is an upper triangular matrix; else,
    !!  set to false if A is a lower triangular matrix.
    !! @param[in] trans Set to true if op(A) = A**H; else, set to false if
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
    !! This routine is based upon the BLAS routine ZTRSM.
    module subroutine solve_tri_mtx_cmplx(lside, upper, trans, nounit, alpha, a, b, err)
        logical, intent(in) :: lside, upper, trans, nounit
        complex(real64), intent(in) :: alpha
        complex(real64), intent(in), dimension(:,:) :: a
        complex(real64), intent(inout), dimension(:,:) :: b
        class(errors), intent(inout), optional, target :: err
    end subroutine

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
    module subroutine solve_tri_vec(upper, trans, nounit, a, x, err)
        logical, intent(in) :: upper, trans, nounit
        real(real64), intent(in), dimension(:,:) :: a
        real(real64), intent(inout), dimension(:) :: x
        class(errors), intent(inout), optional, target :: err
    end subroutine

    !> @brief Solves the system of equations: op(A) * X = B, where A is a
    !!  triangular matrix.
    !!
    !! @param[in] upper Set to true if A is an upper triangular matrix; else,
    !!  set to false if A is a lower triangular matrix.
    !! @param[in] trans Set to true if op(A) = A**H; else, set to false if
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
    !! This routine is based upon the BLAS routine ZTRSV.
    module subroutine solve_tri_vec_cmplx(upper, trans, nounit, a, x, err)
        logical, intent(in) :: upper, trans, nounit
        complex(real64), intent(in), dimension(:,:) :: a
        complex(real64), intent(inout), dimension(:) :: x
        class(errors), intent(inout), optional, target :: err
    end subroutine

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
    module subroutine solve_lu_mtx(a, ipvt, b, err)
        real(real64), intent(in), dimension(:,:) :: a
        integer(int32), intent(in), dimension(:) :: ipvt
        real(real64), intent(inout), dimension(:,:) :: b
        class(errors), intent(inout), optional, target :: err
    end subroutine

    !> @brief Solves a system of complex-valued LU-factored equations.
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
    module subroutine solve_lu_mtx_cmplx(a, ipvt, b, err)
        complex(real64), intent(in), dimension(:,:) :: a
        integer(int32), intent(in), dimension(:) :: ipvt
        complex(real64), intent(inout), dimension(:,:) :: b
        class(errors), intent(inout), optional, target :: err
    end subroutine

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
    module subroutine solve_lu_vec(a, ipvt, b, err)
        real(real64), intent(in), dimension(:,:) :: a
        integer(int32), intent(in), dimension(:) :: ipvt
        real(real64), intent(inout), dimension(:) :: b
        class(errors), intent(inout), optional, target :: err
    end subroutine

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
    module subroutine solve_lu_vec_cmplx(a, ipvt, b, err)
        complex(real64), intent(in), dimension(:,:) :: a
        integer(int32), intent(in), dimension(:) :: ipvt
        complex(real64), intent(inout), dimension(:) :: b
        class(errors), intent(inout), optional, target :: err
    end subroutine

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
    module subroutine solve_qr_no_pivot_mtx(a, tau, b, work, olwork, err)
        real(real64), intent(inout), dimension(:,:) :: a, b
        real(real64), intent(in), dimension(:) :: tau
        real(real64), intent(out), target, optional, dimension(:) :: work
        integer(int32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err
    end subroutine

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
    module subroutine solve_qr_no_pivot_vec(a, tau, b, work, olwork, err)
        real(real64), intent(inout), dimension(:,:) :: a
        real(real64), intent(in), dimension(:) :: tau
        real(real64), intent(inout), dimension(:) :: b
        real(real64), intent(out), target, optional, dimension(:) :: work
        integer(int32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err
    end subroutine

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
    module subroutine solve_qr_pivot_mtx(a, tau, jpvt, b, work, olwork, err)
        real(real64), intent(inout), dimension(:,:) :: a
        real(real64), intent(in), dimension(:) :: tau
        integer(int32), intent(in), dimension(:) :: jpvt
        real(real64), intent(inout), dimension(:,:) :: b
        real(real64), intent(out), target, optional, dimension(:) :: work
        integer(int32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err
    end subroutine

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
    module subroutine solve_qr_pivot_vec(a, tau, jpvt, b, work, olwork, err)
        real(real64), intent(inout), dimension(:,:) :: a
        real(real64), intent(in), dimension(:) :: tau
        integer(int32), intent(in), dimension(:) :: jpvt
        real(real64), intent(inout), dimension(:) :: b
        real(real64), intent(out), target, optional, dimension(:) :: work
        integer(int32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err
    end subroutine

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
    module subroutine solve_cholesky_mtx(upper, a, b, err)
        logical, intent(in) :: upper
        real(real64), intent(in), dimension(:,:) :: a
        real(real64), intent(inout), dimension(:,:) :: b
        class(errors), intent(inout), optional, target :: err
    end subroutine

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
    module subroutine solve_cholesky_vec(upper, a, b, err)
        logical, intent(in) :: upper
        real(real64), intent(in), dimension(:,:) :: a
        real(real64), intent(inout), dimension(:) :: b
        class(errors), intent(inout), optional, target :: err
    end subroutine

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
    module subroutine solve_least_squares_mtx(a, b, work, olwork, err)
        real(real64), intent(inout), dimension(:,:) :: a, b
        real(real64), intent(out), target, optional, dimension(:) :: work
        integer(int32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err
    end subroutine

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
    module subroutine solve_least_squares_vec(a, b, work, olwork, err)
        real(real64), intent(inout), dimension(:,:) :: a
        real(real64), intent(inout), dimension(:) :: b
        real(real64), intent(out), target, optional, dimension(:) :: work
        integer(int32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err
    end subroutine

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
    module subroutine solve_least_squares_mtx_pvt(a, b, ipvt, arnk, work, olwork, err)
        real(real64), intent(inout), dimension(:,:) :: a, b
        integer(int32), intent(inout), target, optional, dimension(:) :: ipvt
        integer(int32), intent(out), optional :: arnk
        real(real64), intent(out), target, optional, dimension(:) :: work
        integer(int32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err
    end subroutine

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
    module subroutine solve_least_squares_vec_pvt(a, b, ipvt, arnk, work, olwork, err)
        real(real64), intent(inout), dimension(:,:) :: a
        real(real64), intent(inout), dimension(:) :: b
        integer(int32), intent(inout), target, optional, dimension(:) :: ipvt
        integer(int32), intent(out), optional :: arnk
        real(real64), intent(out), target, optional, dimension(:) :: work
        integer(int32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err
    end subroutine

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
    !! @param[out] s An optional MIN(M, N)-element array that on output contains
    !!  the singular values of @p a in descending order.  Notice, the condition
    !!  number of @p a can be determined by S(1) / S(MIN(M, N)).
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
    !!  - LA_CONVERGENCE_ERROR: Occurs as a warning if the QR iteration process
    !!      could not converge to a zero value.
    !!
    !! @par Notes
    !! This routine utilizes the LAPACK routine DGELSS.
    module subroutine solve_least_squares_mtx_svd(a, b, s, arnk, work, olwork, err)
        real(real64), intent(inout), dimension(:,:) :: a, b
        integer(int32), intent(out), optional :: arnk
        real(real64), intent(out), target, optional, dimension(:) :: work, s
        integer(int32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err
    end subroutine

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
    !! @param[out] s An optional MIN(M, N)-element array that on output contains
    !!  the singular values of @p a in descending order.  Notice, the condition
    !!  number of @p a can be determined by S(1) / S(MIN(M, N)).
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
    !!  - LA_CONVERGENCE_ERROR: Occurs as a warning if the QR iteration process
    !!      could not converge to a zero value.
    !!
    !! @par Notes
    !! This routine utilizes the LAPACK routine DGELSS.
    module subroutine solve_least_squares_vec_svd(a, b, s, arnk, work, olwork, err)
        real(real64), intent(inout), dimension(:,:) :: a
        real(real64), intent(inout), dimension(:) :: b
        integer(int32), intent(out), optional :: arnk
        real(real64), intent(out), target, optional, dimension(:) :: work, s
        integer(int32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err
    end subroutine

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
    !! @par Notes
    !! This routine utilizes the LAPACK routines DGETRF to perform an LU
    !! factorization of the matrix, and DGETRI to invert the LU factored
    !! matrix.
    !!
    !! @par See Also
    !! - [Wikipedia](https://en.wikipedia.org/wiki/Invertible_matrix)
    !! - [Wolfram MathWorld](http://mathworld.wolfram.com/MatrixInverse.html)
    module subroutine mtx_inverse_dbl(a, iwork, work, olwork, err)
        real(real64), intent(inout), dimension(:,:) :: a
        integer(int32), intent(out), target, optional, dimension(:) :: iwork
        real(real64), intent(out), target, optional, dimension(:) :: work
        integer(int32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err
    end subroutine

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
    !! @par See Also
    !! - [Wikipedia](https://en.wikipedia.org/wiki/Moore%E2%80%93Penrose_pseudoinverse)
    !! - [Wolfram MathWorld](http://mathworld.wolfram.com/Moore-PenroseMatrixInverse.html)
    !! - [MathWorks](http://www.mathworks.com/help/matlab/ref/pinv.html?s_tid=srchtitle)
    module subroutine mtx_pinverse_dbl(a, ainv, tol, work, olwork, err)
        real(real64), intent(inout), dimension(:,:) :: a
        real(real64), intent(out), dimension(:,:) :: ainv
        real(real64), intent(in), optional :: tol
        real(real64), intent(out), target, dimension(:), optional :: work
        integer(int32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err
    end subroutine

end interface

! ******************************************************************************
! LINALG_EIGEN.F90
! ------------------------------------------------------------------------------
interface
    !> @brief Computes the eigenvalues, and optionally the eigenvectors of a
    !! real, symmetric matrix.
    !!
    !! @param[in] vecs Set to true to compute the eigenvectors as well as the
    !!  eigenvalues; else, set to false to just compute the eigenvalues.
    !! @param[in,out] a On input, the N-by-N symmetric matrix on which to
    !!  operate.  On output, and if @p vecs is set to true, the matrix will
    !!  contain the eigenvectors (one per column) corresponding to each
    !!  eigenvalue in @p vals.  If @p vecs is set to false, the lower triangular
    !!  portion of the matrix is overwritten.
    !! @param[out] vals An N-element array that will contain the eigenvalues
    !!  sorted into ascending order.
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
    !!  - LA_CONVERGENCE_ERROR: Occurs if the algorithm failed to converge.
    !!
    !! @par Notes
    !! This routine utilizes the LAPACK routine DSYEV.
    module subroutine eigen_symm(vecs, a, vals, work, olwork, err)
        logical, intent(in) :: vecs
        real(real64), intent(inout), dimension(:,:) :: a
        real(real64), intent(out), dimension(:) :: vals
        real(real64), intent(out), pointer, optional, dimension(:) :: work
        integer(int32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err
    end subroutine

    !> @brief Computes the eigenvalues, and optionally the right eigenvectors of
    !! a square matrix.
    !!
    !! @param[in,out] a On input, the N-by-N matrix on which to operate.  On
    !!  output, the contents of this matrix are overwritten.
    !! @param[out] vals An N-element array containing the eigenvalues of the
    !!  matrix.  The eigenvalues are not sorted.
    !! @param[out] vecs An optional N-by-N matrix, that if supplied, signals to
    !!  compute the right eigenvectors (one per column).  If not provided, only
    !!  the eigenvalues will be computed.
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
    !!  - LA_CONVERGENCE_ERROR: Occurs if the algorithm failed to converge.
    !!
    !! @par Notes
    !! This routine utilizes the LAPACK routine DGEEV.
    module subroutine eigen_asymm(a, vals, vecs, work, olwork, err)
        real(real64), intent(inout), dimension(:,:) :: a
        complex(real64), intent(out), dimension(:) :: vals
        complex(real64), intent(out), optional, dimension(:,:) :: vecs
        real(real64), intent(out), pointer, optional, dimension(:) :: work
        integer(int32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err
    end subroutine

    !> @brief Computes the eigenvalues, and optionally the right eigenvectors of
    !! a square matrix assuming the structure of the eigenvalue problem is
    !! A*X = lambda*B*X.
    !!
    !! @param[in,out] a On input, the N-by-N matrix A.  On output, the contents
    !!  of this matrix are overwritten.
    !! @param[in,out] b On input, the N-by-N matrix B.  On output, the contents
    !!  of this matrix are overwritten.
    !! @param[out] alpha An N-element array that, if @p beta is not supplied,
    !!  contains the eigenvalues.  If @p beta is supplied however, the
    !!  eigenvalues must be computed as ALPHA / BETA.  This however, is not as
    !!  trivial as it seems as it is entirely possible, and likely, that
    !!  ALPHA / BETA can overflow or underflow.  With that said, the values in
    !!  ALPHA will always be less than and usually comparable with the NORM(A).
    !! @param[out] beta An optional N-element array that if provided forces
    !!  @p alpha to return the numerator, and this array contains the
    !!  denominator used to determine the eigenvalues as ALPHA / BETA.  If used,
    !!  the values in this array will always be less than and usually comparable
    !!  with the NORM(B).
    !! @param[out] vecs An optional N-by-N matrix, that if supplied, signals to
    !!  compute the right eigenvectors (one per column).  If not provided, only
    !!  the eigenvalues will be computed.
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
    !!  - LA_CONVERGENCE_ERROR: Occurs if the algorithm failed to converge.
    !!
    !! @par Notes
    !! This routine utilizes the LAPACK routine DGGEV.
    module subroutine eigen_gen(a, b, alpha, beta, vecs, work, olwork, err)
        real(real64), intent(inout), dimension(:,:) :: a, b
        complex(real64), intent(out), dimension(:) :: alpha
        real(real64), intent(out), optional, dimension(:) :: beta
        complex(real64), intent(out), optional, dimension(:,:) :: vecs
        real(real64), intent(out), optional, pointer, dimension(:) :: work
        integer(int32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err
    end subroutine

end interface

! ******************************************************************************
! LINALG_SORTING.F90
! ------------------------------------------------------------------------------
interface
    !> @brief Sorts an array.
    !!
    !! @param[in,out] x On input, the array to sort.  On output, the sorted
    !!  array.
    !! @param[in] ascend An optional input that, if specified, controls if the
    !!  the array is sorted in an ascending order (default), or a descending
    !!  order.
    !!
    !! @par Remarks
    !! The routine utilizes a quick sort algorithm unless the size of the array
    !! is less than or equal to 20.  For such small arrays an insertion sort
    !! algorithm is utilized.
    !!
    !! @par Notes
    !! This routine utilizes the LAPACK routine DLASRT.
    module subroutine sort_dbl_array(x, ascend)
        real(real64), intent(inout), dimension(:) :: x
        logical, intent(in), optional :: ascend
    end subroutine

    !> @brief Sorts an array.
    !!
    !! @param[in,out] x On input, the array to sort.  On output, the sorted
    !!  array.
    !! @param[in,out] ind On input, an integer array.  On output, the contents
    !!  of this array are shifted in the same order as that of @p x as a means
    !!  of tracking the sorting operation.  It is often useful to set this
    !!  array to an ascending group of values (1, 2, ... n) such that this
    !!  array tracks the original positions of the sorted array.  Such an array
    !!  can then be used to align other arrays.  This array must be the same
    !!  size as @p x.
    !! @param[in] ascend An optional input that, if specified, controls if the
    !!  the array is sorted in an ascending order (default), or a descending
    !!  order.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - LA_ARRAY_SIZE_ERROR: Occurs if @p ind is not sized to match @p x.
    !!
    !! @par Remarks
    !! This routine utilizes a quick sort algorithm explained at
    !! http://www.fortran.com/qsort_c.f95.
    module subroutine sort_dbl_array_ind(x, ind, ascend, err)
        real(real64), intent(inout), dimension(:) :: x
        integer(int32), intent(inout), dimension(:) :: ind
        logical, intent(in), optional :: ascend
        class(errors), intent(inout), optional, target :: err
    end subroutine

    !> @brief Sorts an array.
    !!
    !! @param[in,out] x On input, the array to sort.  On output, the sorted
    !!  array.
    !! @param[in] ascend An optional input that, if specified, controls if the
    !!  the array is sorted in an ascending order (default), or a descending
    !!  order.
    !!
    !! @par Remarks
    !! This routine utilizes a quick sort algorithm.  As this routine operates
    !! on complex valued items, the complex values are sorted based upon the
    !! real component of the number.
    !!
    !! @par Notes
    !! This implementation is a slight modification of the code presented at
    !! http://www.fortran.com/qsort_c.f95.
    module subroutine sort_cmplx_array(x, ascend)
        complex(real64), intent(inout), dimension(:) :: x
        logical, intent(in), optional :: ascend
    end subroutine

    !> @brief Sorts an array.
    !!
    !! @param[in,out] x On input, the array to sort.  On output, the sorted
    !!  array.
    !! @param[in,out] ind On input, an integer array.  On output, the contents
    !!  of this array are shifted in the same order as that of @p x as a means
    !!  of tracking the sorting operation.  It is often useful to set this
    !!  array to an ascending group of values (1, 2, ... n) such that this
    !!  array tracks the original positions of the sorted array.  Such an array
    !!  can then be used to align other arrays.  This array must be the same
    !!  size as @p x.
    !! @param[in] ascend An optional input that, if specified, controls if the
    !!  the array is sorted in an ascending order (default), or a descending
    !!  order.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - LA_ARRAY_SIZE_ERROR: Occurs if @p ind is not sized to match @p x.
    !!
    !! @par Remarks
    !! This routine utilizes a quick sort algorithm.  As this routine operates
    !! on complex valued items, the complex values are sorted based upon the
    !! real component of the number.
    !!
    !! @par Notes
    !! This implementation is a slight modification of the code presented at
    !! http://www.fortran.com/qsort_c.f95.
    module subroutine sort_cmplx_array_ind(x, ind, ascend, err)
        complex(real64), intent(inout), dimension(:) :: x
        integer(int32), intent(inout), dimension(:) :: ind
        logical, intent(in), optional :: ascend
        class(errors), intent(inout), optional, target :: err
    end subroutine

    !> @brief A sorting routine specifically tailored for sorting of eigenvalues
    !! and their associated eigenvectors using a quick-sort approach.
    !!
    !! @param[in,out] vals On input, an N-element array containing the
    !!  eigenvalues.  On output, the sorted eigenvalues.
    !! @param[in,out] vecs On input, an N-by-N matrix containing the
    !!  eigenvectors associated with @p vals (one vector per column).  On
    !!  output, the sorted eigenvector matrix.
    !! @param[in] ascend An optional input that, if specified, controls if the
    !!  the array is sorted in an ascending order (default), or a descending
    !!  order.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - LA_ARRAY_SIZE_ERROR: Occurs if @p vecs is not sized to match @p vals.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
    !!      available to comoplete this operation.
    module subroutine sort_eigen_cmplx(vals, vecs, ascend, err)
        complex(real64), intent(inout), dimension(:) :: vals
        complex(real64), intent(inout), dimension(:,:) :: vecs
        logical, intent(in), optional :: ascend
        class(errors), intent(inout), optional, target :: err
    end subroutine

    !> @brief A sorting routine specifically tailored for sorting of eigenvalues
    !! and their associated eigenvectors using a quick-sort approach.
    !!
    !! @param[in,out] vals On input, an N-element array containing the
    !!  eigenvalues.  On output, the sorted eigenvalues.
    !! @param[in,out] vecs On input, an N-by-N matrix containing the
    !!  eigenvectors associated with @p vals (one vector per column).  On
    !!  output, the sorted eigenvector matrix.
    !! @param[in] ascend An optional input that, if specified, controls if the
    !!  the array is sorted in an ascending order (default), or a descending
    !!  order.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - LA_ARRAY_SIZE_ERROR: Occurs if @p vecs is not sized to match @p vals.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
    !!      available to comoplete this operation.
    module subroutine sort_eigen_dbl(vals, vecs, ascend, err)
        real(real64), intent(inout), dimension(:) :: vals
        real(real64), intent(inout), dimension(:,:) :: vecs
        logical, intent(in), optional :: ascend
        class(errors), intent(inout), optional, target :: err
    end subroutine

end interface

end module
