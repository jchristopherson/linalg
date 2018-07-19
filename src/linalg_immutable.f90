! linalg_immutable.f90

!> @brief \b linalg_immutable
!! @par Purpose
!! Provides an immutable interface to many of the core linear algebra routines
!! in this library.  The intent is to allow for ease of use in situations
!! where memory allocation, or absolute speed are of lesser importance to code
!! readability.
!!
!! @par
!! Routines in this module do not provide an error handling interface.  Any
!! errors encountered will result in an error message printed to the prompt,
!! the generation (or appending to) an error file, and termination of the
!! program.  Notice, warning situations will be handled similarily, but without
!! termination of the program.
module linalg_immutable
    use, intrinsic :: iso_fortran_env, only : int32, real64
    use linalg_core
    implicit none
    private
    public :: mat_rank1_update
    public :: mat_mult_diag
    public :: mat_mult_upper_tri
    public :: mat_mult_lower_tri
    public :: mat_det
    public :: mat_lu
    public :: mat_qr
    public :: mat_qr_rank1_update
    public :: mat_svd
    public :: mat_cholesky
    public :: mat_cholesky_rank1_update
    public :: mat_cholesky_rank1_downdate
    public :: mat_inverse
    public :: mat_pinverse
    public :: mat_solve_upper_tri
    public :: mat_solve_lower_tri
    public :: mat_eigen
    public :: lu_results
    public :: lu_results_cmplx
    public :: qr_results
    public :: svd_results
    public :: eigen_results
    public :: identity

! ------------------------------------------------------------------------------
    !> @brief Computes the matrix operation: C = A * B, where A is a
    !! diagonal matrix.
    interface mat_mult_diag
        module procedure :: mat_mult_diag_1
        module procedure :: mat_mult_diag_2
        module procedure :: mat_mult_diag_3
        module procedure :: mat_mult_diag_1_cmplx
        module procedure :: mat_mult_diag_2_cmplx
        module procedure :: mat_mult_diag_3_cmplx
    end interface

! ------------------------------------------------------------------------------
    !> @brief Computes the matrix operation C = A * B, where A is an upper
    !! triangular matrix.
    interface mat_mult_upper_tri
        module procedure :: mat_mult_upper_tri_1
        module procedure :: mat_mult_upper_tri_2
        module procedure :: mat_mult_upper_tri_1_cmplx
        module procedure :: mat_mult_upper_tri_2_cmplx
    end interface

! ------------------------------------------------------------------------------
    !> @brief Computes the matrix operation C = A * B, where A is a lower
    !! triangular matrix.
    interface mat_mult_lower_tri
        module procedure :: mat_mult_lower_tri_1
        module procedure :: mat_mult_lower_tri_2
        module procedure :: mat_mult_lower_tri_1_cmplx
        module procedure :: mat_mult_lower_tri_2_cmplx
    end interface

! ------------------------------------------------------------------------------
    !> @brief Solves the upper triangular system A X = B, where A is an
    !! upper triangular matrix.
    interface mat_solve_upper_tri
        module procedure :: mat_solve_upper_tri_1
        module procedure :: mat_solve_upper_tri_2
        module procedure :: mat_solve_upper_tri_1_cmplx
        module procedure :: mat_solve_upper_tri_2_cmplx
    end interface

! ------------------------------------------------------------------------------
    !> @brief Solves the lower triangular system A X = B, where A is a
    !! lower triangular matrix.
    interface mat_solve_lower_tri
        module procedure :: mat_solve_lower_tri_1
        module procedure :: mat_solve_lower_tri_2
        module procedure :: mat_solve_lower_tri_1_cmplx
        module procedure :: mat_solve_lower_tri_2_cmplx
    end interface

! ------------------------------------------------------------------------------
    !> @brief Computes the LU factorization of a square matrix.  Notice,
    !! partial row pivoting is utilized.
    interface mat_lu
        module procedure :: mat_lu_dbl
        module procedure :: mat_lu_cmplx
    end interface

! ------------------------------------------------------------------------------
    !> @brief Computes the eigenvalues and eigenvectors (right) of a general
    !! N-by-N matrix.
    interface mat_eigen
        module procedure :: mat_eigen_1
        module procedure :: mat_eigen_2
    end interface

! ------------------------------------------------------------------------------
    !> @brief Defines a container for the output of an LU factorization.
    type lu_results
        !> The lower triangular matrix L.
        real(real64), allocatable, dimension(:,:) :: l
        !> The upper triangular matrix U.
        real(real64), allocatable, dimension(:,:) :: u
        !> The row pivot tracking matrix P where P A = L U.
        real(real64), allocatable, dimension(:,:) :: p
    end type

! ------------------------------------------------------------------------------
    !> @brief Defines a container for the output of an LU factorization.
    type lu_results_cmplx
        !> The lower triangular matrix L.
        complex(real64), allocatable, dimension(:,:) :: l
        !> The upper triangular matrix U.
        complex(real64), allocatable, dimension(:,:) :: u
        !> The row pivot tracking matrix P where P A = L U.
        real(real64), allocatable, dimension(:,:) :: p
    end type

! ------------------------------------------------------------------------------
    !> @brief Defines a container for the output of a QR factorization.
    type qr_results
        !> The M-by-M orthogonal matrix Q.
        real(real64), allocatable, dimension(:,:) :: q
        !> The M-by-N upper trapezoidal matrix R.
        real(real64), allocatable, dimension(:,:) :: r
        !> The N-by-N column pivot tracking matrix P where A P = Q R.  If no
        !! column pivoting is utilized, this matrix is left unallocated.
        real(real64), allocatable, dimension(:,:) :: p
    end type

! ------------------------------------------------------------------------------
    !> @brief Defines a container for the output of a singular value
    !! decomposition of a matrix.
    type svd_results
        !> The M-by-M orthogonal matrix U.
        real(real64), allocatable, dimension(:,:) :: u
        !> The M-by-N matrix containing the singular values on its diagonal.
        real(real64), allocatable, dimension(:,:) :: s
        !> The N-by-N transpose of the matrix V.
        real(real64), allocatable, dimension(:,:) :: vt
    end type

! ------------------------------------------------------------------------------
    !> @brief Defines a container for the output of an Eigen analysis of a
    !! square matrix.
    type eigen_results
        !> @brief An N-element array containing the eigenvalues.
        complex(real64), allocatable, dimension(:) :: values
        !> @brief An N-by-N matrix containing the N right eigenvectors (one per
        !! column).
        complex(real64), allocatable, dimension(:,:) :: vectors
    end type

contains
! ------------------------------------------------------------------------------
    !> @brief Performs the rank-1 update to matrix A such that:
    !! B = X * Y**T + A, where A is an M-by-N matrix, X is an M-element array,
    !! and N is an N-element array.
    !!
    !! @param[in] a The M-by-N matrix A.
    !! @param[in] x The M-element array X.
    !! @param[in] y THe N-element array Y.
    !! @return The resulting M-by-N matrix.
    function mat_rank1_update(a, x, y) result(b)
        ! Arguments
        real(real64), intent(in), dimension(:,:) :: a
        real(real64), intent(in), dimension(:) :: x, y
        real(real64), dimension(size(a, 1), size(a, 2)) :: b

        ! Process
        b = a
        call rank1_update(1.0d0, x, y, b)
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the matrix operation: C = A * B, where A is a
    !! diagonal matrix.
    !!
    !! @param[in] a The M-element array containing the diagonal elements of
    !!  the matrix A.
    !! @param[in] b The P-by-N matrix B where P is greater than or equal to M.
    !! @return The resulting M-by-N matrix.
    function mat_mult_diag_1(a, b) result(c)
        ! Arguments
        real(real64), intent(in), dimension(:) :: a
        real(real64), intent(in), dimension(:,:) :: b
        real(real64), dimension(size(a), size(b, 2)) :: c

        ! Process
        if (size(b, 1) > size(a)) then
            call diag_mtx_mult(.true., .false., 1.0d0, a, b(1:size(a),:), &
                0.0d0, c)
        else
            call diag_mtx_mult(.true., .false., 1.0d0, a, b, 0.0d0, c)
        end if
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the matrix operation: C = A * B, where A is a
    !! diagonal matrix.
    !!
    !! @param[in] a The M-element array containing the diagonal elements of
    !!  the matrix A.
    !! @param[in] b The P-element array B where P is greater than or equal to M.
    !! @return The resulting M-element array.
    function mat_mult_diag_2(a, b) result(c)
        ! Arguments
        real(real64), intent(in), dimension(:) :: a, b
        real(real64), dimension(size(a)) :: c

        ! Local Variables
        real(real64), dimension(size(a), 1) :: bc, cc

        ! Process
        bc(:,1) = b(1:min(size(a), size(b)))
        call diag_mtx_mult(.true., .false., 1.0d0, a, bc, 0.0d0, cc)
        c = cc(:,1)
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the matrix operation: C = A * B, where B is a diagonal
    !! matrix.
    !!
    !! @param[in] a The M-by-N matrix A.
    !! @param[in] b The P-element array containing the diagonal matrix B where
    !!  P is at least N.
    !! @return The resulting M-by-P matrix.
    function mat_mult_diag_3(a, b) result(c)
        ! Arguments
        real(real64), intent(in), dimension(:,:) :: a
        real(real64), intent(in), dimension(:) :: b
        real(real64), dimension(size(a, 1), size(b)) :: c

        ! Process
        if (size(a, 2) > size(b)) then
            call diag_mtx_mult(.false., .false., 1.0d0, b, a(:,1:size(b)), &
                0.0d0, c)
        else
            call diag_mtx_mult(.false., .false., 1.0d0, b, a, 0.0d0, c)
        end if
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the matrix operation: C = A * B, where A is a
    !! diagonal matrix.
    !!
    !! @param[in] a The M-element array containing the diagonal elements of
    !!  the matrix A.
    !! @param[in] b The P-by-N matrix B where P is greater than or equal to M.
    !! @return The resulting M-by-N matrix.
    function mat_mult_diag_1_cmplx(a, b) result(c)
        ! Arguments
        complex(real64), intent(in), dimension(:) :: a
        complex(real64), intent(in), dimension(:,:) :: b
        complex(real64), dimension(size(a), size(b, 2)) :: c

        ! Parameters
        complex(real64), parameter :: zero = (0.0d0, 0.0d0)
        complex(real64), parameter :: one = (1.0d0, 0.0d0)

        ! Process
        if (size(b, 1) > size(a)) then
            call diag_mtx_mult(.true., .false., one, a, b(1:size(a),:), &
                zero, c)
        else
            call diag_mtx_mult(.true., .false., one, a, b, zero, c)
        end if
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the matrix operation: C = A * B, where A is a
    !! diagonal matrix.
    !!
    !! @param[in] a The M-element array containing the diagonal elements of
    !!  the matrix A.
    !! @param[in] b The P-element array B where P is greater than or equal to M.
    !! @return The resulting M-element array.
    function mat_mult_diag_2_cmplx(a, b) result(c)
        ! Arguments
        complex(real64), intent(in), dimension(:) :: a, b
        complex(real64), dimension(size(a)) :: c

        ! Parameters
        complex(real64), parameter :: zero = (0.0d0, 0.0d0)
        complex(real64), parameter :: one = (1.0d0, 0.0d0)

        ! Local Variables
        complex(real64), dimension(size(a), 1) :: bc, cc

        ! Process
        bc(:,1) = b(1:min(size(a), size(b)))
        call diag_mtx_mult(.true., .false., one, a, bc, zero, cc)
        c = cc(:,1)
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the matrix operation: C = A * B, where B is a diagonal
    !! matrix.
    !!
    !! @param[in] a The M-by-N matrix A.
    !! @param[in] b The P-element array containing the diagonal matrix B where
    !!  P is at least N.
    !! @return The resulting M-by-P matrix.
    function mat_mult_diag_3_cmplx(a, b) result(c)
        ! Arguments
        complex(real64), intent(in), dimension(:,:) :: a
        complex(real64), intent(in), dimension(:) :: b
        complex(real64), dimension(size(a, 1), size(b)) :: c

        ! Process
        if (size(a, 2) > size(b)) then
            call diag_mtx_mult(.false., .false., 1.0d0, b, a(:,1:size(b)), &
                0.0d0, c)
        else
            call diag_mtx_mult(.false., .false., 1.0d0, b, a, 0.0d0, c)
        end if
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the matrix operation C = A * B, where A is an upper
    !! triangular matrix.
    !!
    !! @param[in] a The M-by-M triangular matrix A.
    !! @param[in] b The M-by-N matrix B.
    !! @return The resulting M-by-N matrix.
    function mat_mult_upper_tri_1(a, b) result(c)
        ! Arguments
        real(real64), intent(in), dimension(:,:) :: a, b
        real(real64), dimension(size(a, 1), size(b, 2)) :: c

        ! Process
        c = b
        call DTRMM('L', 'U', 'N', 'N', size(b, 1), size(b, 2), 1.0d0, &
            a, size(a, 1), c, size(c, 1))
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the matrix operation C = A * B, where A is an upper
    !! triangular matrix.
    !!
    !! @param[in] a The M-by-M triangular matrix A.
    !! @param[in] b The M-element array B.
    !! @return The resulting M-element array.
    function mat_mult_upper_tri_2(a, b) result(c)
        ! Arguments
        real(real64), intent(in), dimension(:,:) :: a
        real(real64), intent(in), dimension(:) :: b
        real(real64), dimension(size(a, 1)) :: c

        ! Process
        c = b
        call DTRMV('U', 'N', 'N', size(a, 1), a, size(a, 1), c, 1)
    end function

    ! ------------------------------------------------------------------------------
        !> @brief Computes the matrix operation C = A * B, where A is a lower
        !! triangular matrix.
        !!
        !! @param[in] a The M-by-M triangular matrix A.
        !! @param[in] b The M-by-N matrix B.
        !! @return The resulting M-by-N matrix.
        function mat_mult_lower_tri_1(a, b) result(c)
            ! Arguments
            real(real64), intent(in), dimension(:,:) :: a, b
            real(real64), dimension(size(a, 1), size(b, 2)) :: c

            ! Process
            c = b
            call DTRMM('L', 'L', 'N', 'N', size(b, 1), size(b, 2), 1.0d0, &
                a, size(a, 1), c, size(c, 1))
        end function

    ! ------------------------------------------------------------------------------
        !> @brief Computes the matrix operation C = A * B, where A is a lower
        !! triangular matrix.
        !!
        !! @param[in] a The M-by-M triangular matrix A.
        !! @param[in] b The M-element array B.
        !! @return The resulting M-element array.
        function mat_mult_lower_tri_2(a, b) result(c)
            ! Arguments
            real(real64), intent(in), dimension(:,:) :: a
            real(real64), intent(in), dimension(:) :: b
            real(real64), dimension(size(a, 1)) :: c

            ! Process
            c = b
            call DTRMV('L', 'N', 'N', size(a, 1), a, size(a, 1), c, 1)
        end function

! ------------------------------------------------------------------------------
    !> @brief Computes the matrix operation C = A * B, where A is an upper
    !! triangular matrix.
    !!
    !! @param[in] a The M-by-M triangular matrix A.
    !! @param[in] b The M-by-N matrix B.
    !! @return The resulting M-by-N matrix.
    function mat_mult_upper_tri_1_cmplx(a, b) result(c)
        ! Arguments
        complex(real64), intent(in), dimension(:,:) :: a, b
        complex(real64), dimension(size(a, 1), size(b, 2)) :: c

        ! Process
        c = b
        call ZTRMM('L', 'U', 'N', 'N', size(b, 1), size(b, 2), 1.0d0, &
            a, size(a, 1), c, size(c, 1))
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the matrix operation C = A * B, where A is an upper
    !! triangular matrix.
    !!
    !! @param[in] a The M-by-M triangular matrix A.
    !! @param[in] b The M-element array B.
    !! @return The resulting M-element array.
    function mat_mult_upper_tri_2_cmplx(a, b) result(c)
        ! Arguments
        complex(real64), intent(in), dimension(:,:) :: a
        complex(real64), intent(in), dimension(:) :: b
        complex(real64), dimension(size(a, 1)) :: c

        ! Process
        c = b
        call ZTRMV('U', 'N', 'N', size(a, 1), a, size(a, 1), c, 1)
    end function

    ! ------------------------------------------------------------------------------
        !> @brief Computes the matrix operation C = A * B, where A is a lower
        !! triangular matrix.
        !!
        !! @param[in] a The M-by-M triangular matrix A.
        !! @param[in] b The M-by-N matrix B.
        !! @return The resulting M-by-N matrix.
        function mat_mult_lower_tri_1_cmplx(a, b) result(c)
            ! Arguments
            complex(real64), intent(in), dimension(:,:) :: a, b
            complex(real64), dimension(size(a, 1), size(b, 2)) :: c

            ! Process
            c = b
            call ZTRMM('L', 'L', 'N', 'N', size(b, 1), size(b, 2), 1.0d0, &
                a, size(a, 1), c, size(c, 1))
        end function

    ! ------------------------------------------------------------------------------
        !> @brief Computes the matrix operation C = A * B, where A is a lower
        !! triangular matrix.
        !!
        !! @param[in] a The M-by-M triangular matrix A.
        !! @param[in] b The M-element array B.
        !! @return The resulting M-element array.
        function mat_mult_lower_tri_2_cmplx(a, b) result(c)
            ! Arguments
            complex(real64), intent(in), dimension(:,:) :: a
            complex(real64), intent(in), dimension(:) :: b
            complex(real64), dimension(size(a, 1)) :: c

            ! Process
            c = b
            call ZTRMV('L', 'N', 'N', size(a, 1), a, size(a, 1), c, 1)
        end function

! ------------------------------------------------------------------------------
    !> @brief Computes the determinant of a square matrix.
    !!
    !! @param[in] a The N-by-N matrix on which to operate.
    !! @return The determinant of the matrix.
    function mat_det(a) result(x)
        ! Arguments
        real(real64), intent(in), dimension(:,:) :: a
        real(real64) :: x

        ! Local Variables
        real(real64), dimension(size(a, 1), size(a, 2)) :: b

        ! Process
        b = a
        x = det(b)
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the LU factorization of a square matrix.  Notice,
    !! partial row pivoting is utilized.
    !!
    !! @param[in] a The N-by-N matrix to factor.
    !! @result The L, U, and P matrices resulting from the factorization.
    function mat_lu_dbl(a) result(x)
        ! Arguments
        real(real64), intent(in), dimension(:,:) :: a
        type(lu_results) :: x

        ! Local Variables
        integer(int32) :: n
        integer(int32), allocatable, dimension(:) :: ipvt

        ! Memory Allocation
        n = size(a, 1)
        allocate(ipvt(n))
        allocate(x%l(n,n))
        allocate(x%u(n,n))
        allocate(x%p(n,n))

        ! Compute the factorization
        x%l = a
        call lu_factor(x%l, ipvt)

        ! Form L, U, and P
        call form_lu(x%l, ipvt, x%u, x%p)
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the LU factorization of a square matrix.  Notice,
    !! partial row pivoting is utilized.
    !!
    !! @param[in] a The N-by-N matrix to factor.
    !! @result The L, U, and P matrices resulting from the factorization.
    function mat_lu_cmplx(a) result(x)
        ! Arguments
        complex(real64), intent(in), dimension(:,:) :: a
        type(lu_results_cmplx) :: x

        ! Local Variables
        integer(int32) :: n
        integer(int32), allocatable, dimension(:) :: ipvt

        ! Memory Allocation
        n = size(a, 1)
        allocate(ipvt(n))
        allocate(x%l(n,n))
        allocate(x%u(n,n))
        allocate(x%p(n,n))

        ! Compute the factorization
        x%l = a
        call lu_factor(x%l, ipvt)

        ! Form L, U, and P
        call form_lu(x%l, ipvt, x%u, x%p)
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the QR factorization of an M-by-N matrix.  column
    !! pivoting can be used by this routine.
    !!
    !! @param[in] a The M-by-N matrix to factor.
    !! @param[in] pvt An optional value that, if supplied, can be used to turn
    !!  on column pivoting.  The default value is false, such that no column
    !!  pivoting is utilized.
    !! @return The Q, R, and optionally P matrices resulting from the
    !!  factorization.
    function mat_qr(a, pvt) result(x)
        ! Arguments
        real(real64), intent(in), dimension(:,:) :: a
        logical, intent(in), optional :: pvt
        type(qr_results) :: x

        ! Local Variables
        logical :: use_pivot
        integer(int32) :: m, n, mn
        integer(int32), allocatable, dimension(:) :: jpvt
        real(real64), allocatable, dimension(:) :: tau

        ! Memory Allocation
        use_pivot = .false.
        if (present(pvt)) use_pivot = pvt
        m = size(a, 1)
        n = size(a, 2)
        mn = min(m, n)
        allocate(tau(mn))
        allocate(x%q(m,m))
        allocate(x%r(m,n))

        ! Compute the factorization, and then form Q, R, and P
        x%r = a
        if (use_pivot) then
            allocate(x%p(n,n))
            allocate(jpvt(n))
            jpvt = 0 ! Ensure all columns are free columns
            call qr_factor(x%r, tau, jpvt)
            call form_qr(x%r, tau, jpvt, x%q, x%p)
        else
            call qr_factor(x%r, tau)
            call form_qr(x%r, tau, x%q)
        end if
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the rank-1 update of a QR-factored system.
    !!
    !! @param[in] q The M-by-M orthogonal matrix Q from the factorization of
    !!  the original system.
    !! @param[in] r The M-by-N upper trapezoidal matrix R from the factorization
    !!  of the original system.
    !! @param[in] x The M-element update vector.
    !! @param[in] y The N-element update vector.
    !! @return The updated Q and R matrices.
    function mat_qr_rank1_update(q, r, x, y) result(rst)
        ! Arguments
        real(real64), intent(in), dimension(:,:) :: q, r
        real(real64), intent(in), dimension(:) :: x, y
        type(qr_results) :: rst

        ! Local Variables
        integer(int32) :: i, m, n
        real(real64), allocatable, dimension(:) :: xc, yc

        ! Memory allocation
        m = size(q, 1)
        n = size(r, 2)
        allocate(xc(m))
        allocate(yc(n))
        allocate(rst%q(m,m))
        allocate(rst%r(m,n))

        ! Process
        do i = 1, m
            xc(i) = x(i)
            rst%q(:,i) = q(:,i)
        end do
        do i = 1, n
            yc(i) = y(i)
            rst%r(:,i) = r(:,i)
        end do
        call qr_rank1_update(rst%q, rst%r, xc, yc)
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the singular value decomposition of an M-by-N matrix.
    !!
    !! @param[in] a The M-by-N matrix to factor.
    !! @result The U, S, and transpose of V matrices resulting from the
    !!  factorization where A = U * S * V**T.
    function mat_svd(a) result(x)
        ! Arguments
        real(real64), intent(in), dimension(:,:) :: a
        type(svd_results) :: x

        ! Local Variables
        integer(int32) :: i, m, n, mn
        real(real64), allocatable, dimension(:) :: s
        real(real64), allocatable, dimension(:,:) :: ac

        ! Memory Allocation
        m = size(a, 1)
        n = size(a, 2)
        mn = min(m, n)
        allocate(s(mn))
        allocate(ac(m,n))
        allocate(x%u(m,m))
        allocate(x%s(m,n))
        allocate(x%vt(n,n))

        ! Process
        ac = a
        call svd(ac, s, x%u, x%vt)

        ! Extract the singular values, and populate the results matrix
        x%s = 0.0d0
        do i = 1, mn
            x%s(i,i) = s(i)
        end do
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the Cholesky factorization of a positive-definite
    !! matrix.
    !!
    !! @param[in] a The M-by-M positive definite matrix to factor.
    !! @param[in] upper An optional input that can be used to determine if the
    !!  upper triangular factorization should be computed such that
    !!  A = R**T * R, or if the lower triangular facotrization should be
    !!  computed such that A = L * L**T.  The default is true such that the
    !!  upper triangular form is computed.
    function mat_cholesky(a, upper) result(r)
        ! Arguments
        real(real64), intent(in), dimension(:,:) :: a
        logical, intent(in), optional :: upper
        real(real64), dimension(size(a, 1), size(a, 2)) :: r

        ! Local Variables
        logical :: compute_upper

        ! Process
        compute_upper = .true.
        if (present(upper)) compute_upper = upper
        r = a
        call cholesky_factor(r, compute_upper)
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the rank 1 update to a Cholesky factored matrix.
    !!
    !! @param[in] a The M-by-M upper triangular Cholesky factored matrix to
    !!  update.
    !! @param[in] x The M-element update array.
    !! @return The updated M-by-M upper triangular matrix.
    function mat_cholesky_rank1_update(a, x) result(r)
        ! Arguments
        real(real64), intent(in), dimension(:,:) :: a
        real(real64), intent(in), dimension(:) :: x
        real(real64), dimension(size(a, 1), size(a, 2)) :: r

        ! Local Variables
        real(real64), dimension(size(x)) :: xc

        ! Process
        r = a
        xc = x
        call cholesky_rank1_update(r, xc)
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the rank 1 downdate to a Cholesky factored matrix.
    !!
    !! @param[in] a The M-by-M upper triangular Cholesky factored matrix to
    !!  downdate.
    !! @param[in] x The M-element downdate array.
    !! @return The downdated M-by-M upper triangular matrix.
    function mat_cholesky_rank1_downdate(a, x) result(r)
        ! Arguments
        real(real64), intent(in), dimension(:,:) :: a
        real(real64), intent(in), dimension(:) :: x
        real(real64), dimension(size(a, 1), size(a, 2)) :: r

        ! Local Variables
        real(real64), dimension(size(x)) :: xc

        ! Process
        r = a
        xc = x
        call cholesky_rank1_downdate(r, xc)
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the inverse of a square matrix.
    !!
    !! @param[in] a The M-by-M matrix to invert.
    !! @return The M-by-M inverted matrix.
    function mat_inverse(a) result(x)
        ! Arguments
        real(real64), intent(in), dimension(:,:) :: a
        real(real64), dimension(size(a, 2), size(a, 1)) :: x

        ! Compute the inverse of A
        x = a
        call mtx_inverse(x)
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the Moore-Penrose pseudo-inverse of a M-by-N matrix.
    !!
    !! @param[in] a The M-by-N matrix to invert.
    !! @return The N-by-M inverted matrix.
    function mat_pinverse(a) result(x)
        ! Arguments
        real(real64), intent(in), dimension(:,:) :: a
        real(real64), dimension(size(a, 2), size(a, 1)) :: x

        ! Local Variables
        real(real64), dimension(size(a, 1), size(a, 2)) :: ac

        ! Compute the inverse of A
        ac = a
        call mtx_pinverse(ac, x)
    end function

! ------------------------------------------------------------------------------
    !> @brief Solves the upper triangular system A X = B, where A is an
    !! upper triangular matrix.
    !!
    !! @param[in] a The M-by-M upper triangluar matrix A.
    !! @param[in] b The M-by-NRHS matrix B.
    !! @return The M-by-NRHS solution matrix X.
    function mat_solve_upper_tri_1(a, b) result(x)
        ! Arguments
        real(real64), intent(in), dimension(:,:) :: a, b
        real(real64), dimension(size(b, 1), size(b, 2)) :: x

        ! Process
        x = b
        call solve_triangular_system(.true., .true., .false., .true., 1.0d0, &
            a, x)
    end function

! ------------------------------------------------------------------------------
    !> @brief Solves the upper triangular system A X = B, where A is an
    !! upper triangular matrix.
    !!
    !! @param[in] a The M-by-M upper triangluar matrix A.
    !! @param[in] b The M-element array B.
    !! @return The M-element solution array X.
    function mat_solve_upper_tri_2(a, b) result(x)
        ! Arguments
        real(real64), intent(in), dimension(:,:) :: a
        real(real64), intent(in), dimension(:) :: b
        real(real64), dimension(size(b)) :: x

        ! Process
        x = b
        call solve_triangular_system(.true., .false., .true., a, x)
    end function

! ------------------------------------------------------------------------------
    !> @brief Solves the upper triangular system A X = B, where A is an
    !! upper triangular matrix.
    !!
    !! @param[in] a The M-by-M upper triangluar matrix A.
    !! @param[in] b The M-by-NRHS matrix B.
    !! @return The M-by-NRHS solution matrix X.
    function mat_solve_upper_tri_1_cmplx(a, b) result(x)
        ! Arguments
        complex(real64), intent(in), dimension(:,:) :: a, b
        complex(real64), dimension(size(b, 1), size(b, 2)) :: x

        ! Process
        x = b
        call solve_triangular_system(.true., .true., .false., .true., &
            (1.0d0, 0.0d0), a, x)
    end function

! ------------------------------------------------------------------------------
    !> @brief Solves the upper triangular system A X = B, where A is an
    !! upper triangular matrix.
    !!
    !! @param[in] a The M-by-M upper triangluar matrix A.
    !! @param[in] b The M-element array B.
    !! @return The M-element solution array X.
    function mat_solve_upper_tri_2_cmplx(a, b) result(x)
        ! Arguments
        complex(real64), intent(in), dimension(:,:) :: a
        complex(real64), intent(in), dimension(:) :: b
        complex(real64), dimension(size(b)) :: x

        ! Process
        x = b
        call solve_triangular_system(.true., .false., .true., a, x)
    end function

! ------------------------------------------------------------------------------
    !> @brief Solves the lower triangular system A X = B, where A is a
    !! lower triangular matrix.
    !!
    !! @param[in] a The M-by-M lower triangluar matrix A.
    !! @param[in] b The M-by-NRHS matrix B.
    !! @return The M-by-NRHS solution matrix X.
    function mat_solve_lower_tri_1(a, b) result(x)
        ! Arguments
        real(real64), intent(in), dimension(:,:) :: a, b
        real(real64), dimension(size(b, 1), size(b, 2)) :: x

        ! Process
        x = b
        call solve_triangular_system(.true., .false., .false., .true., 1.0d0, &
            a, x)
    end function

! ------------------------------------------------------------------------------
    !> @brief Solves the lower triangular system A X = B, where A is a
    !! lower triangular matrix.
    !!
    !! @param[in] a The M-by-M lower triangluar matrix A.
    !! @param[in] b The M-element array B.
    !! @return The M-element solution array X.
    function mat_solve_lower_tri_2(a, b) result(x)
        ! Arguments
        real(real64), intent(in), dimension(:,:) :: a
        real(real64), intent(in), dimension(:) :: b
        real(real64), dimension(size(b)) :: x

        ! Process
        x = b
        call solve_triangular_system(.false., .false., .true., a, x)
    end function

! ------------------------------------------------------------------------------
    !> @brief Solves the lower triangular system A X = B, where A is a
    !! lower triangular matrix.
    !!
    !! @param[in] a The M-by-M lower triangluar matrix A.
    !! @param[in] b The M-by-NRHS matrix B.
    !! @return The M-by-NRHS solution matrix X.
    function mat_solve_lower_tri_1_cmplx(a, b) result(x)
        ! Arguments
        complex(real64), intent(in), dimension(:,:) :: a, b
        complex(real64), dimension(size(b, 1), size(b, 2)) :: x

        ! Process
        x = b
        call solve_triangular_system(.true., .false., .false., .true., &
            (1.0d0, 0.0d0), a, x)
    end function

! ------------------------------------------------------------------------------
    !> @brief Solves the lower triangular system A X = B, where A is a
    !! lower triangular matrix.
    !!
    !! @param[in] a The M-by-M lower triangluar matrix A.
    !! @param[in] b The M-element array B.
    !! @return The M-element solution array X.
    function mat_solve_lower_tri_2_cmplx(a, b) result(x)
        ! Arguments
        complex(real64), intent(in), dimension(:,:) :: a
        complex(real64), intent(in), dimension(:) :: b
        complex(real64), dimension(size(b)) :: x

        ! Process
        x = b
        call solve_triangular_system(.false., .false., .true., a, x)
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the eigenvalues and eigenvectors (right) of a general
    !! N-by-N matrix.
    !!
    !! @param[in] a The N-by-N matrix on which to operate.
    !! @return The eigenvalues and eigenvectors of the matrix.  The results are
    !!  sorted into ascending order.
    function mat_eigen_1(a) result(x)
        ! Arguments
        real(real64), intent(in), dimension(:,:) :: a
        type(eigen_results) :: x

        ! Local Variables
        integer(int32) :: n
        real(real64), dimension(size(a, 1), size(a, 2)) :: ac

        ! Memory Allocation
        n = size(a, 1)
        allocate(x%values(n))
        allocate(x%vectors(n,n))

        ! Process
        ac = a
        call eigen(ac, x%values, x%vectors)

        ! Sort the eigenvalues and eigenvectors.
        call sort(x%values, x%vectors, .true.)
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes eigenvalues and eigenvectors (right) from the eigenvalue
    !! problem: A X = lambda B X.
    !!
    !! @param[in] a The N-by-N matrix A.
    !! @param[in] b The N-by-N matrix B.
    !! @return The eigenvalues and eigenvectors.  The results are sorted into
    !!  ascending order.
    function mat_eigen_2(a, b) result(x)
        ! Arguments
        real(real64), intent(in), dimension(:,:) :: a, b
        type(eigen_results) :: x

        ! Local Variables
        integer(int32) :: i, j, n
        real(real64), dimension(size(a, 1), size(a, 2)) :: ac
        real(real64), dimension(size(b, 1), size(b, 2)) :: bc

        ! Memory Allocation
        n = size(a, 1)
        allocate(x%values(n))
        allocate(x%vectors(n,n))

        ! Process
        do j = 1, n
            do i = 1, n
                ac(i,j) = a(i,j)
                bc(i,j) = b(i,j)
            end do
        end do
        call eigen(ac, bc, x%values, vecs = x%vectors)

        ! Sort the eigenvalues and eigenvectors.
        call sort(x%values, x%vectors, .true.)
    end function

! ------------------------------------------------------------------------------
    !> @brief Creates an N-by-N identity matrix.
    !!
    !! @param[in] n The dimension of the matrix.
    !! @return The N-by-N identity matrix.
    pure function identity(n) result(x)
        integer(int32), intent(in) :: n
        real(real64), dimension(n, n) :: x
        integer(int32) :: i
        x = 0.0d0
        do i = 1, n
            x(i,i) = 1.0d0
        end do
    end function

end module
