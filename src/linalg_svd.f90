! linalg_svd.f90
! SPDX-License-Identifier: MIT
module linalg_svd
    !! Provides singular value decompositions for dense real and complex matrices.
    use iso_fortran_env, only : int32, real64
    use lapack
    use blas
    use linalg_errors
    use linalg_types
    implicit none
    private
    public :: svd_factorization
    public :: svd

    type, extends(matrix_factorization) :: svd_factorization
        !! A container for the singular value decomposition of a matrix and
        !! associated operations.  The decomposition is \(A = U \Sigma V^T\).
        !! The solvers compute the minimum-norm least-squares solution
        !! \(X = A^+ B\) using the Moore-Penrose pseudo-inverse
        !! \(A^+ = V \Sigma^+ U^T\), where singular values below a tolerance
        !! are treated as zero.  Any matrix shape is supported.
        real(real64), private, allocatable, dimension(:,:) :: a
            !! The M-by-N matrix to factor.  Once factored, this is a workspace
            !! matrix as the original matrix is destroyed.
        real(real64), private, allocatable, dimension(:) :: workspace
            !! A workspace array large enough for any of the operations 
            !! undertaken for this type.
        real(real64), private, allocatable, dimension(:,:) :: u
            !! The orthogonal matrix U.  If skinny_u is true, this is 
            !! M-by-MIN(M,N); else, this is the full M-by-M orthogonal matrix
            !! U.
        real(real64), private, allocatable, dimension(:,:) :: vt
            !! The transpose of the orthogonal matrix V.  If skinny_v is true,
            !! this matrix is MIN(M,N)-by-N; else, this is the full N-by-N
            !! transpose of the orthogonal matrix V.
        real(real64), private, allocatable, dimension(:) :: s
            !! The MIN(M,N) element array containing the singular values.
        logical, private :: factored = .false.
            !! True if the matrix is factored; else, false.
        logical, private :: skinny_u = .false.
            !! True if the skinny version of U should be computed; else, false
            !! to compute the full U.
        logical, private :: skinny_v = .false.
            !! True if the skinny version of V should be computed; else, false
            !! to compute the full V.
        integer, private :: numerical_rank = 0
            !! The number of singular values greater than the tolerance used
            !! to define the pseudo-inverse.
    contains
        procedure, public :: is_factored => svdf_is_factored
            !! Checks to see if the matrix has been factored.
        procedure, public :: size => svdf_size
            !! Gets the size of the matrix.
        procedure, public :: factor => svdf_factor
            !! Factors the matrix.
        procedure, public :: solve_multiple_rhs => svdf_solve_multiple_rhs
            !! Solves the factored system with multiple right-hand-side vectors.
        procedure, public :: solve_single_rhs => svdf_solve_single_rhs
            !! Solves the factored system with a single right-hand-side vector.
        procedure, public :: solve_multiple_rhs_in_place => svdf_solve_multiple_rhs_in_place
            !! Solves the factored system with multiple right-hand-side vectors.
        procedure, public :: solve_single_rhs_in_place => svdf_solve_single_rhs_in_place
            !! Solves the factored system with a single right-hand-side vector.
        procedure, public :: get_factored_form => svdf_get_factored_form
            !! Gets the factored matrices \(U\), \(\Sigma\), and \(V^T\).
        procedure, public :: update => svdf_update
            !! Updates the object with a new matrix.
    end type

    interface svd_factorization
        module procedure :: svd_init
    end interface

    interface svd
        module procedure :: svd_dbl
        module procedure :: svd_cmplx
    end interface
contains
! ******************************************************************************
! SVD_FACTORIZATION INTERFACE
! ------------------------------------------------------------------------------
pure function svd_init(a, skinny_u, skinny_v) result(rst)
    !! Initializes a new [[svd_factorization]] object.
    real(real64), intent(in), dimension(:,:) :: a
        !! The M-by-N matrix to factor.
    logical, intent(in), optional :: skinny_u
        !! True if the skinny version of \(U\) should be computed; else, false
        !! to compute the full \(U\).  The default is false such that the full
        !! matrix is computed.
    logical, intent(in), optional :: skinny_v
        !! True if the skinny version of \(V\) should be computed; else, false
        !! to compute the full \(V\).  The default is false such that the full
        !! matrix is computed.
    type(svd_factorization) :: rst
        !! The new [[svd_factorization]] object.

    call rst%update(a, skinny_u = skinny_u, skinny_v = skinny_v)
end function

! ------------------------------------------------------------------------------
pure subroutine svdf_update(this, a, skinny_u, skinny_v)
    !! Updates the [[svd_factorization]] object with a new matrix.  All
    !! memory required for the factorization and the solvers is allocated
    !! here.  The matrix is only stored, not factored; therefore a call to
    !! factor is still required.  Any existing factorization is destroyed by
    !! this routine.
    class(svd_factorization), intent(inout) :: this
        !! The [[svd_factorization]] object.
    real(real64), intent(in), dimension(:,:) :: a
        !! The M-by-N matrix to factor.
    logical, intent(in), optional :: skinny_u
        !! True if the skinny version of \(U\) should be computed; else, false
        !! to compute the full \(U\).  The default is false such that the full
        !! matrix is computed.
    logical, intent(in), optional :: skinny_v
        !! True if the skinny version of \(V\) should be computed; else, false
        !! to compute the full \(V\).  The default is false such that the full
        !! matrix is computed.

    ! Local Variables
    character :: jobu, jobvt
    integer :: m, n, mn
    integer(int32) :: info, lwork
    real(real64), dimension(1) :: temp

    if (allocated(this%a)) deallocate(this%a)
    if (allocated(this%workspace)) deallocate(this%workspace)
    if (allocated(this%u)) deallocate(this%u)
    if (allocated(this%vt)) deallocate(this%vt)
    if (allocated(this%s)) deallocate(this%s)
    allocate(this%a(size(a, 1), size(a, 2)), source = a)
    this%factored = .false.
    this%numerical_rank = 0
    this%skinny_u = .false.
    this%skinny_v = .false.
    if (present(skinny_u)) this%skinny_u = skinny_u
    if (present(skinny_v)) this%skinny_v = skinny_v

    m = size(this%a, 1)
    n = size(this%a, 2)
    mn = min(m, n)
    jobu = merge('S', 'A', this%skinny_u)
    jobvt = merge('S', 'A', this%skinny_v)
    allocate(this%u(m, merge(mn, m, this%skinny_u)), &
        this%vt(merge(mn, n, this%skinny_v), n), this%s(mn))

    ! Workspace query.  The workspace is also large enough to hold the
    ! intermediate results of the solvers.
    call DGESVD(jobu, jobvt, m, n, this%a, m, this%s, this%u, m, this%vt, &
        size(this%vt, 1), temp, -1, info)
    lwork = max(int(temp(1), int32), int(mn, int32), 1_int32)
    allocate(this%workspace(lwork))
end subroutine

! ------------------------------------------------------------------------------
pure function svdf_is_factored(this) result(rst)
    !! Gets a value querying whether the matrix has been factored.
    class(svd_factorization), intent(in) :: this
        !! The [[svd_factorization]] object.
    logical :: rst
        !! True if the matrix has been factored; else, false.
    rst = this%factored
end function

! ------------------------------------------------------------------------------
pure function svdf_size(this, dim) result(rst)
    !! Retrieves the size of the requested matrix dimension.
    class(svd_factorization), intent(in) :: this
        !! The [[svd_factorization]] object.
    integer, intent(in) :: dim
        !! The dimension to query.
    integer :: rst
        !! The requested dimension size.
    if (allocated(this%a)) then
        rst = size(this%a, dim)
    else
        rst = 0
    end if
end function

! ------------------------------------------------------------------------------
subroutine svdf_factor(this, flag)
    !! Factors the matrix using the LAPACK routine DGESVD.  If the matrix is
    !! already factored, this routine does nothing.  Note that the stored
    !! matrix is destroyed by the factorization.
    class(svd_factorization), intent(inout) :: this
        !! The [[svd_factorization]] object.
    integer, intent(out), optional :: flag
        !! An error flag.  LA_NO_ERROR on success; LA_NOT_INITIALIZED_ERROR if
        !! no matrix is stored; LA_CONVERGENCE_ERROR if the algorithm fails to
        !! converge; or LA_INVALID_INPUT_ERROR for an invalid argument.

    ! Local Variables
    character :: jobu, jobvt
    integer :: m, n
    integer(int32) :: info
    real(real64) :: tol

    ! Initialization
    if (present(flag)) flag = LA_NO_ERROR
    if (.not.allocated(this%a)) then
        if (present(flag)) flag = LA_NOT_INITIALIZED_ERROR
        return
    end if
    if (this%factored) return

    m = size(this%a, 1)
    n = size(this%a, 2)
    jobu = merge('S', 'A', this%skinny_u)
    jobvt = merge('S', 'A', this%skinny_v)

    ! Factor the matrix
    call DGESVD(jobu, jobvt, m, n, this%a, m, this%s, this%u, m, this%vt, &
        size(this%vt, 1), this%workspace, size(this%workspace, kind = int32), &
        info)
    if (info > 0) then
        if (present(flag)) flag = LA_CONVERGENCE_ERROR
    else if (info < 0) then
        if (present(flag)) flag = LA_INVALID_INPUT_ERROR
    else
        this%factored = .true.
        tol = real(max(m, n), real64) * epsilon(1.0_real64) * this%s(1)
        this%numerical_rank = count(this%s > tol)
    end if
end subroutine

! ------------------------------------------------------------------------------
function svdf_prepare_solve(this, nrows) result(rst)
    !! Validates the right-hand-side dimension and factors the matrix, if
    !! necessary, prior to a solve.  The right-hand-side must have
    !! MAX(M, N) rows.
    class(svd_factorization), intent(inout) :: this
        !! The [[svd_factorization]] object.
    integer, intent(in) :: nrows
        !! The number of rows in the right-hand-side storage.
    integer :: rst
        !! LA_NO_ERROR on success; otherwise, an LA_* error code.

    rst = LA_NO_ERROR
    if (.not.allocated(this%a)) then
        rst = LA_NOT_INITIALIZED_ERROR
    else if (nrows /= max(size(this%a, 1), size(this%a, 2))) then
        rst = LA_MATRIX_DIMENSION_ERROR
    else if (.not.this%is_factored()) then
        call this%factor(flag = rst)
    end if
end function

! ------------------------------------------------------------------------------
subroutine svdf_solve_vec(this, b)
    !! Overwrites a single right-hand-side vector with the minimum-norm
    !! least-squares solution of the factored system.  The first M entries of
    !! the input vector hold the right-hand-side.  The first N entries of the
    !! output vector hold the solution; the remaining entries are unchanged.
    !! The object's workspace array is used for all intermediate results, so
    !! no allocations occur.
    class(svd_factorization), intent(inout) :: this
        !! The [[svd_factorization]] object.
    real(real64), intent(inout), dimension(:) :: b
        !! On input, the right-hand-side vector.  On output, the solution
        !! vector.

    ! Local Variables
    integer :: m, n, r

    m = size(this%a, 1)
    n = size(this%a, 2)
    r = this%numerical_rank
    if (r == 0) then
        b(1:n) = 0.0_real64
    else
        call DGEMV('T', m, r, 1.0_real64, this%u, m, b, 1, 0.0_real64, &
            this%workspace, 1)
        this%workspace(1:r) = this%workspace(1:r) / this%s(1:r)
        call DGEMV('T', r, n, 1.0_real64, this%vt, size(this%vt, 1), &
            this%workspace, 1, 0.0_real64, b, 1)
    end if
end subroutine

! ------------------------------------------------------------------------------
function svdf_solve_multiple_rhs(this, b, flag) result(rst)
    !! Solves the linear system \(A X = B \) where \(B\) is a matrix defining
    !! multiple right-hand-side vectors.  The solution is the minimum-norm
    !! least-squares solution computed with the Moore-Penrose pseudo-inverse.
    !! The number of rows in \(B\) must equal the number of rows in \(A\).
    class(svd_factorization), intent(inout) :: this
        !! The [[svd_factorization]] object.
    real(real64), intent(in), dimension(:,:) :: b
        !! The right-hand-side matrix.
    integer, intent(out), optional :: flag
        !! An error flag.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The N-by-K solution matrix.

    ! Local Variables
    integer :: info, n
    real(real64), allocatable, dimension(:,:) :: work

    info = LA_NO_ERROR
    if (.not.allocated(this%a)) then
        info = LA_NOT_INITIALIZED_ERROR
    else if (size(b, 1) /= size(this%a, 1)) then
        info = LA_MATRIX_DIMENSION_ERROR
    end if
    if (info /= LA_NO_ERROR) then
        if (present(flag)) flag = info
        allocate(rst(0, 0))
        return
    end if

    n = size(this%a, 2)
    allocate(work(max(size(this%a, 1), n), size(b, 2)))
    work = 0.0_real64
    work(1:size(b, 1), :) = b
    call this%solve_multiple_rhs_in_place(work, flag = info)
    if (present(flag)) flag = info
    allocate(rst(n, size(b, 2)), source = work(1:n, :))
end function

! ------------------------------------------------------------------------------
function svdf_solve_single_rhs(this, b, flag) result(rst)
    !! Solves the linear system \(A x = b \) where \(b\) is a vector defining
    !! a single right-hand-side.  The solution is the minimum-norm
    !! least-squares solution computed with the Moore-Penrose pseudo-inverse.
    !! The number of entries in \(b\) must equal the number of rows in \(A\).
    class(svd_factorization), intent(inout) :: this
        !! The [[svd_factorization]] object.
    real(real64), intent(in), dimension(:) :: b
        !! The right-hand-side vector.
    integer, intent(out), optional :: flag
        !! An error flag.
    real(real64), allocatable, dimension(:) :: rst
        !! The N-element solution vector.

    ! Local Variables
    integer :: info, n
    real(real64), allocatable, dimension(:) :: work

    info = LA_NO_ERROR
    if (.not.allocated(this%a)) then
        info = LA_NOT_INITIALIZED_ERROR
    else if (size(b) /= size(this%a, 1)) then
        info = LA_MATRIX_DIMENSION_ERROR
    end if
    if (info /= LA_NO_ERROR) then
        if (present(flag)) flag = info
        allocate(rst(0))
        return
    end if

    n = size(this%a, 2)
    allocate(work(max(size(this%a, 1), n)))
    work = 0.0_real64
    work(1:size(b)) = b
    call this%solve_single_rhs_in_place(work, flag = info)
    if (present(flag)) flag = info
    allocate(rst(n), source = work(1:n))
end function

! ------------------------------------------------------------------------------
subroutine svdf_solve_multiple_rhs_in_place(this, b, flag)
    !! Solves the linear system \(A X = B \) where \(B\) is a matrix defining
    !! multiple right-hand-side vectors.  The solution is the minimum-norm
    !! least-squares solution computed with the Moore-Penrose pseudo-inverse.
    !! Following the LAPACK DGELS convention, \(B\) must have MAX(M, N) rows.
    !! The first M rows hold the right-hand-side on input, and the first N
    !! rows hold the solution on output.  No allocations occur when the
    !! matrix has already been factored.
    class(svd_factorization), intent(inout) :: this
        !! The [[svd_factorization]] object.
    real(real64), intent(inout), dimension(:,:) :: b
        !! On input, the right-hand-side matrix \(B\).  On output, the solution
        !! matrix \(X\).
    integer, intent(out), optional :: flag
        !! An error flag.

    ! Local Variables
    integer :: j, info

    info = svdf_prepare_solve(this, size(b, 1))
    if (present(flag)) flag = info
    if (info /= LA_NO_ERROR) return
    do j = 1, size(b, 2)
        call svdf_solve_vec(this, b(:,j))
    end do
end subroutine

! ------------------------------------------------------------------------------
subroutine svdf_solve_single_rhs_in_place(this, b, flag)
    !! Solves the linear system \(A x = b \) where \(b\) is a vector defining
    !! a single right-hand-side.  The solution is the minimum-norm
    !! least-squares solution computed with the Moore-Penrose pseudo-inverse.
    !! Following the LAPACK DGELS convention, \(b\) must have MAX(M, N)
    !! entries.  The first M entries hold the right-hand-side on input, and
    !! the first N entries hold the solution on output.  No allocations occur
    !! when the matrix has already been factored.
    class(svd_factorization), intent(inout) :: this
        !! The [[svd_factorization]] object.
    real(real64), intent(inout), dimension(:) :: b
        !! On input, the right-hand-side vector.  On output, the solution
        !! vector.
    integer, intent(out), optional :: flag
        !! An error flag.

    ! Local Variables
    integer :: info

    info = svdf_prepare_solve(this, size(b))
    if (present(flag)) flag = info
    if (info /= LA_NO_ERROR) return
    call svdf_solve_vec(this, b)
end subroutine

! ------------------------------------------------------------------------------
subroutine svdf_get_factored_form(this, u, s, vt, flag)
    !! Gets the factored matrices \(U\), \(\Sigma\), and \(V^T\) such that
    !! \(A = U \Sigma V^T\).  If the matrix has not been factored, this
    !! routine factors it.
    class(svd_factorization), intent(inout) :: this
        !! The [[svd_factorization]] object.
    real(real64), intent(out), allocatable, dimension(:,:) :: u
        !! The orthogonal matrix \(U\).
    real(real64), intent(out), allocatable, dimension(:) :: s
        !! The singular values of the matrix, sorted in descending order.
    real(real64), intent(out), allocatable, dimension(:,:) :: vt
        !! The transpose of the orthogonal matrix \(V\).
    integer, intent(out), optional :: flag
        !! An error flag.

    ! Local Variables
    integer :: info

    if (present(flag)) flag = LA_NO_ERROR
    if (.not.allocated(this%a)) then
        if (present(flag)) flag = LA_NOT_INITIALIZED_ERROR
        return
    end if
    if (.not.this%is_factored()) then
        call this%factor(flag = info)
        if (info /= LA_NO_ERROR) then
            if (present(flag)) flag = info
            return
        end if
    end if
    allocate(u, source = this%u)
    allocate(s, source = this%s)
    allocate(vt, source = this%vt)
end subroutine

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ******************************************************************************
! PROCEDURAL INTERFACE
! ------------------------------------------------------------------------------
pure subroutine svd_dbl(a, s, u, vt)
    !! Computes the singular value decomposition of an M-by-N matrix \(A\) such 
    !! that \(A = U S V^T\) where \(U\) is an M-by-M orthogonal matrix, \(S\)
    !! is an M-by-N diagonal matrix containing the singular values, and \(V\)
    !! is an N-by-N orthogonal matrix.
    real(real64), intent(in), dimension(:,:) :: a
        !! The M-by-N matrix to factor.
    real(real64), intent(out), optional, allocatable, target, dimension(:) :: s
        !! A MIN(M, N)-element array containing the singular values of a sorted 
        !! in descending order.
    real(real64), intent(out), optional, allocatable, dimension(:,:) :: u
        !! An optional argument, that if supplied, is used to contain the 
        !! orthogonal matrix \(U\) from the decomposition.  The matrix \(U\) 
        !! contains the left singular vectors.
    real(real64), intent(out), optional, allocatable, dimension(:,:) :: vt
        !! An optional argument, that if supplied, is used to contain the 
        !! transpose of the N-by-N orthogonal matrix \(V\).  The matrix \(V\) 
        !! contains the right singular vectors.

    ! Local Variables
    character :: jobu, jobvt
    integer(int32) :: m, n, mn, lwork, flag
    real(real64), allocatable, target, dimension(:) :: w, sc
    real(real64), pointer, dimension(:) :: sptr
    real(real64), allocatable, dimension(:,:) :: ac
    real(real64), dimension(1) :: temp

    ! Initialization
    m = size(a, 1)
    n = size(a, 2)
    mn = min(m, n)
    allocate(ac(m, n), source = a)
    if (present(u)) then
        jobu = 'A'
        allocate(u(m,m))
    else
        jobu = 'N'
    end if
    if (present(vt)) then
        jobvt = 'A'
        allocate(vt(n,n))
    else
        jobvt = 'N'
    end if
    if (present(s)) then
        allocate(s(mn))
        sptr => s
    else
        allocate(sc(mn))
        sptr => sc
    end if

    ! Workspace Query
    call DGESVD(jobu, jobvt, m, n, ac, m, temp, temp, m, temp, n, temp, -1, &
        flag)
    lwork = int(temp(1), int32)
    allocate(w(lwork))

    ! Call DGESVD
    if (present(u) .and. present(vt)) then
        call DGESVD(jobu, jobvt, m, n, ac, m, sptr, u, m, vt, n, w, lwork, &
            flag)
    else if (present(u) .and. .not.present(vt)) then
        call DGESVD(jobu, jobvt, m, n, ac, m, sptr, u, m, temp, n, w, &
            lwork, flag)
    else if (.not.present(u) .and. present(vt)) then
        call DGESVD(jobu, jobvt, m, n, ac, m, sptr, temp, m, vt, n, w, &
            lwork, flag)
    else
        call DGESVD(jobu, jobvt, m, n, ac, m, sptr, temp, m, temp, n, w, &
            lwork, flag)
    end if

    ! Check for convergence
    if (flag > 0) then
        error stop LA_CONVERGENCE_ERROR
    end if
end subroutine

! ------------------------------------------------------------------------------
pure subroutine svd_cmplx(a, s, u, vt)
    !! Computes the singular value decomposition of an M-by-N matrix \(A\) such 
    !! that \(A = U S V^H\) where \(U\) is an M-by-M orthogonal matrix, \(S\)
    !! is an M-by-N diagonal matrix containing the singular values, and \(V\)
    !! is an N-by-N orthogonal matrix.
    complex(real64), intent(in), dimension(:,:) :: a
        !! The M-by-N matrix to factor.
    real(real64), intent(out), optional, allocatable, target, dimension(:) :: s
        !! A MIN(M, N)-element array containing the singular values of a sorted 
        !! in descending order.
    complex(real64), intent(out), optional, allocatable, dimension(:,:) :: u
        !! An optional argument, that if supplied, is used to contain the 
        !! orthogonal matrix \(U\) from the decomposition.  The matrix \(U\) 
        !! contains the left singular vectors.
    complex(real64), intent(out), optional, allocatable, dimension(:,:) :: vt
        !! An optional argument, that if supplied, is used to contain the 
        !! conjugate transpose of the N-by-N orthogonal matrix \(V\).  The 
        !! matrix \(V\) contains the right singular vectors.

    ! Local Variables
    character :: jobu, jobvt
    integer(int32) :: m, n, mn, lwork, flag, lrwork
    complex(real64), allocatable, dimension(:) :: w
    complex(real64), allocatable, dimension(:,:) :: ac
    complex(real64), dimension(1) :: temp
    real(real64), dimension(1) :: rtemp
    real(real64), allocatable, target, dimension(:) :: rw, sc
    real(real64), pointer, dimension(:) :: sptr

    ! Initialization
    m = size(a, 1)
    n = size(a, 2)
    mn = min(m, n)
    allocate(ac(m, n), source = a)
    lrwork = 5 * mn
    if (present(u)) then
        jobu = 'A'
        allocate(u(m,m))
    else
        jobu = 'N'
    end if
    if (present(vt)) then
        jobvt = 'A'
        allocate(vt(n,n))
    else
        jobvt = 'N'
    end if
    if (present(s)) then
        allocate(s(mn))
        sptr => s
    else
        allocate(sc(mn))
        sptr => sc
    end if

    ! Workspace Query
    call ZGESVD(jobu, jobvt, m, n, ac, m, rtemp, temp, m, temp, n, temp, -1, &
        rtemp, flag)
    lwork = int(temp(1), int32)
    allocate(w(lwork), rw(lrwork))

    ! Call ZGESVD
    if (present(u) .and. present(vt)) then
        call ZGESVD(jobu, jobvt, m, n, ac, m, sptr, u, m, vt, n, w, lwork, &
            rw, flag)
    else if (present(u) .and. .not.present(vt)) then
        call ZGESVD(jobu, jobvt, m, n, ac, m, sptr, u, m, temp, n, w, &
            lwork, rw, flag)
    else if (.not.present(u) .and. present(vt)) then
        call ZGESVD(jobu, jobvt, m, n, ac, m, sptr, temp, m, vt, n, w, &
            lwork, rw, flag)
    else
        call ZGESVD(jobu, jobvt, m, n, ac, m, sptr, temp, m, temp, n, w, &
            lwork, rw, flag)
    end if

    ! Check for convergence
    if (flag > 0) then
        error stop LA_CONVERGENCE_ERROR
    end if
end subroutine

! ------------------------------------------------------------------------------
end module
