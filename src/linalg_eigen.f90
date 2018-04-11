! linalg_eigen.f90

!> @brief \b linalg_eigen
!!
!! @par Purpose
!! Provides routines for computing the eigenvalues and eigenvectors of matrices.
module linalg_eigen
    use, intrinsic :: iso_fortran_env, only : int32, real64
    use ferror, only : errors
    use linalg_constants
    implicit none
    private
    public :: eigen

! ******************************************************************************
! INTERFACES
! ------------------------------------------------------------------------------
    !> @brief Computes the eigenvalues, and optionally the eigenvectors, of a
    !! matrix.
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


contains
! ******************************************************************************
! EIGENVALUE/EIGENVECTOR ROUTINES
! ------------------------------------------------------------------------------
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
    subroutine eigen_symm(vecs, a, vals, work, olwork, err)
        ! Arguments
        logical, intent(in) :: vecs
        real(real64), intent(inout), dimension(:,:) :: a
        real(real64), intent(out), dimension(:) :: vals
        real(real64), intent(out), pointer, optional, dimension(:) :: work
        integer(int32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        character :: jobz
        integer(int32) :: n, istat, flag, lwork
        real(real64), pointer, dimension(:) :: wptr
        real(real64), allocatable, target, dimension(:) :: wrk
        real(real64), dimension(1) :: temp
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg

        ! Initialization
        n = size(a, 1)
        if (vecs) then
            jobz = 'V'
        else
            jobz = 'N'
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
        else if (size(vals) /= n) then
            flag = 3
        end if
        if (flag /= 0) then
            ! ERROR: One of the input arrays is not sized correctly
            write(errmsg, '(AI0A)') "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("eigen_symm", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Workspace Query
        call DSYEV(jobz, 'L', n, a, n, vals, temp, -1, flag)
        lwork = int(temp(1), int32)
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (size(work) < lwork) then
                ! ERROR: WORK not sized correctly
                call errmgr%report_error("eigen_symm", &
                    "Incorrectly sized input array WORK, argument 5.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lwork)
        else
            allocate(wrk(lwork), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("eigen_symm", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if

        ! Process
        call DSYEV(jobz, 'L', n, a, n, vals, wptr, lwork, flag)
        if (flag > 0) then
            call errmgr%report_error("eigen_symm", &
                "The algorithm failed to converge.", LA_CONVERGENCE_ERROR)
        end if
    end subroutine

! ------------------------------------------------------------------------------
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
    subroutine eigen_asymm(a, vals, vecs, work, olwork, err)
        ! Arguments
        real(real64), intent(inout), dimension(:,:) :: a
        complex(real64), intent(out), dimension(:) :: vals
        complex(real64), intent(out), optional, dimension(:,:) :: vecs
        real(real64), intent(out), pointer, optional, dimension(:) :: work
        integer(int32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        real(real64), parameter :: zero = 0.0d0
        real(real64), parameter :: two = 2.0d0

        ! Local Variables
        character :: jobvl, jobvr
        integer(int32) :: i, j, jp1, n, n1, n2a, n2b, n3a, n3b, istat, flag, &
            lwork, lwork1
        real(real64) :: eps
        real(real64), dimension(1) :: dummy, temp
        real(real64), dimension(1,1) :: dummy_mtx
        real(real64), pointer, dimension(:) :: wr, wi, wptr, w
        real(real64), pointer, dimension(:,:) :: vr
        real(real64), allocatable, target, dimension(:) :: wrk
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg

        ! Initialization
        jobvl = 'N'
        if (present(vecs)) then
            jobvr = 'V'
        else
            jobvr = 'N'
        end if
        n = size(a, 1)
        eps = two * epsilon(eps)
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        flag = 0
        if (size(a, 2) /= n) then
            flag = 1
        else if (size(vals) /= n) then
            flag = 2
        else if (present(vecs)) then
            if (size(vecs, 1) /= n .or. size(vecs, 2) /= n) then
                flag = 3
            end if
        end if
        if (flag /= 0) then
            ! ERROR: One of the input arrays is not sized correctly
            write(errmsg, '(AI0A)') "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("eigen_asymm", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Workspace Query
        call DGEEV(jobvl, jobvr, n, a, n, dummy, dummy, dummy_mtx, n, &
            dummy_mtx, n, temp, -1, flag)
        lwork1 = int(temp(1), int32)
        if (present(vecs)) then
            lwork = lwork1 + 2 * n + n * n
        else
            lwork = lwork1 + 2 * n
        end if
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (size(work) < lwork) then
                ! ERROR: WORK not sized correctly
                call errmgr%report_error("eigen_asymm", &
                    "Incorrectly sized input array WORK, argument 5.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lwork)
        else
            allocate(wrk(lwork), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("eigen_asymm", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if

        ! Locate each array within the workspace array
        n1 = n
        n2a = n1 + 1
        n2b = n2a + n - 1
        n3a = n2b + 1
        n3b = n3a + lwork1 - 1

        ! Assign pointers
        wr => wptr(1:n1)
        wi => wptr(n2a:n2b)
        w => wptr(n3a:n3b)

        ! Process
        if (present(vecs)) then
            ! Assign a pointer to the eigenvector matrix
            vr(1:n,1:n) => wptr(n3b+1:lwork)

            ! Compute the eigenvectors and eigenvalues
            call DGEEV(jobvl, jobvr, n, a, n, wr, wi, dummy_mtx, n, vr, n, &
                w, lwork1, flag)

            ! Check for convergence
            if (flag > 0) then
                call errmgr%report_error("eigen_asymm", &
                    "The algorithm failed to converge.", LA_CONVERGENCE_ERROR)
                return
            end if

            ! Store the eigenvalues and eigenvectors
            j = 1
            do while (j <= n)
                if (abs(wi(j)) < eps) then
                    ! We've got a real-valued eigenvalue
                    vals(j) = cmplx(wr(j), zero, real64)
                    do i = 1, n
                        vecs(i,j) = cmplx(vr(i,j), zero, real64)
                    end do
                else
                    ! We've got a complex cojugate pair of eigenvalues
                    jp1 = j + 1
                    vals(j) = cmplx(wr(j), wi(j), real64)
                    vals(jp1) = conjg(vals(j))
                    do i = 1, n
                        vecs(i,j) = cmplx(vr(i,j), vr(i,jp1), real64)
                        vecs(i,jp1) = conjg(vecs(i,j))
                    end do

                    ! Increment j and continue the loop
                    j = j + 2
                    cycle
                end if

                ! Increment j
                j = j + 1
            end do
        else
            ! Compute just the eigenvalues
            call DGEEV(jobvl, jobvr, n, a, n, wr, wi, dummy_mtx, n, &
                dummy_mtx, n, w, lwork1, flag)

            ! Check for convergence
            if (flag > 0) then
                call errmgr%report_error("eigen_asymm", &
                    "The algorithm failed to converge.", LA_CONVERGENCE_ERROR)
                return
            end if

            ! Store the eigenvalues
            do i = 1, n
                vals(i) = cmplx(wr(i), wi(i), real64)
            end do
        end if
    end subroutine

! ------------------------------------------------------------------------------
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
    !!     use linalg_constants, only : dp, int32
    !!     use linalg_eigen
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
    !!
    !! @par Notes
    !! This routine utilizes the LAPACK routine DGGEV.
    subroutine eigen_gen(a, b, alpha, beta, vecs, work, olwork, err)
        ! Arguments
        real(real64), intent(inout), dimension(:,:) :: a, b
        complex(real64), intent(out), dimension(:) :: alpha
        real(real64), intent(out), optional, dimension(:) :: beta
        complex(real64), intent(out), optional, dimension(:,:) :: vecs
        real(real64), intent(out), optional, pointer, dimension(:) :: work
        integer(int32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        real(real64), parameter :: zero = 0.0d0
        real(real64), parameter :: two = 2.0d0

        ! Local Variables
        character :: jobvl, jobvr
        integer(int32) :: i, j, jp1, n, n1, n2a, n2b, n3a, n3b, n4a, n4b, &
            istat, flag, lwork, lwork1
        real(real64), dimension(1) :: temp
        real(real64), dimension(1,1) :: dummy
        real(real64), pointer, dimension(:) :: wptr, w, alphar, alphai, bptr
        real(real64), pointer, dimension(:,:) :: vr
        real(real64), allocatable, target, dimension(:) :: wrk
        real(real64) :: eps
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg

        ! Initialization
        jobvl = 'N'
        jobvr = 'N'
        if (present(vecs)) then
            jobvr = 'V'
        else
            jobvr = 'N'
        end if
        n = size(a, 1)
        eps = two * epsilon(eps)
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        flag = 0
        if (size(a, 2) /= n) then
            flag = 1
        else if (size(b, 1) /= n .or. size(b, 2) /= n) then
            flag = 2
        else if (size(alpha) /= n) then
            flag = 3
        else if (present(beta)) then
            if (size(beta) /= n) flag = 4
        else if (present(vecs)) then
            if (size(vecs, 1) /= n .or. size(vecs, 2) /= n) flag = 5
        end if
        if (flag /= 0) then
            ! ERROR: One of the input arrays is not sized correctly
            write(errmsg, '(AI0A)') "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("eigen_gen", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Workspace Query
        call DGGEV(jobvl, jobvr, n, a, n, b, n, temp, temp, temp, dummy, n, &
            dummy, n, temp, -1, flag)
        lwork1 = int(temp(1), int32)
        lwork = lwork1 + 2 * n
        if (.not.present(beta)) then
            lwork = lwork + n
        end if
        if (present(vecs)) then
            lwork = lwork + n * n
        end if
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (size(work) < lwork) then
                ! ERROR: WORK not sized correctly
                call errmgr%report_error("eigen_gen", &
                    "Incorrectly sized input array WORK, argument 5.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lwork)
        else
            allocate(wrk(lwork), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("eigen_gen", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if

        ! Locate each array within the workspace array & assign pointers
        n1 = n
        n2a = n1 + 1
        n2b = n2a + n - 1
        n3a = n2b + 1
        n3b = n3a + lwork1 - 1
        n4b = n3b
        alphar => wptr(1:n1)
        alphai => wptr(n2a:n2b)
        w => wptr(n3a:n3b)
        if (.not.present(beta)) then
            n4a = n3b + 1
            n4b = n4a + n - 1
            bptr => wptr(n4a:n4b)
        end if

        ! Process
        if (present(vecs)) then
            ! Assign a pointer to the eigenvector matrix
            vr(1:n,1:n) => wptr(n4b+1:lwork)

            ! Compute the eigenvalues and eigenvectors
            if (present(beta)) then
                call DGGEV(jobvl, jobvr, n, a, n, b, n, alphar, alphai, &
                    beta, dummy, n, vr, n, w, lwork1, flag)
            else
                call DGGEV(jobvl, jobvr, n, a, n, b, n, alphar, alphai, &
                    bptr, dummy, n, vr, n, w, lwork1, flag)
            end if

            ! Check for convergence
            if (flag > 0) then
                call errmgr%report_error("eigen_gen", &
                    "The algorithm failed to converge.", LA_CONVERGENCE_ERROR)
                return
            end if

            ! Store the eigenvalues and eigenvectors
            j = 1
            do while (j <= n)
                if (abs(alphai(j)) < eps) then
                    ! Real-Valued
                    alpha(j) = cmplx(alphar(j), zero, real64)
                    do i = 1, n
                        vecs(i,j) = cmplx(vr(i,j), zero, real64)
                    end do
                else
                    ! Complex-Valued
                    jp1 = j + 1
                    alpha(j) = cmplx(alphar(j), alphai(j), real64)
                    alpha(jp1) = cmplx(alphar(jp1), alphai(jp1), real64)
                    do i = 1, n
                        vecs(i,j) = cmplx(vr(i,j), vr(i,jp1), real64)
                        vecs(i,jp1) = conjg(vecs(i,j))
                    end do

                    ! Increment j and continue
                    j = j + 2
                    cycle
                end if

                ! Increment j
                j = j + 1
            end do
            if (.not.present(beta)) alpha = alpha / bptr
        else
            ! Compute just the eigenvalues
            if (present(beta)) then
                call DGGEV(jobvl, jobvr, n, a, n, b, n, alphar, alphai, &
                    beta, dummy, n, dummy, n, w, lwork1, flag)
            else
                call DGGEV(jobvl, jobvr, n, a, n, b, n, alphar, alphai, &
                    bptr, dummy, n, dummy, n, w, lwork1, flag)
            end if

            ! Check for convergence
            if (flag > 0) then
                call errmgr%report_error("eigen_gen", &
                    "The algorithm failed to converge.", LA_CONVERGENCE_ERROR)
                return
            end if

            ! Store the eigenvalues
            do i = 1, n
                alpha(i) = cmplx(alphar(i), alphai(i), real64)
            end do
            if (.not.present(beta)) alpha = alpha / bptr
        end if
    end subroutine



end module
