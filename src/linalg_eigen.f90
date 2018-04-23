! linalg_eigen.f90

!> @brief \b linalg_eigen
!!
!! @par Purpose
!! Provides routines for computing the eigenvalues and eigenvectors of matrices.
submodule (linalg_core) linalg_eigen
contains
! ------------------------------------------------------------------------------
    module subroutine eigen_symm(vecs, a, vals, work, olwork, err)
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
    module subroutine eigen_asymm(a, vals, vecs, work, olwork, err)
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
    module subroutine eigen_gen(a, b, alpha, beta, vecs, work, olwork, err)
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

end submodule
