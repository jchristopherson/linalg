! linalg_basic.f90

submodule (linalg_core) linalg_basic
contains
! ******************************************************************************
! MATRIX MULTIPLICATION ROUTINES
! ------------------------------------------------------------------------------
    module subroutine mtx_mult_mtx(transa, transb, alpha, a, b, beta, c, err)
        ! Arguments
        logical, intent(in) :: transa, transb
        real(real64), intent(in) :: alpha, beta
        real(real64), intent(in), dimension(:,:) :: a, b
        real(real64), intent(inout), dimension(:,:) :: c
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        real(real64), parameter :: zero = 0.0d0
        real(real64), parameter :: one = 1.0d0

        ! Local Variables
        character :: ta, tb
        integer(int32) :: m, n, k, lda, ldb, flag
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
    module subroutine mtx_mult_vec(trans, alpha, a, b, beta, c, err)
        ! Arguments
        logical, intent(in) :: trans
        real(real64), intent(in) :: alpha, beta
        real(real64), intent(in), dimension(:,:) :: a
        real(real64), intent(in), dimension(:) :: b
        real(real64), intent(inout), dimension(:) :: c
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        character :: t
        integer(int32) :: m, n, flag
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
    module subroutine rank1_update_dbl(alpha, x, y, a, err)
        ! Arguments
        real(real64), intent(in) :: alpha
        real(real64), intent(in), dimension(:) :: x, y
        real(real64), intent(inout), dimension(:,:) :: a
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        real(real64), parameter :: zero = 0.0d0

        ! Local Variables
        integer(int32) :: j, m, n
        real(real64) :: temp
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
    module subroutine diag_mtx_mult_mtx(lside, trans, alpha, a, b, beta, c, err)
        ! Arguments
        logical, intent(in) :: lside, trans
        real(real64) :: alpha, beta
        real(real64), intent(in), dimension(:) :: a
        real(real64), intent(in), dimension(:,:) :: b
        real(real64), intent(inout), dimension(:,:) :: c
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        real(real64), parameter :: zero = 0.0d0
        real(real64), parameter :: one = 1.0d0

        ! Local Variables
        integer(int32) :: i, m, n, k, nrowb, ncolb, flag
        real(real64) :: temp
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
    module subroutine diag_mtx_mult_mtx2(lside, alpha, a, b, err)
        ! Arguments
        logical, intent(in) :: lside
        real(real64), intent(in) :: alpha
        real(real64), intent(in), dimension(:) :: a
        real(real64), intent(inout), dimension(:,:) :: b
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        real(real64), parameter :: zero = 0.0d0
        real(real64), parameter :: one = 1.0d0

        ! Local Variables
        integer(int32) :: i, m, n, k
        real(real64) :: temp
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
            call errmgr%report_error("diag_mtx_mult_mtx2", &
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
    module subroutine diag_mtx_mult_mtx3(lside, trans, alpha, a, b, beta, c, err)
        ! Arguments
        logical, intent(in) :: lside, trans
        real(real64) :: alpha, beta
        complex(real64), intent(in), dimension(:) :: a
        real(real64), intent(in), dimension(:,:) :: b
        complex(real64), intent(inout), dimension(:,:) :: c
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        complex(real64), parameter :: zero = (0.0d0, 0.0d0)
        complex(real64), parameter :: one = (1.0d0, 0.0d0)

        ! Local Variables
        integer(int32) :: i, m, n, k, nrowb, ncolb, flag
        complex(real64) :: temp
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
    module subroutine diag_mtx_mult_mtx4(lside, trans, alpha, a, b, beta, c, err)
        ! Arguments
        logical, intent(in) :: lside, trans
        real(real64) :: alpha, beta
        complex(real64), intent(in), dimension(:) :: a
        complex(real64), intent(in), dimension(:,:) :: b
        complex(real64), intent(inout), dimension(:,:) :: c
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        complex(real64), parameter :: zero = (0.0d0, 0.0d0)
        complex(real64), parameter :: one = (1.0d0, 0.0d0)

        ! Local Variables
        integer(int32) :: i, m, n, k, nrowb, ncolb, flag
        complex(real64) :: temp
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

! ------------------------------------------------------------------------------
    module subroutine diag_mtx_mult_mtx_cmplx(lside, trans, alpha, a, b, beta, c, err)
        ! Arguments
        logical, intent(in) :: lside, trans
        complex(real64) :: alpha, beta
        complex(real64), intent(in), dimension(:) :: a
        complex(real64), intent(in), dimension(:,:) :: b
        complex(real64), intent(inout), dimension(:,:) :: c
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        complex(real64), parameter :: zero = (0.0d0, 0.0d0)
        complex(real64), parameter :: one = (1.0d0, 0.0d0)

        ! Local Variables
        integer(int32) :: i, m, n, k, nrowb, ncolb, flag
        complex(real64) :: temp
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
            call errmgr%report_error("diag_mtx_mult_mtx_cmplx", trim(errmsg), &
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
    module subroutine diag_mtx_mult_mtx2_cmplx(lside, alpha, a, b, err)
        ! Arguments
        logical, intent(in) :: lside
        complex(real64), intent(in) :: alpha
        complex(real64), intent(in), dimension(:) :: a
        complex(real64), intent(inout), dimension(:,:) :: b
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        complex(real64), parameter :: zero = (0.0d0, 0.0d0)
        complex(real64), parameter :: one = (1.0d0, 0.0d0)

        ! Local Variables
        integer(int32) :: i, m, n, k
        complex(real64) :: temp
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
            call errmgr%report_error("diag_mtx_mult_mtx2_cmplx", &
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

! ******************************************************************************
! BASIC OPERATION ROUTINES
! ------------------------------------------------------------------------------
    pure module function trace_dbl(x) result(y)
        ! Arguments
        real(real64), intent(in), dimension(:,:) :: x
        real(real64) :: y

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
    module function mtx_rank_dbl(a, tol, work, olwork, err) result(rnk)
        ! Arguments
        real(real64), intent(inout), dimension(:,:) :: a
        real(real64), intent(in), optional :: tol
        real(real64), intent(out), pointer, optional, dimension(:) :: work
        integer(int32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err
        integer(int32) :: rnk

        ! External Function Interfaces
        interface
            function DLAMCH(cmach) result(x)
                use, intrinsic :: iso_fortran_env, only : real64
                character, intent(in) :: cmach
                real(real64) :: x
            end function
        end interface

        ! Local Variables
        integer(int32) :: i, m, n, mn, istat, lwork, flag
        real(real64), pointer, dimension(:) :: wptr, s, w
        real(real64), allocatable, target, dimension(:) :: wrk
        real(real64) :: t, tref, smlnum
        real(real64), dimension(1) :: dummy, temp
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
        lwork = int(temp(1), int32) + mn
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
    module function det_dbl(a, iwork, err) result(x)
        ! Arguments
        real(real64), intent(inout), dimension(:,:) :: a
        integer(int32), intent(out), pointer, optional, dimension(:) :: iwork
        class(errors), intent(inout), optional, target :: err
        real(real64) :: x

        ! Parameters
        real(real64), parameter :: zero = 0.0d0
        real(real64), parameter :: one = 1.0d0
        real(real64), parameter :: ten = 1.0d1
        real(real64), parameter :: p1 = 1.0d-1

        ! Local Variables
        integer(int32) :: i, ep, n, istat, flag
        integer(int32), pointer, dimension(:) :: ipvt
        integer(int32), allocatable, target, dimension(:) :: iwrk
        real(real64) :: temp
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
    module subroutine swap_dbl(x, y, err)
        ! Arguments
        real(real64), intent(inout), dimension(:) :: x, y
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        integer(int32) :: i, n
        real(real64) :: temp
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
    module subroutine recip_mult_array_dbl(a, x)
        ! Arguments
        real(real64), intent(in) :: a
        real(real64), intent(inout), dimension(:) :: x

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
    module subroutine tri_mtx_mult_dbl(upper, alpha, a, beta, b, err)
        ! Arguments
        logical, intent(in) :: upper
        real(real64), intent(in) :: alpha, beta
        real(real64), intent(in), dimension(:,:) :: a
        real(real64), intent(inout), dimension(:,:) :: b
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        real(real64), parameter :: zero = 0.0d0

        ! Local Variables
        integer(int32) :: i, j, k, n, d1, d2, flag
        real(real64) :: temp
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg

        ! Initialization
        n = size(a, 1)
        d1 = n
        d2 = n
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        flag = 0
        if (size(a, 2) /= n) then
            flag = 3
            d2 = size(a, 2)
        else if (size(b, 1) /= n .or. size(b, 2) /= n) then
            flag = 5
            d1 = size(b, 1)
            d2 = size(b, 2)
        end if
        if (flag /= 0) then
            ! ERROR: Incorrectly sized matrix
            write(errmsg, '(AI0AI0AI0AI0AI0A)') "The matrix at input ", flag, &
                " was not sized appropriately.  A matrix of ", n, "-by-", n, &
                "was expected, but a matrix of ", d1, "-by-", d2, " was found."
            call errmgr%report_error("tri_mtx_mult_dbl", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
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
    module subroutine tri_mtx_mult_cmplx(upper, alpha, a, beta, b, err)
        ! Arguments
        logical, intent(in) :: upper
        complex(real64), intent(in) :: alpha, beta
        complex(real64), intent(in), dimension(:,:) :: a
        complex(real64), intent(inout), dimension(:,:) :: b
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        complex(real64), parameter :: zero = (0.0d0, 0.0d0)

        ! Local Variables
        integer(int32) :: i, j, k, n, d1, d2, flag
        complex(real64) :: temp
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg

        ! Initialization
        n = size(a, 1)
        d1 = n
        d2 = n
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        flag = 0
        if (size(a, 2) /= n) then
            flag = 3
            d2 = size(a, 2)
        else if (size(b, 1) /= n .or. size(b, 2) /= n) then
            flag = 5
            d1 = size(b, 1)
            d2 = size(b, 2)
        end if
        if (flag /= 0) then
            ! ERROR: Incorrectly sized matrix
            write(errmsg, '(AI0AI0AI0AI0AI0A)') "The matrix at input ", flag, &
                " was not sized appropriately.  A matrix of ", n, "-by-", n, &
                "was expected, but a matrix of ", d1, "-by-", d2, " was found."
            call errmgr%report_error("tri_mtx_mult_cmplx", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
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
end submodule
