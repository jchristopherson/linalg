! linalg_basic.f90

submodule (linalg) linalg_basic
    use blas
    use lapack
    implicit none
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
        character(len = :), allocatable :: errmsg

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
            allocate(character(len = 256) :: errmsg)
            write(errmsg, 100) &
                "Matrix dimension mismatch.  Input number ", flag, &
                " was not sized correctly."
            call errmgr%report_error("mtx_mult_mtx", errmsg, &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Call DGEMM
        call DGEMM(ta, tb, m, n, k, alpha, a, lda, b, ldb, beta, c, m)

        ! Formatting
100     format(A, I0, A)
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
        character(len = :), allocatable :: errmsg

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
            allocate(character(len = 256) :: errmsg)
            write(errmsg, 100) &
                "Matrix dimension mismatch.  Input number ", flag, &
                " was not sized correctly."
            call errmgr%report_error("mtx_mult_vec", errmsg, &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Call DGEMV
        call DGEMV(t, m, n, alpha, a, m, b, 1, beta, c, 1)

        ! Formatting
100     format(A, I0, A)
    end subroutine

! xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx !
!                           COMPLEX VALUED VERSIONS                            !
! xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx !
    module subroutine cmtx_mult_mtx(opa, opb, alpha, a, b, beta, c, err)
        ! Arguments
        integer(int32), intent(in) :: opa, opb
        complex(real64), intent(in) :: alpha, beta
        complex(real64), intent(in), dimension(:,:) :: a, b
        complex(real64), intent(inout), dimension(:,:) :: c
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        real(real64), parameter :: zero = 0.0d0
        real(real64), parameter :: one = 1.0d0

        ! Local Variables
        character :: ta, tb
        integer(int32) :: m, n, k, lda, ldb, flag
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = :), allocatable :: errmsg

        ! Initialization
        m = size(c, 1)
        n = size(c, 2)
        if (opa == LA_TRANSPOSE) then ! K = # of columns in op(A) (# of rows in op(B))
            k = size(a, 1)
            ta = 'T'
            lda = k
        else if (opa == LA_HERMITIAN_TRANSPOSE) then
            k = size(a, 1)
            ta = 'C'
            lda = k
        else
            k = size(a, 2)
            ta = 'N'
            lda = m
        end if
        if (opb == LA_TRANSPOSE) then
            tb = 'T'
            ldb = n
        else if (opb ==  LA_HERMITIAN_TRANSPOSE) then
            tb = 'C'
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
        if (opa == LA_TRANSPOSE .or. opa ==  LA_HERMITIAN_TRANSPOSE) then
            if (size(a, 2) /= m) flag = 4
        else
            if (size(a, 1) /= m) flag = 4
        end if
        if (opb == LA_TRANSPOSE .or. opb ==  LA_HERMITIAN_TRANSPOSE) then
            if (size(b, 2) /= k .or. size(b, 1) /= n) flag = 5
        else
            if (size(b, 1) /= k .or. size(b, 2) /= n) flag = 5
        end if
        if (flag /= 0) then
            ! ERROR: Matrix dimensions mismatch
            allocate(character(len = 256) :: errmsg)
            write(errmsg, 100) &
                "Matrix dimension mismatch.  Input number ", flag, &
                " was not sized correctly."
            call errmgr%report_error("cmtx_mult_mtx", errmsg, &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Call ZGEMM
        call ZGEMM(ta, tb, m, n, k, alpha, a, lda, b, ldb, beta, c, m)

        ! Formatting
100     format(A, I0, A)
    end subroutine

! ------------------------------------------------------------------------------
    module subroutine cmtx_mult_vec(opa, alpha, a, b, beta, c, err)
        ! Arguments
        integer(int32), intent(in) :: opa
        complex(real64), intent(in) :: alpha, beta
        complex(real64), intent(in), dimension(:,:) :: a
        complex(real64), intent(in), dimension(:) :: b
        complex(real64), intent(inout), dimension(:) :: c
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        character :: t
        integer(int32) :: m, n, flag
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = :), allocatable :: errmsg

        ! Initialization
        m = size(a, 1)
        n = size(a, 2)
        if (opa == LA_TRANSPOSE) then
            t = 'T'
        else if (opa ==  LA_HERMITIAN_TRANSPOSE) then
            t = 'C'
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
        if (opa == LA_TRANSPOSE .or. opa ==  LA_HERMITIAN_TRANSPOSE) then
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
            allocate(character(len = 256) :: errmsg)
            write(errmsg, 100) &
                "Matrix dimension mismatch.  Input number ", flag, &
                " was not sized correctly."
            call errmgr%report_error("cmtx_mult_vec", errmsg, &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Call ZGEMV
        call ZGEMV(t, m, n, alpha, a, m, b, 1, beta, c, 1)

        ! Formatting
100     format(A, I0, A)
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
            call errmgr%report_error("rank1_update_dbl", &
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

! xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx !
!                           COMPLEX VALUED VERSIONS                            !
! xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx !
    module subroutine rank1_update_cmplx(alpha, x, y, a, err)
        ! Arguments
        complex(real64), intent(in) :: alpha
        complex(real64), intent(in), dimension(:) :: x, y
        complex(real64), intent(inout), dimension(:,:) :: a
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        complex(real64), parameter :: zero = (0.0d0, 0.0d0)

        ! Local Variables
        integer(int32) :: j, m, n
        complex(real64) :: temp
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
            call errmgr%report_error("rank1_update_cmplx", &
                "Matrix dimension mismatch.", LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Process
        do j = 1, n
            if (y(j) /= zero) then
                temp = alpha * conjg(y(j))
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
        character(len = :), allocatable :: errmsg

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
            allocate(character(len = 256) :: errmsg)
            write(errmsg, 100) "Input number ", flag, &
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
                    c(i,:) = c(i,:) + temp * b(:,i)
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
                    c(i,:) = c(i,:) + temp * b(i,:)
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
                    c(:,i) = c(:,i) + temp * b(i,:)
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
                    c(:,i) = c(:,i) + temp * b(:,i)
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

        ! Formatting
100     format(A, I0, A)
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
                b(i,:) = temp * b(i,:)
            end do
            if (m > k) b(k+1:m,:) = zero
        else
            ! Compute B = alpha * B * A
            do i = 1, k
                temp = alpha * a(i)
                b(:,i) = temp * b(:,i)
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
        character(len = :), allocatable :: errmsg

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
            allocate(character(len = 256) :: errmsg)
            write(errmsg, 100) "Input number ", flag, &
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
                    c(i,:) = c(i,:) + temp * b(:,i)
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
                    c(i,:) = c(i,:) + temp * b(i,:)
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
                    c(:,i) = c(:,i) + temp * b(i,:)
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
                    c(:,i) = c(:,i) + temp * b(:,i)
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

        ! Formatting
100     format(A, I0, A)
    end subroutine

! ------------------------------------------------------------------------------
    module subroutine diag_mtx_mult_mtx4(lside, opb, alpha, a, b, beta, c, err)
        ! Arguments
        logical, intent(in) :: lside
        integer(int32), intent(in) :: opb
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
        character(len = :), allocatable :: errmsg

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
                if (opb == LA_TRANSPOSE .or. opb ==  LA_HERMITIAN_TRANSPOSE) then
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
                if (opb == LA_TRANSPOSE .or. opb ==  LA_HERMITIAN_TRANSPOSE) then
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
            allocate(character(len = 256) :: errmsg)
            write(errmsg, 100) "Input number ", flag, &
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
            if (opb == LA_TRANSPOSE) then
                ! Compute C = alpha * A * B**T + beta * C
                do i = 1, k
                    if (beta == zero) then
                        c(i,:) = zero
                    else if (beta /= one) then
                        c(i,:) = beta * c(i,:)
                    end if
                    temp = alpha * a(i)
                    c(i,:) = c(i,:) + temp * b(:,i)
                end do
            else if (opb ==  LA_HERMITIAN_TRANSPOSE) then
                ! Compute C = alpha * A * B**H + beta * C
                do i = 1, k
                    if (beta == zero) then
                        c(i,:) = zero
                    else if (beta /= one) then
                        c(i,:) = beta * c(i,:)
                    end if
                    temp = alpha * a(i)
                    c(i,:) = c(i,:) + temp * conjg(b(:,i))
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
                    c(i,:) = c(i,:) + temp * b(i,:)
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
            if (opb == LA_TRANSPOSE) then
                ! Compute C = alpha * B**T * A + beta * C
                do i = 1, k
                    if (beta == zero) then
                        c(:,i) = zero
                    else if (beta /= one) then
                        c(:,i) = beta * c(:,i)
                    end if
                    temp = alpha * a(i)
                    c(:,i) = c(:,i) + temp * b(i,:)
                end do
            else if (opb ==  LA_HERMITIAN_TRANSPOSE) then
                ! Compute C = alpha * B**H * A + beta * C
                do i = 1, k
                    if (beta == zero) then
                        c(:,i) = zero
                    else if (beta /= one) then
                        c(:,i) = beta * c(:,i)
                    end if
                    temp = alpha * a(i)
                    c(:,i) = c(:,i) + temp * conjg(b(i,:))
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
                    c(:,i) = c(:,i) + temp * b(:,i)
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

        ! Formatting
100     format(A, I0, A)
    end subroutine

! ------------------------------------------------------------------------------
    module subroutine diag_mtx_mult_mtx_cmplx(lside, opb, alpha, a, b, beta, c, err)
        ! Arguments
        logical, intent(in) :: lside
        integer(int32), intent(in) :: opb
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
        character(len = :), allocatable :: errmsg

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
                if (opb == LA_TRANSPOSE .or. opb ==  LA_HERMITIAN_TRANSPOSE) then
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
                if (opb == LA_TRANSPOSE .or. opb ==  LA_HERMITIAN_TRANSPOSE) then
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
            allocate(character(len = 256) :: errmsg)
            write(errmsg, 100) "Input number ", flag, &
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
            if (opb == LA_TRANSPOSE) then
                ! Compute C = alpha * A * B**T + beta * C
                do i = 1, k
                    if (beta == zero) then
                        c(i,:) = zero
                    else if (beta /= one) then
                        c(i,:) = beta * c(i,:)
                    end if
                    temp = alpha * a(i)
                    c(i,:) = c(i,:) + temp * b(:,i)
                end do
            else if (opb ==  LA_HERMITIAN_TRANSPOSE) then
                ! Compute C = alpha * A * B**H + beta * C
                do i = 1, k
                    if (beta == zero) then
                        c(i,:) = zero
                    else if (beta /= one) then
                        c(i,:) = beta * c(i,:)
                    end if
                    temp = alpha * a(i)
                    c(i,:) = c(i,:) + temp * conjg(b(:,i))
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
                    c(i,:) = c(i,:) + temp * b(i,:)
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
            if (opb == LA_TRANSPOSE) then
                ! Compute C = alpha * B**T * A + beta * C
                do i = 1, k
                    if (beta == zero) then
                        c(:,i) = zero
                    else if (beta /= one) then
                        c(:,i) = beta * c(:,i)
                    end if
                    temp = alpha * a(i)
                    c(:,i) = c(:,i) + temp * b(i,:)
                end do
            else if (opb ==  LA_HERMITIAN_TRANSPOSE) then
                ! Compute C = alpha * B**H * A + beta * C
                do i = 1, k
                    if (beta == zero) then
                        c(:,i) = zero
                    else if (beta /= one) then
                        c(:,i) = beta * c(:,i)
                    end if
                    temp = alpha * a(i)
                    c(:,i) = c(:,i) + temp * conjg(b(i,:))
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
                    c(:,i) = c(:,i) + temp * b(:,i)
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

        ! Formatting
100     format(A, I0, A)
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
                b(i,:) = temp * b(i,:)
            end do
            if (m > k) b(k+1:m,:) = zero
        else
            ! Compute B = alpha * B * A
            do i = 1, k
                temp = alpha * a(i)
                b(:,i) = temp * b(:,i)
            end do
            if (n > k) b(:,k+1:n) = zero
        end if
    end subroutine

! ------------------------------------------------------------------------------
    module subroutine diag_mtx_mult_mtx_mix(lside, opb, alpha, a, b, beta, c, err)
        ! Arguments
        logical, intent(in) :: lside
        integer(int32), intent(in) :: opb
        complex(real64) :: alpha, beta
        real(real64), intent(in), dimension(:) :: a
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
        character(len = :), allocatable :: errmsg

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
                if (opb == LA_TRANSPOSE .or. opb ==  LA_HERMITIAN_TRANSPOSE) then
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
                if (opb == LA_TRANSPOSE .or. opb ==  LA_HERMITIAN_TRANSPOSE) then
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
            allocate(character(len = 256) :: errmsg)
            write(errmsg, 100) "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("diag_mtx_mult_mtx_mix", trim(errmsg), &
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
            if (opb == LA_TRANSPOSE) then
                ! Compute C = alpha * A * B**T + beta * C
                do i = 1, k
                    if (beta == zero) then
                        c(i,:) = zero
                    else if (beta /= one) then
                        c(i,:) = beta * c(i,:)
                    end if
                    temp = alpha * a(i)
                    c(i,:) = c(i,:) + temp * b(:,i)
                end do
            else if (opb ==  LA_HERMITIAN_TRANSPOSE) then
                ! Compute C = alpha * A * B**H + beta * C
                do i = 1, k
                    if (beta == zero) then
                        c(i,:) = zero
                    else if (beta /= one) then
                        c(i,:) = beta * c(i,:)
                    end if
                    temp = alpha * a(i)
                    c(i,:) = c(i,:) + temp * conjg(b(:,i))
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
                    c(i,:) = c(i,:) + temp * b(i,:)
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
            if (opb == LA_TRANSPOSE) then
                ! Compute C = alpha * B**T * A + beta * C
                do i = 1, k
                    if (beta == zero) then
                        c(:,i) = zero
                    else if (beta /= one) then
                        c(:,i) = beta * c(:,i)
                    end if
                    temp = alpha * a(i)
                    c(:,i) = c(:,i) + temp * b(i,:)
                end do
            else if (opb ==  LA_HERMITIAN_TRANSPOSE) then
                ! Compute C = alpha * B**H * A + beta * C
                do i = 1, k
                    if (beta == zero) then
                        c(:,i) = zero
                    else if (beta /= one) then
                        c(:,i) = beta * c(:,i)
                    end if
                    temp = alpha * a(i)
                    c(:,i) = c(:,i) + temp * conjg(b(i,:))
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
                    c(:,i) = c(:,i) + temp * b(:,i)
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

        ! Formatting
100     format(A, I0, A)
    end subroutine

! ------------------------------------------------------------------------------
    module subroutine diag_mtx_mult_mtx2_mix(lside, alpha, a, b, err)
        ! Arguments
        logical, intent(in) :: lside
        complex(real64), intent(in) :: alpha
        real(real64), intent(in), dimension(:) :: a
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
                b(i,:) = temp * b(i,:)
            end do
            if (m > k) b(k+1:m,:) = zero
        else
            ! Compute B = alpha * B * A
            do i = 1, k
                temp = alpha * a(i)
                b(:,i) = temp * b(:,i)
            end do
            if (n > k) b(:,k+1:n) = zero
        end if
    end subroutine

! ------------------------------------------------------------------------------
    module subroutine diag_mtx_sparse_mult(lside, alpha, a, b, err)
        ! Arguments
        logical, intent(in) :: lside
        real(real64), intent(in) :: alpha
        real(real64), intent(in), dimension(:) :: a
        type(csr_matrix), intent(inout) :: b
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        integer(int32) :: ii, k, k1, k2, nrow
        real(real64) :: scal
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        
        ! Initialization
        nrow = size(b, 1)
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        if (lside) then
            if (size(a) /= nrow) then
                call errmgr%report_error("diag_mtx_sparse_mult", &
                    "Inner matrix dimension error.", LA_ARRAY_SIZE_ERROR)
                return
            end if
        else
            if (size(a) /= size(b, 2)) then
                call errmgr%report_error("diag_mtx_sparse_mult", &
                    "Inner matrix dimension error.", LA_ARRAY_SIZE_ERROR)
                return
            end if
        end if

        ! Process
        if (lside) then
            ! Compute B = DIAG * B
            do ii = 1, nrow
                k1 = b%ia(ii)
                k2 = b%ia(ii+1) - 1
                if (alpha == 1.0d0) then
                    scal = a(ii)
                else
                    scal = alpha * a(ii)
                end if
                do k = k1, k2
                    b%v(k) = b%v(k) * scal
                end do
            end do
        else
            ! Compute B = B * DIAG
            do ii = 1, nrow
                k1 = b%ia(ii)
                k2 = b%ia(ii+1) - 1
                if (alpha == 1.0d0) then
                    do k = k1, k2
                        b%v(k) = b%v(k) * a(b%ja(k))
                    end do
                else
                    do k = k1, k2
                        b%v(k) = alpha * b%v(k) * a(b%ja(k))
                    end do
                end if
            end do
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
    pure module function trace_cmplx(x) result(y)
        ! Arguments
        complex(real64), intent(in), dimension(:,:) :: x
        complex(real64) :: y

        ! Parameters
        complex(real64), parameter :: zero = (0.0d0, 0.0d0)

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
        real(real64), intent(out), target, optional, dimension(:) :: work
        integer(int32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err
        integer(int32) :: rnk

        ! Local Variables
        integer(int32) :: i, m, n, mn, istat, lwork, flag
        real(real64), pointer, dimension(:) :: wptr, s, w
        real(real64), allocatable, target, dimension(:) :: wrk
        real(real64) :: t, tref, smlnum
        real(real64), dimension(1) :: dummy, temp
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = :), allocatable :: errmsg

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
            allocate(character(len = 256) :: errmsg)
            write(errmsg, 100) flag, " superdiagonals could not " // &
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

        ! Formatting
100     format(I0, A)
    end function

! ------------------------------------------------------------------------------
    module function mtx_rank_cmplx(a, tol, work, olwork, rwork, err) result(rnk)
        ! Arguments
        complex(real64), intent(inout), dimension(:,:) :: a
        real(real64), intent(in), optional :: tol
        complex(real64), intent(out), target, optional, dimension(:) :: work
        integer(int32), intent(out), optional :: olwork
        real(real64), intent(out), target, optional, dimension(:) :: rwork
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
        integer(int32) :: i, m, n, mn, istat, lwork, flag, lrwork
        real(real64), pointer, dimension(:) :: s, rwptr, rw
        real(real64), allocatable, target, dimension(:) :: rwrk
        complex(real64), allocatable, target, dimension(:) :: wrk
        complex(real64), pointer, dimension(:) :: wptr
        real(real64) :: t, tref, smlnum
        real(real64), dimension(1) :: dummy
        complex(real64), dimension(1) :: cdummy, temp
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = :), allocatable :: errmsg

        ! Initialization
        m = size(a, 1)
        n = size(a, 2)
        mn = min(m, n)
        lrwork = 6 * mn
        smlnum = DLAMCH('s')
        rnk = 0
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Workspace Query
        call ZGESVD('N', 'N', m, n, a, m, dummy, cdummy, m, cdummy, n, temp, &
            -1, dummy, flag)
        lwork = int(temp(1), int32)
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (size(work) < lwork) then
                ! ERROR: WORK not sized correctly
                call errmgr%report_error("mtx_rank_cmplx", &
                    "Incorrectly sized input array WORK, argument 5.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lwork)
        else
            allocate(wrk(lwork), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("mtx_rank_cmplx", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if

        if (present(rwork)) then
            if (size(rwork) < lrwork) then
                ! ERROR: RWORK not sized correctly
                call errmgr%report_error("mtx_rank_cmplx", &
                    "Incorrectly sized input array RWORK.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            rwptr => rwork(1:lrwork)
        else
            allocate(rwrk(lrwork), stat = istat)
            if (istat /= 0) then
            end if
            rwptr => rwrk
        end if
        s => rwptr(1:mn)
        rw => rwptr(mn+1:lrwork)

        ! Compute the singular values of A
        call ZGESVD('N', 'N', m, n, a, m, s, cdummy, m, cdummy, n, wptr, &
            lwork - mn, rw, flag)
        if (flag > 0) then
            allocate(character(len = 256) :: errmsg)
            write(errmsg, 100) flag, " superdiagonals could not " // &
                "converge to zero as part of the QR iteration process."
            call errmgr%report_warning("mtx_rank_cmplx", errmsg, LA_CONVERGENCE_ERROR)
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

        ! Formatting
100     format(I0, A)
    end function

! ------------------------------------------------------------------------------
    module function det_dbl(a, iwork, err) result(x)
        ! Arguments
        real(real64), intent(inout), dimension(:,:) :: a
        integer(int32), intent(out), target, optional, dimension(:) :: iwork
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

! ------------------------------------------------------------------------------
    module function det_cmplx(a, iwork, err) result(x)
        ! Arguments
        complex(real64), intent(inout), dimension(:,:) :: a
        integer(int32), intent(out), target, optional, dimension(:) :: iwork
        class(errors), intent(inout), optional, target :: err
        complex(real64) :: x

        ! Parameters
        complex(real64), parameter :: zero = (0.0d0, 0.0d0)
        complex(real64), parameter :: one = (1.0d0, 0.0d0)
        complex(real64), parameter :: ten = (1.0d1, 0.0d0)
        complex(real64), parameter :: p1 = (1.0d-1, 0.0d0)
        real(real64), parameter :: real_one = 1.0d0
        real(real64), parameter :: real_ten = 1.0d1

        ! Local Variables
        integer(int32) :: i, ep, n, istat, flag
        integer(int32), pointer, dimension(:) :: ipvt
        integer(int32), allocatable, target, dimension(:) :: iwrk
        complex(real64) :: temp
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
            call errmgr%report_error("det_cmplx", &
                "The supplied matrix must be square.", LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Local Memory Allocation
        if (present(iwork)) then
            if (size(iwork) < n) then
                ! ERROR: WORK not sized correctly
                call errmgr%report_error("det_cmplx", &
                    "Incorrectly sized input array IWORK, argument 2.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            ipvt => iwork(1:n)
        else
            allocate(iwrk(n), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("det_cmplx", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            ipvt => iwrk
        end if

        ! Compute the LU factorization of A
        call ZGETRF(n, n, a, n, ipvt, flag)
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

            do while (abs(temp) < real_one)
                temp = ten * temp
                ep = ep - 1
            end do

            do while (abs(temp) > real_ten)
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

! ------------------------------------------------------------------------------
    module subroutine swap_cmplx(x, y, err)
        ! Arguments
        complex(real64), intent(inout), dimension(:) :: x, y
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        integer(int32) :: i, n
        complex(real64) :: temp
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
            call errmgr%report_error("swap_cmplx", &
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
        character(len = :), allocatable :: errmsg

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
            allocate(character(len = 256) :: errmsg)
            write(errmsg, 100) "The matrix at input ", flag, &
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

        ! Formatting
100     format(A, I0, A, I0, A, I0, A, I0, A, I0, A)
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
        character(len = :), allocatable :: errmsg

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
            allocate(character(len = 256) :: errmsg)
            write(errmsg, 100) "The matrix at input ", flag, &
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

        ! Formatting
100     format(A, I0, A, I0, A, I0, A, I0, A, I0, A)
    end subroutine

! ******************************************************************************
! BANDED MATRIX MULTIPLICATION ROUTINES
! ------------------------------------------------------------------------------
    module subroutine band_mtx_vec_mult_dbl(trans, kl, ku, alpha, a, x, beta, &
        y, err)
        ! Arguments
        logical, intent(in) :: trans
        integer(int32), intent(in) :: kl, ku
        real(real64), intent(in) :: alpha, beta
        real(real64), intent(in), dimension(:,:) :: a
        real(real64), intent(in), dimension(:) :: x
        real(real64), intent(inout), dimension(:) :: y
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        integer(int32) :: m, n
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        
        ! Initialization
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if
        if (trans) then
            m = size(x)
            n = size(y)
        else
            m = size(y)
            n = size(x)
        end if

        ! Input Checking
        if (kl < 0) go to 10
        if (ku < 0) go to 20
        if (size(a, 1) /= kl + ku + 1) go to 30
        if (size(a, 2) /= n) go to 30

        ! Process
        if (trans) then
            call DGBMV("T", m, n, kl, ku, alpha, a, size(a, 1), x, 1, beta, y, 1)
        else
            call DGBMV("N", m, n, kl, ku, alpha, a, size(a, 1), x, 1, beta, y, 1)
        end if

        ! End
        return

        ! KL < 0
10      continue
        call errmgr%report_error("band_mtx_vec_mult_dbl", &
            "The number of subdiagonals must be at least 0.", &
            LA_INVALID_INPUT_ERROR)
        return

        ! KU < 0
20      continue
        call errmgr%report_error("band_mtx_vec_mult_dbl", &
            "The number of superdiagonals must be at least 0.", &
            LA_INVALID_INPUT_ERROR)
        return

        ! A is incorrectly sized
30      continue
        call errmgr%report_error("band_mtx_vec_mult_dbl", &
            "The size of matrix A is not compatible with the other vectors.", &
            LA_ARRAY_SIZE_ERROR)
        return
    end subroutine

! ------------------------------------------------------------------------------
    module subroutine band_mtx_vec_mult_cmplx(trans, kl, ku, alpha, a, x, &
        beta, y, err)
        ! Arguments
        integer(int32), intent(in) :: trans
        integer(int32), intent(in) :: kl, ku
        complex(real64), intent(in) :: alpha, beta
        complex(real64), intent(in), dimension(:,:) :: a
        complex(real64), intent(in), dimension(:) :: x
        complex(real64), intent(inout), dimension(:) :: y
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        character :: op
        logical :: trns
        integer(int32) :: m, n
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        
        ! Initialization
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if
        if (trans == LA_TRANSPOSE) then
            op = "T"
            trns = .true.
        else if (trans == LA_HERMITIAN_TRANSPOSE) then
            op = "C"
            trns = .true.
        else
            op = "N"
            trns = .false.
        end if
        if (trns) then
            m = size(x)
            n = size(y)
        else
            m = size(y)
            n = size(x)
        end if

        ! Input Checking
        if (kl < 0) go to 10
        if (ku < 0) go to 20
        if (size(a, 1) /= kl + ku + 1) go to 30
        if (size(a, 2) /= n) go to 30

        ! Process
        call ZGBMV(op, m, n, kl, ku, alpha, a, size(a, 1), x, 1, beta, y, 1)

        ! End
        return

        ! KL < 0
10      continue
        call errmgr%report_error("band_mtx_vec_mult_cmplx", &
            "The number of subdiagonals must be at least 0.", &
            LA_INVALID_INPUT_ERROR)
        return

        ! KU < 0
20      continue
        call errmgr%report_error("band_mtx_vec_mult_cmplx", &
            "The number of superdiagonals must be at least 0.", &
            LA_INVALID_INPUT_ERROR)
        return

        ! A is incorrectly sized
30      continue
        call errmgr%report_error("band_mtx_vec_mult_cmplx", &
            "The size of matrix A is not compatible with the other vectors.", &
            LA_ARRAY_SIZE_ERROR)
        return
    end subroutine

! ------------------------------------------------------------------------------
    module subroutine band_to_full_mtx_dbl(kl, ku, b, f, err)
        ! Arguments
        integer(int32), intent(in) :: kl, ku
        real(real64), intent(in), dimension(:,:) :: b
        real(real64), intent(out), dimension(:,:) :: f
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        real(real64), parameter :: zero = 0.0d0

        ! Local Variables
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        integer(int32) :: i, j, k, m, n, i1, i2
        
        ! Initialization
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if
        m = size(f, 1)
        n = size(f, 2)

        ! Input Check
        if (kl < 0) go to 10
        if (ku < 0) go to 20
        if (size(b, 2) /= n) go to 30
        if (size(b, 1) /= kl + ku + 1) go to 40

        ! Process
        do j = 1, n
            k = ku + 1 - j
            i1 = max(1, j - ku)
            i2 = min(m, j + kl)
            do i = 1, i1 - 1
                f(i,j) = zero
            end do
            do i = i1, i2
                f(i,j) = b(k+i,j)
            end do
            do i = i2 + 1, m
                f(i,j) = zero
            end do
        end do

        ! End
        return

        ! KL < 0
10      continue
        call errmgr%report_error("band_to_full_mtx_dbl", &
            "The number of subdiagonals must be at least 0.", &
            LA_INVALID_INPUT_ERROR)
        return

        ! KU < 0
20      continue
        call errmgr%report_error("band_to_full_mtx_dbl", &
            "The number of superdiagonals must be at least 0.", &
            LA_INVALID_INPUT_ERROR)
        return

        ! A is incorrectly sized
30      continue
        call errmgr%report_error("band_to_full_mtx_dbl", &
            "The number of columns in the banded matrix does not match " // &
            "the number of columns in the full matrix.", &
            LA_ARRAY_SIZE_ERROR)
        return

40      continue
        call errmgr%report_error("band_to_full_mtx_dbl", &
            "The number of rows in the banded matrix does not align with " // &
            "the number of sub and super-diagonals specified.", &
            LA_ARRAY_SIZE_ERROR)
        return
    end subroutine

! ------------------------------------------------------------------------------
    module subroutine band_to_full_mtx_cmplx(kl, ku, b, f, err)
        ! Arguments
        integer(int32), intent(in) :: kl, ku
        complex(real64), intent(in), dimension(:,:) :: b
        complex(real64), intent(out), dimension(:,:) :: f
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        complex(real64), parameter :: zero = (0.0d0, 0.0d0)

        ! Local Variables
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        integer(int32) :: i, j, k, m, n, i1, i2
        
        ! Initialization
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if
        m = size(f, 1)
        n = size(f, 2)

        ! Input Check
        if (kl < 0) go to 10
        if (ku < 0) go to 20
        if (size(b, 2) /= n) go to 30
        if (size(b, 1) /= kl + ku + 1) go to 40

        ! Process
        do j = 1, n
            k = ku + 1 - j
            i1 = max(1, j - ku)
            i2 = min(m, j + kl)
            do i = 1, i1 - 1
                f(i,j) = zero
            end do
            do i = i1, i2
                f(i,j) = b(k+i,j)
            end do
            do i = i2 + 1, m
                f(i,j) = zero
            end do
        end do

        ! End
        return

        ! KL < 0
10      continue
        call errmgr%report_error("band_to_full_mtx_cmplx", &
            "The number of subdiagonals must be at least 0.", &
            LA_INVALID_INPUT_ERROR)
        return

        ! KU < 0
20      continue
        call errmgr%report_error("band_to_full_mtx_cmplx", &
            "The number of superdiagonals must be at least 0.", &
            LA_INVALID_INPUT_ERROR)
        return

        ! A is incorrectly sized
30      continue
        call errmgr%report_error("band_to_full_mtx_cmplx", &
            "The number of columns in the banded matrix does not match " // &
            "the number of columns in the full matrix.", &
            LA_ARRAY_SIZE_ERROR)
        return

40      continue
        call errmgr%report_error("band_to_full_mtx_cmplx", &
            "The number of rows in the banded matrix does not align with " // &
            "the number of sub and super-diagonals specified.", &
            LA_ARRAY_SIZE_ERROR)
        return
    end subroutine

! ------------------------------------------------------------------------------
    module subroutine band_diag_mtx_mult_dbl(left, m, kl, ku, alpha, a, b, err)
        ! Arguments
        logical, intent(in) :: left
        integer(int32), intent(in) :: m, kl, ku
        real(real64), intent(in) :: alpha
        real(real64), intent(inout), dimension(:,:) :: a
        real(real64), intent(in), dimension(:) :: b
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        real(real64), parameter :: one = 1.0d0

        ! Local Variables
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        integer(int32) :: i, i1, i2, j, k, n
        real(real64) :: temp
        
        ! Initialization
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if
        n = size(a, 2)

        ! Input Checking
        if (kl < 0) go to 10
        if (ku < 0) go to 20
        if (left) then
            if (size(b) /= n) go to 30
        else
            if (size(b) < m) go to 30
        end if

        ! Process
        if (left) then
            ! Compute A = A * B
            do j = 1, n
                k = ku + 1 - j
                i1 = max(1, j - ku) + k
                i2 = min(m, j + kl) + k
                if (alpha == one) then
                    temp = b(j)
                else
                    temp = alpha * b(j)
                end if
                do i = i1, i2
                    a(i,j) = a(i,j) * temp
                end do
            end do
        else
            ! Compute A = B * A
            do j = 1, n
                k = ku + 1 - j
                i1 = max(1, j - ku)
                i2 = min(m, j + kl)
                if (alpha == 1.0d0) then
                    do i = i1, i2
                        a(i+k,j) = a(i+k,j) * b(i)
                    end do
                else
                    do i = i1, i2
                        a(i+k,j) = alpha * a(i+k,j) * b(i)
                    end do
                end if
            end do
        end if


        ! End
        return

        ! KL < 0
10      continue
        call errmgr%report_error("band_diag_mtx_mult_dbl", &
            "The number of subdiagonals must be at least 0.", &
            LA_INVALID_INPUT_ERROR)
        return

        ! KU < 0
20      continue
        call errmgr%report_error("band_diag_mtx_mult_dbl", &
            "The number of superdiagonals must be at least 0.", &
            LA_INVALID_INPUT_ERROR)
        return

        ! B is not sized correctly
30      continue
        call errmgr%report_error("band_diag_mtx_mult_dbl", &
            "Inner matrix dimensions do not agree.", &
            LA_ARRAY_SIZE_ERROR)
        return
    end subroutine

! ------------------------------------------------------------------------------
    module subroutine band_diag_mtx_mult_cmplx(left, m, kl, ku, alpha, a, b, err)
        ! Arguments
        logical, intent(in) :: left
        integer(int32), intent(in) :: m, kl, ku
        complex(real64), intent(in) :: alpha
        complex(real64), intent(inout), dimension(:,:) :: a
        complex(real64), intent(in), dimension(:) :: b
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        complex(real64), parameter :: one = (1.0d0, 0.0d0)

        ! Local Variables
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        integer(int32) :: i, i1, i2, j, k, n
        complex(real64) :: temp
        
        ! Initialization
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if
        n = size(a, 2)

        ! Input Checking
        if (kl < 0) go to 10
        if (ku < 0) go to 20
        if (left) then
            if (size(b) /= n) go to 30
        else
            if (size(b) < m) go to 30
        end if

        ! Process
        if (left) then
            ! Compute A = A * B
            do j = 1, n
                k = ku + 1 - j
                i1 = max(1, j - ku) + k
                i2 = min(m, j + kl) + k
                if (alpha == one) then
                    temp = b(j)
                else
                    temp = alpha * b(j)
                end if
                do i = i1, i2
                    a(i,j) = a(i,j) * temp
                end do
            end do
        else
            ! Compute A = B * A
            do j = 1, n
                k = ku + 1 - j
                i1 = max(1, j - ku)
                i2 = min(m, j + kl)
                if (alpha == 1.0d0) then
                    do i = i1, i2
                        a(i+k,j) = a(i+k,j) * b(i)
                    end do
                else
                    do i = i1, i2
                        a(i+k,j) = alpha * a(i+k,j) * b(i)
                    end do
                end if
            end do
        end if


        ! End
        return

        ! KL < 0
10      continue
        call errmgr%report_error("band_diag_mtx_mult_cmplx", &
            "The number of subdiagonals must be at least 0.", &
            LA_INVALID_INPUT_ERROR)
        return

        ! KU < 0
20      continue
        call errmgr%report_error("band_diag_mtx_mult_cmplx", &
            "The number of superdiagonals must be at least 0.", &
            LA_INVALID_INPUT_ERROR)
        return

        ! B is not sized correctly
30      continue
        call errmgr%report_error("band_diag_mtx_mult_cmplx", &
            "Inner matrix dimensions do not agree.", &
            LA_ARRAY_SIZE_ERROR)
        return
    end subroutine

! ------------------------------------------------------------------------------
    module subroutine banded_to_dense_dbl(m, kl, ku, a, x, err)
        ! Arguments
        integer(int32), intent(in) :: m, kl, ku
        real(real64), intent(in), dimension(:,:) :: a
        real(real64), intent(out), dimension(:,:) :: x
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        real(real64), parameter :: zero = 0.0d0

        ! Local Variables
        integer(int32) :: i, j, k, n, i1, i2
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        
        ! Initialization
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if
        n = size(a, 2)

        ! Input Checking
        if (kl < 0 .or. ku < 0) then
            call errmgr%report_error("banded_to_dense_dbl", &
                "The bandwidth dimensions must not be negative-valued.", &
                LA_INVALID_INPUT_ERROR)
            return
        end if
        if (size(a, 1) /= kl + ku + 1) then
            call errmgr%report_error("banded_to_dense_dbl", "The size of " // &
                "the input matrix does not match the specified bandwidth.", &
                LA_MATRIX_FORMAT_ERROR)
            return
        end if
        if (size(x, 1) /= m .or. size(x, 2) /= n) then
            call errmgr%report_error("banded_to_dense_dbl", &
                "The output matrix dimensions are not correct.", &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Process
        do j = 1, n
            k = ku + 1 - j
            i1 = max(1, j - ku)
            i2 = min(m, j + kl)
            x(:i1-1,j) = zero
            do i = i1, i2
                x(i, j) = a(k + i, j)
            end do
            x(i2+1:,j) = zero
        end do
    end subroutine

! ------------------------------------------------------------------------------
    module subroutine banded_to_dense_cmplx(m, kl, ku, a, x, err)
        ! Arguments
        integer(int32), intent(in) :: m, kl, ku
        complex(real64), intent(in), dimension(:,:) :: a
        complex(real64), intent(out), dimension(:,:) :: x
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        complex(real64), parameter :: zero = (0.0d0, 0.0d0)

        ! Local Variables
        integer(int32) :: i, j, k, n, i1, i2
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        
        ! Initialization
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if
        n = size(a, 2)

        ! Input Checking
        if (kl < 0 .or. ku < 0) then
            call errmgr%report_error("banded_to_dense_cmplx", &
                "The bandwidth dimensions must not be negative-valued.", &
                LA_INVALID_INPUT_ERROR)
            return
        end if
        if (size(a, 1) /= kl + ku + 1) then
            call errmgr%report_error("banded_to_dense_cmplx", "The size of " // &
                "the input matrix does not match the specified bandwidth.", &
                LA_MATRIX_FORMAT_ERROR)
            return
        end if
        if (size(x, 1) /= m .or. size(x, 2) /= n) then
            call errmgr%report_error("banded_to_dense_cmplx", &
                "The output matrix dimensions are not correct.", &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Process
        do j = 1, n
            k = ku + 1 - j
            i1 = max(1, j - ku)
            i2 = min(m, j + kl)
            x(:i1-1,j) = zero
            do i = i1, i2
                x(i, j) = a(k + i, j)
            end do
            x(i2+1:,j) = zero
        end do
    end subroutine

! ------------------------------------------------------------------------------
    module subroutine dense_to_banded_dbl(a, kl, ku, x, err)
        ! Arguments
        real(real64), intent(in), dimension(:,:) :: a
        integer(int32), intent(in) :: kl, ku
        real(real64), intent(out), dimension(:,:) :: x
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        integer(int32) :: i, j, k, m, n, mm, flag
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        
        ! Initialization
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if
        m = size(a, 1)
        n = size(a, 2)
        mm = kl + ku + 1

        ! Input Check
        if (kl < 0 .or. ku < 0) then
            call errmgr%report_error("dense_to_banded_dbl", &
                "The bandwidth dimensions must not be negative-valued.", &
                LA_INVALID_INPUT_ERROR)
            return
        end if
        if (size(x, 1) /= mm .or. size(x, 2) /= n) then
            call errmgr%report_error("dense_to_banded_dbl", &
                "The output matrix dimensions are not correct.", &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Process
        do j = 1, n
            k = ku + 1 - j
            do i = max(1, j - ku), min(m, j + kl)
                x(k + i, j) = a(i,j)
            end do
        end do
    end subroutine

! ------------------------------------------------------------------------------
    module subroutine dense_to_banded_cmplx(a, kl, ku, x, err)
        ! Arguments
        complex(real64), intent(in), dimension(:,:) :: a
        integer(int32), intent(in) :: kl, ku
        complex(real64), intent(out), dimension(:,:) :: x
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        integer(int32) :: i, j, k, m, n, mm, flag
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        
        ! Initialization
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if
        m = size(a, 1)
        n = size(a, 2)
        mm = kl + ku + 1

        ! Input Check
        if (kl < 0 .or. ku < 0) then
            call errmgr%report_error("dense_to_banded_cmplx", &
                "The bandwidth dimensions must not be negative-valued.", &
                LA_INVALID_INPUT_ERROR)
            return
        end if
        if (size(x, 1) /= mm .or. size(x, 2) /= n) then
            call errmgr%report_error("dense_to_banded_cmplx", &
                "The output matrix dimensions are not correct.", &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Process
        do j = 1, n
            k = ku + 1 - j
            do i = max(1, j - ku), min(m, j + kl)
                x(k + i, j) = a(i,j)
            end do
        end do
    end subroutine

! ------------------------------------------------------------------------------
    module subroutine extract_diagonal_dbl(a, diag, err)
        ! Arguments
        real(real64), intent(in), dimension(:,:) :: a
        real(real64), intent(out), dimension(:) :: diag
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        integer(int32) :: i, m, n, mn
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        
        ! Initialization
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if
        m = size(a, 1)
        n = size(a, 2)
        mn = min(m, n)

        ! Input Checking
        if (size(diag) /= mn) then
            call errmgr%report_error("extract_diagonal_dbl", &
                "The is expected to have MIN(M, N) elements.", &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Process
        do i = 1, mn
            diag(i) = a(i,i)
        end do
    end subroutine

! ------------------------------------------------------------------------------
    module subroutine extract_diagonal_cmplx(a, diag, err)
        ! Arguments
        complex(real64), intent(in), dimension(:,:) :: a
        complex(real64), intent(out), dimension(:) :: diag
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        integer(int32) :: i, m, n, mn
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        
        ! Initialization
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if
        m = size(a, 1)
        n = size(a, 2)
        mn = min(m, n)

        ! Input Checking
        if (size(diag) /= mn) then
            call errmgr%report_error("extract_diagonal_cmplx", &
                "The is expected to have MIN(M, N) elements.", &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Process
        do i = 1, mn
            diag(i) = a(i,i)
        end do
    end subroutine

! ------------------------------------------------------------------------------
end submodule
