! linalg_c_api.f90

!> @brief Provides a C-friendly API to the LINALG library.  Notice, all C-API 
!! LINALG routines begin with the prefix "la_".
module linalg_c_api
    use iso_c_binding
    use linalg_core
    use linalg_constants
    use ferror
    implicit none

contains
! ------------------------------------------------------------------------------
    !> @brief Performs the rank-1 update to matrix A such that:
    !! A = alpha * X * Y**T + A, where A is an M-by-N matrix, alpha is a scalar,
    !! X is an M-element array, and N is an N-element array.
    !!
    !! @param[in] m The number of rows in the matrix.
    !! @param[in] n The number of columns in the matrix.
    !! @param[in] alpha The scalar multiplier.
    !! @param[in] x An M-element array.
    !! @param[in] y An N-element array.
    !! @param[in,out] a On input, the M-by-N matrix to update.  On output, the
    !!  updated M-by-N matrix.
    !! @param[in] lda The leading dimension of matrix A.
    !!
    !! @return An error code.  The following codes are possible.
    !!  - LA_NO_ERROR: No error occurred.  Successful operation.
    !!  - LA_INVALID_INPUT_ERROR: Occurs if @p lda is not correct.
    function la_rank1_update(m, n, alpha, x, y, a, lda) &
            bind(C, name = "la_rank1_update") result(flag)
        ! Arguments
        integer(c_int), intent(in), value :: m, n, lda
        real(c_double), intent(in), value :: alpha
        real(c_double), intent(in) :: x(*), y(*)
        real(c_double), intent(inout) :: a(lda,*)
        integer(c_int) :: flag

        ! Initialization
        flag = LA_NO_ERROR

        ! Input Checking
        if (lda < m) then
            flag = LA_INVALID_INPUT_ERROR
            return
        end if

        ! Process
        call rank1_update(alpha, x(1:m), y(1:n), a(1:m,1:n))
    end function

! ------------------------------------------------------------------------------
    !> @brief Performs the rank-1 update to matrix A such that:
    !! A = alpha * X * Y**T + A, where A is an M-by-N matrix, alpha is a scalar,
    !! X is an M-element array, and N is an N-element array.
    !!
    !! @param[in] m The number of rows in the matrix.
    !! @param[in] n The number of columns in the matrix.
    !! @param[in] alpha The scalar multiplier.
    !! @param[in] x An M-element array.
    !! @param[in] y An N-element array.
    !! @param[in,out] a On input, the M-by-N matrix to update.  On output, the
    !!  updated M-by-N matrix.
    !! @param[in] lda The leading dimension of matrix A.
    !!
    !! @return An error code.  The following codes are possible.
    !!  - LA_NO_ERROR: No error occurred.  Successful operation.
    !!  - LA_INVALID_INPUT_ERROR: Occurs if @p lda is not correct.
    function la_rank1_update_cmplx(m, n, alpha, x, y, a, lda) &
            bind(C, name = "la_rank1_update_cmplx") result(flag)
        ! Arguments
        integer(c_int), intent(in), value :: m, n, lda
        complex(c_double), intent(in), value :: alpha
        complex(c_double), intent(in) :: x(*), y(*)
        complex(c_double), intent(inout) :: a(lda,*)
        integer(c_int) :: flag

        ! Initialization
        flag = LA_NO_ERROR

        ! Input Checking
        if (lda < m) then
            flag = LA_INVALID_INPUT_ERROR
            return
        end if

        ! Process
        call rank1_update(alpha, x(1:m), y(1:n), a(1:m,1:n))
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the trace of a matrix (the sum of the main diagonal
    !! elements).
    !!
    !! @param[in] m The number of rows in the matrix.
    !! @param[in] n The number of columns in the matrix.
    !! @param[in] a The M-by-N matrix on which to operate.
    !! @param[in] lda The leading dimension of the matrix.
    !! @param[out] rst The results of the operation.
    !!
    !! @return An error code.  The following codes are possible.
    !!  - LA_NO_ERROR: No error occurred.  Successful operation.
    !!  - LA_INVALID_INPUT_ERROR: Occurs if @p lda is not correct.
    function la_trace(m, n, a, lda, rst) bind(C, name = "la_trace") &
            result(flag)
        ! Arguments
        integer(c_int), intent(in), value :: m, n, lda
        real(c_double), intent(in) :: a(lda,*)
        real(c_double), intent(out) :: rst
        integer(c_int) :: flag

        ! Initialization
        flag = LA_NO_ERROR

        ! Input Checking
        if (lda < m) then
            flag = LA_INVALID_INPUT_ERROR
            return
        end if

        ! Process
        rst = trace(a(1:m,1:n))
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the trace of a matrix (the sum of the main diagonal
    !! elements).
    !!
    !! @param[in] m The number of rows in the matrix.
    !! @param[in] n The number of columns in the matrix.
    !! @param[in] a The M-by-N matrix on which to operate.
    !! @param[in] lda The leading dimension of the matrix.
    !! @param[out] rst The results of the operation.
    !!
    !! @return An error code.  The following codes are possible.
    !!  - LA_NO_ERROR: No error occurred.  Successful operation.
    !!  - LA_INVALID_INPUT_ERROR: Occurs if @p lda is not correct.
    function la_trace_cmplx(m, n, a, lda, rst) &
            bind(C, name = "la_trace_cmplx") result(flag)
        ! Arguments
        integer(c_int), intent(in), value :: m, n, lda
        complex(c_double), intent(in) :: a(lda,*)
        complex(c_double), intent(out) :: rst
        integer(c_int) :: flag

        ! Initialization
        flag = LA_NO_ERROR

        ! Input Checking
        if (lda < m) then
            flag = LA_INVALID_INPUT_ERROR
            return
        end if

        ! Process
        rst = trace(a(1:m,1:n))
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the matrix operation C = alpha * op(A) * op(B) + beta * C.
    !!
    !! @param transa Set to true to compute op(A) as the transpose of A; else,
    !!  set to false to compute op(A) as A.
    !! @param transb Set to true to compute op(B) as the transpose of B; else,
    !!  set to false to compute op(B) as B.
    !! @param m The number of rows in @p c.
    !! @param n The number of columns in @p c.
    !! @param k The interior dimension of the product @p a and @p b.
    !! @param alpha A scalar multiplier.
    !! @param a If @p transa is true, this matrix must be @p k by @p m; else,
    !!  if @p transa is false, this matrix must be @p m by @p k.
    !! @param lda The leading dimension of matrix @p a.
    !! @param b If @p transb is true, this matrix must be @p n by @p k; else,
    !!  if @p transb is false, this matrix must be @p k by @p n.
    !! @param ldb The leading dimension of matrix @p b.
    !! @param beta A scalar multiplier.
    !! @param c The @p m by @p n matrix C.
    !! @param ldc The leading dimension of matrix @p c.
    !!
    !! @return An error code.  The following codes are possible.
    !!  - LA_NO_ERROR: No error occurred.  Successful operation.
    !!  - LA_INVALID_INPUT_ERROR: Occurs if @p lda, @p ldb, or @p ldc are not
    !!      correct.
    function la_mtx_mult(transa, transb, m, n, k, alpha, a, lda, b, ldb, &
            beta, c, ldc) bind(C, name="la_mtx_mult") result(flag)
        ! Arugments
        logical(c_bool), intent(in), value :: transa, transb
        integer(c_int), intent(in), value :: m, n, k, lda, ldb, ldc
        real(c_double), intent(in), value :: alpha, beta
        real(c_double), intent(in) :: a(lda,*), b(ldb,*)
        real(c_double), intent(inout) :: c(ldc,*)
        integer(c_int) :: flag

        ! Local Variables
        character :: ta, tb
        integer(c_int) :: nrowa, nrowb

        ! Initialization
        flag = LA_NO_ERROR
        ta = "N"
        if (transa) ta = "T"

        tb = "N"
        if (transb) tb = "T"

        if (transa) then
            nrowa = k
        else
            nrowa = m
        end if

        if (transb) then
            nrowb = n
        else
            nrowb = k
        end if

        ! Input Checking
        if (lda < nrowa .or. ldb < nrowb .or. ldc < m) then
            flag = LA_INVALID_INPUT_ERROR
            return
        end if

        ! Call DGEMM directly
        call DGEMM(ta, tb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the matrix operation C = alpha * op(A) * op(B) + beta * C.
    !!
    !! @param opa Set to TRANSPOSE to compute op(A) as a direct transpose of A,
    !!  set to HERMITIAN_TRANSPOSE to compute op(A) as the Hermitian transpose
    !!  of A, otherwise, set to NO_OPERATION to compute op(A) as A.
    !! @param opb Set to TRANSPOSE to compute op(B) as a direct transpose of B,
    !!  set to HERMITIAN_TRANSPOSE to compute op(B) as the Hermitian transpose
    !!  of B, otherwise, set to NO_OPERATION to compute op(B) as B.
    !! @param mThe number of rows in @p c.
    !! @param n The number of columns in @p c.
    !! @param k The interior dimension of the product @p a and @p b.
    !! @param alpha A scalar multiplier.
    !! @param a If @p opa is TRANSPOSE or HERMITIAN_TRANSPOSE, this matrix must
    !!  be @p k by @p m; else, this matrix must be @p m by @p k.
    !! @param lda The leading dimension of matrix @p a.
    !! @param b If @p opb is TRANSPOSE or HERMITIAN_TRANSPOSE, this matrix must
    !!  be @p n by @p k; else, this matrix must be @p k by @p n.
    !! @param ldb The leading dimension of matrix @p b.
    !! @param beta A scalar multiplier.
    !! @param c The @p m by @p n matrix C.
    !! @param ldc The leading dimension of matrix @p c.
    !!
    !! @return An error code.  The following codes are possible.
    !!  - LA_NO_ERROR: No error occurred.  Successful operation.
    !!  - LA_INVALID_INPUT_ERROR: Occurs if @p lda, @p ldb, or @p ldc are not
    !!      correct.
    function la_mtx_mult_cmplx(opa, opb, m, n, k, alpha, a, lda, b, ldb, &
            beta, c, ldc) bind(C, name="la_mtx_mult_cmplx") result(flag)
        ! Arguments
        integer(c_int), intent(in), value :: opa, opb, m, n, k, lda, ldb, ldc
        complex(c_double), intent(in), value :: alpha, beta
        complex(c_double), intent(in) :: a(lda,*), b(ldb,*)
        complex(c_double), intent(inout) :: c(ldc,*)
        integer(c_int) :: flag

        ! Local Variables
        character :: ta, tb
        integer(c_int) :: nrowa, nrowb

        ! Initialization
        flag = LA_NO_ERROR
        if (opa == TRANSPOSE) then
            ta = "T"
        else if (opa == HERMITIAN_TRANSPOSE) then
            ta = "H"
        else
            ta = "N"
        end if

        if (opb == TRANSPOSE) then
            tb = "T"
        else if (opb == HERMITIAN_TRANSPOSE) then
            tb = "H"
        else
            tb = "N"
        end if

        if (opa == TRANSPOSE .or. opa == HERMITIAN_TRANSPOSE) then
            nrowa = k
        else
            nrowa = m
        end if

        if (opb == TRANSPOSE .or. opb == HERMITIAN_TRANSPOSE) then
            nrowb = n
        else
            nrowb = k
        end if

        ! Input Checking
        if (lda < nrowa .or. ldb < nrowb .or. ldc < m) then
            flag = LA_INVALID_INPUT_ERROR
            return
        end if
        
        ! Call ZGEMM directly
        call ZGEMM(ta, tb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the matrix operation: C = alpha * A * op(B) + beta * C,
    !! or C = alpha * op(B) * A + beta * C.
    !!
    !! @param lside Set to true to apply matrix A from the left; else, set
    !!  to false to apply matrix A from the left.
    !! @param trans Set to true if op(B) == B**T; else, set to false if
    !!  op(B) == B.
    !! @param m The number of rows in the matrix C.
    !! @param n The number of columns in the matrix C.
    !! @param k The inner dimension of the matrix product A * op(B).
    !! @param alpha A scalar multiplier.
    !! @param a A P-element array containing the diagonal elements of matrix A
    !!  where P = MIN(@p m, @p k) if @p lside is true; else, P = MIN(@p n, @p k)
    !!  if @p lside is false.
    !! @param b The LDB-by-TDB matrix B where (LDB = leading dimension of B,
    !!  and TDB = trailing dimension of B):
    !!  - @p lside == true & @p trans == true: LDB = @p n, TDB = @p k
    !!  - @p lside == true & @p trans == false: LDB = @p k, TDB = @p n
    !!  - @p lside == false & @p trans == true: LDB = @p k, TDB = @p m
    !!  - @p lside == false & @p trans == false: LDB = @p m, TDB = @p k
    !! @param ldb The leading dimension of matrix B.
    !! @param beta A scalar multiplier.
    !! @param c The @p m by @p n matrix C.
    !! @param ldc The leading dimension of matrix C.
    !!
    !! @return An error code.  The following codes are possible.
    !!  - LA_NO_ERROR: No error occurred.  Successful operation.
    !!  - LA_INVALID_INPUT_ERROR: Occurs if @p ldb, or @p ldc are not
    !!      correct.
    !!  - LA_ARRAY_SIZE_ERROR: Occurs if any of the input array sizes are
    !!      incorrect.
    function la_diag_mtx_mult(lside, transb, m, n, k, alpha, a, b, ldb, &
            beta, c, ldc) bind(C, name="la_diag_mtx_mult") result(flag)
        ! Arguments
        logical(c_bool), intent(in), value :: lside, transb
        integer(c_int), intent(in), value :: m, n, k, ldb, ldc
        real(c_double), intent(in), value :: alpha, beta
        real(c_double), intent(in) :: a(*), b(ldb,*)
        real(c_double), intent(inout) :: c(ldc,*)
        integer(c_int) :: flag

        ! Local Variabes
        integer(c_int) :: nrows, ncols, p
        logical :: ls, tb
        type(errors) :: err

        ! Initialization
        call err%set_exit_on_error(.false.)
        flag = LA_NO_ERROR
        if (lside .and. transb) then
            nrows = n
            ncols = k
            p = min(k, m)
            ls = .true.
            tb = .true.
        else if (lside .and. .not. transb) then
            nrows = k
            ncols = n
            p = min(k, m)
            ls = .true.
            tb = .false.
        else if (.not. lside .and. transb) then
            nrows = k
            ncols = m
            p = min(k, n)
            ls = .false.
            tb = .true.
        else
            nrows = m
            ncols = k
            p = min(k, n)
            ls = .false.
            tb = .false.
        end if

        ! Error Checking
        if (ldb < nrows .or. ldc < m) then
            flag = LA_INVALID_INPUT_ERROR
            return
        end if

        ! Process
        call diag_mtx_mult(ls, tb, alpha, a(1:p), b(1:nrows,1:ncols), &
            beta, c(1:m,1:n), err)
        if (err%has_error_occurred()) flag = err%get_error_flag()
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the matrix operation: C = alpha * A * op(B) + beta * C,
    !! or C = alpha * op(B) * A + beta * C.
    !!
    !! @param lside Set to true to apply matrix A from the left; else, set
    !!  to false to apply matrix A from the left.
    !! @param opb Set to TRANSPOSE to compute op(B) as a direct transpose of B,
    !!  set to HERMITIAN_TRANSPOSE to compute op(B) as the Hermitian transpose
    !!  of B, otherwise, set to NO_OPERATION to compute op(B) as B.
    !! @param m The number of rows in the matrix C.
    !! @param n The number of columns in the matrix C.
    !! @param k The inner dimension of the matrix product A * op(B).
    !! @param alpha A scalar multiplier.
    !! @param a A P-element array containing the diagonal elements of matrix A
    !!  where P = MIN(@p m, @p k) if @p lside is true; else, P = MIN(@p n, @p k)
    !!  if @p lside is false.
    !! @param b The LDB-by-TDB matrix B where (LDB = leading dimension of B,
    !!  and TDB = trailing dimension of B):
    !!  - @p lside == true & @p trans == true: LDB = @p n, TDB = @p k
    !!  - @p lside == true & @p trans == false: LDB = @p k, TDB = @p n
    !!  - @p lside == false & @p trans == true: LDB = @p k, TDB = @p m
    !!  - @p lside == false & @p trans == false: LDB = @p m, TDB = @p k
    !! @param ldb The leading dimension of matrix B.
    !! @param beta A scalar multiplier.
    !! @param c The @p m by @p n matrix C.
    !! @param ldc The leading dimension of matrix C.
    !!
    !! @return An error code.  The following codes are possible.
    !!  - LA_NO_ERROR: No error occurred.  Successful operation.
    !!  - LA_INVALID_INPUT_ERROR: Occurs if @p ldb, or @p ldc are not
    !!      correct.
    !!  - LA_ARRAY_SIZE_ERROR: Occurs if any of the input array sizes are
    !!      incorrect.
    function la_diag_mtx_mult_cmplx(lside, opb, m, n, k, alpha, a, b, &
            ldb, beta, c, ldc) bind(C, name="la_diag_mtx_mult_cmplx") &
            result(flag)
        ! Arguments
        logical(c_bool), intent(in), value :: lside
        integer(c_int), intent(in), value :: opb, m, n, k, ldb, ldc
        complex(c_double), intent(in), value :: alpha, beta
        complex(c_double), intent(in) :: a(*), b(ldb,*)
        complex(c_double), intent(inout) :: c(ldc,*)
        integer(c_int) :: flag

        ! Local Variabes
        integer(c_int) :: nrows, ncols, p
        logical :: ls, tb
        type(errors) :: err

        ! Initialization
        call err%set_exit_on_error(.false.)
        flag = LA_NO_ERROR
        tb = .false.
        if (opb == TRANSPOSE .or. opb == HERMITIAN_TRANSPOSE) tb = .true.
        if (lside .and. tb) then
            nrows = n
            ncols = k
            p = min(k, m)
            ls = .true.
        else if (lside .and. .not. tb) then
            nrows = k
            ncols = n
            p = min(k, m)
            ls = .true.
        else if (.not. lside .and. tb) then
            nrows = k
            ncols = m
            p = min(k, n)
            ls = .false.
        else
            nrows = m
            ncols = k
            p = min(k, n)
            ls = .false.
        end if

        ! Error Checking
        if (ldb < nrows .or. ldc < m) then
            flag = LA_INVALID_INPUT_ERROR
            return
        end if

        ! Process
        call diag_mtx_mult(ls, opb, alpha, a(1:p), b(1:nrows,1:ncols), &
            beta, c(1:m,1:n))
        if (err%has_error_occurred()) flag = err%get_error_flag()
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the rank of a matrix.
    !!
    !! @param m The number of rows in the matrix.
    !! @param n The number of columns in the matrix.
    !! @param a The M-by-N matrix.  The matrix is overwritten as part of this
    !!  operation.
    !! @param lda The leading dimension of matrix A.
    !! @param[out] rnk The rank of @p a.
    !!
    !! @return An error code.  The following codes are possible.
    !!  - LA_NO_ERROR: No error occurred.  Successful operation.
    !!  - LA_INVALID_INPUT_ERROR: Occurs if @p lda is not correct.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
    !!      there is insufficient memory available.
    !!  - LA_CONVERGENCE_ERROR: Occurs as a warning if the QR iteration process
    !!      could not converge to a zero value.
    function la_rank(m, n, a, lda, rnk) bind(C, name="la_rank") result(flag)
        ! Arguments
        integer(c_int), intent(in), value :: m, n, lda
        real(c_double), intent(inout) :: a(lda,*)
        integer(c_int), intent(out) :: rnk
        integer(c_int) :: flag

        ! Local Variables
        type(errors) :: err

        ! Input Check
        call err%set_exit_on_error(.false.)
        flag = LA_NO_ERROR
        if (lda < m) then
            flag = LA_INVALID_INPUT_ERROR
            return
        end if

        ! Process
        rnk = mtx_rank(a(1:m,1:n), err =err)
        if (err%has_error_occurred()) flag = err%get_error_flag()
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the rank of a matrix.
    !!
    !! @param m The number of rows in the matrix.
    !! @param n The number of columns in the matrix.
    !! @param a The M-by-N matrix.  The matrix is overwritten as part of this
    !!  operation.
    !! @param lda The leading dimension of matrix A.
    !! @param[out] rnk The rank of @p a.
    !!
    !! @return An error code.  The following codes are possible.
    !!  - LA_NO_ERROR: No error occurred.  Successful operation.
    !!  - LA_INVALID_INPUT_ERROR: Occurs if @p lda is not correct.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
    !!      there is insufficient memory available.
    !!  - LA_CONVERGENCE_ERROR: Occurs as a warning if the QR iteration process
    !!      could not converge to a zero value.
    function la_rank_cmplx(m, n, a, lda, rnk) bind(C, name="la_rank_cmplx") &
            result(flag)
        ! Arguments
        integer(c_int), intent(in), value :: m, n, lda
        complex(c_double), intent(inout) :: a(lda,*)
        integer(c_int), intent(out) :: rnk
        integer(c_int) :: flag

        ! Local Variables
        type(errors) :: err

        ! Input Check
        call err%set_exit_on_error(.false.)
        flag = LA_NO_ERROR
        if (lda < m) then
            flag = LA_INVALID_INPUT_ERROR
            return
        end if

        ! Process
        rnk = mtx_rank(a(1:m,1:n), err = err)
        if (err%has_error_occurred()) flag = err%get_error_flag()
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the determinant of a square matrix.
    !!
    !! @param n The dimension of the matrix.
    !! @param a The N-by-N matrix.  The matrix is overwritten on output.
    !! @param lda The leading dimension of the matrix.
    !! @param[out] d The determinant of @p a.
    !!
    !! @return An error code.  The following codes are possible.
    !!  - LA_NO_ERROR: No error occurred.  Successful operation.
    !!  - LA_ARRAY_SIZE_ERROR: Occurs if any of the input arrays are not sized
    !!      appropriately.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory 
    !!      available.
    function la_det(n, a, lda, d) bind(C, name="la_det") result(flag)
        ! Arguments
        integer(c_int), intent(in), value :: n, lda
        real(c_double), intent(inout) :: a(lda,*)
        real(c_double), intent(out) :: d
        integer(c_int) :: flag

        ! Local Variables
        type(errors) :: err

        ! Error Checking
        call err%set_exit_on_error(.false.)
        flag = LA_NO_ERROR
        if (lda < n) then
            flag = LA_INVALID_INPUT_ERROR
            return
        end if

        ! Process
        d = det(a(1:n,1:n), err = err)
        if (err%has_error_occurred()) flag = err%get_error_flag()
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the determinant of a square matrix.
    !!
    !! @param n The dimension of the matrix.
    !! @param a The N-by-N matrix.  The matrix is overwritten on output.
    !! @param lda The leading dimension of the matrix.
    !! @param[out] d The determinant of @p a.
    !!
    !! @return An error code.  The following codes are possible.
    !!  - LA_NO_ERROR: No error occurred.  Successful operation.
    !!  - LA_ARRAY_SIZE_ERROR: Occurs if any of the input arrays are not sized
    !!      appropriately.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory 
    !!      available.
    function la_det_cmplx(n, a, lda, d) bind(C, name="la_det_cmplx") result(flag)
        ! Arguments
        integer(c_int), intent(in), value :: n, lda
        complex(c_double), intent(inout) :: a(lda,*)
        complex(c_double), intent(out) :: d
        integer(c_int) :: flag

        ! Local Variables
        type(errors) :: err

        ! Error Checking
        call err%set_exit_on_error(.false.)
        flag = LA_NO_ERROR
        if (lda < n) then
            flag = LA_INVALID_INPUT_ERROR
            return
        end if

        ! Process
        d = det(a(1:n,1:n), err = err)
        if (err%has_error_occurred()) flag = err%get_error_flag()
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the triangular matrix operation:
    !! B = alpha * A**T * A + beta * B, or B = alpha * A * A**T + beta * B,
    !! where A is a triangular matrix.
    !!
    !! @param upper Set to true if matrix A is upper triangular, and
    !!  B = alpha * A**T * A + beta * B is to be calculated; else, set to false
    !!  if A is lower triangular, and B = alpha * A * A**T + beta * B is to
    !!  be computed.
    !! @param alpha A scalar multiplier.
    !! @param n The dimension of the matrix.
    !! @param a The @p n by @p n triangular matrix A.  Notice, if @p upper is
    !!  true, only the upper triangular portion of this matrix is referenced;
    !!  else, if @p upper is false, only the lower triangular portion of this
    !!  matrix is referenced.
    !! @param lda The leading dimension of matrix A.
    !! @param beta A scalar multiplier.
    !! @param b On input, the @p n by @p n matrix B.  On output, the @p n by
    !!  @p n resulting matrix.
    !! @param ldb The leading dimension of matrix B.
    !!
    !! @return An error code.  The following codes are possible.
    !!  - LA_NO_ERROR: No error occurred.  Successful operation.
    !!  - LA_INVALID_INPUT_ERROR: Occurs if @p lda or @p ldb are not correct.
    function la_tri_mtx_mult(upper, alpha, n, a, lda, beta, b, ldb) &
            bind(C, name = "la_tri_mtx_mult") result(flag)
        ! Arguments
        logical(c_bool), intent(in), value :: upper
        integer(c_int), intent(in), value :: n, lda, ldb
        real(c_double), intent(in), value :: alpha, beta
        real(c_double), intent(in) :: a(lda,*)
        real(c_double), intent(inout) :: b(ldb,*)
        integer(c_int) :: flag

        ! Error Checking
        flag = LA_NO_ERROR
        if (lda < n .or. ldb < n) then
            flag = LA_INVALID_INPUT_ERROR
            return
        end if

        ! Process
        call tri_mtx_mult(logical(upper), alpha, a(1:n,1:n), beta, b(1:n,1:n))
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the triangular matrix operation:
    !! B = alpha * A**T * A + beta * B, or B = alpha * A * A**T + beta * B,
    !! where A is a triangular matrix.
    !!
    !! @param upper Set to true if matrix A is upper triangular, and
    !!  B = alpha * A**T * A + beta * B is to be calculated; else, set to false
    !!  if A is lower triangular, and B = alpha * A * A**T + beta * B is to
    !!  be computed.
    !! @param alpha A scalar multiplier.
    !! @param n The dimension of the matrix.
    !! @param a The @p n by @p n triangular matrix A.  Notice, if @p upper is
    !!  true, only the upper triangular portion of this matrix is referenced;
    !!  else, if @p upper is false, only the lower triangular portion of this
    !!  matrix is referenced.
    !! @param lda The leading dimension of matrix A.
    !! @param beta A scalar multiplier.
    !! @param b On input, the @p n by @p n matrix B.  On output, the @p n by
    !!  @p n resulting matrix.
    !! @param ldb The leading dimension of matrix B.
    !!
    !! @return An error code.  The following codes are possible.
    !!  - LA_NO_ERROR: No error occurred.  Successful operation.
    !!  - LA_INVALID_INPUT_ERROR: Occurs if @p lda or @p ldb are not correct.
    function la_tri_mtx_mult_cmplx(upper, alpha, n, a, lda, beta, b, ldb) &
            bind(C, name = "la_tri_mtx_mult_cmplx") result(flag)
        ! Arguments
        logical(c_bool), intent(in), value :: upper
        integer(c_int), intent(in), value :: n, lda, ldb
        complex(c_double), intent(in), value :: alpha, beta
        complex(c_double), intent(in) :: a(lda,*)
        complex(c_double), intent(inout) :: b(ldb,*)
        integer(c_int) :: flag

        ! Error Checking
        flag = LA_NO_ERROR
        if (lda < n .or. ldb < n) then
            flag = LA_INVALID_INPUT_ERROR
            return
        end if

        ! Process
        call tri_mtx_mult(logical(upper), alpha, a(1:n,1:n), beta, b(1:n,1:n))
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the LU factorization of an M-by-N matrix.
    !!
    !! @param m The number of rows in the matrix.
    !! @param n The number of columns in the matrix.
    !! @param[in,out] a On input, the M-by-N matrix on which to operate.  On
    !! output, the LU factored matrix in the form [L\\U] where the unit diagonal
    !! elements of L are not stored.
    !! @param lda The leading dimension of matrix A.
    !! @param[out] ipvt An MIN(M, N)-element array used to track row-pivot
    !!  operations.  The array stored pivot information such that row I is
    !!  interchanged with row IPVT(I).
    !!
    !! @return An error code.  The following codes are possible.
    !!  - LA_NO_ERROR: No error occurred.  Successful operation.
    !!  - LA_INVALID_INPUT_ERROR: Occurs if @p lda is not correct.
    !!  - LA_SINGULAR_MATRIX_ERROR: Occurs as a warning if @p a is found to be
    !!      singular.
    function la_lu_factor(m, n, a, lda, ipvt) bind(C, name = "la_lu_factor") &
            result(flag)
        ! Arguments
        integer(c_int), intent(in), value :: m, n, lda
        real(c_double), intent(inout) :: a(lda,*)
        integer(c_int), intent(out) :: ipvt(*)
        integer(c_int) :: flag

        ! Local Variables
        type(errors) :: err
        integer(c_int) :: mn

        ! Error Checking
        call err%set_exit_on_error(.false.)
        flag = LA_NO_ERROR
        if (lda < m) then
            flag = LA_INVALID_INPUT_ERROR
            return
        end if

        ! Process
        mn = min(m, n)
        call lu_factor(a(1:m,1:n), ipvt(1:mn), err)
        if (err%has_error_occurred()) then
            flag = err%get_error_flag()
            return
        end if
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the LU factorization of an M-by-N matrix.
    !!
    !! @param m The number of rows in the matrix.
    !! @param n The number of columns in the matrix.
    !! @param[in,out] a On input, the M-by-N matrix on which to operate.  On
    !! output, the LU factored matrix in the form [L\\U] where the unit diagonal
    !! elements of L are not stored.
    !! @param lda The leading dimension of matrix A.
    !! @param[out] ipvt An MIN(M, N)-element array used to track row-pivot
    !!  operations.  The array stored pivot information such that row I is
    !!  interchanged with row IPVT(I).
    !!
    !! @return An error code.  The following codes are possible.
    !!  - LA_NO_ERROR: No error occurred.  Successful operation.
    !!  - LA_INVALID_INPUT_ERROR: Occurs if @p lda is not correct.
    !!  - LA_SINGULAR_MATRIX_ERROR: Occurs as a warning if @p a is found to be
    !!      singular.
    function la_lu_factor_cmplx(m, n, a, lda, ipvt) &
            bind(C, name = "la_lu_factor_cmplx") result(flag)
        ! Arguments
        integer(c_int), intent(in), value :: m, n, lda
        complex(c_double), intent(inout) :: a(lda,*)
        integer(c_int), intent(out) :: ipvt(*)
        integer(c_int) :: flag

        ! Local Variables
        type(errors) :: err
        integer(c_int) :: mn

        ! Error Checking
        call err%set_exit_on_error(.false.)
        flag = LA_NO_ERROR
        if (lda < m) then
            flag = LA_INVALID_INPUT_ERROR
            return
        end if

        ! Process
        mn = min(m, n)
        call lu_factor(a(1:m,1:n), ipvt(1:mn), err)
        if (err%has_error_occurred()) then
            flag = err%get_error_flag()
            return
        end if
    end function

! ------------------------------------------------------------------------------
    !> @brief Extracts the L, U, and P matrices from the LU factorization
    !! output from la_lu_factor.
    !!
    !! @param n The dimension of the input matrix.
    !! @param[in,out] a On input, the N-by-N matrix as output by
    !!  @ref la_lu_factor.  On output, the N-by-N lower triangular matrix L.
    !! @param lda The leading dimension of @p a.
    !! @param ipvt The N-element pivot array as output by
    !!  @ref la_lu_factor.
    !! @param[out] u An N-by-N matrix where the U matrix will be written.
    !! @param ldu The leading dimension of @p u.
    !! @param[out] p An N-by-N matrix where the row permutation matrix will be
    !!  written.
    !! @param ldp The leading dimension of @p p.
    !!
    !! @return An error code.  The following codes are possible.
    !!  - LA_NO_ERROR: No error occurred.  Successful operation.
    !!  - LA_INVALID_INPUT_ERROR: Occurs if @p lda, @p ldu, or @p ldp is not 
    !!      correct.
    function la_form_lu(n, a, lda, ipvt, u, ldu, p, ldp) &
            bind(C, name = "la_form_lu") result(flag)
        ! Arguments
        integer(c_int), intent(in), value :: n, lda, ldu, ldp
        real(c_double), intent(inout) :: a(lda,*)
        real(c_double), intent(out) :: u(ldu,*), p(ldp,*)
        integer(c_int), intent(in) :: ipvt(*)
        integer(c_int) :: flag

        ! Input Checking
        flag = LA_NO_ERROR
        if (lda < n .or. ldu < n .or. ldp < n) then
            flag = LA_INVALID_INPUT_ERROR
            return
        end if

        ! Process
        call form_lu(a(1:n,1:n), ipvt(1:n), u(1:n,1:n), p(1:n,1:n))
    end function

! ------------------------------------------------------------------------------
    !> @brief Extracts the L, U, and P matrices from the LU factorization
    !! output from la_lu_factor.
    !!
    !! @param n The dimension of the input matrix.
    !! @param[in,out] a On input, the N-by-N matrix as output by
    !!  @ref la_lu_factor.  On output, the N-by-N lower triangular matrix L.
    !! @param lda The leading dimension of @p a.
    !! @param ipvt The N-element pivot array as output by
    !!  @ref la_lu_factor.
    !! @param[out] u An N-by-N matrix where the U matrix will be written.
    !! @param ldu The leading dimension of @p u.
    !! @param[out] p An N-by-N matrix where the row permutation matrix will be
    !!  written.
    !! @param ldp The leading dimension of @p p.
    !!
    !! @return An error code.  The following codes are possible.
    !!  - LA_NO_ERROR: No error occurred.  Successful operation.
    !!  - LA_INVALID_INPUT_ERROR: Occurs if @p lda, @p ldu, or @p ldp is not 
    !!      correct.
    function la_form_lu_cmplx(n, a, lda, ipvt, u, ldu, p, ldp) &
            bind(C, name = "la_form_lu_cmplx") result(flag)
        ! Arguments
        integer(c_int), intent(in), value :: n, lda, ldu, ldp
        complex(c_double), intent(inout) :: a(lda,*)
        complex(c_double), intent(out) :: u(ldu,*)
        real(c_double), intent(out) :: p(ldp,*)
        integer(c_int), intent(in) :: ipvt(*)
        integer(c_int) :: flag

        ! Input Checking
        flag = LA_NO_ERROR
        if (lda < n .or. ldu < n .or. ldp < n) then
            flag = LA_INVALID_INPUT_ERROR
            return
        end if

        ! Process
        call form_lu(a(1:n,1:n), ipvt(1:n), u(1:n,1:n), p(1:n,1:n))
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the QR factorization of an M-by-N matrix without
    !! pivoting.
    !!
    !! @param m The number of rows in the matrix.
    !! @param n The number of columns in the matrix.
    !! @param[in,out] a  On input, the M-by-N matrix to factor.  On output, the
    !!  elements on and above the diagonal contain the MIN(M, N)-by-N upper
    !!  trapezoidal matrix R (R is upper triangular if M >= N).  The elements
    !!  below the diagonal, along with the array @p tau, represent the
    !!  orthogonal matrix Q as a product of elementary reflectors.
    !! @param lda The leading dimension of matrix A.
    !! @param[out] tau A MIN(M, N)-element array used to store the scalar
    !!  factors of the elementary reflectors.
    !!
    !! @return An error code.  The following codes are possible.
    !!  - LA_NO_ERROR: No error occurred.  Successful operation.
    !!  - LA_INVALID_INPUT_ERROR: Occurs if @p lda is not correct.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory 
    !!      available.
    function la_qr_factor(m, n, a, lda, tau) bind(C, name = "la_qr_factor") &
            result(flag)
        ! Arguments
        integer(c_int), intent(in), value :: m, n, lda
        real(c_double), intent(inout) :: a(lda,*)
        real(c_double), intent(out) :: tau(*)
        integer(c_int) :: flag

        ! Local Variables
        type(errors) :: err
        integer(c_int) :: mn

        ! Error Checking
        call err%set_exit_on_error(.false.)
        flag = LA_NO_ERROR
        if (lda < m) then
            flag = LA_INVALID_INPUT_ERROR
            return
        end if

        ! Process
        mn = min(m, n)
        call qr_factor(a(1:m,1:n), tau(1:mn), err = err)
        if (err%has_error_occurred()) then
            flag = err%get_error_flag()
            return
        end if
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the QR factorization of an M-by-N matrix without
    !! pivoting.
    !!
    !! @param m The number of rows in the matrix.
    !! @param n The number of columns in the matrix.
    !! @param[in,out] a  On input, the M-by-N matrix to factor.  On output, the
    !!  elements on and above the diagonal contain the MIN(M, N)-by-N upper
    !!  trapezoidal matrix R (R is upper triangular if M >= N).  The elements
    !!  below the diagonal, along with the array @p tau, represent the
    !!  orthogonal matrix Q as a product of elementary reflectors.
    !! @param lda The leading dimension of matrix A.
    !! @param[out] tau A MIN(M, N)-element array used to store the scalar
    !!  factors of the elementary reflectors.
    !!
    !! @return An error code.  The following codes are possible.
    !!  - LA_NO_ERROR: No error occurred.  Successful operation.
    !!  - LA_INVALID_INPUT_ERROR: Occurs if @p lda is not correct.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory 
    !!      available.
    function la_qr_factor_cmplx(m, n, a, lda, tau) &
            bind(C, name = "la_qr_factor_cmplx") result(flag)
        ! Arguments
        integer(c_int), intent(in), value :: m, n, lda
        complex(c_double), intent(inout) :: a(lda,*)
        complex(c_double), intent(out) :: tau(*)
        integer(c_int) :: flag

        ! Local Variables
        type(errors) :: err
        integer(c_int) :: mn

        ! Error Checking
        call err%set_exit_on_error(.false.)
        flag = LA_NO_ERROR
        if (lda < m) then
            flag = LA_INVALID_INPUT_ERROR
            return
        end if

        ! Process
        mn = min(m, n)
        call qr_factor(a(1:m,1:n), tau(1:mn), err = err)
        if (err%has_error_occurred()) then
            flag = err%get_error_flag()
            return
        end if
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the QR factorization of an M-by-N matrix with column
    !! pivoting.
    !!
    !! @param m The number of rows in the matrix.
    !! @param n The number of columns in the matrix.
    !! @param[in,out] a  On input, the M-by-N matrix to factor.  On output, the
    !!  elements on and above the diagonal contain the MIN(M, N)-by-N upper
    !!  trapezoidal matrix R (R is upper triangular if M >= N).  The elements
    !!  below the diagonal, along with the array @p tau, represent the
    !!  orthogonal matrix Q as a product of elementary reflectors.
    !! @param lda The leading dimension of matrix A.
    !! @param[out] tau A MIN(M, N)-element array used to store the scalar
    !!  factors of the elementary reflectors.
    !! @param[in,out] jpvt On input, an N-element array that if JPVT(I) .ne. 0,
    !!  the I-th column of A is permuted to the front of A * P; if JPVT(I) = 0,
    !!  the I-th column of A is a free column.  On output, if JPVT(I) = K, then
    !!  the I-th column of A * P was the K-th column of A.
    !!
    !! @return An error code.  The following codes are possible.
    !!  - LA_NO_ERROR: No error occurred.  Successful operation.
    !!  - LA_INVALID_INPUT_ERROR: Occurs if @p lda is not correct.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory 
    !!      available.
    function la_qr_factor_pvt(m, n, a, lda, tau, jpvt) &
            bind(C, name = "la_qr_factor_pvt") result(flag)
        ! Arguments
        integer(c_int), intent(in), value :: m, n, lda
        real(c_double), intent(inout) :: a(lda,*)
        real(c_double), intent(out) :: tau(*)
        integer(c_int), intent(inout) :: jpvt(*)
        integer(c_int) :: flag

        ! Local Variables
        type(errors) :: err
        integer(c_int) :: mn

        ! Error Checking
        call err%set_exit_on_error(.false.)
        flag = LA_NO_ERROR
        if (lda < m) then
            flag = LA_INVALID_INPUT_ERROR
            return
        end if

        ! Process
        mn = min(m, n)
        call qr_factor(a(1:m,1:n), tau(1:mn), jpvt(1:n), err = err)
        if (err%has_error_occurred()) then
            flag = err%get_error_flag()
            return
        end if
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the QR factorization of an M-by-N matrix with column
    !! pivoting.
    !!
    !! @param m The number of rows in the matrix.
    !! @param n The number of columns in the matrix.
    !! @param[in,out] a  On input, the M-by-N matrix to factor.  On output, the
    !!  elements on and above the diagonal contain the MIN(M, N)-by-N upper
    !!  trapezoidal matrix R (R is upper triangular if M >= N).  The elements
    !!  below the diagonal, along with the array @p tau, represent the
    !!  orthogonal matrix Q as a product of elementary reflectors.
    !! @param lda The leading dimension of matrix A.
    !! @param[out] tau A MIN(M, N)-element array used to store the scalar
    !!  factors of the elementary reflectors.
    !! @param[in,out] jpvt On input, an N-element array that if JPVT(I) .ne. 0,
    !!  the I-th column of A is permuted to the front of A * P; if JPVT(I) = 0,
    !!  the I-th column of A is a free column.  On output, if JPVT(I) = K, then
    !!  the I-th column of A * P was the K-th column of A.
    !!
    !! @return An error code.  The following codes are possible.
    !!  - LA_NO_ERROR: No error occurred.  Successful operation.
    !!  - LA_INVALID_INPUT_ERROR: Occurs if @p lda is not correct.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory 
    !!      available.
    function la_qr_factor_cmplx_pvt(m, n, a, lda, tau, jpvt) &
            bind(C, name = "la_qr_factor_cmplx_pvt") result(flag)
        ! Arguments
        integer(c_int), intent(in), value :: m, n, lda
        complex(c_double), intent(inout) :: a(lda,*)
        complex(c_double), intent(out) :: tau(*)
        integer(c_int), intent(inout) :: jpvt(*)
        integer(c_int) :: flag

        ! Local Variables
        type(errors) :: err
        integer(c_int) :: mn

        ! Error Checking
        call err%set_exit_on_error(.false.)
        flag = LA_NO_ERROR
        if (lda < m) then
            flag = LA_INVALID_INPUT_ERROR
            return
        end if

        ! Process
        mn = min(m, n)
        call qr_factor(a(1:m,1:n), tau(1:mn), jpvt(1:n), err = err)
        if (err%has_error_occurred()) then
            flag = err%get_error_flag()
            return
        end if
    end function

! ------------------------------------------------------------------------------
    !> @brief Forms the full M-by-M orthogonal matrix Q from the elementary
    !! reflectors returned by the base QR factorization algorithm.
    !!
    !! @param[in] fullq Set to true to always return the full Q matrix; else,
    !!  set to false, and in the event that M > N, Q may be supplied as M-by-N,
    !!  and therefore only return the useful submatrix Q1 (Q = [Q1, Q2]) as the
    !!  factorization can be written as Q * R = [Q1, Q2] * [R1; 0].
    !! @param[in] m The number of rows in R.
    !! @param[in] n The number of columns in R.
    !! @param[in,out] r On input, the M-by-N factored matrix as returned by the
    !!  QR factorization routine.  On output, the upper triangular matrix R.
    !! @param[in] ldr The leading dimension of matrix R.
    !! @param[in] tau A MIN(M, N)-element array containing the scalar factors of
    !!  each elementary reflector defined in @p r.
    !! @param[out] q An M-by-M matrix where the full Q matrix will be written.
    !!  In the event that @p fullq is set to false, and M > N, this matrix need
    !!  only by M-by-N.
    !! @param[in] ldq The leading dimension of matrix Q.
    !!
    !! @return An error code.  The following codes are possible.
    !!  - LA_NO_ERROR: No error occurred.  Successful operation.
    !!  - LA_INVALID_INPUT_ERROR: Occurs if @p lda is not correct.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory 
    !!      available.
    function la_form_qr(fullq, m, n, r, ldr, tau, q, ldq) &
            bind(C, name = "la_form_qr") result(flag)
        ! Arguments
        logical(c_bool), intent(in), value :: fullq
        integer(c_int), intent(in), value :: m, n, ldr, ldq
        real(c_double), intent(inout) :: r(ldr,*)
        real(c_double), intent(in) :: tau(*)
        real(c_double), intent(out) :: q(ldq,*)
        integer(c_int) :: flag

        ! Local Variables
        type(errors) :: err
        integer(c_int) :: mn, nq

        ! Error Checking
        call err%set_exit_on_error(.false.)
        flag = LA_NO_ERROR
        if (ldr < m .or. ldq < m) then
            flag = LA_INVALID_INPUT_ERROR
            return
        end if

        ! Process
        mn = min(m, n)
        nq = m
        if (.not.fullq) nq = n
        call form_qr(r(1:m,1:n), tau(1:mn), q(1:m,1:nq), err = err)
        if (err%has_error_occurred()) then
            flag = err%get_error_flag()
            return
        end if
    end function

! ------------------------------------------------------------------------------
    !> @brief Forms the full M-by-M orthogonal matrix Q from the elementary
    !! reflectors returned by the base QR factorization algorithm.
    !!
    !! @param[in] fullq Set to true to always return the full Q matrix; else,
    !!  set to false, and in the event that M > N, Q may be supplied as M-by-N,
    !!  and therefore only return the useful submatrix Q1 (Q = [Q1, Q2]) as the
    !!  factorization can be written as Q * R = [Q1, Q2] * [R1; 0].
    !! @param[in] m The number of rows in R.
    !! @param[in] n The number of columns in R.
    !! @param[in,out] r On input, the M-by-N factored matrix as returned by the
    !!  QR factorization routine.  On output, the upper triangular matrix R.
    !! @param[in] ldr The leading dimension of matrix R.
    !! @param[in] tau A MIN(M, N)-element array containing the scalar factors of
    !!  each elementary reflector defined in @p r.
    !! @param[out] q An M-by-M matrix where the full Q matrix will be written.
    !!  In the event that @p fullq is set to false, and M > N, this matrix need
    !!  only by M-by-N.
    !! @param[in] ldq The leading dimension of matrix Q.
    !!
    !! @return An error code.  The following codes are possible.
    !!  - LA_NO_ERROR: No error occurred.  Successful operation.
    !!  - LA_INVALID_INPUT_ERROR: Occurs if @p lda is not correct.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory 
    !!      available.
    function la_form_qr_cmplx(fullq, m, n, r, ldr, tau, q, ldq) &
            bind(C, name = "la_form_qr_cmplx") result(flag)
        ! Arguments
        logical(c_bool), intent(in), value :: fullq
        integer(c_int), intent(in), value :: m, n, ldr, ldq
        complex(c_double), intent(inout) :: r(ldr,*)
        complex(c_double), intent(in) :: tau(*)
        complex(c_double), intent(out) :: q(ldq,*)
        integer(c_int) :: flag

        ! Local Variables
        type(errors) :: err
        integer(c_int) :: mn, nq

        ! Error Checking
        call err%set_exit_on_error(.false.)
        flag = LA_NO_ERROR
        if (ldr < m .or. ldq < m) then
            flag = LA_INVALID_INPUT_ERROR
            return
        end if

        ! Process
        mn = min(m, n)
        nq = m
        if (.not.fullq) nq = n
        call form_qr(r(1:m,1:n), tau(1:mn), q(1:m,1:nq), err = err)
        if (err%has_error_occurred()) then
            flag = err%get_error_flag()
            return
        end if
    end function

! ------------------------------------------------------------------------------
    !> @brief Forms the full M-by-M orthogonal matrix Q from the elementary
    !! reflectors returned by the base QR factorization algorithm.  This
    !! routine also inflates the pivot array into an N-by-N matrix P such
    !! that A * P = Q * R.
    !!
    !! @param[in] fullq Set to true to always return the full Q matrix; else,
    !!  set to false, and in the event that M > N, Q may be supplied as M-by-N,
    !!  and therefore only return the useful submatrix Q1 (Q = [Q1, Q2]) as the
    !!  factorization can be written as Q * R = [Q1, Q2] * [R1; 0].
    !! @param[in] m The number of rows in R.
    !! @param[in] n The number of columns in R.
    !! @param[in,out] r On input, the M-by-N factored matrix as returned by the
    !!  QR factorization routine.  On output, the upper triangular matrix R.
    !! @param[in] ldr The leading dimension of matrix R.
    !! @param[in] tau A MIN(M, N)-element array containing the scalar factors of
    !!  each elementary reflector defined in @p r.
    !! @param[in] pvt An N-element array containing the pivot information from
    !!  the QR factorization.
    !! @param[out] q An M-by-M matrix where the full Q matrix will be written.
    !!  In the event that @p fullq is set to false, and M > N, this matrix need
    !!  only by M-by-N.
    !! @param[in] ldq The leading dimension of matrix Q.
    !! @param[out] p An N-by-N matrix where the pivot matrix P will be written.
    !! @param[in] ldp The leading dimension of matrix P.
    !!
    !! @return An error code.  The following codes are possible.
    !!  - LA_NO_ERROR: No error occurred.  Successful operation.
    !!  - LA_INVALID_INPUT_ERROR: Occurs if @p lda is not correct.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory 
    !!      available.
    function la_form_qr_pvt(fullq, m, n, r, ldr, tau, pvt, q, ldq, p, ldp) &
            bind(C, name = "la_form_qr_pvt") result(flag)
        ! Arguments
        logical(c_bool), intent(in), value :: fullq
        integer(c_int), intent(in), value :: m, n, ldr, ldq, ldp
        real(c_double), intent(inout) :: r(ldr,*)
        real(c_double), intent(in) :: tau(*)
        integer(c_int), intent(in) :: pvt(*)
        real(c_double), intent(out) :: q(ldq,*), p(ldp,*)
        integer(c_int) :: flag

        ! Local Variables
        type(errors) :: err
        integer(c_int) :: mn, nq

        ! Error Checking
        call err%set_exit_on_error(.false.)
        flag = LA_NO_ERROR
        if (ldr < m .or. ldq < m .or. ldp < n) then
            flag = LA_INVALID_INPUT_ERROR
            return
        end if

        ! Process
        mn = min(m, n)
        nq = m
        if (.not.fullq) nq = n
        call form_qr(r(1:m,1:n), tau(1:mn), pvt(1:n), q(1:m,1:nq), p(1:n,1:n), &
            err = err)
        if (err%has_error_occurred()) then
            flag = err%get_error_flag()
            return
        end if
    end function

! ------------------------------------------------------------------------------
    !> @brief Forms the full M-by-M orthogonal matrix Q from the elementary
    !! reflectors returned by the base QR factorization algorithm.  This
    !! routine also inflates the pivot array into an N-by-N matrix P such
    !! that A * P = Q * R.
    !!
    !! @param[in] fullq Set to true to always return the full Q matrix; else,
    !!  set to false, and in the event that M > N, Q may be supplied as M-by-N,
    !!  and therefore only return the useful submatrix Q1 (Q = [Q1, Q2]) as the
    !!  factorization can be written as Q * R = [Q1, Q2] * [R1; 0].
    !! @param[in] m The number of rows in R.
    !! @param[in] n The number of columns in R.
    !! @param[in,out] r On input, the M-by-N factored matrix as returned by the
    !!  QR factorization routine.  On output, the upper triangular matrix R.
    !! @param[in] ldr The leading dimension of matrix R.
    !! @param[in] tau A MIN(M, N)-element array containing the scalar factors of
    !!  each elementary reflector defined in @p r.
    !! @param[in] pvt An N-element array containing the pivot information from
    !!  the QR factorization.
    !! @param[out] q An M-by-M matrix where the full Q matrix will be written.
    !!  In the event that @p fullq is set to false, and M > N, this matrix need
    !!  only by M-by-N.
    !! @param[in] ldq The leading dimension of matrix Q.
    !! @param[out] p An N-by-N matrix where the pivot matrix P will be written.
    !! @param[in] ldp The leading dimension of matrix P.
    !!
    !! @return An error code.  The following codes are possible.
    !!  - LA_NO_ERROR: No error occurred.  Successful operation.
    !!  - LA_INVALID_INPUT_ERROR: Occurs if @p lda is not correct.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory 
    !!      available.
    function la_form_qr_cmplx_pvt(fullq, m, n, r, ldr, tau, pvt, q, ldq, p, &
            ldp) bind(C, name = "la_form_qr_cmplx_pvt") result(flag)
        ! Arguments
        logical(c_bool), intent(in), value :: fullq
        integer(c_int), intent(in), value :: m, n, ldr, ldq, ldp
        complex(c_double), intent(inout) :: r(ldr,*)
        complex(c_double), intent(in) :: tau(*)
        integer(c_int), intent(in) :: pvt(*)
        complex(c_double), intent(out) :: q(ldq,*), p(ldp,*)
        integer(c_int) :: flag

        ! Local Variables
        type(errors) :: err
        integer(c_int) :: mn, nq

        ! Error Checking
        call err%set_exit_on_error(.false.)
        flag = LA_NO_ERROR
        if (ldr < m .or. ldq < m .or. ldp < n) then
            flag = LA_INVALID_INPUT_ERROR
            return
        end if

        ! Process
        mn = min(m, n)
        nq = m
        if (.not.fullq) nq = n
        call form_qr(r(1:m,1:n), tau(1:mn), pvt(1:n), q(1:m,1:nq), p(1:n,1:n), &
            err = err)
        if (err%has_error_occurred()) then
            flag = err%get_error_flag()
            return
        end if
    end function

! ------------------------------------------------------------------------------
    !> @brief Multiplies a general matrix by the orthogonal matrix Q from a QR
    !! factorization such that: C = op(Q) * C, or C = C * op(Q).
    !!
    !! @param[in] lside Set to true to apply Q or Q**T from the left; else, set
    !!  to false to apply Q or Q**T from the right.
    !! @param[in] trans Set to true to apply Q**T; else, set to false.
    !! @param[in] m The number of rows in matrix C.
    !! @param[in] n The number of columns in matrix C.
    !! @param[in] k The number of elementary reflectors whose product defines 
    !!  the matrix Q.
    !! @param[in] a On input, an LDA-by-K matrix containing the elementary
    !!  reflectors output from the QR factorization.  If @p lside is set to
    !!  true, LDA = M, and M >= K >= 0; else, if @p lside is set to false,
    !!  LDA = N, and N >= K >= 0.  Notice, the contents of this matrix are
    !!  restored on exit.
    !! @param[in] lda The leading dimension of matrix A.
    !! @param[in] tau A K-element array containing the scalar factors of each
    !!  elementary reflector defined in @p a.
    !! @param[in,out] c On input, the M-by-N matrix C.  On output, the product
    !!  of the orthogonal matrix Q and the original matrix C.
    !! @param[in] ldc THe leading dimension of matrix C.
    !!
    !! @return An error code.  The following codes are possible.
    !!  - LA_NO_ERROR: No error occurred.  Successful operation.
    !!  - LA_INVALID_INPUT_ERROR: Occurs if @p lda is not correct.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory 
    !!      available.
    function la_mult_qr(lside, trans, m, n, k, a, lda, tau, c, ldc) &
            bind(C, name = "la_mult_qr") result(flag)
        ! Local Variables
        logical(c_bool), intent(in), value :: lside, trans
        integer(c_int), intent(in), value :: m, n, k, lda, ldc
        real(c_double), intent(inout) :: a(lda,*), c(ldc,*)
        real(c_double), intent(in) :: tau(*)
        integer(c_int) :: flag

        ! Local Variables
        type(errors) :: err
        integer(c_int) :: ma, na

        ! Initialization
        if (lside) then
            ma = m
            na = m
        else
            ma = n
            na = n
        end if

        ! Error Checking
        call err%set_exit_on_error(.false.)
        flag = LA_NO_ERROR
        if (lda < ma .or. ldc < m) then
            flag = LA_INVALID_INPUT_ERROR
            return
        end if
        if (k > na .or. k < 0) then
            flag = LA_INVALID_INPUT_ERROR
            return
        end if

        ! Process
        call mult_qr(logical(lside), logical(trans), a(1:ma,1:k), tau(1:k), &
            c(1:m,1:n), err = err)
        if (err%has_error_occurred()) then
            flag = err%get_error_flag()
            return
        end if
    end function

! ------------------------------------------------------------------------------
    !> @brief Multiplies a general matrix by the orthogonal matrix Q from a QR
    !! factorization such that: C = op(Q) * C, or C = C * op(Q).
    !!
    !! @param[in] lside Set to true to apply Q or Q**H from the left; else, set
    !!  to false to apply Q or Q**H from the right.
    !! @param[in] trans Set to true to apply Q**H; else, set to false.
    !! @param[in] m The number of rows in matrix C.
    !! @param[in] n The number of columns in matrix C.
    !! @param[in] k The number of elementary reflectors whose product defines 
    !!  the matrix Q.
    !! @param[in] a On input, an LDA-by-K matrix containing the elementary
    !!  reflectors output from the QR factorization.  If @p lside is set to
    !!  true, LDA = M, and M >= K >= 0; else, if @p lside is set to false,
    !!  LDA = N, and N >= K >= 0.  Notice, the contents of this matrix are
    !!  restored on exit.
    !! @param[in] lda The leading dimension of matrix A.
    !! @param[in] tau A K-element array containing the scalar factors of each
    !!  elementary reflector defined in @p a.
    !! @param[in,out] c On input, the M-by-N matrix C.  On output, the product
    !!  of the orthogonal matrix Q and the original matrix C.
    !! @param[in] ldc THe leading dimension of matrix C.
    !!
    !! @return An error code.  The following codes are possible.
    !!  - LA_NO_ERROR: No error occurred.  Successful operation.
    !!  - LA_INVALID_INPUT_ERROR: Occurs if @p lda is not correct.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory 
    !!      available.
    function la_mult_qr_cmplx(lside, trans, m, n, k, a, lda, tau, c, ldc) &
            bind(C, name = "la_mult_qr_cmplx") result(flag)
        ! Local Variables
        logical(c_bool), intent(in), value :: lside, trans
        integer(c_int), intent(in), value :: m, n, k, lda, ldc
        complex(c_double), intent(inout) :: a(lda,*), c(ldc,*)
        complex(c_double), intent(in) :: tau(*)
        integer(c_int) :: flag

        ! Local Variables
        type(errors) :: err
        integer(c_int) :: ma, na

        ! Initialization
        if (lside) then
            ma = m
            na = m
        else
            ma = n
            na = n
        end if

        ! Error Checking
        call err%set_exit_on_error(.false.)
        flag = LA_NO_ERROR
        if (lda < ma .or. ldc < m) then
            flag = LA_INVALID_INPUT_ERROR
            return
        end if
        if (k > na .or. k < 0) then
            flag = LA_INVALID_INPUT_ERROR
            return
        end if

        ! Process
        call mult_qr(logical(lside), logical(trans), a(1:ma,1:k), tau(1:k), &
            c(1:m,1:n), err = err)
        if (err%has_error_occurred()) then
            flag = err%get_error_flag()
            return
        end if
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the rank 1 update to an M-by-N QR factored matrix A
    !! (M >= N) where A = Q * R, and A1 = A + U * V**T such that A1 = Q1 * R1.
    !!
    !! @param[in] m The number of rows in R.
    !! @param[in] n The number of columns in R.
    !! @param[in,out] q On input, the original M-by-K orthogonal matrix Q.  On
    !!  output, the updated matrix Q1.
    !! @param[in] ldq The leading dimension of matrix Q.
    !! @param[in,out] r On input, the M-by-N matrix R.  On output, the updated
    !!  matrix R1.
    !! @param[in] ldr The leading dimension of matrix R.
    !! @param[in,out] u On input, the M-element U update vector.  On output,
    !!  the original content of the array is overwritten.
    !! @param[in,out] v On input, the N-element V update vector.  On output,
    !!  the original content of the array is overwritten.
    !!
    !! @return An error code.  The following codes are possible.
    !!  - LA_NO_ERROR: No error occurred.  Successful operation.
    !!  - LA_INVALID_INPUT_ERROR: Occurs if @p ldq or @p ldr is not correct.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory 
    !!      available.
    function la_qr_rank1_update(m, n, q, ldq, r, ldr, u, v) &
            bind(C, name = "la_qr_rank1_update") result(flag)
        ! Arguments
        integer(c_int), intent(in), value :: m, n, ldq, ldr
        real(c_double), intent(inout) :: q(ldq,*), r(ldr,*), u(*), v(*)
        integer(c_int) :: flag

        ! Local Variables
        type(errors) :: err
        integer(c_int) :: mn

        ! Error Checking
        call err%set_exit_on_error(.false.)
        flag = LA_NO_ERROR
        if (ldq < m .or. ldr < m) then
            flag = LA_INVALID_INPUT_ERROR
            return
        end if

        ! Process
        mn = min(m, n)
        call qr_rank1_update(q(1:m,1:mn), r(1:m,1:n), u(1:m), v(1:n), err = err)
        if (err%has_error_occurred()) then
            flag = err%get_error_flag()
            return
        end if
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the rank 1 update to an M-by-N QR factored matrix A
    !! (M >= N) where A = Q * R, and A1 = A + U * V**T such that A1 = Q1 * R1.
    !!
    !! @param[in] m The number of rows in R.
    !! @param[in] n The number of columns in R.
    !! @param[in,out] q On input, the original M-by-K orthogonal matrix Q.  On
    !!  output, the updated matrix Q1.
    !! @param[in] ldq The leading dimension of matrix Q.
    !! @param[in,out] r On input, the M-by-N matrix R.  On output, the updated
    !!  matrix R1.
    !! @param[in] ldr The leading dimension of matrix R.
    !! @param[in,out] u On input, the M-element U update vector.  On output,
    !!  the original content of the array is overwritten.
    !! @param[in,out] v On input, the N-element V update vector.  On output,
    !!  the original content of the array is overwritten.
    !!
    !! @return An error code.  The following codes are possible.
    !!  - LA_NO_ERROR: No error occurred.  Successful operation.
    !!  - LA_INVALID_INPUT_ERROR: Occurs if @p ldq or @p ldr is not correct.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory 
    !!      available.
    function la_qr_rank1_update_cmplx(m, n, q, ldq, r, ldr, u, v) &
            bind(C, name = "la_qr_rank1_update_cmplx") result(flag)
        ! Arguments
        integer(c_int), intent(in), value :: m, n, ldq, ldr
        complex(c_double), intent(inout) :: q(ldq,*), r(ldr,*), u(*), v(*)
        integer(c_int) :: flag

        ! Local Variables
        type(errors) :: err
        integer(c_int) :: mn

        ! Error Checking
        call err%set_exit_on_error(.false.)
        flag = LA_NO_ERROR
        if (ldq < m .or. ldr < m) then
            flag = LA_INVALID_INPUT_ERROR
            return
        end if

        ! Process
        mn = min(m, n)
        call qr_rank1_update(q(1:m,1:mn), r(1:m,1:n), u(1:m), v(1:n), err = err)
        if (err%has_error_occurred()) then
            flag = err%get_error_flag()
            return
        end if
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the Cholesky factorization of a symmetric, positive
    !! definite matrix.
    !!
    !! @param[in] upper Set to true to compute the upper triangular factoriztion
    !!  A = U**T * U; else, set to false to compute the lower triangular
    !!  factorzation A = L * L**T.
    !! @param[in] n The dimension of matrix A.
    !! @param[in,out] a On input, the N-by-N matrix to factor.  On output, the
    !!  factored matrix is returned in either the upper or lower triangular
    !!  portion of the matrix, dependent upon the value of @p upper.
    !! @param[in] lda The leading dimension of matrix A.
    !!
    !! @return An error code.  The following codes are possible.
    !!  - LA_NO_ERROR: No error occurred.  Successful operation.
    !!  - LA_INVALID_INPUT_ERROR: Occurs if @p lda is not correct.
    !!  - LA_MATRIX_FORMAT_ERROR: Occurs if @p a is not positive definite.
    function la_cholesky_factor(upper, n, a, lda) &
            bind(C, name = "la_cholesky_factor") result(flag)
        ! Arguments
        logical(c_bool), intent(in), value :: upper
        integer(c_int), intent(in), value :: n, lda
        real(c_double), intent(inout) :: a(lda,*)
        integer(c_int) :: flag

        ! Local Variables
        type(errors) :: err

        ! Error Checking
        call err%set_exit_on_error(.false.)
        flag = LA_NO_ERROR
        if (lda < n) then
            flag = LA_INVALID_INPUT_ERROR
            return
        end if

        ! Process
        call cholesky_factor(a(1:n,1:n), logical(upper), err = err)
        if (err%has_error_occurred()) then
            flag = err%get_error_flag()
            return
        end if
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the Cholesky factorization of a symmetric, positive
    !! definite matrix.
    !!
    !! @param[in] upper Set to true to compute the upper triangular factoriztion
    !!  A = U**H * U; else, set to false to compute the lower triangular
    !!  factorzation A = L * L**H.
    !! @param[in] n The dimension of matrix A.
    !! @param[in,out] a On input, the N-by-N matrix to factor.  On output, the
    !!  factored matrix is returned in either the upper or lower triangular
    !!  portion of the matrix, dependent upon the value of @p upper.
    !! @param[in] lda The leading dimension of matrix A.
    !!
    !! @return An error code.  The following codes are possible.
    !!  - LA_NO_ERROR: No error occurred.  Successful operation.
    !!  - LA_INVALID_INPUT_ERROR: Occurs if @p lda is not correct.
    !!  - LA_MATRIX_FORMAT_ERROR: Occurs if @p a is not positive definite.
    function la_cholesky_factor_cmplx(upper, n, a, lda) &
            bind(C, name = "la_cholesky_factor_cmplx") result(flag)
        ! Arguments
        logical(c_bool), intent(in), value :: upper
        integer(c_int), intent(in), value :: n, lda
        complex(c_double), intent(inout) :: a(lda,*)
        integer(c_int) :: flag

        ! Local Variables
        type(errors) :: err

        ! Error Checking
        call err%set_exit_on_error(.false.)
        flag = LA_NO_ERROR
        if (lda < n) then
            flag = LA_INVALID_INPUT_ERROR
            return
        end if

        ! Process
        call cholesky_factor(a(1:n,1:n), logical(upper), err = err)
        if (err%has_error_occurred()) then
            flag = err%get_error_flag()
            return
        end if
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the rank 1 update to a Cholesky factored matrix (upper
    !! triangular).
    !!
    !! @param[in] n The dimension of the matrix.
    !! @param[in,out] r On input, the N-by-N upper triangular matrix R.  On
    !!  output, the updated matrix R1.
    !! @param[in] ldr The leading dimension of matrix R.
    !! @param[in,out] u On input, the N-element update vector U.  On output,
    !!  the rotation sines used to transform R to R1.
    !!
    !! @return An error code.  The following codes are possible.
    !!  - LA_NO_ERROR: No error occurred.  Successful operation.
    !!  - LA_INVALID_INPUT_ERROR: Occurs if @p ldr is not correct.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
    !!      there is insufficient memory available.
    function la_cholesky_rank1_update(n, r, ldr, u) &
            bind(C, name = "la_cholesky_rank1_update") result(flag)
        ! Arguments
        integer(c_int), intent(in), value :: n, ldr
        real(c_double), intent(inout) :: r(ldr,*), u(*)
        integer(c_int) :: flag

        ! Local Variables
        type(errors) :: err

        ! Error Checking
        call err%set_exit_on_error(.false.)
        flag = LA_NO_ERROR
        if (ldr < n) then
            flag = LA_INVALID_INPUT_ERROR
            return
        end if

        ! Process
        call cholesky_rank1_update(r(1:n,1:n), u(1:n), err = err)
        if (err%has_error_occurred()) then
            flag = err%get_error_flag()
            return
        end if
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the rank 1 update to a Cholesky factored matrix (upper
    !! triangular).
    !!
    !! @param[in] n The dimension of the matrix.
    !! @param[in,out] r On input, the N-by-N upper triangular matrix R.  On
    !!  output, the updated matrix R1.
    !! @param[in] ldr The leading dimension of matrix R.
    !! @param[in,out] u On input, the N-element update vector U.  On output,
    !!  the rotation sines used to transform R to R1.
    !!
    !! @return An error code.  The following codes are possible.
    !!  - LA_NO_ERROR: No error occurred.  Successful operation.
    !!  - LA_INVALID_INPUT_ERROR: Occurs if @p ldr is not correct.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
    !!      there is insufficient memory available.
    function la_cholesky_rank1_update_cmplx(n, r, ldr, u) &
            bind(C, name = "la_cholesky_rank1_update_cmplx") result(flag)
        ! Arguments
        integer(c_int), intent(in), value :: n, ldr
        complex(c_double), intent(inout) :: r(ldr,*), u(*)
        integer(c_int) :: flag

        ! Local Variables
        type(errors) :: err

        ! Error Checking
        call err%set_exit_on_error(.false.)
        flag = LA_NO_ERROR
        if (ldr < n) then
            flag = LA_INVALID_INPUT_ERROR
            return
        end if

        ! Process
        call cholesky_rank1_update(r(1:n,1:n), u(1:n), err = err)
        if (err%has_error_occurred()) then
            flag = err%get_error_flag()
            return
        end if
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the rank 1 downdate to a Cholesky factored matrix (upper
    !! triangular).
    !!
    !! @param[in] n The dimension of the matrix.
    !! @param[in,out] r On input, the N-by-N upper triangular matrix R.  On
    !!  output, the updated matrix R1.
    !! @param[in] ldr The leading dimension of matrix R.
    !! @param[in,out] u On input, the N-element update vector U.  On output,
    !!  the rotation sines used to transform R to R1.
    !!
    !! @return An error code.  The following codes are possible.
    !!  - LA_NO_ERROR: No error occurred.  Successful operation.
    !!  - LA_INVALID_INPUT_ERROR: Occurs if @p ldr is not correct.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
    !!      there is insufficient memory available.
    !!  - LA_MATRIX_FORMAT_ERROR: Occurs if the downdated matrix is not
    !!      positive definite.
    function la_cholesky_rank1_downdate(n, r, ldr, u) &
            bind(C, name = "la_cholesky_rank1_downdate") result(flag)
        ! Arguments
        integer(c_int), intent(in), value :: n, ldr
        real(c_double), intent(inout) :: r(ldr,*), u(*)
        integer(c_int) :: flag

        ! Local Variables
        type(errors) :: err

        ! Error Checking
        call err%set_exit_on_error(.false.)
        flag = LA_NO_ERROR
        if (ldr < n) then
            flag = LA_INVALID_INPUT_ERROR
            return
        end if

        ! Process
        call cholesky_rank1_downdate(r(1:n,1:n), u(1:n), err = err)
        if (err%has_error_occurred()) then
            flag = err%get_error_flag()
            return
        end if
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the rank 1 downdate to a Cholesky factored matrix (upper
    !! triangular).
    !!
    !! @param[in] n The dimension of the matrix.
    !! @param[in,out] r On input, the N-by-N upper triangular matrix R.  On
    !!  output, the updated matrix R1.
    !! @param[in] ldr The leading dimension of matrix R.
    !! @param[in,out] u On input, the N-element update vector U.  On output,
    !!  the rotation sines used to transform R to R1.
    !!
    !! @return An error code.  The following codes are possible.
    !!  - LA_NO_ERROR: No error occurred.  Successful operation.
    !!  - LA_INVALID_INPUT_ERROR: Occurs if @p ldr is not correct.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
    !!      there is insufficient memory available.
    !!  - LA_MATRIX_FORMAT_ERROR: Occurs if the downdated matrix is not
    !!      positive definite.
    function la_cholesky_rank1_downdate_cmplx(n, r, ldr, u) &
            bind(C, name = "la_cholesky_rank1_downdate_cmplx") result(flag)
        ! Arguments
        integer(c_int), intent(in), value :: n, ldr
        complex(c_double), intent(inout) :: r(ldr,*), u(*)
        integer(c_int) :: flag

        ! Local Variables
        type(errors) :: err

        ! Error Checking
        call err%set_exit_on_error(.false.)
        flag = LA_NO_ERROR
        if (ldr < n) then
            flag = LA_INVALID_INPUT_ERROR
            return
        end if

        ! Process
        call cholesky_rank1_downdate(r(1:n,1:n), u(1:n), err = err)
        if (err%has_error_occurred()) then
            flag = err%get_error_flag()
            return
        end if
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the singular value decomposition of a matrix A.  The
    !!  SVD is defined as: A = U * S * V**T, where U is an M-by-M orthogonal
    !!  matrix, S is an M-by-N diagonal matrix, and V is an N-by-N orthogonal
    !!  matrix.
    !!
    !! @param[in] m The number of rows in the matrix.
    !! @param[in] n The number of columns in the matrix.
    !! @param[in,out] a On input, the M-by-N matrix to factor.  The matrix is
    !!  overwritten on output.
    !! @param[in] lda The leading dimension of matrix A.
    !! @param[out] s A MIN(M, N)-element array containing the singular values
    !!  of @p a sorted in descending order.
    !! @param[out] u An M-by-M matrix where the orthogonal U matrix will be
    !!  written.
    !! @param[in] ldu The leading dimension of matrix U.
    !! @param[out] vt An N-by-N matrix where the transpose of the right 
    !!  singular vector matrix V.
    !! @param[in] ldv The leading dimension of matrix V.
    !!
    !! @return An error code.  The following codes are possible.
    !!  - LA_NO_ERROR: No error occurred.  Successful operation.
    !!  - LA_INVALID_INPUT_ERROR: Occurs if @p lda, @p ldu, or @p ldv is not 
    !!      correct.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
    !!      there is insufficient memory available.
    !!  - LA_CONVERGENCE_ERROR: Occurs as a warning if the QR iteration process
    !!      could not converge to a zero value.
    function la_svd(m, n, a, lda, s, u, ldu, vt, ldv) &
            bind(C, name = "la_svd") result(flag)
        ! Arguments
        integer(c_int), intent(in), value :: m, n, lda, ldu, ldv
        real(c_double), intent(inout) :: a(lda,*)
        real(c_double), intent(out) :: s(*), u(ldu,*), vt(ldv,*)
        integer(c_int) :: flag

        ! Local Variables
        type(errors) :: err
        integer(c_int) :: mn

        ! Error Checking
        call err%set_exit_on_error(.false.)
        flag = LA_NO_ERROR
        if (lda < m .or. ldu < m .or. ldv < n) then
            flag = LA_INVALID_INPUT_ERROR
            return
        end if

        ! Process
        mn = min(m, n)
        call svd(a(1:m,1:n), s(1:mn), u(1:m,1:m), vt(1:n,1:n), err = err)
        if (err%has_error_occurred()) then
            flag = err%get_error_flag()
            return
        end if
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the singular value decomposition of a matrix A.  The
    !!  SVD is defined as: A = U * S * V**T, where U is an M-by-M orthogonal
    !!  matrix, S is an M-by-N diagonal matrix, and V is an N-by-N orthogonal
    !!  matrix.
    !!
    !! @param[in] m The number of rows in the matrix.
    !! @param[in] n The number of columns in the matrix.
    !! @param[in,out] a On input, the M-by-N matrix to factor.  The matrix is
    !!  overwritten on output.
    !! @param[in] lda The leading dimension of matrix A.
    !! @param[out] s A MIN(M, N)-element array containing the singular values
    !!  of @p a sorted in descending order.
    !! @param[out] u An M-by-M matrix where the orthogonal U matrix will be
    !!  written.
    !! @param[in] ldu The leading dimension of matrix U.
    !! @param[out] vt An N-by-N matrix where the transpose of the right 
    !!  singular vector matrix V.
    !! @param[in] ldv The leading dimension of matrix V.
    !!
    !! @return An error code.  The following codes are possible.
    !!  - LA_NO_ERROR: No error occurred.  Successful operation.
    !!  - LA_INVALID_INPUT_ERROR: Occurs if @p lda, @p ldu, or @p ldv is not 
    !!      correct.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
    !!      there is insufficient memory available.
    !!  - LA_CONVERGENCE_ERROR: Occurs as a warning if the QR iteration process
    !!      could not converge to a zero value.
    function la_svd_cmplx(m, n, a, lda, s, u, ldu, vt, ldv) &
            bind(C, name = "la_svd_cmplx") result(flag)
        ! Arguments
        integer(c_int), intent(in), value :: m, n, lda, ldu, ldv
        complex(c_double), intent(inout) :: a(lda,*)
        real(c_double), intent(out) :: s(*)
        complex(c_double), intent(out) :: u(ldu,*), vt(ldv,*)
        integer(c_int) :: flag

        ! Local Variables
        type(errors) :: err
        integer(c_int) :: mn

        ! Error Checking
        call err%set_exit_on_error(.false.)
        flag = LA_NO_ERROR
        if (lda < m .or. ldu < m .or. ldv < n) then
            flag = LA_INVALID_INPUT_ERROR
            return
        end if

        ! Process
        mn = min(m, n)
        call svd(a(1:m,1:n), s(1:mn), u(1:m,1:m), vt(1:n,1:n), err = err)
        if (err%has_error_occurred()) then
            flag = err%get_error_flag()
            return
        end if
    end function

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
    !! @param[in] m The number of rows in matrix B.
    !! @param[in] n The number of columns in matrix B.
    !! @param[in] alpha The scalar multiplier to B.
    !! @param[in] a If @p lside is true, the M-by-M triangular matrix on which
    !!  to operate; else, if @p lside is false, the N-by-N triangular matrix on
    !!  which to operate.
    !! @param[in] lda The leading dimension of matrix A.
    !! @param[in,out] b On input, the M-by-N right-hand-side.  On output, the
    !!  M-by-N solution.
    !! @param[in] ldb The leading dimension of matrix B.
    !!
    !! @return An error code.  The following codes are possible.
    !!  - LA_NO_ERROR: No error occurred.  Successful operation.
    !!  - LA_INVALID_INPUT_ERROR: Occurs if @p lda, or @p ldb is not correct.
    function la_solve_tri_mtx(lside, upper, trans, nounit, m, n, alpha, a, &
            lda, b, ldb) bind(C, name = "la_solve_tri_mtx") result(flag)
        ! Arguments
        logical(c_bool), intent(in), value :: lside, upper, trans, nounit
        integer(c_int), intent(in), value :: m, n, lda, ldb
        real(c_double), intent(in), value :: alpha
        real(c_double), intent(in) :: a(lda,*)
        real(c_double), intent(inout) :: b(ldb,*)
        integer(c_int) :: flag

        ! Local Variables
        type(errors) :: err
        integer(c_int) :: ma

        ! Initialization
        if (lside) then
            ma = m
        else
            ma = n
        end if

        ! Error Checking
        call err%set_exit_on_error(.false.)
        flag = LA_NO_ERROR
        if (lda < ma .or. ldb < m) then
            flag = LA_INVALID_INPUT_ERROR
            return
        end if

        ! Process
        call solve_triangular_system(logical(lside), logical(upper), &
            logical(trans), logical(nounit), alpha, a(1:ma,1:ma), b(1:m,1:n))
    end function

! ------------------------------------------------------------------------------
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
    !! @param[in] m The number of rows in matrix B.
    !! @param[in] n The number of columns in matrix B.
    !! @param[in] alpha The scalar multiplier to B.
    !! @param[in] a If @p lside is true, the M-by-M triangular matrix on which
    !!  to operate; else, if @p lside is false, the N-by-N triangular matrix on
    !!  which to operate.
    !! @param[in] lda The leading dimension of matrix A.
    !! @param[in,out] b On input, the M-by-N right-hand-side.  On output, the
    !!  M-by-N solution.
    !! @param[in] ldb The leading dimension of matrix B.
    !!
    !! @return An error code.  The following codes are possible.
    !!  - LA_NO_ERROR: No error occurred.  Successful operation.
    !!  - LA_INVALID_INPUT_ERROR: Occurs if @p lda, or @p ldb is not correct.
    function la_solve_tri_mtx_cmplx(lside, upper, trans, nounit, m, n, &
            alpha, a, lda, b, ldb) &
            bind(C, name = "la_solve_tri_mtx_cmplx") result(flag)
        ! Arguments
        logical(c_bool), intent(in), value :: lside, upper, trans, nounit
        integer(c_int), intent(in), value :: m, n, lda, ldb
        complex(c_double), intent(in), value :: alpha
        complex(c_double), intent(in) :: a(lda,*)
        complex(c_double), intent(inout) :: b(ldb,*)
        integer(c_int) :: flag

        ! Local Variables
        type(errors) :: err
        integer(c_int) :: ma

        ! Initialization
        if (lside) then
            ma = m
        else
            ma = n
        end if

        ! Error Checking
        call err%set_exit_on_error(.false.)
        flag = LA_NO_ERROR
        if (lda < ma .or. ldb < m) then
            flag = LA_INVALID_INPUT_ERROR
            return
        end if

        ! Process
        call solve_triangular_system(logical(lside), logical(upper), &
            logical(trans), logical(nounit), alpha, a(1:ma,1:ma), b(1:m,1:n))
    end function

! ------------------------------------------------------------------------------
    !> @brief Solves a system of LU-factored equations.
    !!
    !! @param[in] m The number of rows in matrix B.
    !! @param[in] n The number of columns in matrix B.
    !! @param[in] a The M-by-M LU factored matrix.
    !! @param[in] lda The leading dimension of matrix A.
    !! @param[in] ipvt The M-element pivot array from the LU factorization.
    !! @param[in,out] b On input, the M-by-N right-hand-side.  On output, the
    !!  M-by-N solution.
    !! @param[in] ldb The leading dimension of matrix B.
    !!
    !! @return An error code.  The following codes are possible.
    !!  - LA_NO_ERROR: No error occurred.  Successful operation.
    !!  - LA_INVALID_INPUT_ERROR: Occurs if @p lda, or @p ldb is not correct.
    function la_solve_lu(m, n, a, lda, ipvt, b, ldb) &
            bind(C, name = "la_solve_lu") result(flag)
        ! Arguments
        integer(c_int), intent(in), value :: m, n, lda, ldb
        real(c_double), intent(in) :: a(lda,*)
        integer(c_int), intent(in) :: ipvt(*)
        real(c_double), intent(inout) :: b(ldb,*)
        integer(c_int) :: flag

        ! Local Variables
        type(errors) :: err

        ! Error Checking
        call err%set_exit_on_error(.false.)
        flag = LA_NO_ERROR
        if (lda < m .or. ldb < m) then
            flag = LA_INVALID_INPUT_ERROR
            return
        end if

        ! Process
        call solve_lu(a(1:m,1:m), ipvt(1:m), b(1:m,1:n))
    end function

! ------------------------------------------------------------------------------
    !> @brief Solves a system of LU-factored equations.
    !!
    !! @param[in] m The number of rows in matrix B.
    !! @param[in] n The number of columns in matrix B.
    !! @param[in] a The M-by-M LU factored matrix.
    !! @param[in] lda The leading dimension of matrix A.
    !! @param[in] ipvt The M-element pivot array from the LU factorization.
    !! @param[in,out] b On input, the M-by-N right-hand-side.  On output, the
    !!  M-by-N solution.
    !! @param[in] ldb The leading dimension of matrix B.
    !!
    !! @return An error code.  The following codes are possible.
    !!  - LA_NO_ERROR: No error occurred.  Successful operation.
    !!  - LA_INVALID_INPUT_ERROR: Occurs if @p lda, or @p ldb is not correct.
    function la_solve_lu_cmplx(m, n, a, lda, ipvt, b, ldb) &
            bind(C, name = "la_solve_lu_cmplx") result(flag)
        ! Arguments
        integer(c_int), intent(in), value :: m, n, lda, ldb
        complex(c_double), intent(in) :: a(lda,*)
        integer(c_int), intent(in) :: ipvt(*)
        complex(c_double), intent(inout) :: b(ldb,*)
        integer(c_int) :: flag

        ! Local Variables
        type(errors) :: err

        ! Error Checking
        call err%set_exit_on_error(.false.)
        flag = LA_NO_ERROR
        if (lda < m .or. ldb < m) then
            flag = LA_INVALID_INPUT_ERROR
            return
        end if

        ! Process
        call solve_lu(a(1:m,1:m), ipvt(1:m), b(1:m,1:n))
    end function

! ------------------------------------------------------------------------------
    !> @brief Solves a system of M QR-factored equations of N unknowns where
    !! M >= N.
    !!
    !! @param[in] m The number of equations (rows in matrix A).
    !! @param[in] n The number of unknowns (columns in matrix A).
    !! @param[in] k The number of columns in the right-hand-side matrix.
    !! @param[in,out] a On input, the M-by-N QR factored matrix as returned by
    !!  qr_factor.  On output, the contents of this matrix are restored.
    !! @param[in] lda The leading dimension of matrix A.
    !! @param[in] tau A MIN(M, N)-element array containing the scalar factors of
    !!  the elementary reflectors as returned by qr_factor.
    !! @param[in,out] b On input, the M-by-K right-hand-side matrix.  On output,
    !!  the first N rows are overwritten by the solution matrix X.
    !! @param[in] ldb The leading dimension of matrix B.
    !!
    !! @return An error code.  The following codes are possible.
    !!  - LA_NO_ERROR: No error occurred.  Successful operation.
    !!  - LA_INVALID_INPUT_ERROR: Occurs if @p lda, or @p ldb is not correct, or
    !!      if @p m is less than @p n.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
    !!      there is insufficient memory available.
    function la_solve_qr(m, n, k, a, lda, tau, b, ldb) &
            bind(C, name = "la_solve_qr") result(flag)
        ! Arguments
        integer(c_int), intent(in), value :: m, n, k, lda, ldb
        real(c_double), intent(inout) :: a(lda,*), b(ldb,*)
        real(c_double), intent(in) :: tau(*)
        integer(c_int) :: flag

        ! Local Variables
        type(errors) :: err
        integer(c_int) :: minmn

        ! Error Checking
        call err%set_exit_on_error(.false.)
        flag = LA_NO_ERROR
        if (lda < m .or. ldb < m .or. m < n) then
            flag = LA_INVALID_INPUT_ERROR
            return
        end if

        ! Process
        minmn = min(m, n)
        call solve_qr(a(1:m,1:n), tau(1:minmn), b(1:m,1:k), err = err)
        if (err%has_error_occurred()) then
            flag = err%get_error_flag()
            return
        end if
    end function

! ------------------------------------------------------------------------------
    !> @brief Solves a system of M QR-factored equations of N unknowns where
    !! M >= N.
    !!
    !! @param[in] m The number of equations (rows in matrix A).
    !! @param[in] n The number of unknowns (columns in matrix A).
    !! @param[in] k The number of columns in the right-hand-side matrix.
    !! @param[in,out] a On input, the M-by-N QR factored matrix as returned by
    !!  qr_factor.  On output, the contents of this matrix are restored.
    !! @param[in] lda The leading dimension of matrix A.
    !! @param[in] tau A MIN(M, N)-element array containing the scalar factors of
    !!  the elementary reflectors as returned by qr_factor.
    !! @param[in,out] b On input, the M-by-K right-hand-side matrix.  On output,
    !!  the first N rows are overwritten by the solution matrix X.
    !! @param[in] ldb The leading dimension of matrix B.
    !!
    !! @return An error code.  The following codes are possible.
    !!  - LA_NO_ERROR: No error occurred.  Successful operation.
    !!  - LA_INVALID_INPUT_ERROR: Occurs if @p lda, or @p ldb is not correct, or
    !!      if @p m is less than @p n.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
    !!      there is insufficient memory available.
    function la_solve_qr_cmplx(m, n, k, a, lda, tau, b, ldb) &
            bind(C, name = "la_solve_qr_cmplx") result(flag)
        ! Arguments
        integer(c_int), intent(in), value :: m, n, k, lda, ldb
        complex(c_double), intent(inout) :: a(lda,*), b(ldb,*)
        complex(c_double), intent(in) :: tau(*)
        integer(c_int) :: flag

        ! Local Variables
        type(errors) :: err
        integer(c_int) :: minmn

        ! Error Checking
        call err%set_exit_on_error(.false.)
        flag = LA_NO_ERROR
        if (lda < m .or. ldb < m .or. m < n) then
            flag = LA_INVALID_INPUT_ERROR
            return
        end if

        ! Process
        minmn = min(m, n)
        call solve_qr(a(1:m,1:n), tau(1:minmn), b(1:m,1:k), err = err)
        if (err%has_error_occurred()) then
            flag = err%get_error_flag()
            return
        end if
    end function

! ------------------------------------------------------------------------------
    !> @brief Solves a system of M QR-factored equations of N unknowns.
    !!
    !! @param[in] m The number of equations (rows in matrix A).
    !! @param[in] n The number of unknowns (columns in matrix A).
    !! @param[in] k The number of columns in the right-hand-side matrix.
    !! @param[in,out] a On input, the M-by-N QR factored matrix as returned by
    !!  qr_factor.  On output, the contents of this matrix are restored.
    !! @param[in] lda The leading dimension of matrix A.
    !! @param[in] tau A MIN(M, N)-element array containing the scalar factors of
    !!  the elementary reflectors as returned by qr_factor.
    !! @param[in] jpvt The N-element array that was used to track the column
    !!  pivoting operations in the QR factorization.
    !! @param[in,out] b On input, the M-by-K right-hand-side matrix.  On output,
    !!  the first N rows are overwritten by the solution matrix X.
    !! @param[in] ldb The leading dimension of matrix B.
    !!
    !! @return An error code.  The following codes are possible.
    !!  - LA_NO_ERROR: No error occurred.  Successful operation.
    !!  - LA_INVALID_INPUT_ERROR: Occurs if @p lda, or @p ldb is not correct.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
    !!      there is insufficient memory available.
    function la_solve_qr_pvt(m, n, k, a, lda, tau, jpvt, b, ldb) &
            bind(C, name = "la_solve_qr_pvt") result(flag)
        ! Arguments
        integer(c_int), intent(in), value :: m, n, k, lda, ldb
        real(c_double), intent(inout) :: a(lda,*), b(ldb,*)
        real(c_double), intent(in) :: tau(*)
        integer(c_int), intent(in) :: jpvt(*)
        integer(c_int) :: flag

        ! Local Variables
        type(errors) :: err
        integer(c_int) :: minmn, maxmn

        ! Error Checking
        minmn = min(m, n)
        maxmn = max(m, n)
        call err%set_exit_on_error(.false.)
        flag = LA_NO_ERROR
        if (lda < m .or. ldb < maxmn) then
            flag = LA_INVALID_INPUT_ERROR
            return
        end if

        ! Process
        call solve_qr(a(1:m,1:n), tau(1:minmn), jpvt(1:n), b(1:maxmn,1:k), &
            err = err)
        if (err%has_error_occurred()) then
            flag = err%get_error_flag()
            return
        end if
    end function

! ------------------------------------------------------------------------------
    !> @brief Solves a system of M QR-factored equations of N unknowns.
    !!
    !! @param[in] m The number of equations (rows in matrix A).
    !! @param[in] n The number of unknowns (columns in matrix A).
    !! @param[in] k The number of columns in the right-hand-side matrix.
    !! @param[in,out] a On input, the M-by-N QR factored matrix as returned by
    !!  qr_factor.  On output, the contents of this matrix are restored.
    !! @param[in] lda The leading dimension of matrix A.
    !! @param[in] tau A MIN(M, N)-element array containing the scalar factors of
    !!  the elementary reflectors as returned by qr_factor.
    !! @param[in] jpvt The N-element array that was used to track the column
    !!  pivoting operations in the QR factorization.
    !! @param[in,out] b On input, the M-by-K right-hand-side matrix.  On output,
    !!  the first N rows are overwritten by the solution matrix X.
    !! @param[in] ldb The leading dimension of matrix B.
    !!
    !! @return An error code.  The following codes are possible.
    !!  - LA_NO_ERROR: No error occurred.  Successful operation.
    !!  - LA_INVALID_INPUT_ERROR: Occurs if @p lda, or @p ldb is not correct.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
    !!      there is insufficient memory available.
    function la_solve_qr_cmplx_pvt(m, n, k, a, lda, tau, jpvt, b, ldb) &
            bind(C, name = "la_solve_qr_cmplx_pvt") result(flag)
        ! Arguments
        integer(c_int), intent(in), value :: m, n, k, lda, ldb
        complex(c_double), intent(inout) :: a(lda,*), b(ldb,*)
        complex(c_double), intent(in) :: tau(*)
        integer(c_int), intent(in) :: jpvt(*)
        integer(c_int) :: flag

        ! Local Variables
        type(errors) :: err
        integer(c_int) :: minmn, maxmn

        ! Error Checking
        minmn = min(m, n)
        maxmn = max(m, n)
        call err%set_exit_on_error(.false.)
        flag = LA_NO_ERROR
        if (lda < m .or. ldb < maxmn) then
            flag = LA_INVALID_INPUT_ERROR
            return
        end if

        ! Process
        call solve_qr(a(1:m,1:n), tau(1:minmn), jpvt(1:n), b(1:maxmn,1:k), &
            err = err)
        if (err%has_error_occurred()) then
            flag = err%get_error_flag()
            return
        end if
    end function

! ------------------------------------------------------------------------------
    !> @brief Solves a system of Cholesky factored equations.
    !!
    !! @param[in] upper Set to true if the original matrix A was factored such
    !!  that A = U**T * U; else, set to false if the factorization of A was
    !!  A = L**T * L.
    !! @param[in] m The number of rows in matrix B.
    !! @param[in] n The number of columns in matrix B.
    !! @param[in] a The M-by-M Cholesky factored matrix.
    !! @param[in] lda The leading dimension of matrix A.
    !! @param[in,out] b On input, the M-by-N right-hand-side matrix B.  On
    !!  output, the M-by-N solution matrix X.
    !! @param[in] ldb The leading dimension of matrix B.
    !!
    !! @return An error code.  The following codes are possible.
    !!  - LA_NO_ERROR: No error occurred.  Successful operation.
    !!  - LA_INVALID_INPUT_ERROR: Occurs if @p lda, or @p ldb is not correct.
    function la_solve_cholesky(upper, m, n, a, lda, b, ldb) &
            bind(C, name = "la_solve_cholesky") result(flag)
        ! Arguments
        logical(c_bool), intent(in), value :: upper
        integer(c_int), intent(in), value :: m, n, lda, ldb
        real(c_double), intent(in) :: a(lda,*)
        real(c_double), intent(inout) :: b(ldb,*)
        integer(c_int) :: flag

        ! Local Variables
        type(errors) :: err

        ! Error Checking
        call err%set_exit_on_error(.false.)
        flag = LA_NO_ERROR
        if (lda < m .or. ldb < m) then
            flag = LA_INVALID_INPUT_ERROR
            return
        end if

        ! Process
        call solve_cholesky(logical(upper), a(1:m,1:m), b(1:m,1:n))
    end function

! ------------------------------------------------------------------------------
    !> @brief Solves a system of Cholesky factored equations.
    !!
    !! @param[in] upper Set to true if the original matrix A was factored such
    !!  that A = U**T * U; else, set to false if the factorization of A was
    !!  A = L**T * L.
    !! @param[in] m The number of rows in matrix B.
    !! @param[in] n The number of columns in matrix B.
    !! @param[in] a The M-by-M Cholesky factored matrix.
    !! @param[in] lda The leading dimension of matrix A.
    !! @param[in,out] b On input, the M-by-N right-hand-side matrix B.  On
    !!  output, the M-by-N solution matrix X.
    !! @param[in] ldb The leading dimension of matrix B.
    !!
    !! @return An error code.  The following codes are possible.
    !!  - LA_NO_ERROR: No error occurred.  Successful operation.
    !!  - LA_INVALID_INPUT_ERROR: Occurs if @p lda, or @p ldb is not correct.
    function la_solve_cholesky_cmplx(upper, m, n, a, lda, b, ldb) &
            bind(C, name = "la_solve_cholesky_cmplx") result(flag)
        ! Arguments
        logical(c_bool), intent(in), value :: upper
        integer(c_int), intent(in), value :: m, n, lda, ldb
        complex(c_double), intent(in) :: a(lda,*)
        complex(c_double), intent(inout) :: b(ldb,*)
        integer(c_int) :: flag

        ! Local Variables
        type(errors) :: err

        ! Error Checking
        call err%set_exit_on_error(.false.)
        flag = LA_NO_ERROR
        if (lda < m .or. ldb < m) then
            flag = LA_INVALID_INPUT_ERROR
            return
        end if

        ! Process
        call solve_cholesky(logical(upper), a(1:m,1:m), b(1:m,1:n))
    end function

! ------------------------------------------------------------------------------
    !> @brief Solves the overdetermined or underdetermined system (A*X = B) of
    !! M equations of N unknowns using a QR or LQ factorization of the matrix A.
    !! Notice, it is assumed that matrix A has full rank.
    !!
    !! @param[in] m The number of equations (rows in matrix A).
    !! @param[in] n The number of unknowns (columns in matrix A).
    !! @param[in] k The number of columns in the right-hand-side matrix.
    !! @param[in,out] a On input, the M-by-N matrix A.  On output, if M >= N,
    !!  the QR factorization of A in the form as output by qr_factor; else,
    !!  if M < N, the LQ factorization of A.
    !! @param[in] lda The leading dimension of matrix A.
    !! @param[in,out] b If M >= N, the M-by-NRHS matrix B.  On output, the first
    !!  N rows contain the N-by-NRHS solution matrix X.  If M < N, an
    !!  N-by-NRHS matrix with the first M rows containing the matrix B.  On
    !!  output, the N-by-NRHS solution matrix X.
    !! @param[in] ldb The leading dimension of matrix B.
    !!
    !! @return An error code.  The following codes are possible.
    !!  - LA_NO_ERROR: No error occurred.  Successful operation.
    !!  - LA_INVALID_INPUT_ERROR: Occurs if @p lda, or @p ldb is not correct.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
    !!      there is insufficient memory available.
    !!  - LA_INVALID_OPERATION_ERROR: Occurs if @p a is not of full rank.
    function la_solve_least_squares(m, n, k, a, lda, b, ldb) &
            bind(C, name = "la_solve_least_squares") result(flag)
        ! Arguments
        integer(c_int), intent(in), value :: m, n, k, lda, ldb
        real(c_double), intent(inout) :: a(lda,*), b(ldb,*)
        integer(c_int) :: flag

        ! Local Variables
        type(errors) :: err
        integer(c_int) :: maxmn

        ! Error Checking
        maxmn = max(m, n)
        call err%set_exit_on_error(.false.)
        flag = LA_NO_ERROR
        if (lda < m .or. ldb < maxmn) then
            flag = LA_INVALID_INPUT_ERROR
            return
        end if

        ! Process
        call solve_least_squares(a(1:m,1:n), b(1:maxmn,1:k), err = err)
        if (err%has_error_occurred()) then
            flag = err%get_error_flag()
            return
        end if
    end function

! ------------------------------------------------------------------------------
    !> @brief Solves the overdetermined or underdetermined system (A*X = B) of
    !! M equations of N unknowns using a QR or LQ factorization of the matrix A.
    !! Notice, it is assumed that matrix A has full rank.
    !!
    !! @param[in] m The number of equations (rows in matrix A).
    !! @param[in] n The number of unknowns (columns in matrix A).
    !! @param[in] k The number of columns in the right-hand-side matrix.
    !! @param[in,out] a On input, the M-by-N matrix A.  On output, if M >= N,
    !!  the QR factorization of A in the form as output by qr_factor; else,
    !!  if M < N, the LQ factorization of A.
    !! @param[in] lda The leading dimension of matrix A.
    !! @param[in,out] b If M >= N, the M-by-NRHS matrix B.  On output, the first
    !!  N rows contain the N-by-NRHS solution matrix X.  If M < N, an
    !!  N-by-NRHS matrix with the first M rows containing the matrix B.  On
    !!  output, the N-by-NRHS solution matrix X.
    !! @param[in] ldb The leading dimension of matrix B.
    !!
    !! @return An error code.  The following codes are possible.
    !!  - LA_NO_ERROR: No error occurred.  Successful operation.
    !!  - LA_INVALID_INPUT_ERROR: Occurs if @p lda, or @p ldb is not correct.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
    !!      there is insufficient memory available.
    !!  - LA_INVALID_OPERATION_ERROR: Occurs if @p a is not of full rank.
    function la_solve_least_squares_cmplx(m, n, k, a, lda, b, ldb) &
            bind(C, name = "la_solve_least_squares_cmplx") result(flag)
        ! Arguments
        integer(c_int), intent(in), value :: m, n, k, lda, ldb
        complex(c_double), intent(inout) :: a(lda,*), b(ldb,*)
        integer(c_int) :: flag

        ! Local Variables
        type(errors) :: err
        integer(c_int) :: maxmn

        ! Error Checking
        maxmn = max(m, n)
        call err%set_exit_on_error(.false.)
        flag = LA_NO_ERROR
        if (lda < m .or. ldb < maxmn) then
            flag = LA_INVALID_INPUT_ERROR
            return
        end if

        ! Process
        call solve_least_squares(a(1:m,1:n), b(1:maxmn,1:k), err = err)
        if (err%has_error_occurred()) then
            flag = err%get_error_flag()
            return
        end if
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the inverse of a square matrix.
    !!
    !! @param[in] n The dimension of matrix A.
    !! @param[in,out] a On input, the N-by-N matrix to invert.  On output, the
    !!  inverted matrix.
    !! @param[in] lda The leading dimension of matrix A.
    !!
    !! @return An error code.  The following codes are possible.
    !!  - LA_NO_ERROR: No error occurred.  Successful operation.
    !!  - LA_INVALID_INPUT_ERROR: Occurs if @p lda is not correct.
    !!  - LA_SINGULAR_MATRIX_ERROR: Occurs if the input matrix is singular.
    function la_inverse(n, a, lda) bind(C, name = "la_inverse") result(flag)
        ! Arguments
        integer(c_int), intent(in), value :: n, lda
        real(c_double), intent(inout) :: a(lda,*)
        integer(c_int) :: flag

        ! Local Variables
        type(errors) :: err

        ! Error Checking
        call err%set_exit_on_error(.false.)
        flag = LA_NO_ERROR
        if (lda < n) then
            flag = LA_INVALID_INPUT_ERROR
            return
        end if

        ! Process
        call mtx_inverse(a(1:n,1:n), err = err)
        if (err%has_error_occurred()) then
            flag = err%get_error_flag()
            return
        end if
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the inverse of a square matrix.
    !!
    !! @param[in] n The dimension of matrix A.
    !! @param[in,out] a On input, the N-by-N matrix to invert.  On output, the
    !!  inverted matrix.
    !! @param[in] lda The leading dimension of matrix A.
    !!
    !! @return An error code.  The following codes are possible.
    !!  - LA_NO_ERROR: No error occurred.  Successful operation.
    !!  - LA_INVALID_INPUT_ERROR: Occurs if @p lda is not correct.
    !!  - LA_SINGULAR_MATRIX_ERROR: Occurs if the input matrix is singular.
    function la_inverse_cmplx(n, a, lda) bind(C, name = "la_inverse_cmplx") &
            result(flag)
        ! Arguments
        integer(c_int), intent(in), value :: n, lda
        complex(c_double), intent(inout) :: a(lda,*)
        integer(c_int) :: flag

        ! Local Variables
        type(errors) :: err

        ! Error Checking
        call err%set_exit_on_error(.false.)
        flag = LA_NO_ERROR
        if (lda < n) then
            flag = LA_INVALID_INPUT_ERROR
            return
        end if

        ! Process
        call mtx_inverse(a(1:n,1:n), err = err)
        if (err%has_error_occurred()) then
            flag = err%get_error_flag()
            return
        end if
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the Moore-Penrose pseudo-inverse of an M-by-N matrix by
    !! means of singular value decomposition.
    !!
    !! @param[in] m The number of rows in the matrix.
    !! @parma[in] n The number of columns in the matrix.
    !! @param[in,out] a On input, the M-by-N matrix to invert.  The matrix is
    !!  overwritten on output.
    !! @param[in] lda The leading dimension of matrix A.
    !! @param[out] ainv The N-by-M matrix where the pseudo-inverse of @p a
    !!  will be written.
    !! @param[in] ldai The leading dimension of matrix AINV.
    !!
    !! @return An error code.  The following codes are possible.
    !!  - LA_NO_ERROR: No error occurred.  Successful operation.
    !!  - LA_INVALID_INPUT_ERROR: Occurs if @p lda, or @p ldai is not correct.
    function la_pinverse(m, n, a, lda, ainv, ldai) &
            bind(C, name = "la_pinverse") result(flag)
        ! Arguments
        integer(c_int), intent(in), value :: m, n, lda, ldai
        real(c_double), intent(inout) :: a(lda,*)
        real(c_double), intent(out) :: ainv(ldai,*)
        integer(c_int) :: flag

        ! Local Variables
        type(errors) :: err

        ! Error Checking
        call err%set_exit_on_error(.false.)
        flag = LA_NO_ERROR
        if (lda < m .or. ldai < n) then
            flag = LA_INVALID_INPUT_ERROR
            return
        end if

        ! Process
        call mtx_pinverse(a(1:m,1:n), ainv(1:n,1:m), err = err)
        if (err%has_error_occurred()) then
            flag = err%get_error_flag()
            return
        end if
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the Moore-Penrose pseudo-inverse of an M-by-N matrix by
    !! means of singular value decomposition.
    !!
    !! @param[in] m The number of rows in the matrix.
    !! @parma[in] n The number of columns in the matrix.
    !! @param[in,out] a On input, the M-by-N matrix to invert.  The matrix is
    !!  overwritten on output.
    !! @param[in] lda The leading dimension of matrix A.
    !! @param[out] ainv The N-by-M matrix where the pseudo-inverse of @p a
    !!  will be written.
    !! @param[in] ldai The leading dimension of matrix AINV.
    !!
    !! @return An error code.  The following codes are possible.
    !!  - LA_NO_ERROR: No error occurred.  Successful operation.
    !!  - LA_INVALID_INPUT_ERROR: Occurs if @p lda, or @p ldai is not correct.
    function la_pinverse_cmplx(m, n, a, lda, ainv, ldai) &
            bind(C, name = "la_pinverse_cmplx") result(flag)
        ! Arguments
        integer(c_int), intent(in), value :: m, n, lda, ldai
        complex(c_double), intent(inout) :: a(lda,*)
        complex(c_double), intent(out) :: ainv(ldai,*)
        integer(c_int) :: flag

        ! Local Variables
        type(errors) :: err

        ! Error Checking
        call err%set_exit_on_error(.false.)
        flag = LA_NO_ERROR
        if (lda < m .or. ldai < n) then
            flag = LA_INVALID_INPUT_ERROR
            return
        end if

        ! Process
        call mtx_pinverse(a(1:m,1:n), ainv(1:n,1:m), err = err)
        if (err%has_error_occurred()) then
            flag = err%get_error_flag()
            return
        end if
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the eigenvalues, and optionally the eigenvectors of a
    !! real, symmetric matrix.
    !!
    !! @param[in] vecs Set to true to compute the eigenvectors as well as the
    !!  eigenvalues; else, set to false to just compute the eigenvalues.
    !! @param[in] n The dimension of the matrix.
    !! @param[in,out] a On input, the N-by-N symmetric matrix on which to
    !!  operate.  On output, and if @p vecs is set to true, the matrix will
    !!  contain the eigenvectors (one per column) corresponding to each
    !!  eigenvalue in @p vals.  If @p vecs is set to false, the lower triangular
    !!  portion of the matrix is overwritten.
    !! @param[in] lda The leading dimension of matrix A.
    !! @param[out] vals An N-element array that will contain the eigenvalues
    !!  sorted into ascending order.
    !!
    !! @return An error code.  The following codes are possible.
    !!  - LA_NO_ERROR: No error occurred.  Successful operation.
    !!  - LA_INVALID_INPUT_ERROR: Occurs if @p lda is not correct.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
    !!      there is insufficient memory available.
    !!  - LA_CONVERGENCE_ERROR: Occurs if the algorithm failed to converge.
    function la_eigen_symm(vecs, n, a, lda, vals) &
            bind(C, name = "la_eigen_symm") result(flag)
        ! Arguments
        logical(c_bool), intent(in), value :: vecs
        integer(c_int), intent(in), value :: n, lda
        real(c_double), intent(inout) :: a(lda,*)
        real(c_double), intent(out) :: vals(*)
        integer(c_int) :: flag

        ! Local Variables
        type(errors) :: err

        ! Error Checking
        call err%set_exit_on_error(.false.)
        flag = LA_NO_ERROR
        if (lda < n) then
            flag = LA_INVALID_INPUT_ERROR
            return
        end if

        ! Process
        call eigen(logical(vecs), a(1:n,1:n), vals(1:n), err = err)
        if (err%has_error_occurred()) then
            flag = err%get_error_flag()
            return
        end if
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the eigenvalues, and optionally the right eigenvectors of
    !! a square matrix.
    !!
    !! @param[in] vecs Set to true to compute the eigenvectors as well as the
    !!  eigenvalues; else, set to false to just compute the eigenvalues.
    !! @param[in] n The dimension of the matrix.
    !! @param[in,out] a On input, the N-by-N matrix on which to operate.  On
    !!  output, the contents of this matrix are overwritten.
    !! @param[in] lda The leading dimension of matrix A.
    !! @param[out] vals An N-element array containing the eigenvalues of the
    !!  matrix.  The eigenvalues are not sorted.
    !! @param[out] v An N-by-N matrix where the right eigenvectors will be
    !!  written (one per column).
    !!
    !! @return An error code.  The following codes are possible.
    !!  - LA_NO_ERROR: No error occurred.  Successful operation.
    !!  - LA_INVALID_INPUT_ERROR: Occurs if @p lda or @p ldv is not correct.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
    !!      there is insufficient memory available.
    !!  - LA_CONVERGENCE_ERROR: Occurs if the algorithm failed to converge.
    function la_eigen_asymm(vecs, n, a, lda, vals, v, ldv) &
            bind(C, name = "la_eigen_asymm") result(flag)
        ! Arguments
        logical(c_bool), intent(in), value :: vecs
        integer(c_int), intent(in), value :: n, lda, ldv
        real(c_double), intent(inout) :: a(lda,*)
        complex(c_double), intent(out) :: vals(*), v(ldv,*)
        integer(c_int) :: flag

        ! Local Variables
        type(errors) :: err

        ! Error Checking
        call err%set_exit_on_error(.false.)
        flag = LA_NO_ERROR
        if (vecs) then
            if (lda < n .or. ldv < n) then
                flag = LA_INVALID_INPUT_ERROR
                return
            end if
        else
            if (lda < n) then
                flag = LA_INVALID_INPUT_ERROR
                return
            end if
        end if

        ! Process
        if (vecs) then
            call eigen(a(1:n,1:n), vals(1:n), v(1:n,1:n), err = err)
        else
            call eigen(a(1:n,1:n), vals(1:n))
        end if
        if (err%has_error_occurred()) then
            flag = err%get_error_flag()
            return
        end if
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the eigenvalues, and optionally the right eigenvectors of
    !! a square matrix assuming the structure of the eigenvalue problem is
    !! A*X = lambda*B*X.
    !!
    !! @param[in] vecs Set to true to compute the eigenvectors as well as the
    !!  eigenvalues; else, set to false to just compute the eigenvalues.
    !! @param[in] n The dimension of the matrix.
    !! @param[in,out] a On input, the N-by-N matrix A.  On output, the contents
    !!  of this matrix are overwritten.
    !! @param[in] lda The leading dimension of matrix A.
    !! @param[in,out] b On input, the N-by-N matrix B.  On output, the contents
    !!  of this matrix are overwritten.
    !! @param[in] ldb The leading dimension of matrix B.
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
    !! @param[out] v An N-by-N matrix where the right eigenvectors will be
    !!  written (one per column).
    !!
    !! @return An error code.  The following codes are possible.
    !!  - LA_NO_ERROR: No error occurred.  Successful operation.
    !!  - LA_INVALID_INPUT_ERROR: Occurs if @p lda or @p ldv is not correct.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
    !!      there is insufficient memory available.
    !!  - LA_CONVERGENCE_ERROR: Occurs if the algorithm failed to converge.
    function la_eigen_gen(vecs, n, a, lda, b, ldb, alpha, beta, v, ldv) &
            bind(C, name = "la_eigen_gen") result(flag)
        ! Arguments
        logical(c_bool), intent(in), value :: vecs
        integer(c_int), intent(in), value :: n, lda, ldb, ldv
        real(c_double), intent(inout) :: a(lda,*), b(ldb,*)
        real(c_double), intent(out) :: beta(*)
        complex(c_double), intent(out) :: alpha(*), v(ldv,*)
        integer(c_int) :: flag

        ! Local Variables
        type(errors) :: err

        ! Error Checking
        call err%set_exit_on_error(.false.)
        flag = LA_NO_ERROR
        if (vecs) then
            if (lda < n .or. ldb < n .or. ldv < n) then
                flag = LA_INVALID_INPUT_ERROR
                return
            end if
        else
            if (lda < n .or. ldb < n) then
                flag = LA_INVALID_INPUT_ERROR
                return
            end if
        end if

        ! Process
        if (vecs) then
            call eigen(a(1:n,1:n), b(1:n,1:n), alpha(1:n), beta(1:n), &
                v(1:n,1:n), err = err)
        else
            call eigen(a(1:n,1:n), b(1:n,1:n), alpha(1:n), beta(1:n), err = err)
        end if
        if (err%has_error_occurred()) then
            flag = err%get_error_flag()
            return
        end if
    end function
    
! ------------------------------------------------------------------------------
    !> @brief Computes the eigenvalues, and optionally the right eigenvectors of
    !! a square matrix.
    !!
    !! @param[in] vecs Set to true to compute the eigenvectors as well as the
    !!  eigenvalues; else, set to false to just compute the eigenvalues.
    !! @param[in] n The dimension of the matrix.
    !! @param[in,out] a On input, the N-by-N matrix on which to operate.  On
    !!  output, the contents of this matrix are overwritten.
    !! @param[in] lda The leading dimension of matrix A.
    !! @param[out] vals An N-element array containing the eigenvalues of the
    !!  matrix.  The eigenvalues are not sorted.
    !! @param[out] v An N-by-N matrix where the right eigenvectors will be
    !!  written (one per column).
    !!
    !! @return An error code.  The following codes are possible.
    !!  - LA_NO_ERROR: No error occurred.  Successful operation.
    !!  - LA_INVALID_INPUT_ERROR: Occurs if @p lda or @p ldv is not correct.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
    !!      there is insufficient memory available.
    !!  - LA_CONVERGENCE_ERROR: Occurs if the algorithm failed to converge.
    function la_eigen_cmplx(vecs, n, a, lda, vals, v, ldv) &
            bind(C, name = "la_eigen_cmplx") result(flag)
        ! Arguments
        logical(c_bool), intent(in), value :: vecs
        integer(c_int), intent(in), value :: n, lda, ldv
        complex(c_double), intent(inout) :: a(lda,*)
        complex(c_double), intent(out) :: vals(*), v(ldv,*)
        integer(c_int) :: flag

        ! Local Variables
        type(errors) :: err

        ! Error Checking
        call err%set_exit_on_error(.false.)
        flag = LA_NO_ERROR
        if (vecs) then
            if (lda < n .or. ldv < n) then
                flag = LA_INVALID_INPUT_ERROR
                return
            end if
        else
            if (lda < n) then
                flag = LA_INVALID_INPUT_ERROR
                return
            end if
        end if

        ! Process
        if (vecs) then
            call eigen(a(1:n,1:n), vals(1:n), v(1:n,1:n), err = err)
        else
            call eigen(a(1:n,1:n), vals(1:n))
        end if
        if (err%has_error_occurred()) then
            flag = err%get_error_flag()
            return
        end if
    end function

! ------------------------------------------------------------------------------
    !> @brief A sorting routine specifically tailored for sorting of eigenvalues
    !! and their associated eigenvectors using a quick-sort approach.
    !!
    !! @param[in] ascend
    !! @param[in] n The number of eigenvalues.
    !! @param[in,out] vals On input, an N-element array containing the
    !!  eigenvalues.  On output, the sorted eigenvalues.
    !! @param[in,out] vecs On input, an N-by-N matrix containing the
    !!  eigenvectors associated with @p vals (one vector per column).  On
    !!  output, the sorted eigenvector matrix.
    !! @param[in] ldv The leading dimension of @p vecs.
    !!
    !! @return An error code.  The following codes are possible.
    !!  - LA_NO_ERROR: No error occurred.  Successful operation.
    !!  - LA_INVALID_INPUT_ERROR: Occurs if @p ldv is not correct.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
    !!      there is insufficient memory available.
    function la_sort_eigen(ascend, n, vals, vecs, ldv) &
            bind(C, name = "la_sort_eigen") result(flag)
        ! Arguments
        logical(c_bool), intent(in), value :: ascend
        integer(c_int), intent(in), value :: n, ldv
        real(c_double), intent(inout) :: vals(*), vecs(ldv,*)
        integer(c_int) :: flag

        ! Local Variables
        type(errors) :: err

        ! Error Checking
        call err%set_exit_on_error(.false.)
        flag = LA_NO_ERROR
        if (ldv < n) then
            flag = LA_INVALID_INPUT_ERROR
            return
        end if

        ! Process
        call sort(vals(1:n), vecs(1:n,1:n), logical(ascend), err = err)
        if (err%has_error_occurred()) then
            flag = err%get_error_flag()
            return
        end if
    end function

! ------------------------------------------------------------------------------
    !> @brief A sorting routine specifically tailored for sorting of eigenvalues
    !! and their associated eigenvectors using a quick-sort approach.
    !!
    !! @param[in] ascend
    !! @param[in] n The number of eigenvalues.
    !! @param[in,out] vals On input, an N-element array containing the
    !!  eigenvalues.  On output, the sorted eigenvalues.
    !! @param[in,out] vecs On input, an N-by-N matrix containing the
    !!  eigenvectors associated with @p vals (one vector per column).  On
    !!  output, the sorted eigenvector matrix.
    !! @param[in] ldv The leading dimension of @p vecs.
    !!
    !! @return An error code.  The following codes are possible.
    !!  - LA_NO_ERROR: No error occurred.  Successful operation.
    !!  - LA_INVALID_INPUT_ERROR: Occurs if @p ldv is not correct.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
    !!      there is insufficient memory available.
    function la_sort_eigen_cmplx(ascend, n, vals, vecs, ldv) &
            bind(C, name = "la_sort_eigen_cmplx") result(flag)
        ! Arguments
        logical(c_bool), intent(in), value :: ascend
        integer(c_int), intent(in), value :: n, ldv
        complex(c_double), intent(inout) :: vals(*), vecs(ldv,*)
        integer(c_int) :: flag

        ! Local Variables
        type(errors) :: err

        ! Error Checking
        call err%set_exit_on_error(.false.)
        flag = LA_NO_ERROR
        if (ldv < n) then
            flag = LA_INVALID_INPUT_ERROR
            return
        end if

        ! Process
        call sort(vals(1:n), vecs(1:n,1:n), logical(ascend), err = err)
        if (err%has_error_occurred()) then
            flag = err%get_error_flag()
            return
        end if
    end function

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end module
