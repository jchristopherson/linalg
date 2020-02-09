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
    !!@brief Computes the matrix operation C = alpha * op(A) * op(B) + beta * C.
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
        if (transb) tb = "N"

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
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
    !!      there is insufficient memory available.
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
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
    !!      there is insufficient memory available.
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

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end module
