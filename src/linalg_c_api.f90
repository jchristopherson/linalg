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
        real(c_double), intent(inout) :: a(n,*)
        real(c_double), intent(out) :: u(n,*), p(n,*)
        integer(c_int), intent(in) :: ipvt(n)
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
        complex(c_double), intent(inout) :: a(n,*)
        complex(c_double), intent(out) :: u(n,*)
        real(c_double), intent(out) :: p(n,*)
        integer(c_int), intent(in) :: ipvt(n)
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
        integer(c_int), intent(inout) :: jpvt(n)
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
        integer(c_int), intent(inout) :: jpvt(n)
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

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end module
