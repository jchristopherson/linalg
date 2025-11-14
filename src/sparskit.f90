!> @brief An interface to the SPARSKIT library available at 
!! https://www-users.cse.umn.edu/~saad/software/SPARSKIT/.
module sparskit
    implicit none

    ! BLASSM.F
    interface
        !> @brief Computes the matrix product: C = A * B.
        !!
        !! @param[in] nrow The row dimension of matrices A & C.
        !! @param[in] ncol The column dimension of matrices B & C.
        !! @param[in] job Set to 0 to compute only the structure (JC & IC);
        !!  else, set to any non-zero value.
        !! @param[in] a The non-zero elements of matrix A.
        !! @param[in] ja The column indices of matrix A.
        !! @param[in] ia The index in A where the requested row starts.
        !! @param[in] b The non-zero elements of matrix B.
        !! @param[in] jb The column indices of matrix B.
        !! @param[in] ib The index in B where the requested row starts.
        !! @param[out] c The non-zero elements of matrix C.
        !! @param[out] jc The column indices of matrix C.
        !! @param[out] ic The index in C where the requested row starts.
        !! @param[in] nzmax The length of arrays C & JC.  The routine will stop
        !!  if the results matrix C has a number of elements that exceeds NZMAX.
        !! @param[out] iw A workspace array with a length equal to the number of
        !!  of columns in matrix C.
        !! @param[out] ierr An error message indicator.
        !!  * 0: Normal return
        !!  * .gt. 0: Routine failed in row I with IERR = I because the number
        !!     of elements in C exceeds NZMAX.
        pure subroutine amub(nrow, ncol, job, a, ja, ia, b, jb, ib, c, jc, ic, &
            nzmax, iw, ierr)
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: nrow, ncol, job, nzmax
            integer(int32), intent(in) :: ja(*), ia(nrow+1), jb(*), ib(*)
            integer(int32), intent(out) :: jc(*), ic(*), iw(ncol), ierr
            real(real64), intent(in) :: a(*), b(*)
            real(real64), intent(out) :: c(*)
        end subroutine

        !> @brief Computes the matrix sum: C = A + B, where the matrices are
        !! given in CSR format.
        !!
        !! @param[in] nrow The number of rows in the matrices.
        !! @param[in] ncol The number of columns in the matrices.
        !! @param[in] a The non-zero elements of matrix A.
        !! @param[in] ja The column indices of matrix A.
        !! @param[in] ia The index in A where the requested row starts.
        !! @param[in] b The non-zero elements of matrix B.
        !! @param[in] jb The column indices of matrix B.
        !! @param[in] ib The index in B where the requested row starts.
        !! @param[out] c The non-zero elements of matrix C.
        !! @param[out] jc The column indices of matrix C.
        !! @param[out] ic The index in C where the requested row starts.
        !! @param[in] nzmax The length of arrays C & JC.  The routine will stop
        !!  if the results matrix C has a number of elements that exceeds NZMAX.
        !! @param[out] iw A workspace array with a length equal to the number of
        !!  of columns in matrix A.
        !! @param[out] ierr An error message indicator.
        !!  * 0: Normal return
        !!  * .gt. 0: Routine failed in row I with IERR = I because the number
        !!     of elements in C exceeds NZMAX.
        pure subroutine aplb(nrow, ncol, job, a, ja, ia, b, jb, ib, c, jc, ic, &
            nzmax, iw, ierr)
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: nrow, ncol, job, nzmax
            integer(int32), intent(in) :: ja(*), ia(nrow+1), jb(*), ib(nrow+1)
            real(real64), intent(in) :: a(*), b(*)
            real(real64), intent(out) :: c(*)
            integer(int32), intent(out) :: jc(*), ic(nrow+1), iw(ncol), ierr
        end subroutine

        !> @brief Computes the matrix sum: C = A + s * B, where the matrices
        !! are given in CSR format.
        !!
        !! @param[in] nrow The number of rows in the matrices.
        !! @param[in] ncol The number of columns in the matrices.
        !! @param[in] a The non-zero elements of matrix A.
        !! @param[in] ja The column indices of matrix A.
        !! @param[in] ia The index in A where the requested row starts.
        !! @param[in] s The scalar multiplier.
        !! @param[in] b The non-zero elements of matrix B.
        !! @param[in] jb The column indices of matrix B.
        !! @param[in] ib The index in B where the requested row starts.
        !! @param[out] c The non-zero elements of matrix C.
        !! @param[out] jc The column indices of matrix C.
        !! @param[out] ic The index in C where the requested row starts.
        !! @param[in] nzmax The length of arrays C & JC.  The routine will stop
        !!  if the results matrix C has a number of elements that exceeds NZMAX.
        !! @param[out] iw A workspace array with a length equal to the number of
        !!  of columns in matrix A.
        !! @param[out] ierr An error message indicator.
        !!  * 0: Normal return
        !!  * .gt. 0: Routine failed in row I with IERR = I because the number
        !!     of elements in C exceeds NZMAX.
        pure subroutine aplsb(nrow, ncol, a, ja, ia, s, b, jb, ib, c, jc, ic, &
            nzmax, iw, ierr)
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: nrow, ncol, nzmax
            integer(int32), intent(in) :: ja(*), ia(nrow+1), jb(*), ib(nrow+1)
            integer(int32), intent(out) :: ierr
            integer(int32), intent(out) :: jc(*), ic(nrow+1), iw(ncol)
            real(real64), intent(in) :: s
            real(real64), intent(in) :: a(*), b(*)
            real(real64), intent(out) :: c(*)
        end subroutine
    end interface

    ! FORMATS.F
    interface
        !> @brief Converts a CSR matrix into a CSC matrix (transposition).
        !!
        !! @param[in] n The number of rows in the CSR matrix.
        !! @param[in] n2 The number of columns in the CSC matrix.
        !! @param[in] job Fill the values (job == 1) or only the pattern 
        !!  (job /= 1).
        !! @param[in] ipos Starting position of A0 in JA0.
        !! @param[in] a The non-zero elements of matrix A.
        !! @param[in] ja The column indices of matrix A.
        !! @param[in] ia The index in A where the requested row starts.
        !! @param[out] a0 The non-zero elements of the transposed array.
        !! @param[out] ja0 The size NNZ array containing the column indices.
        !! @param[out] ia0 The N+1 size array containing the column starts.
        pure subroutine csrcsc2(n, n2, job, ipos, a, ja, ia, a0, ja0, ia0)
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: n, n2, job, ipos, ja(*), ia(n+1)
            integer(int32), intent(out) :: ja0(*), ia0(n2+1)
            real(real64), intent(in) :: a(*)
            real(real64), intent(out) :: a0(*)
        end subroutine

        !> @brief Converts the LINPACK, BLAS, LAPACK banded matrix format into
        !! a CSR format.
        !!
        !! @param[in] n The row dimension of the matrix.
        !! @param[in] abd The banded matrix.
        !! @param[in] nabd The leading dimension of @p abd.
        !! @param[in] lowd The row index where the lowest diagonal (leftmost) of
        !!  A is located.  LINPACK uses LOWD = 2 * ML + MU + 1.
        !! @param[in] ml The bandwidth of the strict lower part of A.
        !! @param[in] mu The bandwidth of the strict upper part of A.
        !! @param[out] a The non-zero elements of matrix A.
        !! @param[out] ja The column indices of matrix A.
        !! @param[out] ia The index in A where the requested row starts.
        !! @param[in] len The length of @p a and @p ja.
        !! @param[out] ierr Error message output.
        !!  * 0: Normal return.
        !!  * -1: Invalid @p lowd value.
        !!  * Positive Valued: Not enough storage in @p a and @p ja.
        pure subroutine bndcsr(n, abd, nabd, lowd, ml, mu, a, ja, ia, len, ierr)
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: n, nabd, lowd, ml, mu, len
            real(real64), intent(in) :: abd(nabd,*)
            real(real64), intent(out) :: a(*)
            integer(int32), intent(out) :: ia(n+1), ja(*), ierr
        end subroutine

        !> @brief Converts a CSR matrix to an MSR matrix.
        !!
        !! @param[in] n The row dimension of the matrix.
        !! @param[in,out] a The non-zero elements of matrix A.
        !! @param[in,out] ja The column indices of matrix A.
        !! @param[in,out] ia The index in A where the requested row starts.
        !! @param[out] ao An NNZ-element array containing the non-zero elements
        !!  for the MSR matrix.
        !! @param[out] jao An NNZ-element index tracking array for the MSR
        !!  matrix.
        !! @param[out] wk An N-element workspace array.
        !! @param[out] iwk An N+1 element workspace array.
        pure subroutine csrmsr(n, a, ja, ia, ao, jao, wk, iwk)
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: n
            integer(int32), intent(inout) :: ja(*), ia(n+1)
            integer(int32), intent(out) :: jao(*), iwk(n+1)
            real(real64), intent(inout) :: a(*)
            real(real64), intent(out) :: ao(*), wk(n)
        end subroutine

        !> @brief Converts and MSR matrix to a CSR matrix.
        !!
        !! @param[in] n The row dimension of the matrix.
        !! @param[in] a An NNZ-element array containing the non-zero elements
        !!  for the MSR matrix.
        !! @param[in] ja An NNZ-element index tracking array for the MSR
        !!  matrix.
        !! @param[out] ao The non-zero elements of matrix A.
        !! @param[out] jao The column indices of matrix A.
        !! @param[out] iao The index in A where the requested row starts.
        !! @param[out] wk An N-element workspace array.
        !! @param[out] iwk An N+1 element workspace array.
        pure subroutine msrcsr(n, a, ja, ao, jao, iao, wk, iwk)
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: n, ja(*)
            integer(int32), intent(out) :: jao(*), iao(n+1), iwk(n+1)
            real(real64), intent(in) :: a(*)
            real(real64), intent(out) :: ao(*), wk(n)
        end subroutine

        !> @brief Converte a matrix stored in coordinate format to CSR format.
        !!
        !! @param[in] nrow The number of rows in the matrix.
        !! @param[in] nnz The number of non-zero elements in the matrix.
        !! @param[in] a An NNZ-element array containing the non-zero elements
        !!  of the matrix.
        !! @param[in,out] ir An NNZ-element array containing the row indices of
        !!  each non-zero element.
        !! @param[in] jc An NNZ-element array containing the column indices of
        !!  each non-zero element.
        !! @param[out] ao The non-zero elements of matrix A.
        !! @param[out] jao The column indices of matrix A.
        !! @param[out] iao The index in A where the requested row starts.
        pure subroutine coocsr(nrow, nnz, a, ir, jc, ao, jao, iao)
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: nrow, nnz, jc(*)
            integer(int32), intent(inout) :: ir(*)
            real(real64), intent(in) :: a(*)
            integer(int32), intent(out) :: jao(*), iao(*)
            real(real64), intent(out) :: ao(*)
        end subroutine
    end interface

    ! UNARY.F
    interface
        !> @brief Gets element A(i,j) of matrix A for any pair (i,j).
        !!
        !! @param[in] i The row index.
        !! @param[in] j The column index.
        !! @param[in] a The non-zero elements of matrix A.
        !! @param[in] ja The column indices of matrix A.
        !! @param[in] ia The index in A where the requested row starts.
        !! @param[out] iadd The address of element A(i,j) in arrays A & JA, if 
        !!  found; else, zero if not found.
        !! @param[in] sorted Indicates whether the matrix is known to be sorted.
        !!
        !! @return The requested value.
        function getelm(i, j, a, ja, ia, iadd, sorted) result(rst)
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: i, j, ia(*), ja(*)
            real(real64), intent(in) :: a(*)
            integer(int32), intent(out) :: iadd
            logical, intent(in) :: sorted
            real(real64) :: rst
        end function

        !> @brief Extracts the diagonal from a matrix.
        !!
        !! @param[in] nrow The number of rows.
        !! @param[in] ncol The number of columns.
        !! @param[in] job Set to 0 to not alter @p a, @p ja, and @p ia; else,
        !! set to a non-zero value to perform this as an in-place operation.
        !! @param[out] len The number of non-zero elements found in @p diag.
        !! @param[out] idiag An array of length @p len containing the original
        !!  positions in the original arrays @p a and @p ja of the diagonal
        !!  elements collected in diagl.
        !! @param[in] ioff The offset of the wanted diagonal.
        pure subroutine getdia(nrow, ncol, job, a, ja, ia, len, diag, idiag, ioff)
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: nrow, ncol, job, ja(*), ia(*), ioff
            integer(int32), intent(out) :: len, idiag(*)
            real(real64), intent(in) :: a(*)
            real(real64), intent(out) :: diag(*)
        end subroutine

        !> @brief Sorces the elements of a CSR matrix in increasing order of 
        !! their column indices within each row.
        !!
        !! @param[in] n The number of rows in the matrix.
        !! @param[in,out] a The non-zero values.
        !! @param[in,out] ja An array of column indices of the elements in A.
        !! @param[in] ia An array of pointers to the rows.
        !! @param[in] values Idicates whether A must also be permuted.  If
        !!  false, A can be a dummy array.
        pure subroutine csort(n, a, ja, ia, values)
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: n
            real(real64), intent(inout) :: a(*)
            integer(int32), intent(inout) :: ja(*)
            integer(int32), intent(in) :: ia(*)
            logical, intent(in) :: values
        end subroutine

        !> @breif Cleans up a CSR matrix.
        !!
        !! @param[in] job The job to perform.
        !!  - 0: Nothing is done
        !!  - 1: Eliminate duplicate entries and zero entries.
        !!  - 2: Eliminate duplicate entries and perform partial ordering.
        !!  - 3: Eliminate duplicate entries and sort the entries in increasing
        !!      order of column indices.
        !! @param[in] value2 0 if the matrix is pattern only (A is not touched),
        !!  or 1 if the matrix has values.
        !! @param[in] nrow The number of rows in the matrix.
        !! @param[in,out] a The non-zero values.
        !! @param[in,out] ja An array of column indices of the elements in A.
        !! @param[in,out] ia An array of pointers to the rows.
        !! @param[out] indu An NROW array containing pointers to the beginning 
        !!  of the upper triangular portion if job > 1.
        !! @param[out] iwk An NROW+1 element workspace array.

        pure subroutine clncsr(job, value2, nrow, a, ja, ia, indu, iwk)
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: job, value2, nrow
            real(real64), intent(inout) :: a(*)
            integer(int32), intent(inout) :: ja(*), ia(*)
            integer(int32), intent(inout) :: indu(*), iwk(*)
        end subroutine
    end interface

    ! ILUT.F
    interface
        !> @brief Computes the incomplete LU factorization of a sparse matrix
        !! in CSR format using a dual truncation mechanism.
        !!
        !! @param[in] n The row dimension of the matrix.
        !! @param[in] a The non-zero elements of matrix A.
        !! @param[in] ja The column indices of matrix A.
        !! @param[in] ia The index in A where the requested row starts.
        !! @param[in] lfil The fill-in parameter.  Each row of L and each row
        !!  of U will have a maximum of @p lfil elements, excluding the 
        !!  diagonal element.  @p lfil must be greater than or equal to zero.
        !! @param[in] droptol The threshold for dropping small terms in the
        !!  factorization.
        !! @param[out] alu The factored matrix stored in Modified Sparse Row
        !!  (MSR) format containing the L and U factors together.  The diagonal,
        !!  stored in ALU(1:N), is inverted.  Each i-th row of the ALU, JLU
        !!  matrix contains the i-th row of L, excluding the diagonal entry,
        !!  followed by the i-th row of U.
        !! @param[out] jlu The column indices for the factored matrix.
        !! @param[out] ju An N-element array containing the pointers to the
        !!  beginning of each row of U in the factored matrix.
        !! @param[in] iwk The lengths of @p alu and @p jlu.
        !! @param[out] w An N+1 element workspace array.
        !! @param[out] jw A 2*N element workspace array.
        !! @param[out] ierr Error flag:
        !!  * 0: Successful return
        !!  * .gt. 0: Zero pivot encountered at step number IERR.
        !!  * -1: Input matrix is incorrect.  The elimination process generated
        !!      a row in L or U whose length is greater than N.
        !!  * -2: The matrix L overflows the output array.
        !!  * -3: The matrix U overflows the output array.
        !!  * -4: Illegal value for @P lfil.
        !!  * -5: Zero-valued row encountered.
        pure subroutine ilut(n, a, ja, ia, lfil, droptol, alu, jlu, ju, iwk, w, &
            jw, ierr)
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: n, ja(*), ia(n+1), lfil, iwk
            integer(int32), intent(out) :: jlu(*), ju(n), jw(2*n), ierr
            real(real64), intent(in) :: a(*), droptol
            real(real64), intent(out) :: alu(*), w(n+1)
        end subroutine

        !> @brief Computes the incomplete LU factorization of a sparse matrix
        !! in CSR format using a dual truncation mechanism and pivoting.
        !!
        !! @param[in] n The row dimension of the matrix.
        !! @param[in,out] a The non-zero elements of matrix A.  On output, the
        !!  columns are permuted.
        !! @param[in,out] ja The column indices of matrix A.  On output, the
        !!  columns are permuted.
        !! @param[in,out] ia The index in A where the requested row starts.  On 
        !!  output, the columns are permuted.
        !! @param[in] lfil The fill-in parameter.  Each row of L and each row
        !!  of U will have a maximum of @p lfil elements, excluding the 
        !!  diagonal element.  @p lfil must be greater than or equal to zero.
        !! @param[in] droptol The threshold for dropping small terms in the
        !!  factorization.
        !! @param[in] permtol A tolerance ratio used to determine whether or
        !!  not to permute two columns.  At step I, columns I and J are 
        !!  permuted when ABS(A(I,J)) * PERMTOL .GT. ABS(A(I,I)).  Good values
        !!  are typically between 0.1 to 0.01.
        !! @param[in] mbloc If desired, permuting can be done only within the
        !!  diagonal blocks of size MBLOC.  Useful for PDE problems with many
        !!  degrees of freedom.  If this feature is not required, simply set
        !!  MBLOC equal to N.
        !! @param[out] alu The factored matrix stored in Modified Sparse Row
        !!  (MSR) format containing the L and U factors together.  The diagonal,
        !!  stored in ALU(1:N), is inverted.  Each i-th row of the ALU, JLU
        !!  matrix contains the i-th row of L, excluding the diagonal entry,
        !!  followed by the i-th row of U.
        !! @param[out] jlu The column indices for the factored matrix.
        !! @param[out] ju An N-element array containing the pointers to the
        !!  beginning of each row of U in the factored matrix.
        !! @param[in] iwk The lengths of @p alu and @p jlu.
        !! @param[out] w An N+1 element workspace array.
        !! @param[out] jw A 2*N element workspace array.
        !! @param[out] iperm A 2*N element array containing the permutation
        !!  arrays.  IPERM(1:N) contains the old numbers of unknowns, and 
        !!  IPERM(N+1:) contains the new unknowns.
        !! @param[out] ierr Error flag:
        !!  * 0: Successful return
        !!  * .gt. 0: Zero pivot encountered at step number IERR.
        !!  * -1: Input matrix is incorrect.  The elimination process generated
        !!      a row in L or U whose length is greater than N.
        !!  * -2: The matrix L overflows the output array.
        !!  * -3: The matrix U overflows the output array.
        !!  * -4: Illegal value for @P lfil.
        !!  * -5: Zero-valued row encountered.
        !!
        !! @par Remarks
        !! To avoid permuting the solution vector arrays for each LU-solve, the
        !! matrix A is permuted on return.  Similarly for the U matrix.  To
        !! permute the matrix back to its original state, use the following
        !! code.
        !! @code{.f90}
        !! do k = ia(1), ia(n+1) - 1
        !!  ja(k) = iperm(ja(k))
        !! end do
        !! @endcode
        pure subroutine ilutp(n, a, ja, ia, lfil, droptol, permtol, mbloc, alu, &
            jlu, ju, iwk, w, jw, iperm, ierr)
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: n, lfil, iwk, mbloc
            integer(int32), intent(inout) :: ja(*), ia(n+1)
            integer(int32), intent(out) :: jlu(*), ju(n), jw(2*n), &
                iperm(2*n), ierr
            real(real64), intent(in) :: droptol, permtol
            real(real64), intent(inout) :: a(*)
            real(real64), intent(out) :: alu(*), w(n+1)
        end subroutine

        !> @brief Computes the incomplete LU factorization of a sparse matrix
        !! in CSR format with standard dropping strategy.
        !!
        !! @param[in] n The row dimension of the matrix.
        !! @param[in] a The non-zero elements of matrix A.
        !! @param[in] ja The column indices of matrix A.
        !! @param[in] ia The index in A where the requested row starts.
        !! @param[in] alph The diagonal compensation parameter.  If ALPH = 0,
        !!  the process is approximately equivalent to ILU with threshold; else,
        !!  if ALPH = 1, the process is approximately equivalent to MILU with
        !!  threshold.
        !! @param[in] tol The threshold parameter for dropping small terms in
        !!  the factorization.
        !! @param[out] alu The factored matrix stored in Modified Sparse Row
        !!  (MSR) format containing the L and U factors together.  The diagonal,
        !!  stored in ALU(1:N), is inverted.  Each i-th row of the ALU, JLU
        !!  matrix contains the i-th row of L, excluding the diagonal entry,
        !!  followed by the i-th row of U.
        !! @param[out] jlu The column indices for the factored matrix.
        !! @param[out] ju An N-element array containing the pointers to the
        !!  beginning of each row of U in the factored matrix.
        !! @param[in] iwk The lengths of @p alu and @p jlu.
        !! @param[out] w An N+1 element workspace array.
        !! @param[out] jw A 2*N element workspace array.
        !! @param[out] ierr Error flag:
        !!  * 0: Successful return
        !!  * .gt. 0: Zero pivot encountered at step number IERR.
        !!  * -1: Input matrix is incorrect.  The elimination process generated
        !!      a row in L or U whose length is greater than N.
        !!  * -2: Insufficient storage for the LU factors.
        !!  * -3: Zero-valued row encountered.
        pure subroutine ilud(n, a, ja, ia, alph, tol, alu, jlu, ju, iwk, w, jw, ierr)
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: n, iwk, ja(*), ia(n+1)
            integer(int32), intent(out) :: jlu(*), ju(n), jw(2*n), ierr
            real(real64), intent(in) :: a(*), alph, tol
            real(real64), intent(out) :: alu(*), w(2*n)
        end subroutine

        !> @brief Computes the incomplete LU factorization of a sparse matrix
        !! in CSR format with standard dropping strategy.
        !!
        !! @param[in] n The row dimension of the matrix.
        !! @param[in] a The non-zero elements of matrix A.
        !! @param[in] ja The column indices of matrix A.
        !! @param[in] ia The index in A where the requested row starts.
        !! @param[in] alph The diagonal compensation parameter.  If ALPH = 0,
        !!  the process is approximately equivalent to ILU with threshold; else,
        !!  if ALPH = 1, the process is approximately equivalent to MILU with
        !!  threshold.
        !! @param[in] droptol The threshold for dropping small terms in the
        !!  factorization.
        !! @param[in] permtol A tolerance ratio used to determine whether or
        !!  not to permute two columns.  At step I, columns I and J are 
        !!  permuted when ABS(A(I,J)) * PERMTOL .GT. ABS(A(I,I)).  Good values
        !!  are typically between 0.1 to 0.01.
        !! @param[in] mbloc If desired, permuting can be done only within the
        !!  diagonal blocks of size MBLOC.  Useful for PDE problems with many
        !!  degrees of freedom.  If this feature is not required, simply set
        !!  MBLOC equal to N.
        !! @param[out] alu The factored matrix stored in Modified Sparse Row
        !!  (MSR) format containing the L and U factors together.  The diagonal,
        !!  stored in ALU(1:N), is inverted.  Each i-th row of the ALU, JLU
        !!  matrix contains the i-th row of L, excluding the diagonal entry,
        !!  followed by the i-th row of U.
        !! @param[out] jlu The column indices for the factored matrix.
        !! @param[out] ju An N-element array containing the pointers to the
        !!  beginning of each row of U in the factored matrix.
        !! @param[in] iwk The lengths of @p alu and @p jlu.
        !! @param[out] w An N+1 element workspace array.
        !! @param[out] jw A 2*N element workspace array.
        !! @param[out] ierr Error flag:
        !!  * 0: Successful return
        !!  * .gt. 0: Zero pivot encountered at step number IERR.
        !!  * -1: Input matrix is incorrect.  The elimination process generated
        !!      a row in L or U whose length is greater than N.
        !!  * -2: Insufficient storage for the LU factors.
        !!  * -3: Zero-valued row encountered.
        pure subroutine iludp(n, a, ja, ia, alph, droptol, permtol, mbloc, alu, &
            jlu, ju, iwk, w, jw, iperm, ierr)
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: n, iwk, mbloc
            integer(int32), intent(inout) :: ja(*), ia(n+1)
            integer(int32), intent(out) :: jlu(*), ju(n), jw(2*n), iperm(2*n), &
                ierr
            real(real64), intent(in) :: alph, droptol, permtol
            real(real64), intent(inout) :: a(*)
            real(real64), intent(out) :: alu(*), w(2*n)
        end subroutine

        !> @brief An ILUT preconditioned GMRES algorithm.  This routine utilizes
        !! the L and U matrices generated by the ILUT routine to precondition
        !! the GMRES algorithm.  The stopping criteria utilized is based simply 
        !! on reducing the residual norm to the requested tolerance.
        !!
        !! @param[in] n The row dimension of the matrix.
        !! @param[in] im The size of the Krylov subspace.  This value should
        !!  not exceed 50.
        !! @param[in,out] rhs The N-element right-hand-side vector.  On output,
        !!  the contents of this array are overwritten.
        !! @param[in,out] sol On input, the N-element solution estimate.  On
        !!  output, the computed solution.
        !! @param[out] vv An N-by-IM+1 workspace matrix.
        !! @param[in] eps The convergence tolerance against which the norm of
        !!  the residual is checked.
        !! @param[in] maxits The maximum number of iterations to allow.
        !! @param[in] iout The device output number for printing intermediate
        !!  results.  Set to a value less than or equal to zero to suppress
        !!  printing.
        !! @param[in] aa The non-zero elements of matrix A.
        !! @param[in] ja The column indices of matrix A.
        !! @param[in] ia The index in A where the requested row starts.
        !! @param[in] alu The LU-factored matrix from ILUT.
        !! @param[in] jlu The LU-factored matrix from ILUT.
        !! @param[in] ju The LU-factored matrix from ILUT.
        !! @param[out] ierr Error flag:
        !!  * 0: Successful return
        !!  * 1: Convergence not achieved.
        !!  * -1: The initial guess seems to be the exact solution.
        pure subroutine pgmres(n, im, rhs, sol, vv, eps, maxits, iout, aa, ja, ia, &
            alu, jlu, ju, ierr)
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: n, im, maxits, iout, ja(*), ia(n+1), &
                jlu(*), ju(n)
            integer(int32), intent(out) :: ierr
            real(real64), intent(in) :: aa(*), eps, alu(*)
            real(real64), intent(inout) :: rhs(n), sol(n)
            real(real64), intent(out) :: vv(n,*)
        end subroutine

        !> @brief Solves the LU-factored system (LU) x = y.
        !!
        !! @param[in] n The dimension of the system.
        !! @param[in] y The N-element right-hand-side vector.
        !! @param[out] x The N-element solution vector.
        !! @param[in] alu The LU-factored matrix.
        !! @param[in] jlu The LU-factored matrix.
        !! @param[in] ju The LU-factored matrix.
        pure subroutine lusol(n, y, x, alu, jlu, ju)
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: n, jlu(*), ju(*)
            real(real64), intent(in) :: y(n), alu(*)
            real(real64), intent(out) :: x(n)
        end subroutine
    end interface
end module