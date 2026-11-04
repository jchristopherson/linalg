! sparskit.f90
! SPDX-License-Identifier: MIT
module sparskit
    !! An interface to the SPARSKIT library available at
    !! https://www-users.cse.umn.edu/~saad/software/SPARSKIT/.
    implicit none

    interface
        subroutine gmres(n, rhs, sol, ipar, fpar, w)
            !! Restarted GMRES with reverse communication from ITERS.F.
            !! Saved backend state makes complete solves nonreentrant: do not
            !! interleave concurrent or nested invocations of this routine.
            use iso_fortran_env, only : real64
            integer, intent(in) :: n
                !! Positive system order.
            real(real64), intent(in) :: rhs(n)
                !! Right-hand side, preserved across communication requests.
            real(real64), intent(inout) :: sol(n)
                !! Initial guess on entry; current solution on return.
            integer, intent(inout) :: ipar(16)
                !! Protocol and counters. Start with IPAR(1)=0. Returns 1 for
                !! A*X, 5 for right-preconditioner application, or a nonpositive
                !! termination status. IPAR(2) selects preconditioning (0 none,
                !! 2 right); IPAR(3)=2 selects RHS-relative residual stopping;
                !! IPAR(4) is W length; IPAR(5) is restart length (>1, else
                !! backend uses 15); IPAR(6) is the matvec budget. IPAR(8:9)
                !! point to input/output vectors in W. Other entries retain
                !! backend state and must be preserved between requests.
            real(real64), intent(inout) :: fpar(16)
                !! Relative and absolute tolerances in entries 1 and 2, with
                !! 0 <= FPAR(1) < 1 and FPAR(2) >= 0, not both zero. Other
                !! entries hold residual estimates, counters, and backend state.
            real(real64), intent(inout) :: w(*)
                !! Caller-owned workspace, retained between requests. For
                !! restart length M, allocate at least
                !! (N+3)*(M+2) + M*(M+1)/2 real elements. Requested vectors
                !! occupy W(IPAR(8):IPAR(8)+N-1) and the analogous IPAR(9) slice.
        end subroutine
    end interface

    ! BLASSM.F
    interface
        pure subroutine amub(nrow, ncol, job, a, ja, ia, b, jb, ib, c, jc, ic, &
            nzmax, iw, ierr)
            !! Computes the matrix product C = A * B in sparse matrix form.
            !! The routine multiplies two matrices stored in compressed sparse row
            !! format and returns the product in a new sparse representation.
            use iso_fortran_env, only : real64

            integer, intent(in) :: nrow
                !! Number of rows of A.
            integer, intent(in) :: ncol
                !! Number of columns of B.
            integer, intent(in) :: job
                !! Controls the product computation mode.
            integer, intent(in) :: nzmax
                !! Maximum number of entries in the output matrix.
            integer, intent(in) :: ja(*)
                !! Column indices of A.
            integer, intent(in) :: ia(nrow+1)
                !! Row pointers for A.
            integer, intent(in) :: jb(*)
                !! Column indices of B.
            integer, intent(in) :: ib(*)
                !! Row pointers for B.
            integer, intent(out) :: jc(*)
                !! Column indices of the product matrix.
            integer, intent(out) :: ic(*)
                !! Row pointers for the product matrix.
            integer, intent(out) :: iw(ncol)
                !! Workspace array.
            integer, intent(out) :: ierr
                !! Error status indicator.
            real(real64), intent(in) :: a(*)
                !! Values of matrix A.
            real(real64), intent(in) :: b(*)
                !! Values of matrix B.
            real(real64), intent(out) :: c(*)
                !! Values of the product matrix.
        end subroutine

        pure subroutine aplb(nrow, ncol, job, a, ja, ia, b, jb, ib, c, jc, ic, &
            nzmax, iw, ierr)
            !! Computes the matrix sum C = A + B for sparse matrices stored in CSR format.
            !! The result is returned in a sparse matrix representation suitable for
            !! subsequent linear algebra operations.
            use iso_fortran_env, only : real64

            integer, intent(in) :: nrow
                !! Number of rows in the matrices.
            integer, intent(in) :: ncol
                !! Number of columns in the matrices.
            integer, intent(in) :: job
                !! Controls the sparse-matrix summation mode.
            integer, intent(in) :: nzmax
                !! Maximum number of entries in the output matrix.
            integer, intent(in) :: ja(*)
                !! Column indices of A.
            integer, intent(in) :: ia(nrow+1)
                !! Row pointers for A.
            integer, intent(in) :: jb(*)
                !! Column indices of B.
            integer, intent(in) :: ib(nrow+1)
                !! Row pointers for B.
            real(real64), intent(in) :: a(*)
                !! Values of matrix A.
            real(real64), intent(in) :: b(*)
                !! Values of matrix B.
            real(real64), intent(out) :: c(*)
                !! Values of the matrix sum.
            integer, intent(out) :: jc(*)
                !! Column indices of the matrix sum.
            integer, intent(out) :: ic(nrow+1)
                !! Row pointers for the matrix sum.
            integer, intent(out) :: iw(ncol)
                !! Workspace array.
            integer, intent(out) :: ierr
                !! Error status indicator.
        end subroutine

        pure subroutine aplsb(nrow, ncol, a, ja, ia, s, b, jb, ib, c, jc, ic, &
            nzmax, iw, ierr)
            !! Computes C = A + s * B for matrices stored in CSR format.
            !! The routine scales one sparse matrix by a scalar and adds it to another.
            use iso_fortran_env, only : real64

            integer, intent(in) :: nrow
                !! Number of rows in the matrices.
            integer, intent(in) :: ncol
                !! Number of columns in the matrices.
            integer, intent(in) :: nzmax
                !! Maximum number of entries in the output matrix.
            integer, intent(in) :: ja(*)
                !! Column indices of A.
            integer, intent(in) :: ia(nrow+1)
                !! Row pointers for A.
            integer, intent(in) :: jb(*)
                !! Column indices of B.
            integer, intent(in) :: ib(nrow+1)
                !! Row pointers for B.
            real(real64), intent(in) :: s
                !! Scalar multiplier applied to B.
            real(real64), intent(in) :: a(*)
                !! Values of matrix A.
            real(real64), intent(in) :: b(*)
                !! Values of matrix B.
            real(real64), intent(out) :: c(*)
                !! Values of the result matrix.
            integer, intent(out) :: jc(*)
                !! Column indices of the result matrix.
            integer, intent(out) :: ic(nrow+1)
                !! Row pointers of the result matrix.
            integer, intent(out) :: iw(ncol)
                !! Workspace array.
            integer, intent(out) :: ierr
                !! Error status indicator.
        end subroutine
    end interface

    ! FORMATS.F
    interface
        pure subroutine csrcsc2(n, n2, job, ipos, a, ja, ia, a0, ja0, ia0)
            !! Converts a CSR matrix into a CSC matrix by transposing the sparse pattern.
            !! The routine can fill values or only the structure depending on the job flag.
            use iso_fortran_env, only : real64

            integer, intent(in) :: n
                !! Number of rows in the input matrix.
            integer, intent(in) :: n2
                !! Number of columns in the output matrix.
            integer, intent(in) :: job
                !! Controls whether values or only structure are generated.
            integer, intent(in) :: ipos
                !! Starting position in the output arrays.
            integer, intent(in) :: ja(*)
                !! Column indices of the input matrix.
            integer, intent(in) :: ia(n+1)
                !! Row pointers of the input matrix.
            real(real64), intent(in) :: a(*)
                !! Values of the input matrix.
            real(real64), intent(out) :: a0(*)
                !! Values of the transposed matrix.
            integer, intent(out) :: ja0(*)
                !! Column indices of the transposed matrix.
            integer, intent(out) :: ia0(n2+1)
                !! Row pointers of the transposed matrix.
        end subroutine

        pure subroutine bndcsr(n, abd, nabd, lowd, ml, mu, a, ja, ia, len, ierr)
            !! Converts a banded matrix in LINPACK/BLAS/LAPACK format into CSR format.
            !! This is useful when a banded sparse matrix needs to be treated as a general CSR matrix.
            use iso_fortran_env, only : real64

            integer, intent(in) :: n
                !! Matrix dimension.
            integer, intent(in) :: nabd
                !! Leading dimension of the band storage.
            integer, intent(in) :: lowd
                !! Position of the lowest diagonal in the band storage.
            integer, intent(in) :: ml
                !! Lower bandwidth.
            integer, intent(in) :: mu
                !! Upper bandwidth.
            integer, intent(in) :: len
                !! Available storage length in the sparse output arrays.
            real(real64), intent(in) :: abd(nabd,*)
                !! Band matrix stored in banded form.
            real(real64), intent(out) :: a(*)
                !! Values of the CSR matrix.
            integer, intent(out) :: ia(n+1)
                !! Row pointers of the CSR matrix.
            integer, intent(out) :: ja(*)
                !! Column indices of the CSR matrix.
            integer, intent(out) :: ierr
                !! Error flag returned by the conversion.
        end subroutine

        pure subroutine csrmsr(n, a, ja, ia, ao, jao, wk, iwk)
            !! Converts a CSR matrix to an MSR matrix.
            !! The routine rewrites the sparse structure and values into the modified sparse row format.
            use iso_fortran_env, only : real64
            integer, intent(in) :: n
                !! Matrix dimension.
            integer, intent(inout) :: ja(*)
                !! Column indices in the CSR input, overwritten by the MSR representation.
            integer, intent(inout) :: ia(n+1)
                !! Row pointers in the CSR input, overwritten by the MSR representation.
            integer, intent(out) :: jao(*)
                !! Column indices of the MSR matrix.
            integer, intent(out) :: iwk(n+1)
                !! Workspace for the conversion.
            real(real64), intent(inout) :: a(*)
                !! Matrix values converted in place.
            real(real64), intent(out) :: ao(*)
                !! Values in the MSR output.
            real(real64), intent(out) :: wk(n)
                !! Workspace for the conversion.
        end subroutine

        pure subroutine msrcsr(n, a, ja, ao, jao, iao, wk, iwk)
            !! Converts an MSR matrix to CSR format.
            !! This routine reconstructs the row pointer array and column indices for a general CSR representation.
            use iso_fortran_env, only : real64
            integer, intent(in) :: n
                !! Matrix dimension.
            integer, intent(in) :: ja(*)
                !! Column index array of the MSR matrix.
            integer, intent(out) :: jao(*)
                !! Column indices of the CSR matrix.
            integer, intent(out) :: iao(n+1)
                !! Row pointers for the CSR matrix.
            integer, intent(out) :: iwk(n+1)
                !! Workspace array.
            real(real64), intent(in) :: a(*)
                !! Values of the MSR matrix.
            real(real64), intent(out) :: ao(*)
                !! Values of the CSR matrix.
            real(real64), intent(out) :: wk(n)
                !! Workspace array.
        end subroutine

        pure subroutine coocsr(nrow, nnz, a, ir, jc, ao, jao, iao)
            !! Converts a matrix stored in coordinate format to CSR format.
            !! The routine compresses the nonzero entries into row-wise sparse storage.
            use iso_fortran_env, only : real64
            integer, intent(in) :: nrow
                !! Number of rows in the matrix.
            integer, intent(in) :: nnz
                !! Number of nonzero entries.
            integer, intent(in) :: jc(*)
                !! Column indices in the coordinate representation.
            integer, intent(inout) :: ir(*)
                !! Row indices rewritten into CSR order.
            real(real64), intent(in) :: a(*)
                !! Values in the coordinate array.
            integer, intent(out) :: jao(*)
                !! Column indices of the CSR matrix.
            integer, intent(out) :: iao(*)
                !! Row pointers of the CSR matrix.
            real(real64), intent(out) :: ao(*)
                !! Values in the CSR matrix.
        end subroutine
    end interface

    ! UNARY.F
    interface
        function getelm(i, j, a, ja, ia, iadd, sorted) result(rst)
            !! Returns the value A(i, j) from a sparse matrix stored in CSR format.
            !! If the entry is not present, the routine reports its location as zero.
            use iso_fortran_env, only : real64
            integer, intent(in) :: i
                !! Row index to retrieve.
            integer, intent(in) :: j
                !! Column index to retrieve.
            integer, intent(in) :: ia(*)
                !! Row pointers of the CSR matrix.
            integer, intent(in) :: ja(*)
                !! Column indices of the CSR matrix.
            real(real64), intent(in) :: a(*)
                !! Values of the CSR matrix.
            integer, intent(out) :: iadd
                !! Pointer to the located entry, or zero if absent.
            logical, intent(in) :: sorted
                !! Whether the rows are already sorted by column index.
            real(real64) :: rst
                !! Value A(i, j), or zero if the entry is not present.
        end function

        pure subroutine getdia(nrow, ncol, job, a, ja, ia, len, diag, idiag, ioff)
            !! Extracts a diagonal from a sparse matrix.
            !! This routine can operate on a specified diagonal offset and return the extracted values.
            use iso_fortran_env, only : real64
            integer, intent(in) :: nrow
                !! Number of rows in the matrix.
            integer, intent(in) :: ncol
                !! Number of columns in the matrix.
            integer, intent(in) :: job
                !! Selects which diagonal is extracted.
            integer, intent(in) :: ja(*)
                !! Column indices of the CSR matrix.
            integer, intent(in) :: ia(*)
                !! Row pointers of the CSR matrix.
            integer, intent(in) :: ioff
                !! Offset of the diagonal to extract.
            integer, intent(out) :: len
                !! Number of diagonal entries returned.
            integer, intent(out) :: idiag(*)
                !! Indices of the extracted diagonal entries.
            real(real64), intent(in) :: a(*)
                !! Values of the CSR matrix.
            real(real64), intent(out) :: diag(*)
                !! Extracted diagonal values.
        end subroutine

        pure subroutine csort(n, a, ja, ia, values)
            !! Sorts the entries in each row of a CSR matrix by column index.
            !! Optionally, the values array is permuted in the same way as the indices.
            use iso_fortran_env, only : real64
            integer, intent(in) :: n
                !! Matrix dimension.
            real(real64), intent(inout) :: a(*)
                !! Sparse values reordered in place.
            integer, intent(inout) :: ja(*)
                !! Column indices reordered in place.
            integer, intent(in) :: ia(*)
                !! Row pointers describing the sparse structure.
            logical, intent(in) :: values
                !! Whether the values array should be reordered with the indices.
        end subroutine

        pure subroutine clncsr(job, value2, nrow, a, ja, ia, indu, iwk)
            !! Cleans up a CSR matrix by removing duplicates and sorting entries.
            !! The routine can optionally perform partial ordering of the structure.
            use iso_fortran_env, only : real64
            integer, intent(in) :: job
                !! Controls cleaning and ordering behavior.
            integer, intent(in) :: value2
                !! Additional option controlling the cleanup strategy.
            integer, intent(in) :: nrow
                !! Number of rows in the matrix.
            real(real64), intent(inout) :: a(*)
                !! Matrix values, cleaned and reordered in place.
            integer, intent(inout) :: ja(*)
                !! Column indices, cleaned and reordered in place.
            integer, intent(inout) :: ia(*)
                !! Row pointers, cleaned and reordered in place.
            integer, intent(inout) :: indu(*)
                !! Indirection array used during the cleanup process.
            integer, intent(inout) :: iwk(*)
                !! Workspace array.
        end subroutine
    end interface

    ! ILUT.F
    interface
        pure subroutine ilut(n, a, ja, ia, lfil, droptol, alu, jlu, ju, iwk, w, &
            jw, ierr)
            !! Computes an incomplete LU factorization of a sparse matrix in CSR format.
            !! The factorization uses a dual truncation strategy to limit fill-in and drop small values.
            use iso_fortran_env, only : real64
            integer, intent(in) :: n
                !! Matrix order.
            integer, intent(in) :: ja(*)
                !! Column indices of the matrix.
            integer, intent(in) :: ia(n+1)
                !! Row pointers of the matrix.
            integer, intent(in) :: lfil
                !! Maximum fill level per row.
            integer, intent(in) :: iwk
                !! Workspace size.
            integer, intent(out) :: jlu(*)
                !! Column indices of the incomplete LU factors.
            integer, intent(out) :: ju(n)
                !! Diagonal position array for the factors.
            integer, intent(out) :: jw(2*n)
                !! Working array for factorization bookkeeping.
            integer, intent(out) :: ierr
                !! Error status indicator.
            real(real64), intent(in) :: a(*)
                !! Original sparse matrix entries.
            real(real64), intent(in) :: droptol
                !! Drop tolerance used during factorization.
            real(real64), intent(out) :: alu(*)
                !! Nonzero entries of the factors.
            real(real64), intent(out) :: w(n+1)
                !! Workspace array.
        end subroutine

        pure subroutine ilutp(n, a, ja, ia, lfil, droptol, permtol, mbloc, alu, &
            jlu, ju, iwk, w, jw, iperm, ierr)
            !! Computes an incomplete LU factorization with pivoting and threshold dropping.
            !! The routine may permute rows and columns to improve stability while preserving the sparse structure.
            use iso_fortran_env, only : real64
            integer, intent(in) :: n
                !! Matrix order.
            integer, intent(in) :: lfil
                !! Maximum fill level per row.
            integer, intent(in) :: iwk
                !! Workspace size.
            integer, intent(in) :: mbloc
                !! Block size for pivoting and ordering.
            integer, intent(inout) :: ja(*)
                !! Column indices of the matrix, possibly permuted in place.
            integer, intent(inout) :: ia(n+1)
                !! Row pointers of the matrix, possibly permuted in place.
            integer, intent(out) :: jlu(*)
                !! Column indices of the incomplete LU factors.
            integer, intent(out) :: ju(n)
                !! Diagonal position array for the factors.
            integer, intent(out) :: jw(2*n)
                !! Working array for factorization bookkeeping.
            integer, intent(out) :: iperm(2*n)
                !! Row and column permutation indices.
            integer, intent(out) :: ierr
                !! Error status indicator.
            real(real64), intent(in) :: droptol
                !! Drop tolerance used during factorization.
            real(real64), intent(in) :: permtol
                !! Permutation threshold for pivoting.
            real(real64), intent(inout) :: a(*)
                !! Matrix values updated during factorization.
            real(real64), intent(out) :: alu(*)
                !! Nonzero entries of the factors.
            real(real64), intent(out) :: w(n+1)
                !! Workspace array.
        end subroutine

        pure subroutine ilud(n, a, ja, ia, alph, tol, alu, jlu, ju, iwk, w, jw, ierr)
            !! Computes the incomplete LU factorization of a sparse matrix using a standard dropping rule.
            use iso_fortran_env, only : real64
            integer, intent(in) :: n
                !! Matrix order.
            integer, intent(in) :: iwk
                !! Workspace size.
            integer, intent(in) :: ja(*)
                !! Column indices of the matrix.
            integer, intent(in) :: ia(n+1)
                !! Row pointers of the matrix.
            integer, intent(out) :: jlu(*)
                !! Column indices of the incomplete LU factors.
            integer, intent(out) :: ju(n)
                !! Diagonal position array for the factors.
            integer, intent(out) :: jw(2*n)
                !! Working array for factorization bookkeeping.
            integer, intent(out) :: ierr
                !! Error status indicator.
            real(real64), intent(in) :: a(*)
                !! Original sparse matrix entries.
            real(real64), intent(in) :: alph
                !! Fill parameter used by the dropping rule.
            real(real64), intent(in) :: tol
                !! Drop tolerance.
            real(real64), intent(out) :: alu(*)
                !! Nonzero entries of the factors.
            real(real64), intent(out) :: w(2*n)
                !! Workspace array.
        end subroutine

        pure subroutine iludp(n, a, ja, ia, alph, droptol, permtol, mbloc, alu, &
            jlu, ju, iwk, w, jw, iperm, ierr)
            !! Computes a pivoted incomplete LU factorization with standard dropping and tolerance controls.
            use iso_fortran_env, only : real64
            integer, intent(in) :: n
                !! Matrix order.
            integer, intent(in) :: iwk
                !! Workspace size.
            integer, intent(in) :: mbloc
                !! Block size for pivoting and ordering.
            integer, intent(inout) :: ja(*)
                !! Column indices of the matrix, potentially permuted in place.
            integer, intent(inout) :: ia(n+1)
                !! Row pointers of the matrix, potentially permuted in place.
            integer, intent(out) :: jlu(*)
                !! Column indices of the incomplete LU factors.
            integer, intent(out) :: ju(n)
                !! Diagonal position array for the factors.
            integer, intent(out) :: jw(2*n)
                !! Working array for factorization bookkeeping.
            integer, intent(out) :: iperm(2*n)
                !! Row and column permutation indices.
            integer, intent(out) :: ierr
                !! Error status indicator.
            real(real64), intent(in) :: alph
                !! Fill parameter used by the dropping rule.
            real(real64), intent(in) :: droptol
                !! Drop tolerance.
            real(real64), intent(in) :: permtol
                !! Permutation threshold.
            real(real64), intent(inout) :: a(*)
                !! Matrix values updated during the pivoted factorization.
            real(real64), intent(out) :: alu(*)
                !! Nonzero entries of the factors.
            real(real64), intent(out) :: w(2*n)
                !! Workspace array.
        end subroutine

        pure subroutine pgmres(n, im, rhs, sol, vv, eps, maxits, iout, aa, ja, ia, &
            alu, jlu, ju, ierr)
            !! Solves a sparse linear system using ILUT-preconditioned GMRES.
            !! This routine uses the incomplete LU factors as a preconditioner and performs a Krylov iteration.
            use iso_fortran_env, only : real64
            integer, intent(in) :: n
                !! Matrix order.
            integer, intent(in) :: im
                !! Restart length for GMRES.
            integer, intent(in) :: maxits
                !! Maximum number of Krylov iterations.
            integer, intent(in) :: iout
                !! Output level for convergence information.
            integer, intent(in) :: ja(*)
                !! Column indices of the matrix.
            integer, intent(in) :: ia(n+1)
                !! Row pointers of the matrix.
            integer, intent(in) :: jlu(*)
                !! Column indices of the ILU preconditioner.
            integer, intent(in) :: ju(n)
                !! Diagonal positions of the ILU preconditioner.
            integer, intent(out) :: ierr
                !! Error status indicator.
            real(real64), intent(in) :: aa(*)
                !! Matrix values for the linear system.
            real(real64), intent(in) :: eps
                !! Convergence tolerance.
            real(real64), intent(in) :: alu(*)
                !! Incomplete LU factors used as a preconditioner.
            real(real64), intent(inout) :: rhs(n)
                !! Right-hand side vector, overwritten by the solution estimate.
            real(real64), intent(inout) :: sol(n)
                !! Initial guess, overwritten by the final iterate.
            real(real64), intent(out) :: vv(n,*)
                !! Krylov subspace basis vectors.
        end subroutine

        pure subroutine lusol(n, y, x, alu, jlu, ju)
            !! Solves a system given its LU factors from the incomplete LU decomposition.
            use iso_fortran_env, only : real64
            integer, intent(in) :: n
                !! Matrix order.
            integer, intent(in) :: jlu(*)
                !! Column indices of the LU factors.
            integer, intent(in) :: ju(*)
                !! Diagonal positions of the LU factors.
            real(real64), intent(in) :: y(n)
                !! Right-hand side vector.
            real(real64), intent(in) :: alu(*)
                !! LU factor entries.
            real(real64), intent(out) :: x(n)
                !! Solution vector.
        end subroutine
    end interface
end module
