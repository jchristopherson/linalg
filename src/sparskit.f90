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
        subroutine amub(nrow, ncol, job, a, ja, ia, b, jb, ib, c, jc, ic, &
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
        subroutine aplb(nrow, ncol, job, a, ja, ia, b, jb, ib, c, jc, ic, &
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
        subroutine aplsb(nrow, ncol, a, ja, ia, s, b, jb, ib, c, jc, ic, &
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
        subroutine csrcsc2(n, n2, job, ipos, a, ja, ia, a0, ja0, ia0)
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: n, n2, job, ipos, ja(*), ia(n+1)
            integer(int32), intent(out) :: ja0(*), ia0(n2+1)
            real(real64), intent(in) :: a(*)
            real(real64), intent(out) :: a0(*)
        end subroutine
    end interface

    ! MATVEC.F
    interface
    ! AMUX
    ! ATMUXR
    end interface
end module