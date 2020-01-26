! linalg_c_api.f90

!> @brief Provides a C-friendly API to the LINALG library.  Notice, all C-API 
!! LINALG routines begin with the prefix "la_".
module linalg_c_api
    use iso_c_binding
    use linalg_core
    use linalg_constants
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
    subroutine la_mtx_mult(transa, transb, m, n, k, alpha, a, lda, b, ldb, &
            beta, c, ldc) bind(C, name="la_mtx_mult")
        ! Arugments
        logical(c_bool), intent(in), value :: transa, transb
        integer(c_int), intent(in), value :: m, n, k, lda, ldb, ldc
        real(c_double), intent(in), value :: alpha, beta
        real(c_double), intent(in) :: a(lda,*), b(ldb,*)
        real(c_double), intent(out) :: c(ldc,*)

        ! Local Variables
        character :: ta, tb

        ! Initialization
        ta = "N"
        if (transa) ta = "T"

        tb = "N"
        if (transb) tb = "N"

        ! Call DGEMM directly
        call DGEMM(ta, tb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    end subroutine

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
    subroutine la_mtx_mult_cmplx(opa, opb, m, n, k, alpha, a, lda, b, ldb, &
            beta, c, ldc) bind(C, name="la_mtx_mult_cmplx")
        ! Arguments
        integer(c_int), intent(in), value :: opa, opb, m, n, k, lda, ldb, ldc
        complex(c_double), intent(in), value :: alpha, beta
        complex(c_double), intent(in) :: a(lda,*), b(ldb,*)
        complex(c_double), intent(out) :: c(ldc,*)

        ! Local Variables
        character :: ta, tb

        ! Initialization
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
        
        ! Call ZGEMM directly
        call ZGEMM(ta, tb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    end subroutine

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end module
