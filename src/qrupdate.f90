module qrupdate
    !! A module providing explicit interfaces for the QRUPDATE library.
    !! The routines support rank-one updates and downdates for QR
    !! factorizations in both real and complex arithmetic.
    implicit none

    interface
        pure subroutine DQR1UP(m, n, k, q, ldq, r, ldr, u, v, w)
            !! Performs a real QR update by appending one column and one row to the factorization.
            use iso_fortran_env, only : int32, real64

            integer(int32), intent(in) :: m
                !! Number of rows in the factorization.
            integer(int32), intent(in) :: n
                !! Number of columns in the factorization.
            integer(int32), intent(in) :: k
                !! Rank-one update dimension.
            integer(int32), intent(in) :: ldq
                !! Leading dimension of Q.
            integer(int32), intent(in) :: ldr
                !! Leading dimension of R.
            real(real64), intent(inout) :: q(ldq,*)
                !! Orthogonal factor updated in place.
            real(real64), intent(inout) :: r(ldr,*)
                !! Upper-triangular factor updated in place.
            real(real64), intent(inout) :: u(*)
                !! First update vector.
            real(real64), intent(inout) :: v(*)
                !! Second update vector.
            real(real64), intent(out) :: w(*)
                !! Working array for the update.
        end subroutine

        pure subroutine ZQR1UP(m, n, k, q, ldq, r, ldr, u, v, w, rw)
            !! Performs a complex QR update by appending one column and one row to the factorization.
            use iso_fortran_env, only : int32, real64

            integer(int32), intent(in) :: m
                !! Number of rows in the factorization.
            integer(int32), intent(in) :: n
                !! Number of columns in the factorization.
            integer(int32), intent(in) :: k
                !! Rank-one update dimension.
            integer(int32), intent(in) :: ldq
                !! Leading dimension of Q.
            integer(int32), intent(in) :: ldr
                !! Leading dimension of R.
            complex(real64), intent(inout) :: q(ldq,*)
                !! Orthogonal factor updated in place.
            complex(real64), intent(inout) :: r(ldr,*)
                !! Upper-triangular factor updated in place.
            complex(real64), intent(inout) :: u(*)
                !! First update vector.
            complex(real64), intent(inout) :: v(*)
                !! Second update vector.
            complex(real64), intent(out) :: w(*)
                !! Working array for the update.
            real(real64), intent(out) :: rw(*)
                !! Real workspace for the complex update.
        end subroutine

        pure subroutine DCH1UP(n, r, ldr, u, w)
            !! Updates a real Cholesky factorization by a rank-one matrix.
            use iso_fortran_env, only : int32, real64

            integer(int32), intent(in) :: n
                !! Order of the Cholesky factor.
            integer(int32), intent(in) :: ldr
                !! Leading dimension of R.
            real(real64), intent(inout) :: r(ldr,*)
                !! Cholesky factor updated in place.
            real(real64), intent(inout) :: u(*)
                !! Rank-one update vector.
            real(real64), intent(out) :: w(*)
                !! Work array for the update.
        end subroutine

        pure subroutine ZCH1UP(n, r, ldr, u, w)
            !! Updates a complex Cholesky factorization by a rank-one matrix.
            use iso_fortran_env, only : int32, real64

            integer(int32), intent(in) :: n
                !! Order of the Cholesky factor.
            integer(int32), intent(in) :: ldr
                !! Leading dimension of R.
            complex(real64), intent(inout) :: r(ldr,*)
                !! Cholesky factor updated in place.
            complex(real64), intent(inout) :: u(*)
                !! Rank-one update vector.
            real(real64), intent(out) :: w(*)
                !! Work array for the update.
        end subroutine

        pure subroutine DCH1DN(n, r, ldr, u, w, info)
            !! Downdates a real Cholesky factorization by a rank-one matrix.
            use iso_fortran_env, only : int32, real64

            integer(int32), intent(in) :: n
                !! Order of the Cholesky factor.
            integer(int32), intent(in) :: ldr
                !! Leading dimension of R.
            real(real64), intent(inout) :: r(ldr,*)
                !! Cholesky factor updated in place.
            real(real64), intent(inout) :: u(*)
                !! Rank-one downdate vector.
            real(real64), intent(out) :: w(*)
                !! Work array for the downdate.
            integer(int32), intent(out) :: info
                !! Status flag for the downdate.
        end subroutine

        pure subroutine ZCH1DN(n, r, ldr, u, rw, info)
            !! Downdates a complex Cholesky factorization by a rank-one matrix.
            use iso_fortran_env, only : int32, real64

            integer(int32), intent(in) :: n
                !! Order of the Cholesky factor.
            integer(int32), intent(in) :: ldr
                !! Leading dimension of R.
            complex(real64), intent(inout) :: r(ldr,*)
                !! Cholesky factor updated in place.
            complex(real64), intent(inout) :: u(*)
                !! Rank-one downdate vector.
            real(real64), intent(out) :: rw(*)
                !! Real workspace for the downdate.
            integer(int32), intent(out) :: info
                !! Status flag for the downdate.
        end subroutine
    end interface
end module