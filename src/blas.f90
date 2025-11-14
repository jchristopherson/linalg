!> @brief A module providing explicit interfaces to BLAS routines.
module blas
    implicit none

    interface
        pure subroutine DGEMM(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, &
            c, ldc)
            use iso_fortran_env, only : int32, real64
            character, intent(in) :: transa, transb
            integer(int32), intent(in) :: m, n, k, lda, ldb, ldc
            real(real64), intent(in) :: alpha, beta, a(lda,*), b(ldb,*)
            real(real64), intent(inout) :: c(ldc,*)
        end subroutine

        pure subroutine DGEMV(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
            use iso_fortran_env, only : int32, real64
            character, intent(in) :: trans
            integer(int32), intent(in) :: m, n, lda, incx, incy
            real(real64), intent(in) :: alpha, beta, a(lda,*), x(*)
            real(real64), intent(inout) :: y(*)
        end subroutine

        pure subroutine ZGEMM(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, &
            c, ldc)
            use iso_fortran_env, only : int32, real64
            character, intent(in) :: transa, transb
            integer(int32), intent(in) :: m, n, k, lda, ldb, ldc
            complex(real64), intent(in) :: alpha, beta, a(lda,*), b(ldb,*)
            complex(real64), intent(inout) :: c(ldc,*)
        end subroutine

        pure subroutine ZGEMV(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
            use iso_fortran_env, only : int32, real64
            character, intent(in) :: trans
            integer(int32), intent(in) :: m, n, lda, incx, incy
            complex(real64), intent(in) :: alpha, beta, a(lda,*), x(*)
            complex(real64), intent(inout) :: y(*)
        end subroutine

        pure subroutine DTRSM(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
            use iso_fortran_env, only : int32, real64
            character, intent(in) :: side, uplo, transa, diag
            integer(int32), intent(in) :: m, n, lda, ldb
            real(real64), intent(in) :: alpha, a(lda,*)
            real(real64), intent(inout) :: b(ldb,*)
        end subroutine

        pure subroutine ZTRSM(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
            use iso_fortran_env, only : int32, real64
            character, intent(in) :: side, uplo, transa, diag
            integer(int32), intent(in) :: m, n, lda, ldb
            complex(real64), intent(in) :: alpha, a(lda,*)
            complex(real64), intent(inout) :: b(ldb,*)
        end subroutine

        pure subroutine DTRSV(uplo, trans, diag, n, a, lda, x, incx)
            use iso_fortran_env, only : int32, real64
            character, intent(in) :: uplo, trans, diag
            integer(int32), intent(in) :: n, lda, incx
            real(real64), intent(in) :: a(lda,*)
            real(real64), intent(inout) :: x(*)
        end subroutine

        pure subroutine ZTRSV(uplo, trans, diag, n, a, lda, x, incx)
            use iso_fortran_env, only : int32, real64
            character, intent(in) :: uplo, trans, diag
            integer(int32), intent(in) :: n, lda, incx
            complex(real64), intent(in) :: a(lda,*)
            complex(real64), intent(inout) :: x(*)
        end subroutine

        pure subroutine DSCAL(n, da, dx, incx)
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: n, incx
            real(real64), intent(in) :: da
            real(real64), intent(inout) :: dx(*)
        end subroutine

        pure subroutine ZSCAL(n, za, zx, incx)
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: n, incx
            complex(real64), intent(in) :: za
            complex(real64), intent(inout) :: zx(*)
        end subroutine

        pure subroutine ZDSCAL(n, da, zx, incx)
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: n, incx
            real(real64), intent(in) :: da
            complex(real64), intent(inout) :: zx(*)
        end subroutine

        pure subroutine DGBMV(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, &
            incy)
            use iso_fortran_env, only : int32, real64
            character, intent(in) :: trans
            integer(int32), intent(in) :: m, n, kl, ku, lda, incx, incy
            real(real64), intent(in) :: alpha, beta, a(lda,*), x(*)
            real(real64), intent(inout) :: y(*)
        end subroutine

        pure subroutine ZGBMV(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, &
            incy)
            use iso_fortran_env, only : int32, real64
            character, intent(in) :: trans
            integer(int32), intent(in) :: m, n, kl, ku, lda, incx, incy
            complex(real64), intent(in) :: alpha, beta, a(lda,*), x(*)
            complex(real64), intent(inout) :: y(*)
        end subroutine

        pure function DDOT(n, dx, incx, dy, incy) result(rst)
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: n, incx, incy
            real(real64), intent(in) :: dx(*), dy(*)
            real(real64) :: rst
        end function

        pure subroutine DSWAP(n, dx, incx, dy, incy)
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: n, incx, incy
            real(real64), intent(inout) :: dx(*), dy(*)
        end subroutine
    end interface
end module