module blas
    !! A module providing explicit interfaces to BLAS routines.
    !! These declarations expose the core dense linear algebra kernels used by
    !! the higher-level matrix routines in the project, including matrix-matrix
    !! and matrix-vector operations, triangular solves, scaling, and dot-product
    !! helpers.
    implicit none

    interface
        pure subroutine DGEMM(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
            !! Computes C = alpha * op(A) * op(B) + beta * C for real matrices.
            use iso_fortran_env, only : int32, real64
            character, intent(in) :: transa
                !! Operation on the first matrix factor.
            character, intent(in) :: transb
                !! Operation on the second matrix factor.
            integer(int32), intent(in) :: m
                !! Number of rows of the product matrix.
            integer(int32), intent(in) :: n
                !! Number of columns of the product matrix.
            integer(int32), intent(in) :: k
                !! Number of columns of A and rows of B.
            integer(int32), intent(in) :: lda
                !! Leading dimension of A.
            integer(int32), intent(in) :: ldb
                !! Leading dimension of B.
            integer(int32), intent(in) :: ldc
                !! Leading dimension of C.
            real(real64), intent(in) :: alpha
                !! Scalar multiplier for the matrix product.
            real(real64), intent(in) :: beta
                !! Scalar multiplier for the existing C matrix.
            real(real64), intent(in) :: a(lda,*)
                !! First matrix operand.
            real(real64), intent(in) :: b(ldb,*)
                !! Second matrix operand.
            real(real64), intent(inout) :: c(ldc,*)
                !! Result matrix updated in place.
        end subroutine

        pure subroutine DGEMV(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
            !! Computes y = alpha * op(A) * x + beta * y for real matrices.
            use iso_fortran_env, only : int32, real64
            character, intent(in) :: trans
                !! Matrix operation to apply to A.
            integer(int32), intent(in) :: m
                !! Number of rows in the matrix operand.
            integer(int32), intent(in) :: n
                !! Number of columns in the matrix operand.
            integer(int32), intent(in) :: lda
                !! Leading dimension of A.
            integer(int32), intent(in) :: incx
                !! Stride for the x vector.
            integer(int32), intent(in) :: incy
                !! Stride for the y vector.
            real(real64), intent(in) :: alpha
                !! Scalar multiplier for A*x.
            real(real64), intent(in) :: beta
                !! Scalar multiplier for the existing y vector.
            real(real64), intent(in) :: a(lda,*)
                !! Matrix operand.
            real(real64), intent(in) :: x(*)
                !! Input vector.
            real(real64), intent(inout) :: y(*)
                !! Updated vector.
        end subroutine

        pure subroutine ZGEMM(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
            !! Computes C = alpha * op(A) * op(B) + beta * C for complex matrices.
            use iso_fortran_env, only : int32, real64
            character, intent(in) :: transa
                !! Operation on the first matrix factor.
            character, intent(in) :: transb
                !! Operation on the second matrix factor.
            integer(int32), intent(in) :: m
                !! Number of rows of the product matrix.
            integer(int32), intent(in) :: n
                !! Number of columns of the product matrix.
            integer(int32), intent(in) :: k
                !! Number of columns of A and rows of B.
            integer(int32), intent(in) :: lda
                !! Leading dimension of A.
            integer(int32), intent(in) :: ldb
                !! Leading dimension of B.
            integer(int32), intent(in) :: ldc
                !! Leading dimension of C.
            complex(real64), intent(in) :: alpha
                !! Scalar multiplier for the matrix product.
            complex(real64), intent(in) :: beta
                !! Scalar multiplier for the existing C matrix.
            complex(real64), intent(in) :: a(lda,*)
                !! First matrix operand.
            complex(real64), intent(in) :: b(ldb,*)
                !! Second matrix operand.
            complex(real64), intent(inout) :: c(ldc,*)
                !! Result matrix updated in place.
        end subroutine

        pure subroutine ZGEMV(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
            !! Computes y = alpha * op(A) * x + beta * y for complex matrices.
            use iso_fortran_env, only : int32, real64
            character, intent(in) :: trans
                !! Matrix operation to apply to A.
            integer(int32), intent(in) :: m
                !! Number of rows in the matrix operand.
            integer(int32), intent(in) :: n
                !! Number of columns in the matrix operand.
            integer(int32), intent(in) :: lda
                !! Leading dimension of A.
            integer(int32), intent(in) :: incx
                !! Stride for the x vector.
            integer(int32), intent(in) :: incy
                !! Stride for the y vector.
            complex(real64), intent(in) :: alpha
                !! Scalar multiplier for A*x.
            complex(real64), intent(in) :: beta
                !! Scalar multiplier for the existing y vector.
            complex(real64), intent(in) :: a(lda,*)
                !! Matrix operand.
            complex(real64), intent(in) :: x(*)
                !! Input vector.
            complex(real64), intent(inout) :: y(*)
                !! Updated vector.
        end subroutine

        pure subroutine DTRSM(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
            !! Solves a triangular system with one right-hand side in real arithmetic.
            use iso_fortran_env, only : int32, real64
            character, intent(in) :: side
                !! Whether the triangular matrix is on the left or right.
            character, intent(in) :: uplo
                !! Which triangle of A is referenced.
            character, intent(in) :: transa
                !! Whether A is transposed or conjugate-transposed.
            character, intent(in) :: diag
                !! Whether A is unit triangular.
            integer(int32), intent(in) :: m
                !! Number of rows of B.
            integer(int32), intent(in) :: n
                !! Number of columns of B.
            integer(int32), intent(in) :: lda
                !! Leading dimension of A.
            integer(int32), intent(in) :: ldb
                !! Leading dimension of B.
            real(real64), intent(in) :: alpha
                !! Scalar multiplier for the solution.
            real(real64), intent(in) :: a(lda,*)
                !! Triangular matrix factor.
            real(real64), intent(inout) :: b(ldb,*)
                !! Right-hand side matrix to be overwritten by the solution.
        end subroutine

        pure subroutine ZTRSM(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
            !! Solves a triangular system with one right-hand side in complex arithmetic.
            use iso_fortran_env, only : int32, real64
            character, intent(in) :: side
                !! Whether the triangular matrix is on the left or right.
            character, intent(in) :: uplo
                !! Which triangle of A is referenced.
            character, intent(in) :: transa
                !! Whether A is transposed or conjugate-transposed.
            character, intent(in) :: diag
                !! Whether A is unit triangular.
            integer(int32), intent(in) :: m
                !! Number of rows of B.
            integer(int32), intent(in) :: n
                !! Number of columns of B.
            integer(int32), intent(in) :: lda
                !! Leading dimension of A.
            integer(int32), intent(in) :: ldb
                !! Leading dimension of B.
            complex(real64), intent(in) :: alpha
                !! Scalar multiplier for the solution.
            complex(real64), intent(in) :: a(lda,*)
                !! Triangular matrix factor.
            complex(real64), intent(inout) :: b(ldb,*)
                !! Right-hand side matrix to be overwritten by the solution.
        end subroutine

        pure subroutine DTRSV(uplo, trans, diag, n, a, lda, x, incx)
            !! Solves a triangular system with a single real vector.
            use iso_fortran_env, only : int32, real64
            character, intent(in) :: uplo
                !! Which triangle of A is referenced.
            character, intent(in) :: trans
                !! Whether A is transposed.
            character, intent(in) :: diag
                !! Whether A is unit triangular.
            integer(int32), intent(in) :: n
                !! Length of the vector.
            integer(int32), intent(in) :: lda
                !! Leading dimension of A.
            integer(int32), intent(in) :: incx
                !! Stride in the x vector.
            real(real64), intent(in) :: a(lda,*)
                !! Triangular matrix factor.
            real(real64), intent(inout) :: x(*)
                !! Solution vector updated in place.
        end subroutine

        pure subroutine ZTRSV(uplo, trans, diag, n, a, lda, x, incx)
            !! Solves a triangular system with a single complex vector.
            use iso_fortran_env, only : int32, real64
            character, intent(in) :: uplo
                !! Which triangle of A is referenced.
            character, intent(in) :: trans
                !! Whether A is transposed.
            character, intent(in) :: diag
                !! Whether A is unit triangular.
            integer(int32), intent(in) :: n
                !! Length of the vector.
            integer(int32), intent(in) :: lda
                !! Leading dimension of A.
            integer(int32), intent(in) :: incx
                !! Stride in the x vector.
            complex(real64), intent(in) :: a(lda,*)
                !! Triangular matrix factor.
            complex(real64), intent(inout) :: x(*)
                !! Solution vector updated in place.
        end subroutine

        pure subroutine DSCAL(n, da, dx, incx)
            !! Scales a real vector by a real scalar.
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: n
                !! Length of the vector.
            integer(int32), intent(in) :: incx
                !! Stride in the vector.
            real(real64), intent(in) :: da
                !! Scalar used to scale the vector.
            real(real64), intent(inout) :: dx(*)
                !! Vector to scale in place.
        end subroutine

        pure subroutine ZSCAL(n, za, zx, incx)
            !! Scales a complex vector by a complex scalar.
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: n
                !! Length of the vector.
            integer(int32), intent(in) :: incx
                !! Stride in the vector.
            complex(real64), intent(in) :: za
                !! Complex scalar used to scale the vector.
            complex(real64), intent(inout) :: zx(*)
                !! Vector to scale in place.
        end subroutine

        pure subroutine ZDSCAL(n, da, zx, incx)
            !! Scales a complex vector by a real scalar.
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: n
                !! Length of the vector.
            integer(int32), intent(in) :: incx
                !! Stride in the vector.
            real(real64), intent(in) :: da
                !! Real scalar used to scale the entries.
            complex(real64), intent(inout) :: zx(*)
                !! Complex vector to scale in place.
        end subroutine

        pure subroutine DGBMV(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy)
            !! Computes a banded matrix-vector product for real matrices.
            use iso_fortran_env, only : int32, real64
            character, intent(in) :: trans
                !! Matrix operation to apply to A.
            integer(int32), intent(in) :: m
                !! Number of rows in the result.
            integer(int32), intent(in) :: n
                !! Number of columns in the matrix.
            integer(int32), intent(in) :: kl
                !! Number of sub-diagonals in the band.
            integer(int32), intent(in) :: ku
                !! Number of super-diagonals in the band.
            integer(int32), intent(in) :: lda
                !! Leading dimension of A.
            integer(int32), intent(in) :: incx
                !! Stride in the x vector.
            integer(int32), intent(in) :: incy
                !! Stride in the y vector.
            real(real64), intent(in) :: alpha
                !! Scalar multiplier for the matrix product.
            real(real64), intent(in) :: beta
                !! Scalar multiplier for the existing y vector.
            real(real64), intent(in) :: a(lda,*)
                !! Banded matrix operand.
            real(real64), intent(in) :: x(*)
                !! Input vector.
            real(real64), intent(inout) :: y(*)
                !! Updated result vector.
        end subroutine

        pure subroutine ZGBMV(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy)
            !! Computes a banded matrix-vector product for complex matrices.
            use iso_fortran_env, only : int32, real64
            character, intent(in) :: trans
                !! Matrix operation to apply to A.
            integer(int32), intent(in) :: m
                !! Number of rows in the result.
            integer(int32), intent(in) :: n
                !! Number of columns in the matrix.
            integer(int32), intent(in) :: kl
                !! Number of sub-diagonals in the band.
            integer(int32), intent(in) :: ku
                !! Number of super-diagonals in the band.
            integer(int32), intent(in) :: lda
                !! Leading dimension of A.
            integer(int32), intent(in) :: incx
                !! Stride in the x vector.
            integer(int32), intent(in) :: incy
                !! Stride in the y vector.
            complex(real64), intent(in) :: alpha
                !! Scalar multiplier for the matrix product.
            complex(real64), intent(in) :: beta
                !! Scalar multiplier for the existing y vector.
            complex(real64), intent(in) :: a(lda,*)
                !! Banded matrix operand.
            complex(real64), intent(in) :: x(*)
                !! Input vector.
            complex(real64), intent(inout) :: y(*)
                !! Updated result vector.
        end subroutine

        pure function DDOT(n, dx, incx, dy, incy) result(rst)
            !! Computes the dot product of two real vectors.
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: n
                !! Length of the vectors.
            integer(int32), intent(in) :: incx
                !! Stride in dx.
            integer(int32), intent(in) :: incy
                !! Stride in dy.
            real(real64), intent(in) :: dx(*)
                !! First vector operand.
            real(real64), intent(in) :: dy(*)
                !! Second vector operand.
            real(real64) :: rst
                !! Dot product result.
        end function

        pure subroutine DSWAP(n, dx, incx, dy, incy)
            !! Swaps two real vectors.
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: n
                !! Length of the vectors.
            integer(int32), intent(in) :: incx
                !! Stride in dx.
            integer(int32), intent(in) :: incy
                !! Stride in dy.
            real(real64), intent(inout) :: dx(*)
                !! First vector to swap.
            real(real64), intent(inout) :: dy(*)
                !! Second vector to swap.
        end subroutine
    end interface
end module