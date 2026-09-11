module lapack
    !! A module providing explicit interfaces to LAPACK routines.
    !! These declarations expose the factorization, eigensolver, and
    !! orthogonalization routines used throughout the library for decomposition
    !! and solving tasks.
    implicit none

    interface
        pure function DLAMCH(cmach) result(x)
            !! Returns a machine parameter value for the selected floating-point characteristic.
            use iso_fortran_env, only : real64
            character, intent(in) :: cmach
                !! Machine parameter identifier.
            real(real64) :: x
                !! Requested machine constant value.
        end function

        pure subroutine DGESVD(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, &
            work, lwork, info)
            !! Computes the singular value decomposition of a real matrix.
            use iso_fortran_env, only : int32, real64
            character, intent(in) :: jobu
                !! Controls the computation of U.
            character, intent(in) :: jobvt
                !! Controls the computation of V^T.
            integer(int32), intent(in) :: m
                !! Number of rows of A.
            integer(int32), intent(in) :: n
                !! Number of columns of A.
            integer(int32), intent(in) :: lda
                !! Leading dimension of A.
            integer(int32), intent(in) :: ldu
                !! Leading dimension of U.
            integer(int32), intent(in) :: ldvt
                !! Leading dimension of VT.
            integer(int32), intent(in) :: lwork
                !! Dimension of work.
            integer(int32), intent(out) :: info
                !! Status flag on exit.
            real(real64), intent(inout) :: a(lda,*)
                !! Matrix to decompose.
            real(real64), intent(out) :: s(*)
                !! Singular values.
            real(real64), intent(out) :: u(ldu,*)
                !! Left singular vectors.
            real(real64), intent(out) :: vt(ldvt,*)
                !! Right singular vectors in row form.
            real(real64), intent(out) :: work(*)
                !! Workspace array.
        end subroutine

        pure subroutine ZGESVD(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, &
            work, lwork, rwork, info)
            !! Computes the singular value decomposition of a complex matrix.
            use iso_fortran_env, only : int32, real64
            character, intent(in) :: jobu
                !! Controls the computation of U.
            character, intent(in) :: jobvt
                !! Controls the computation of V^H.
            integer(int32), intent(in) :: m
                !! Number of rows of A.
            integer(int32), intent(in) :: n
                !! Number of columns of A.
            integer(int32), intent(in) :: lda
                !! Leading dimension of A.
            integer(int32), intent(in) :: ldu
                !! Leading dimension of U.
            integer(int32), intent(in) :: ldvt
                !! Leading dimension of VT.
            integer(int32), intent(in) :: lwork
                !! Dimension of work.
            integer(int32), intent(out) :: info
                !! Status flag on exit.
            complex(real64), intent(inout) :: a(lda,*)
                !! Matrix to decompose.
            real(real64), intent(out) :: s(*)
                !! Singular values.
            complex(real64), intent(out) :: u(ldu,*)
                !! Left singular vectors.
            complex(real64), intent(out) :: vt(ldvt,*)
                !! Right singular vectors in row form.
            complex(real64), intent(out) :: work(*)
                !! Workspace array.
            real(real64), intent(out) :: rwork(*)
                !! Real workspace array.
        end subroutine

        pure subroutine DGETRF(m, n, a, lda, ipiv, info)
            !! Computes an LU factorization of a real matrix with partial pivoting.
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: m
                !! Number of rows of A.
            integer(int32), intent(in) :: n
                !! Number of columns of A.
            integer(int32), intent(in) :: lda
                !! Leading dimension of A.
            real(real64), intent(inout) :: a(lda,*)
                !! Matrix to factorize.
            integer(int32), intent(out) :: ipiv(*)
                !! Pivot indices from the factorization.
            integer(int32), intent(out) :: info
                !! Status flag on exit.
        end subroutine

        pure subroutine ZGETRF(m, n, a, lda, ipiv, info)
            !! Computes an LU factorization of a complex matrix with partial pivoting.
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: m
                !! Number of rows of A.
            integer(int32), intent(in) :: n
                !! Number of columns of A.
            integer(int32), intent(in) :: lda
                !! Leading dimension of A.
            complex(real64), intent(inout) :: a(lda,*)
                !! Matrix to factorize.
            integer(int32), intent(out) :: ipiv(*)
                !! Pivot indices from the factorization.
            integer(int32), intent(out) :: info
                !! Status flag on exit.
        end subroutine

        pure subroutine DSYEV(jobz, uplo, n, a, lda, w, work, lwork, info)
            !! Computes the eigenvalues and, optionally, eigenvectors of a real symmetric matrix.
            use iso_fortran_env, only : int32, real64
            character, intent(in) :: jobz
                !! Whether eigenvectors are required.
            character, intent(in) :: uplo
                !! Which triangle of A contains the matrix.
            integer(int32), intent(in) :: n
                !! Order of the matrix.
            integer(int32), intent(in) :: lda
                !! Leading dimension of A.
            integer(int32), intent(in) :: lwork
                !! Workspace size.
            real(real64), intent(inout) :: a(lda,*)
                !! Symmetric matrix on input, overwritten by eigenvectors.
            real(real64), intent(out) :: w(*)
                !! Eigenvalues of A.
            real(real64), intent(out) :: work(*)
                !! Workspace array.
            integer(int32), intent(out) :: info
                !! Status flag on exit.
        end subroutine

        pure subroutine DGEEV(jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, &
            work, lwork, info)
            !! Computes the eigenvalues and, optionally, eigenvectors of a real matrix.
            use iso_fortran_env, only : int32, real64
            character, intent(in) :: jobvl
                !! Whether left eigenvectors are computed.
            character, intent(in) :: jobvr
                !! Whether right eigenvectors are computed.
            integer(int32), intent(in) :: n
                !! Order of the matrix.
            integer(int32), intent(in) :: lda
                !! Leading dimension of A.
            integer(int32), intent(in) :: ldvl
                !! Leading dimension of VL.
            integer(int32), intent(in) :: ldvr
                !! Leading dimension of VR.
            integer(int32), intent(in) :: lwork
                !! Workspace size.
            real(real64), intent(inout) :: a(lda,*)
                !! Matrix whose eigenvalues are requested.
            real(real64), intent(out) :: wr(*)
                !! Real parts of eigenvalues.
            real(real64), intent(out) :: wi(*)
                !! Imaginary parts of eigenvalues.
            real(real64), intent(out) :: vl(ldvl,*)
                !! Left eigenvectors if requested.
            real(real64), intent(out) :: vr(ldvr,*)
                !! Right eigenvectors if requested.
            real(real64), intent(out) :: work(*)
                !! Workspace array.
            integer(int32), intent(out) :: info
                !! Status flag on exit.
        end subroutine

        pure subroutine DGGEV(jobvl, jobvr, n, a, lda, b, ldb, alphar, alphai, &
            beta, vl, ldvl, vr, ldvr, work, lwork, info)
            !! Computes the generalized eigenvalues of a real matrix pencil.
            use iso_fortran_env, only : int32, real64
            character, intent(in) :: jobvl
                !! Whether left eigenvectors are computed.
            character, intent(in) :: jobvr
                !! Whether right eigenvectors are computed.
            integer(int32), intent(in) :: n
                !! Order of the pencil matrices.
            integer(int32), intent(in) :: lda
                !! Leading dimension of A.
            integer(int32), intent(in) :: ldb
                !! Leading dimension of B.
            integer(int32), intent(in) :: ldvl
                !! Leading dimension of VL.
            integer(int32), intent(in) :: ldvr
                !! Leading dimension of VR.
            integer(int32), intent(in) :: lwork
                !! Workspace size.
            real(real64), intent(inout) :: a(lda,*)
                !! First matrix in the pencil.
            real(real64), intent(inout) :: b(ldb,*)
                !! Second matrix in the pencil.
            real(real64), intent(out) :: alphar(*)
                !! Real parts of generalized eigenvalues.
            real(real64), intent(out) :: alphai(*)
                !! Imaginary parts of generalized eigenvalues.
            real(real64), intent(out) :: beta(*)
                !! Scaling factors for generalized eigenvalues.
            real(real64), intent(out) :: vl(ldvl,*)
                !! Left eigenvectors if requested.
            real(real64), intent(out) :: vr(ldvr,*)
                !! Right eigenvectors if requested.
            real(real64), intent(out) :: work(*)
                !! Workspace array.
            integer(int32), intent(out) :: info
                !! Status flag on exit.
        end subroutine

        pure subroutine ZGEEV(jobvl, jobvr, n, a, lda, w, vl, ldvl, vr, ldvr, work, &
            lwork, rwork, info)
            !! Computes the eigenvalues and, optionally, eigenvectors of a complex matrix.
            use iso_fortran_env, only : int32, real64
            character, intent(in) :: jobvl
                !! Whether left eigenvectors are computed.
            character, intent(in) :: jobvr
                !! Whether right eigenvectors are computed.
            integer(int32), intent(in) :: n
                !! Order of the matrix.
            integer(int32), intent(in) :: lda
                !! Leading dimension of A.
            integer(int32), intent(in) :: ldvl
                !! Leading dimension of VL.
            integer(int32), intent(in) :: ldvr
                !! Leading dimension of VR.
            integer(int32), intent(in) :: lwork
                !! Workspace size.
            complex(real64), intent(inout) :: a(lda,*)
                !! Matrix whose eigenvalues are requested.
            complex(real64), intent(out) :: w(*)
                !! Eigenvalues.
            complex(real64), intent(out) :: vl(ldvl,*)
                !! Left eigenvectors if requested.
            complex(real64), intent(out) :: vr(ldvr,*)
                !! Right eigenvectors if requested.
            complex(real64), intent(out) :: work(*)
                !! Workspace array.
            real(real64), intent(out) :: rwork(*)
                !! Real workspace array.
            integer(int32), intent(out) :: info
                !! Status flag on exit.
        end subroutine

        pure subroutine DLASET(uplo, m, n, alpha, beta, a, lda)
            !! Fills a real matrix with a constant value on the chosen triangular part.
            use iso_fortran_env, only : int32, real64
            character, intent(in) :: uplo
                !! Which triangle of A is set.
            integer(int32), intent(in) :: m
                !! Number of rows of A.
            integer(int32), intent(in) :: n
                !! Number of columns of A.
            integer(int32), intent(in) :: lda
                !! Leading dimension of A.
            real(real64), intent(in) :: alpha
                !! Value assigned to the selected triangle.
            real(real64), intent(in) :: beta
                !! Value assigned to the complement of the selected triangle.
            real(real64), intent(out) :: a(lda,*)
                !! Matrix filled as requested.
        end subroutine

        pure subroutine DGEQRF(m, n, a, lda, tau, work, lwork, info)
            !! Computes the QR factorization of a real matrix.
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: m
                !! Number of rows of A.
            integer(int32), intent(in) :: n
                !! Number of columns of A.
            integer(int32), intent(in) :: lda
                !! Leading dimension of A.
            integer(int32), intent(in) :: lwork
                !! Workspace size.
            real(real64), intent(inout) :: a(lda,*)
                !! Matrix to factorize.
            real(real64), intent(out) :: tau(*)
                !! Scalar factors of the Householder reflectors.
            real(real64), intent(out) :: work(*)
                !! Workspace array.
            integer(int32), intent(out) :: info
                !! Status flag on exit.
        end subroutine

        pure subroutine ZGEQRF(m, n, a, lda, tau, work, lwork, info)
            !! Computes the QR factorization of a complex matrix.
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: m
                !! Number of rows of A.
            integer(int32), intent(in) :: n
                !! Number of columns of A.
            integer(int32), intent(in) :: lda
                !! Leading dimension of A.
            integer(int32), intent(in) :: lwork
                !! Workspace size.
            complex(real64), intent(inout) :: a(lda,*)
                !! Matrix to factorize.
            complex(real64), intent(out) :: tau(*)
                !! Scalar factors of the Householder reflectors.
            complex(real64), intent(out) :: work(*)
                !! Workspace array.
            integer(int32), intent(out) :: info
                !! Status flag on exit.
        end subroutine

        pure subroutine DGEQP3(m, n, a, lda, jpvt, tau, work, lwork, info)
            !! Computes the QR factorization with column pivoting for a real matrix.
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: m
                !! Number of rows of A.
            integer(int32), intent(in) :: n
                !! Number of columns of A.
            integer(int32), intent(in) :: lda
                !! Leading dimension of A.
            integer(int32), intent(in) :: lwork
                !! Workspace size.
            real(real64), intent(inout) :: a(lda,*)
                !! Matrix to factorize.
            integer(int32), intent(inout) :: jpvt(*)
                !! Column pivot indices.
            real(real64), intent(out) :: tau(*)
                !! Scalar factors of the Householder reflectors.
            real(real64), intent(out) :: work(*)
                !! Workspace array.
            integer(int32), intent(out) :: info
                !! Status flag on exit.
        end subroutine

        pure subroutine ZGEQP3(m, n, a, lda, jpvt, tau, work, lwork, rwork, info)
            !! Computes the QR factorization with column pivoting for a complex matrix.
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: m
                !! Number of rows of A.
            integer(int32), intent(in) :: n
                !! Number of columns of A.
            integer(int32), intent(in) :: lda
                !! Leading dimension of A.
            integer(int32), intent(in) :: lwork
                !! Workspace size.
            complex(real64), intent(inout) :: a(lda,*)
                !! Matrix to factorize.
            integer(int32), intent(inout) :: jpvt(*)
                !! Column pivot indices.
            complex(real64), intent(out) :: tau(*)
                !! Scalar factors of the Householder reflectors.
            complex(real64), intent(out) :: work(*)
                !! Workspace array.
            real(real64), intent(out) :: rwork(*)
                !! Real workspace array.
            integer(int32), intent(out) :: info
                !! Status flag on exit.
        end subroutine

        pure subroutine DORGQR(m, n, k, a, lda, tau, work, lwork, info)
            !! Generates the real orthogonal matrix Q from a QR factorization.
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: m
                !! Number of rows of the generated Q.
            integer(int32), intent(in) :: n
                !! Number of columns of the generated Q.
            integer(int32), intent(in) :: k
                !! Number of elementary reflectors.
            integer(int32), intent(in) :: lda
                !! Leading dimension of A.
            integer(int32), intent(in) :: lwork
                !! Workspace size.
            real(real64), intent(inout) :: a(lda,*)
                !! QR factor data overwritten by Q.
            real(real64), intent(in) :: tau(*)
                !! Householder scalars from the QR factorization.
            real(real64), intent(out) :: work(*)
                !! Workspace array.
            integer(int32), intent(out) :: info
                !! Status flag on exit.
        end subroutine

        pure subroutine ZUNGQR(m, n, k, a, lda, tau, work, lwork, info)
            !! Generates the complex unitary matrix Q from a QR factorization.
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: m
                !! Number of rows of the generated Q.
            integer(int32), intent(in) :: n
                !! Number of columns of the generated Q.
            integer(int32), intent(in) :: k
                !! Number of elementary reflectors.
            integer(int32), intent(in) :: lda
                !! Leading dimension of A.
            integer(int32), intent(in) :: lwork
                !! Workspace size.
            complex(real64), intent(inout) :: a(lda,*)
                !! QR factor data overwritten by Q.
            complex(real64), intent(in) :: tau(*)
                !! Householder scalars from the QR factorization.
            complex(real64), intent(out) :: work(*)
                !! Workspace array.
            integer(int32), intent(out) :: info
                !! Status flag on exit.
        end subroutine

        pure subroutine DORMQR(side, trans, m, n, k, a, lda, tau, c, ldc, work, &
            lwork, info)
            !! Applies a real orthogonal matrix Q to a matrix from a QR factorization.
            use iso_fortran_env, only : int32, real64
            character, intent(in) :: side
                !! Whether Q multiplies from the left or the right.
            character, intent(in) :: trans
                !! Whether Q is transposed.
            integer(int32), intent(in) :: m
                !! Number of rows of C.
            integer(int32), intent(in) :: n
                !! Number of columns of C.
            integer(int32), intent(in) :: k
                !! Number of elementary reflectors.
            integer(int32), intent(in) :: lda
                !! Leading dimension of A.
            integer(int32), intent(in) :: ldc
                !! Leading dimension of C.
            integer(int32), intent(in) :: lwork
                !! Workspace size.
            real(real64), intent(in) :: a(lda,*)
                !! Householder data from the QR factorization.
            real(real64), intent(in) :: tau(*)
                !! Householder scalars.
            real(real64), intent(inout) :: c(ldc,*)
                !! Matrix to which Q is applied.
            real(real64), intent(out) :: work(*)
                !! Workspace array.
            integer(int32), intent(out) :: info
                !! Status flag on exit.
        end subroutine

        pure subroutine ZUNMQR(side, trans, m, n, k, a, lda, tau, c, ldc, work, &
            lwork, info)
            !! Applies a complex unitary matrix Q to a matrix from a QR factorization.
            use iso_fortran_env, only : int32, real64
            character, intent(in) :: side
                !! Whether Q multiplies from the left or the right.
            character, intent(in) :: trans
                !! Whether Q is transposed.
            integer(int32), intent(in) :: m
                !! Number of rows of C.
            integer(int32), intent(in) :: n
                !! Number of columns of C.
            integer(int32), intent(in) :: k
                !! Number of elementary reflectors.
            integer(int32), intent(in) :: lda
                !! Leading dimension of A.
            integer(int32), intent(in) :: ldc
                !! Leading dimension of C.
            integer(int32), intent(in) :: lwork
                !! Workspace size.
            complex(real64), intent(in) :: a(lda,*)
                !! Householder data from the QR factorization.
            complex(real64), intent(in) :: tau(*)
                !! Householder scalars.
            complex(real64), intent(inout) :: c(ldc,*)
                !! Matrix to which Q is applied.
            complex(real64), intent(out) :: work(*)
                !! Workspace array.
            integer(int32), intent(out) :: info
                !! Status flag on exit.
        end subroutine

        pure subroutine DPOTRF(uplo, n, a, lda, info)
            !! Computes the Cholesky factorization of a real symmetric positive-definite matrix.
            use iso_fortran_env, only : int32, real64
            character, intent(in) :: uplo
                !! Which triangle of A contains the matrix.
            integer(int32), intent(in) :: n
                !! Order of the matrix.
            integer(int32), intent(in) :: lda
                !! Leading dimension of A.
            real(real64), intent(inout) :: a(lda,*)
                !! Symmetric positive-definite matrix to factorize.
            integer(int32), intent(out) :: info
                !! Status flag on exit.
        end subroutine

        pure subroutine ZPOTRF(uplo, n, a, lda, info)
            !! Computes the Cholesky factorization of a complex Hermitian positive-definite matrix.
            use iso_fortran_env, only : int32, real64
            character, intent(in) :: uplo
                !! Which triangle of A contains the matrix.
            integer(int32), intent(in) :: n
                !! Order of the matrix.
            integer(int32), intent(in) :: lda
                !! Leading dimension of A.
            complex(real64), intent(inout) :: a(lda,*)
                !! Hermitian positive-definite matrix to factorize.
            integer(int32), intent(out) :: info
                !! Status flag on exit.
        end subroutine

        pure subroutine DTZRZF(m, n, a, lda, tau, work, lwork, info)
            !! Computes the RZ factorization of a real matrix using Householder reflectors.
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: m
                !! Number of rows of A.
            integer(int32), intent(in) :: n
                !! Number of columns of A.
            integer(int32), intent(in) :: lda
                !! Leading dimension of A.
            integer(int32), intent(in) :: lwork
                !! Workspace size.
            real(real64), intent(inout) :: a(lda,*)
                !! Matrix to factorize.
            real(real64), intent(out) :: tau(*)
                !! Householder scalars.
            real(real64), intent(out) :: work(*)
                !! Workspace array.
            integer(int32), intent(out) :: info
                !! Status flag on exit.
        end subroutine

        pure subroutine ZTZRZF(m, n, a, lda, tau, work, lwork, info)
            !! Computes the RZ factorization of a complex matrix using Householder reflectors.
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: m
                !! Number of rows of A.
            integer(int32), intent(in) :: n
                !! Number of columns of A.
            integer(int32), intent(in) :: lda
                !! Leading dimension of A.
            integer(int32), intent(in) :: lwork
                !! Workspace size.
            complex(real64), intent(inout) :: a(lda,*)
                !! Matrix to factorize.
            complex(real64), intent(out) :: tau(*)
                !! Householder scalars.
            complex(real64), intent(out) :: work(*)
                !! Workspace array.
            integer(int32), intent(out) :: info
                !! Status flag on exit.
        end subroutine

        pure subroutine DORMRZ(side, trans, m, n, k, l, a, lda, tau, c, ldc, work, &
            lwork, info)
            !! Applies the orthogonal matrix from an RZ factorization to a real matrix.
            use iso_fortran_env, only : int32, real64
            character, intent(in) :: side
                !! Whether Q multiplies from the left or the right.
            character, intent(in) :: trans
                !! Whether Q is transposed.
            integer(int32), intent(in) :: m
                !! Number of rows of C.
            integer(int32), intent(in) :: n
                !! Number of columns of C.
            integer(int32), intent(in) :: k
                !! Number of rows of the reflected matrix.
            integer(int32), intent(in) :: l
                !! Number of columns of the reflected matrix.
            integer(int32), intent(in) :: lda
                !! Leading dimension of A.
            integer(int32), intent(in) :: ldc
                !! Leading dimension of C.
            integer(int32), intent(in) :: lwork
                !! Workspace size.
            real(real64), intent(in) :: a(lda,*)
                !! Householder data from the RZ factorization.
            real(real64), intent(in) :: tau(*)
                !! Householder scalars.
            real(real64), intent(inout) :: c(ldc,*)
                !! Matrix to which the orthogonal factor is applied.
            real(real64), intent(out) :: work(*)
                !! Workspace array.
            integer(int32), intent(out) :: info
                !! Status flag on exit.
        end subroutine

        pure subroutine zunmrz(side, trans, m, n, k, l, a, lda, tau, c, ldc, work, &
            lwork, info)
            !! Applies the unitary matrix from a complex RZ factorization to a matrix.
            use iso_fortran_env, only : int32, real64
            character, intent(in) :: side
                !! Whether Q multiplies from the left or the right.
            character, intent(in) :: trans
                !! Whether Q is transposed.
            integer(int32), intent(in) :: m
                !! Number of rows of C.
            integer(int32), intent(in) :: n
                !! Number of columns of C.
            integer(int32), intent(in) :: k
                !! Number of rows of the reflected matrix.
            integer(int32), intent(in) :: l
                !! Number of columns of the reflected matrix.
            integer(int32), intent(in) :: lda
                !! Leading dimension of A.
            integer(int32), intent(in) :: ldc
                !! Leading dimension of C.
            integer(int32), intent(in) :: lwork
                !! Workspace size.
            complex(real64), intent(in) :: a(lda,*)
                !! Householder data from the RZ factorization.
            complex(real64), intent(in) :: tau(*)
                !! Householder scalars.
            complex(real64), intent(inout) :: c(ldc,*)
                !! Matrix to which the unitary factor is applied.
            complex(real64), intent(out) :: work(*)
                !! Workspace array.
            integer(int32), intent(out) :: info
                !! Status flag on exit.
        end subroutine

        pure subroutine DGELQF(m, n, a, lda, tau, work, lwork, info)
            !! Computes the LQ factorization of a real matrix.
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: m
                !! Number of rows of A.
            integer(int32), intent(in) :: n
                !! Number of columns of A.
            integer(int32), intent(in) :: lda
                !! Leading dimension of A.
            integer(int32), intent(in) :: lwork
                !! Workspace size.
            real(real64), intent(inout) :: a(lda,*)
                !! Matrix to factorize.
            real(real64), intent(out) :: tau(*)
                !! Householder scalars.
            real(real64), intent(out) :: work(*)
                !! Workspace array.
            integer(int32), intent(out) :: info
                !! Status flag on exit.
        end subroutine

        pure subroutine ZGELQF(m, n, a, lda, tau, work, lwork, info)
            !! Computes the LQ factorization of a complex matrix.
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: m
                !! Number of rows of A.
            integer(int32), intent(in) :: n
                !! Number of columns of A.
            integer(int32), intent(in) :: lda
                !! Leading dimension of A.
            integer(int32), intent(in) :: lwork
                !! Workspace size.
            complex(real64), intent(inout) :: a(lda,*)
                !! Matrix to factorize.
            complex(real64), intent(out) :: tau(*)
                !! Householder scalars.
            complex(real64), intent(out) :: work(*)
                !! Workspace array.
            integer(int32), intent(out) :: info
                !! Status flag on exit.
        end subroutine

        pure subroutine DORGLQ(m, n, k, a, lda, tau, work, lwork, info)
            !! Generates the orthogonal matrix Q from an LQ factorization of a real matrix.
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: m
                !! Number of rows of the generated Q.
            integer(int32), intent(in) :: n
                !! Number of columns of the generated Q.
            integer(int32), intent(in) :: k
                !! Number of elementary reflectors.
            integer(int32), intent(in) :: lda
                !! Leading dimension of A.
            integer(int32), intent(in) :: lwork
                !! Workspace size.
            real(real64), intent(inout) :: a(lda,*)
                !! LQ factor data overwritten by Q.
            real(real64), intent(in) :: tau(*)
                !! Householder scalars from the LQ factorization.
            real(real64), intent(out) :: work(*)
                !! Workspace array.
            integer(int32), intent(out) :: info
                !! Status flag on exit.
        end subroutine

        pure subroutine ZUNGLQ(m, n, k, a, lda, tau, work, lwork, info)
            !! Generates the unitary matrix Q from an LQ factorization of a complex matrix.
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: m
                !! Number of rows of the generated Q.
            integer(int32), intent(in) :: n
                !! Number of columns of the generated Q.
            integer(int32), intent(in) :: k
                !! Number of elementary reflectors.
            integer(int32), intent(in) :: lda
                !! Leading dimension of A.
            integer(int32), intent(in) :: lwork
                !! Workspace size.
            complex(real64), intent(inout) :: a(lda,*)
                !! LQ factor data overwritten by Q.
            complex(real64), intent(in) :: tau(*)
                !! Householder scalars from the LQ factorization.
            complex(real64), intent(out) :: work(*)
                !! Workspace array.
            integer(int32), intent(out) :: info
                !! Status flag on exit.
        end subroutine

        pure subroutine DORMLQ(side, trans, m, n, k, a, lda, tau, c, ldc, work, &
            lwork, info)
            !! Applies the orthogonal matrix from an LQ factorization to a real matrix.
            use iso_fortran_env, only : int32, real64
            character, intent(in) :: side
                !! Whether Q multiplies from the left or the right.
            character, intent(in) :: trans
                !! Whether the orthogonal matrix is transposed.
            integer(int32), intent(in) :: m
                !! Number of rows of C.
            integer(int32), intent(in) :: n
                !! Number of columns of C.
            integer(int32), intent(in) :: k
                !! Number of elementary reflectors.
            integer(int32), intent(in) :: lda
                !! Leading dimension of A.
            integer(int32), intent(in) :: ldc
                !! Leading dimension of C.
            integer(int32), intent(in) :: lwork
                !! Workspace size.
            real(real64), intent(in) :: a(lda,*)
                !! Householder data from the LQ factorization.
            real(real64), intent(in) :: tau(*)
                !! Householder scalars.
            real(real64), intent(inout) :: c(ldc,*)
                !! Matrix to which Q is applied.
            real(real64), intent(out) :: work(*)
                !! Workspace array.
            integer(int32), intent(out) :: info
                !! Status flag on exit.
        end subroutine

        pure subroutine ZUNMLQ(side, trans, m, n, k, a, lda, tau, c, ldc, work, &
            lwork, info)
            !! Applies the unitary matrix from an LQ factorization to a complex matrix.
            use iso_fortran_env, only : int32, real64
            character, intent(in) :: side
                !! Whether Q multiplies from the left or the right.
            character, intent(in) :: trans
                !! Whether the unitary matrix is transposed.
            integer(int32), intent(in) :: m
                !! Number of rows of C.
            integer(int32), intent(in) :: n
                !! Number of columns of C.
            integer(int32), intent(in) :: k
                !! Number of elementary reflectors.
            integer(int32), intent(in) :: lda
                !! Leading dimension of A.
            integer(int32), intent(in) :: ldc
                !! Leading dimension of C.
            integer(int32), intent(in) :: lwork
                !! Workspace size.
            complex(real64), intent(in) :: a(lda,*)
                !! Householder data from the LQ factorization.
            complex(real64), intent(in) :: tau(*)
                !! Householder scalars.
            complex(real64), intent(inout) :: c(ldc,*)
                !! Matrix to which Q is applied.
            complex(real64), intent(out) :: work(*)
                !! Workspace array.
            integer(int32), intent(out) :: info
                !! Status flag on exit.
        end subroutine

        pure subroutine DGETRS(trans, n, nrhs, a, lda, ipiv, b, ldb, info)
            !! Solves a real system of linear equations using an LU factorization.
            use iso_fortran_env, only : int32, real64
            character, intent(in) :: trans
                !! Whether the matrix is transposed or not.
            integer(int32), intent(in) :: n
                !! Order of the system.
            integer(int32), intent(in) :: nrhs
                !! Number of right-hand sides.
            integer(int32), intent(in) :: lda
                !! Leading dimension of A.
            integer(int32), intent(in) :: ldb
                !! Leading dimension of B.
            real(real64), intent(in) :: a(lda,*)
                !! LU factorization of the coefficient matrix.
            integer(int32), intent(in) :: ipiv(*)
                !! Pivot indices from the factorization.
            real(real64), intent(inout) :: b(ldb,*)
                !! Right-hand side matrix overwritten by the solution.
            integer(int32), intent(out) :: info
                !! Status flag on exit.
        end subroutine

        pure subroutine ZGETRS(trans, n, nrhs, a, lda, ipiv, b, ldb, info)
            !! Solves a complex system of linear equations using an LU factorization.
            use iso_fortran_env, only : int32, real64
            character, intent(in) :: trans
                !! Whether the matrix is transposed or not.
            integer(int32), intent(in) :: n
                !! Order of the system.
            integer(int32), intent(in) :: nrhs
                !! Number of right-hand sides.
            integer(int32), intent(in) :: lda
                !! Leading dimension of A.
            integer(int32), intent(in) :: ldb
                !! Leading dimension of B.
            complex(real64), intent(in) :: a(lda,*)
                !! LU factorization of the coefficient matrix.
            integer(int32), intent(in) :: ipiv(*)
                !! Pivot indices from the factorization.
            complex(real64), intent(inout) :: b(ldb,*)
                !! Right-hand side matrix overwritten by the solution.
            integer(int32), intent(out) :: info
                !! Status flag on exit.
        end subroutine

        pure subroutine DLAIC1(job, j, x, sest, w, gamma, sestpr, s, c)
            !! Solves a small scalar update associated with a real rank-one modification.
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: job
                !! Selects the update mode.
            integer(int32), intent(in) :: j
                !! Index of the modified element.
            real(real64), intent(in) :: x(j)
                !! Vector participating in the rank-one correction.
            real(real64), intent(in) :: w(j)
                !! Auxiliary vector in the scalar update.
            real(real64), intent(in) :: sest
                !! Current estimate of the singular value.
            real(real64), intent(in) :: gamma
                !! Scalar defining the rank-one perturbation.
            real(real64), intent(out) :: sestpr
                !! Updated scalar estimate.
            real(real64), intent(out) :: s
                !! First scalar output of the update.
            real(real64), intent(out) :: c
                !! Second scalar output of the update.
        end subroutine

        pure subroutine ZLAIC1(job, j, x, sest, w, gamma, sestpr, s, c)
            !! Solves a small scalar update associated with a complex rank-one modification.
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: job
                !! Selects the update mode.
            integer(int32), intent(in) :: j
                !! Index of the modified element.
            complex(real64), intent(in) :: x(j)
                !! Vector participating in the rank-one correction.
            complex(real64), intent(in) :: w(j)
                !! Auxiliary vector in the scalar update.
            complex(real64), intent(in) :: gamma
                !! Scalar defining the rank-one perturbation.
            real(real64), intent(in) :: sest
                !! Current estimate of the singular value.
            real(real64), intent(out) :: sestpr
                !! Updated scalar estimate.
            complex(real64), intent(out) :: s
                !! First scalar output of the update.
            complex(real64), intent(out) :: c
                !! Second scalar output of the update.
        end subroutine

        pure subroutine DPOTRS(uplo, n, nrhs, a, lda, b, ldb, info)
            !! Solves a real symmetric positive-definite system using a Cholesky factorization.
            use iso_fortran_env, only : int32, real64
            character, intent(in) :: uplo
                !! Which triangle of A contains the Cholesky factor.
            integer(int32), intent(in) :: n
                !! Order of the system.
            integer(int32), intent(in) :: nrhs
                !! Number of right-hand sides.
            integer(int32), intent(in) :: lda
                !! Leading dimension of A.
            integer(int32), intent(in) :: ldb
                !! Leading dimension of B.
            real(real64), intent(in) :: a(lda,*)
                !! Cholesky factorization of the coefficient matrix.
            real(real64), intent(inout) :: b(ldb,*)
                !! Right-hand side matrix overwritten by the solution.
            integer(int32), intent(out) :: info
                !! Status flag on exit.
        end subroutine

        pure subroutine ZPOTRS(uplo, n, nrhs, a, lda, b, ldb, info)
            !! Solves a complex Hermitian positive-definite system using a Cholesky factorization.
            use iso_fortran_env, only : int32, real64
            character, intent(in) :: uplo
                !! Which triangle of A contains the Cholesky factor.
            integer(int32), intent(in) :: n
                !! Order of the system.
            integer(int32), intent(in) :: nrhs
                !! Number of right-hand sides.
            integer(int32), intent(in) :: lda
                !! Leading dimension of A.
            integer(int32), intent(in) :: ldb
                !! Leading dimension of B.
            complex(real64), intent(in) :: a(lda,*)
                !! Cholesky factorization of the coefficient matrix.
            complex(real64), intent(inout) :: b(ldb,*)
                !! Right-hand side matrix overwritten by the solution.
            integer(int32), intent(out) :: info
                !! Status flag on exit.
        end subroutine

        pure subroutine DGETRI(n, a, lda, ipiv, work, lwork, info)
            !! Computes the inverse of a real matrix from its LU factorization.
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: n
                !! Order of the matrix.
            integer(int32), intent(in) :: lda
                !! Leading dimension of A.
            integer(int32), intent(in) :: lwork
                !! Workspace size.
            integer(int32), intent(in) :: ipiv(*)
                !! Pivot indices from the LU factorization.
            real(real64), intent(inout) :: a(lda,*)
                !! Matrix to invert, overwritten by its inverse.
            real(real64), intent(out) :: work(*)
                !! Workspace array.
            integer(int32), intent(out) :: info
                !! Status flag on exit.
        end subroutine

        pure subroutine ZGETRI(n, a, lda, ipiv, work, lwork, info)
            !! Computes the inverse of a complex matrix from its LU factorization.
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: n
                !! Order of the matrix.
            integer(int32), intent(in) :: lda
                !! Leading dimension of A.
            integer(int32), intent(in) :: lwork
                !! Workspace size.
            integer(int32), intent(in) :: ipiv(*)
                !! Pivot indices from the LU factorization.
            complex(real64), intent(inout) :: a(lda,*)
                !! Matrix to invert, overwritten by its inverse.
            complex(real64), intent(out) :: work(*)
                !! Workspace array.
            integer(int32), intent(out) :: info
                !! Status flag on exit.
        end subroutine

        pure subroutine DGELS(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info)
            !! Solves a real least-squares problem using QR or LQ factorization.
            use iso_fortran_env, only : int32, real64
            character, intent(in) :: trans
                !! Whether A is transposed or not.
            integer(int32), intent(in) :: m
                !! Number of rows of A.
            integer(int32), intent(in) :: n
                !! Number of columns of A.
            integer(int32), intent(in) :: nrhs
                !! Number of right-hand sides.
            integer(int32), intent(in) :: lda
                !! Leading dimension of A.
            integer(int32), intent(in) :: ldb
                !! Leading dimension of B.
            integer(int32), intent(in) :: lwork
                !! Workspace size.
            real(real64), intent(inout) :: a(lda,*)
                !! Matrix in the least-squares problem.
            real(real64), intent(inout) :: b(ldb,*)
                !! Right-hand side matrix overwritten by the solution.
            real(real64), intent(out) :: work(*)
                !! Workspace array.
            integer(int32), intent(out) :: info
                !! Status flag on exit.
        end subroutine

        pure subroutine ZGELS(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info)
            !! Solves a complex least-squares problem using QR or LQ factorization.
            use iso_fortran_env, only : int32, real64
            character, intent(in) :: trans
                !! Whether A is transposed or not.
            integer(int32), intent(in) :: m
                !! Number of rows of A.
            integer(int32), intent(in) :: n
                !! Number of columns of A.
            integer(int32), intent(in) :: nrhs
                !! Number of right-hand sides.
            integer(int32), intent(in) :: lda
                !! Leading dimension of A.
            integer(int32), intent(in) :: ldb
                !! Leading dimension of B.
            integer(int32), intent(in) :: lwork
                !! Workspace size.
            complex(real64), intent(inout) :: a(lda,*)
                !! Matrix in the least-squares problem.
            complex(real64), intent(inout) :: b(ldb,*)
                !! Right-hand side matrix overwritten by the solution.
            complex(real64), intent(out) :: work(*)
                !! Workspace array.
            integer(int32), intent(out) :: info
                !! Status flag on exit.
        end subroutine

        pure subroutine DGELSY(m, n, nrhs, a, lda, b, ldb, jpvt, rcond, rank, work, &
            lwork, info)
            !! Solves a real rank-deficient least-squares problem using a complete orthogonal factorization.
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: m
                !! Number of rows of A.
            integer(int32), intent(in) :: n
                !! Number of columns of A.
            integer(int32), intent(in) :: nrhs
                !! Number of right-hand sides.
            integer(int32), intent(in) :: lda
                !! Leading dimension of A.
            integer(int32), intent(in) :: ldb
                !! Leading dimension of B.
            integer(int32), intent(in) :: lwork
                !! Workspace size.
            real(real64), intent(inout) :: a(lda,*)
                !! Matrix in the least-squares problem.
            real(real64), intent(inout) :: b(ldb,*)
                !! Right-hand side matrix overwritten by the solution.
            integer(int32), intent(inout) :: jpvt(*)
                !! Column pivoting information.
            real(real64), intent(in) :: rcond
                !! Relative rank tolerance.
            integer(int32), intent(out) :: rank
                !! Effective rank of the matrix.
            real(real64), intent(out) :: work(*)
                !! Workspace array.
            integer(int32), intent(out) :: info
                !! Status flag on exit.
        end subroutine

        pure subroutine ZGELSY(m, n, nrhs, a, lda, b, ldb, jpvt, rcond, rank, work, &
            lwork, rwork, info)
            !! Solves a complex rank-deficient least-squares problem using a complete orthogonal factorization.
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: m
                !! Number of rows of A.
            integer(int32), intent(in) :: n
                !! Number of columns of A.
            integer(int32), intent(in) :: nrhs
                !! Number of right-hand sides.
            integer(int32), intent(in) :: lda
                !! Leading dimension of A.
            integer(int32), intent(in) :: ldb
                !! Leading dimension of B.
            integer(int32), intent(in) :: lwork
                !! Workspace size.
            complex(real64), intent(inout) :: a(lda,*)
                !! Matrix in the least-squares problem.
            complex(real64), intent(inout) :: b(ldb,*)
                !! Right-hand side matrix overwritten by the solution.
            integer(int32), intent(inout) :: jpvt(*)
                !! Column pivoting information.
            real(real64), intent(in) :: rcond
                !! Relative rank tolerance.
            integer(int32), intent(out) :: rank
                !! Effective rank of the matrix.
            complex(real64), intent(out) :: work(*)
                !! Workspace array.
            real(real64), intent(out) :: rwork(*)
                !! Real workspace array.
            integer(int32), intent(out) :: info
                !! Status flag on exit.
        end subroutine

        pure subroutine DGELSS(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, &
            lwork, info)
            !! Solves a real least-squares problem using the SVD-based method.
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: m
                !! Number of rows of A.
            integer(int32), intent(in) :: n
                !! Number of columns of A.
            integer(int32), intent(in) :: nrhs
                !! Number of right-hand sides.
            integer(int32), intent(in) :: lda
                !! Leading dimension of A.
            integer(int32), intent(in) :: ldb
                !! Leading dimension of B.
            integer(int32), intent(in) :: lwork
                !! Workspace size.
            real(real64), intent(inout) :: a(lda,*)
                !! Matrix in the least-squares problem.
            real(real64), intent(inout) :: b(ldb,*)
                !! Right-hand side matrix overwritten by the solution.
            real(real64), intent(out) :: s(*)
                !! Singular values of A.
            real(real64), intent(in) :: rcond
                !! Relative tolerance for rank determination.
            integer(int32), intent(out) :: rank
                !! Effective rank.
            real(real64), intent(out) :: work(*)
                !! Workspace array.
            integer(int32), intent(out) :: info
                !! Status flag on exit.
        end subroutine

        pure subroutine ZGELSS(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, &
            lwork, rwork, info)
            !! Solves a complex least-squares problem using the SVD-based method.
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: m
                !! Number of rows of A.
            integer(int32), intent(in) :: n
                !! Number of columns of A.
            integer(int32), intent(in) :: nrhs
                !! Number of right-hand sides.
            integer(int32), intent(in) :: lda
                !! Leading dimension of A.
            integer(int32), intent(in) :: ldb
                !! Leading dimension of B.
            integer(int32), intent(in) :: lwork
                !! Workspace size.
            complex(real64), intent(inout) :: a(lda,*)
                !! Matrix in the least-squares problem.
            complex(real64), intent(inout) :: b(ldb,*)
                !! Right-hand side matrix overwritten by the solution.
            real(real64), intent(out) :: s(*)
                !! Singular values of A.
            real(real64), intent(in) :: rcond
                !! Relative tolerance for rank determination.
            integer(int32), intent(out) :: rank
                !! Effective rank.
            complex(real64), intent(out) :: work(*)
                !! Workspace array.
            real(real64), intent(out) :: rwork(*)
                !! Real workspace array.
            integer(int32), intent(out) :: info
                !! Status flag on exit.
        end subroutine

        pure subroutine DLASRT(id, n, d, info)
            !! Sorts a real array of values into ascending or descending order.
            use iso_fortran_env, only : int32, real64
            character, intent(in) :: id
                !! Sort order indicator.
            integer(int32), intent(in) :: n
                !! Length of the array.
            real(real64), intent(inout) :: d(*)
                !! Array to sort in place.
            integer(int32), intent(out) :: info
                !! Status flag on exit.
        end subroutine

        pure subroutine DGESV(n, nrhs, a, lda, ipiv, b, ldb, info)
            !! Solves a real system of linear equations with partial pivoting.
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: n
                !! Order of the matrix.
            integer(int32), intent(in) :: nrhs
                !! Number of right-hand sides.
            integer(int32), intent(in) :: lda
                !! Leading dimension of A.
            integer(int32), intent(in) :: ldb
                !! Leading dimension of B.
            real(real64), intent(inout) :: a(lda,*)
                !! Coefficient matrix to factorize in place.
            integer(int32), intent(out) :: ipiv(*)
                !! Pivot indices.
            real(real64), intent(inout) :: b(ldb,*)
                !! Right-hand side matrix overwritten by the solution.
            integer(int32), intent(out) :: info
                !! Status flag on exit.
        end subroutine

        pure subroutine DGGEV3(jobvl, jobvr, n, a, lda, b, ldb, alphar, alphai, &
            beta, vl, ldvl, vr, ldvr, work, lwork, info)
            !! Computes the generalized eigenvalues of a real matrix pencil with balancing.
            use iso_fortran_env, only : int32, real64
            character, intent(in) :: jobvl
                !! Whether left eigenvectors are computed.
            character, intent(in) :: jobvr
                !! Whether right eigenvectors are computed.
            integer(int32), intent(in) :: n
                !! Order of the pencil matrices.
            integer(int32), intent(in) :: lda
                !! Leading dimension of A.
            integer(int32), intent(in) :: ldb
                !! Leading dimension of B.
            integer(int32), intent(in) :: ldvl
                !! Leading dimension of VL.
            integer(int32), intent(in) :: ldvr
                !! Leading dimension of VR.
            integer(int32), intent(in) :: lwork
                !! Workspace size.
            real(real64), intent(inout) :: a(lda,*)
                !! First matrix in the pencil.
            real(real64), intent(inout) :: b(ldb,*)
                !! Second matrix in the pencil.
            real(real64), intent(out) :: alphar(*)
                !! Real parts of generalized eigenvalues.
            real(real64), intent(out) :: alphai(*)
                !! Imaginary parts of generalized eigenvalues.
            real(real64), intent(out) :: beta(*)
                !! Scaling factors for generalized eigenvalues.
            real(real64), intent(out) :: vl(ldvl,*)
                !! Left eigenvectors if requested.
            real(real64), intent(out) :: vr(ldvr,*)
                !! Right eigenvectors if requested.
            real(real64), intent(out) :: work(*)
                !! Workspace array.
            integer(int32), intent(out) :: info
                !! Status flag on exit.
        end subroutine

        pure subroutine DLADIV(a, b, c, d, p, q)
            !! Computes the quotient of complex numbers represented as real pairs.
            use iso_fortran_env, only : real64
            real(real64), intent(in) :: a
                !! Real part of the first complex number.
            real(real64), intent(in) :: b
                !! Imaginary part of the first complex number.
            real(real64), intent(in) :: c
                !! Real part of the second complex number.
            real(real64), intent(in) :: d
                !! Imaginary part of the second complex number.
            real(real64), intent(out) :: p
                !! Real part of the quotient.
            real(real64), intent(out) :: q
                !! Imaginary part of the quotient.
        end subroutine
    end interface
end module