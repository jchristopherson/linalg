module lapack
    implicit none

    interface
        pure function DLAMCH(cmach) result(x)
            use iso_fortran_env, only : real64
            character, intent(in) :: cmach
            real(real64) :: x
        end function

        pure subroutine DGESVD(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, &
            work, lwork, info)
            use iso_fortran_env, only : int32, real64
            character, intent(in) :: jobu, jobvt
            integer(int32), intent(in) :: m, n, lda, ldu, ldvt, lwork
            integer(int32), intent(out) :: info
            real(real64), intent(inout) :: a(lda,*)
            real(real64), intent(out) :: s(*), u(ldu,*), vt(ldvt,*), work(*)
        end subroutine

        pure subroutine ZGESVD(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, &
            work, lwork, rwork, info)
            use iso_fortran_env, only : int32, real64
            character, intent(in) :: jobu, jobvt
            integer(int32), intent(in) :: m, n, lda, ldu, ldvt, lwork
            integer(int32), intent(out) :: info
            complex(real64), intent(inout) :: a(lda,*)
            real(real64), intent(out) :: s(*), rwork(*)
            complex(real64), intent(out) :: u(ldu,*), vt(ldvt,*), work(*)
        end subroutine

        pure subroutine DGETRF(m, n, a, lda, ipiv, info)
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: m, n, lda
            real(real64), intent(inout) :: a(lda,*)
            integer(int32), intent(out) :: ipiv(*), info
        end subroutine

        pure subroutine ZGETRF(m, n, a, lda, ipiv, info)
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: m, n, lda
            complex(real64), intent(inout) :: a(lda,*)
            integer(int32), intent(out) :: ipiv(*), info
        end subroutine

        pure subroutine DSYEV(jobz, uplo, n, a, lda, w, work, lwork, info)
            use iso_fortran_env, only : int32, real64
            character, intent(in) :: jobz, uplo
            integer(int32), intent(in) :: n, lda, lwork
            real(real64), intent(inout) :: a(lda,*)
            real(real64), intent(out) :: w(*), work(*)
            integer(int32), intent(out) :: info
        end subroutine

        pure subroutine DGEEV(jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, &
            work, lwork, info)
            use iso_fortran_env, only : int32, real64
            character, intent(in) :: jobvl, jobvr
            integer(int32), intent(in) :: n, lda, ldvl, ldvr, lwork
            real(real64), intent(inout) :: a(lda,*)
            real(real64), intent(out) :: wr(*), wi(*), vl(ldvl,*), vr(ldvr,*), &
                work(*)
            integer(int32), intent(out) :: info
        end subroutine

        pure subroutine DGGEV(jobvl, jobvr, n, a, lda, b, ldb, alphar, alphai, &
            beta, vl, ldvl, vr, ldvr, work, lwork, info)
            use iso_fortran_env, only : int32, real64
            character, intent(in) :: jobvl, jobvr
            integer(int32), intent(in) :: n, lda, ldb, ldvl, ldvr, lwork
            real(real64), intent(inout) :: a(lda,*), b(ldb,*)
            real(real64), intent(out) :: alphar(*), alphai(*), beta(*), &
                vl(ldvl,*), vr(ldvr,*), work(*)
            integer(int32), intent(out) :: info
        end subroutine

        pure subroutine ZGEEV(jobvl, jobvr, n, a, lda, w, vl, ldvl, vr, ldvr, work, &
            lwork, rwork, info)
            use iso_fortran_env, only : int32, real64
            character, intent(in) :: jobvl, jobvr
            integer(int32), intent(in) :: n, lda, ldvl, ldvr, lwork
            complex(real64), intent(inout) :: a(lda,*)
            complex(real64), intent(out) :: w(*), vl(ldvl,*), vr(ldvr,*), &
                work(*)
            real(real64), intent(out) :: rwork(*)
            integer(int32), intent(out) :: info
        end subroutine

        pure subroutine DLASET(uplo, m, n, alpha, beta, a, lda)
            use iso_fortran_env, only : int32, real64
            character, intent(in) :: uplo
            integer(int32), intent(in) :: m, n, lda
            real(real64), intent(in) :: alpha, beta
            real(real64), intent(out) :: a(lda,*)
        end subroutine

        pure subroutine DGEQRF(m, n, a, lda, tau, work, lwork, info)
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: m, n, lda, lwork
            real(real64), intent(inout) :: a(lda,*)
            real(real64), intent(out) :: tau(*), work(*)
            integer(int32), intent(out) :: info
        end subroutine

        pure subroutine ZGEQRF(m, n, a, lda, tau, work, lwork, info)
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: m, n, lda, lwork
            complex(real64), intent(inout) :: a(lda,*)
            complex(real64), intent(out) :: tau(*), work(*)
            integer(int32), intent(out) :: info
        end subroutine

        pure subroutine DGEQP3(m, n, a, lda, jpvt, tau, work, lwork, info)
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: m, n, lda, lwork
            real(real64), intent(inout) :: a(lda,*)
            integer(int32), intent(inout) :: jpvt(*)
            real(real64), intent(out) :: tau(*), work(*)
            integer(int32), intent(out) :: info
        end subroutine

        pure subroutine ZGEQP3(m, n, a, lda, jpvt, tau, work, lwork, rwork, info)
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: m, n, lda, lwork
            complex(real64), intent(inout) :: a(lda,*)
            integer(int32), intent(inout) :: jpvt(*)
            complex(real64), intent(out) :: tau(*), work(*)
            real(real64), intent(out) :: rwork(*)
            integer(int32), intent(out) :: info
        end subroutine

        pure subroutine DORGQR(m, n, k, a, lda, tau, work, lwork, info)
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: m, n, k, lda, lwork
            real(real64), intent(inout) :: a(lda,*)
            real(real64), intent(in) :: tau(*)
            real(real64), intent(out) :: work(*)
            integer(int32), intent(out) :: info
        end subroutine

        pure subroutine ZUNGQR(m, n, k, a, lda, tau, work, lwork, info)
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: m, n, k, lda, lwork
            complex(real64), intent(inout) :: a(lda,*)
            complex(real64), intent(in) :: tau(*)
            complex(real64), intent(out) :: work(*)
            integer(int32), intent(out) :: info
        end subroutine

        pure subroutine DORMQR(side, trans, m, n, k, a, lda, tau, c, ldc, work, &
            lwork, info)
            use iso_fortran_env, only : int32, real64
            character, intent(in) :: side, trans
            integer(int32), intent(in) :: m, n, k, lda, ldc, lwork
            real(real64), intent(in) :: a(lda,*), tau(*)
            real(real64), intent(inout) :: c(ldc,*)
            real(real64), intent(out) :: work(*)
            integer(int32), intent(out) :: info
        end subroutine

        pure subroutine ZUNMQR(side, trans, m, n, k, a, lda, tau, c, ldc, work, &
            lwork, info)
            use iso_fortran_env, only : int32, real64
            character, intent(in) :: side, trans
            integer(int32), intent(in) :: m, n, k, lda, ldc, lwork
            complex(real64), intent(in) :: a(lda,*), tau(*)
            complex(real64), intent(inout) :: c(ldc,*)
            complex(real64), intent(out) :: work(*)
            integer(int32), intent(out) :: info
        end subroutine

        pure subroutine DPOTRF(uplo, n, a, lda, info)
            use iso_fortran_env, only : int32, real64
            character, intent(in) :: uplo
            integer(int32), intent(in) :: n, lda
            real(real64), intent(inout) :: a(lda,*)
            integer(int32), intent(out) :: info
        end subroutine

        pure subroutine ZPOTRF(uplo, n, a, lda, info)
            use iso_fortran_env, only : int32, real64
            character, intent(in) :: uplo
            integer(int32), intent(in) :: n, lda
            complex(real64), intent(inout) :: a(lda,*)
            integer(int32), intent(out) :: info
        end subroutine

        pure subroutine DTZRZF(m, n, a, lda, tau, work, lwork, info)
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: m, n, lda, lwork
            real(real64), intent(inout) :: a(lda,*)
            real(real64), intent(out) :: tau(*), work(*)
            integer(int32), intent(out) :: info
        end subroutine

        pure subroutine ZTZRZF(m, n, a, lda, tau, work, lwork, info)
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: m, n, lda, lwork
            complex(real64), intent(inout) :: a(lda,*)
            complex(real64), intent(out) :: tau(*), work(*)
            integer(int32), intent(out) :: info
        end subroutine

        pure subroutine DORMRZ(side, trans, m, n, k, l, a, lda, tau, c, ldc, work, &
            lwork, info)
            use iso_fortran_env, only : int32, real64
            character, intent(in) :: side, trans
            integer(int32), intent(in) :: m, n, k, l, lda, ldc, lwork
            real(real64), intent(in) :: a(lda,*), tau(*)
            real(real64), intent(inout) :: c(ldc,*)
            real(real64), intent(out) :: work(*)
            integer(int32), intent(out) :: info
        end subroutine

        pure subroutine zunmrz(side, trans, m, n, k, l, a, lda, tau, c, ldc, work, &
            lwork, info)
            use iso_fortran_env, only : int32, real64
            character, intent(in) :: side, trans
            integer(int32), intent(in) :: m, n, k, l, lda, ldc, lwork
            complex(real64), intent(in) :: a(lda,*), tau(*)
            complex(real64), intent(inout) :: c(ldc,*)
            complex(real64), intent(out) :: work(*)
            integer(int32), intent(out) :: info
        end subroutine

        pure subroutine DGELQF(m, n, a, lda, tau, work, lwork, info)
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: m, n, lda, lwork
            real(real64), intent(inout) :: a(lda,*)
            real(real64), intent(out) :: tau(*), work(*)
            integer(int32), intent(out) :: info
        end subroutine

        pure subroutine ZGELQF(m, n, a, lda, tau, work, lwork, info)
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: m, n, lda, lwork
            complex(real64), intent(inout) :: a(lda,*)
            complex(real64), intent(out) :: tau(*), work(*)
            integer(int32), intent(out) :: info
        end subroutine

        pure subroutine DORGLQ(m, n, k, a, lda, tau, work, lwork, info)
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: m, n, k, lda, lwork
            real(real64), intent(inout) :: a(lda,*)
            real(real64), intent(in) :: tau(*)
            real(real64), intent(out) :: work(*)
            integer(int32), intent(out) :: info
        end subroutine

        pure subroutine ZUNGLQ(m, n, k, a, lda, tau, work, lwork, info)
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: m, n, k, lda, lwork
            complex(real64), intent(inout) :: a(lda,*)
            complex(real64), intent(in) :: tau(*)
            complex(real64), intent(out) :: work(*)
            integer(int32), intent(out) :: info
        end subroutine

        pure subroutine DORMLQ(side, trans, m, n, k, a, lda, tau, c, ldc, work, &
            lwork, info)
            use iso_fortran_env, only : int32, real64
            character, intent(in) :: side, trans
            integer(int32), intent(in) :: m, n, k, lda, ldc, lwork
            real(real64), intent(in) :: a(lda,*), tau(*)
            real(real64), intent(inout) :: c(ldc,*)
            real(real64), intent(out) :: work(*)
            integer(int32), intent(out) :: info
        end subroutine

        pure subroutine ZUNMLQ(side, trans, m, n, k, a, lda, tau, c, ldc, work, &
            lwork, info)
            use iso_fortran_env, only : int32, real64
            character, intent(in) :: side, trans
            integer(int32), intent(in) :: m, n, k, lda, ldc, lwork
            complex(real64), intent(in) :: a(lda,*), tau(*)
            complex(real64), intent(inout) :: c(ldc,*)
            complex(real64), intent(out) :: work(*)
            integer(int32), intent(out) :: info
        end subroutine

        pure subroutine DGETRS(trans, n, nrhs, a, lda, ipiv, b, ldb, info)
            use iso_fortran_env, only : int32, real64
            character, intent(in) :: trans
            integer(int32), intent(in) :: n, nrhs, lda, ldb
            real(real64), intent(in) :: a(lda,*)
            integer(int32), intent(in) :: ipiv(*)
            real(real64), intent(inout) :: b(ldb,*)
            integer(int32), intent(out) :: info
        end subroutine

        pure subroutine ZGETRS(trans, n, nrhs, a, lda, ipiv, b, ldb, info)
            use iso_fortran_env, only : int32, real64
            character, intent(in) :: trans
            integer(int32), intent(in) :: n, nrhs, lda, ldb
            complex(real64), intent(in) :: a(lda,*)
            integer(int32), intent(in) :: ipiv(*)
            complex(real64), intent(inout) :: b(ldb,*)
            integer(int32), intent(out) :: info
        end subroutine

        pure subroutine DLAIC1(job, j, x, sest, w, gamma, sestpr, s, c)
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: job, j
            real(real64), intent(in) :: x(j), w(j), sest, gamma
            real(real64), intent(out) :: sestpr, s, c
        end subroutine

        pure subroutine ZLAIC1(job, j, x, sest, w, gamma, sestpr, s, c)
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: job, j
            complex(real64), intent(in) :: x(j), w(j), gamma
            real(real64), intent(in) :: sest
            real(real64), intent(out) :: sestpr
            complex(real64), intent(out) :: s, c
        end subroutine

        pure subroutine DPOTRS(uplo, n, nrhs, a, lda, b, ldb, info)
            use iso_fortran_env, only : int32, real64
            character, intent(in) :: uplo
            integer(int32), intent(in) :: n, nrhs, lda, ldb
            real(real64), intent(in) :: a(lda,*)
            real(real64), intent(inout) :: b(ldb,*)
            integer(int32), intent(out) :: info
        end subroutine

        pure subroutine ZPOTRS(uplo, n, nrhs, a, lda, b, ldb, info)
            use iso_fortran_env, only : int32, real64
            character, intent(in) :: uplo
            integer(int32), intent(in) :: n, nrhs, lda, ldb
            complex(real64), intent(in) :: a(lda,*)
            complex(real64), intent(inout) :: b(ldb,*)
            integer(int32), intent(out) :: info
        end subroutine

        pure subroutine DGETRI(n, a, lda, ipiv, work, lwork, info)
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: n, lda, lwork, ipiv(*)
            real(real64), intent(inout) :: a(lda,*)
            real(real64), intent(out) :: work(*)
            integer(int32), intent(out) :: info
        end subroutine

        pure subroutine ZGETRI(n, a, lda, ipiv, work, lwork, info)
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: n, lda, lwork, ipiv(*)
            complex(real64), intent(inout) :: a(lda,*)
            complex(real64), intent(out) :: work(*)
            integer(int32), intent(out) :: info
        end subroutine

        pure subroutine DGELS(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info)
            use iso_fortran_env, only : int32, real64
            character, intent(in) :: trans
            integer(int32), intent(in) :: m, n, nrhs, lda, ldb, lwork
            real(real64), intent(inout) :: a(lda,*), b(ldb,*)
            real(real64), intent(out) :: work(*)
            integer(int32), intent(out) :: info
        end subroutine

        pure subroutine ZGELS(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info)
            use iso_fortran_env, only : int32, real64
            character, intent(in) :: trans
            integer(int32), intent(in) :: m, n, nrhs, lda, ldb, lwork
            complex(real64), intent(inout) :: a(lda,*), b(ldb,*)
            complex(real64), intent(out) :: work(*)
            integer(int32), intent(out) :: info
        end subroutine

        pure subroutine DGELSY(m, n, nrhs, a, lda, b, ldb, jpvt, rcond, rank, work, &
            lwork, info)
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: m, n, nrhs, lda, ldb, lwork
            real(real64), intent(inout) :: a(lda,*), b(ldb,*)
            integer(int32), intent(inout) :: jpvt(*)
            real(real64), intent(in) :: rcond
            integer(int32), intent(out) :: rank, info
            real(real64), intent(out) :: work(*)
        end subroutine

        pure subroutine ZGELSY(m, n, nrhs, a, lda, b, ldb, jpvt, rcond, rank, work, &
            lwork, rwork, info)
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: m, n, nrhs, lda, ldb, lwork
            complex(real64), intent(inout) :: a(lda,*), b(ldb,*)
            integer(int32), intent(inout) :: jpvt(*)
            real(real64), intent(in) :: rcond
            integer(int32), intent(out) :: rank, info
            complex(real64), intent(out) :: work(*)
            real(real64), intent(out) :: rwork(*)
        end subroutine

        pure subroutine DGELSS(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, &
            lwork, info)
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: m, n, nrhs, lda, ldb, lwork
            real(real64), intent(inout) :: a(lda,*), b(ldb,*)
            real(real64), intent(out) :: s(*)
            real(real64), intent(in) :: rcond
            integer(int32), intent(out) :: rank, info
            real(real64), intent(out) :: work(*)
        end subroutine

        pure subroutine ZGELSS(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, &
            lwork, rwork, info)
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: m, n, nrhs, lda, ldb, lwork
            complex(real64), intent(inout) :: a(lda,*), b(ldb,*)
            real(real64), intent(out) :: s(*)
            real(real64), intent(in) :: rcond
            integer(int32), intent(out) :: rank, info
            complex(real64), intent(out) :: work(*)
            real(real64), intent(out) :: rwork(*)
        end subroutine

        pure subroutine DLASRT(id, n, d, info)
            use iso_fortran_env, only : int32, real64
            character, intent(in) :: id
            integer(int32), intent(in) :: n
            real(real64), intent(inout) :: d(*)
            integer(int32), intent(out) :: info
        end subroutine

        pure subroutine DGESV(n, nrhs, a, lda, ipiv, b, ldb, info)
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: n, nrhs, lda, ldb
            integer(int32), intent(out) :: ipiv(*)
            real(real64), intent(inout) :: a(lda,*), b(ldb,*)
            integer(int32), intent(out) :: info
        end subroutine

        pure subroutine DGGEV3(jobvl, jobvr, n, a, lda, b, ldb, alphar, alphai, &
            beta, vl, ldvl, vr, ldvr, work, lwork, info)
            use iso_fortran_env, only : int32, real64
            character, intent(in) :: jobvl, jobvr
            integer(int32), intent(in) :: n, lda, ldb, ldvl, ldvr, lwork
            real(real64), intent(inout) :: a(lda,*), b(ldb,*)
            real(real64), intent(out) :: alphar(*), alphai(*), beta(*), &
                vl(ldvl,*), vr(ldvr,*), work(*)
            integer(int32), intent(out) :: info
        end subroutine
    end interface
end module