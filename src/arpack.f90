module arpack
    !! An interface to the ARPACK library available at
    !! https://github.com/opencollab/arpack-ng.
    implicit none

    interface
        subroutine dsaupd(ido, bmat, n, which, nev, tol, resid, ncv, v, ldv, &
            iparam, ipntr, workd, workl, lworkl, info)
            !! Implements the implicitly restarted Arnoldi iteration (Lanczos
            !! for symmetric problems) by means of a reverse communication
            !! interface.
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(inout) :: ido
                !! The reverse communication flag.
            character(len = 1), intent(in) :: bmat
                !! 'I' for a standard problem, or 'G' for a generalized problem.
            integer(int32), intent(in) :: n
                !! The dimension of the problem.
            character(len = 2), intent(in) :: which
                !! Identifies which eigenvalues to compute ('LA', 'SA', 'LM',
                !! 'SM', or 'BE').
            integer(int32), intent(inout) :: nev
                !! The number of eigenvalues to compute.
            real(real64), intent(inout) :: tol
                !! The convergence tolerance.  If zero, machine precision is
                !! used.
            real(real64), intent(inout) :: resid(*)
                !! An N-element residual vector.
            integer(int32), intent(in) :: ncv
                !! The number of Lanczos basis vectors to use.
            real(real64), intent(inout) :: v(ldv,*)
                !! An LDV-by-NCV matrix containing the Lanczos basis vectors.
            integer(int32), intent(in) :: ldv
                !! The leading dimension of V.
            integer(int32), intent(inout) :: iparam(*)
                !! An 11-element array containing the iteration parameters.
            integer(int32), intent(inout) :: ipntr(*)
                !! An 11-element array containing pointers into WORKD and WORKL.
            real(real64), intent(inout) :: workd(*)
                !! A 3*N-element reverse communication workspace array.
            real(real64), intent(inout) :: workl(*)
                !! An LWORKL-element workspace array.
            integer(int32), intent(in) :: lworkl
                !! The length of WORKL.  Must be at least NCV**2 + 8 * NCV.
            integer(int32), intent(inout) :: info
                !! The status flag.
        end subroutine

        subroutine dseupd(rvec, howmny, sel, d, z, ldz, sigma, bmat, n, which, &
            nev, tol, resid, ncv, v, ldv, iparam, ipntr, workd, workl, &
            lworkl, info)
            !! Extracts the eigenvalues and, optionally, the eigenvectors, from
            !! the results of DSAUPD.
            use iso_fortran_env, only : int32, real64
            logical, intent(in) :: rvec
                !! Set to true to compute the eigenvectors; else, false.
            character(len = 1), intent(in) :: howmny
                !! 'A' to compute all NEV eigenvectors, or 'S' to compute those
                !! flagged in SEL.
            logical, intent(inout) :: sel(*)
                !! An NCV-element array used when HOWMNY is 'S'.  The array is
                !! also used as workspace when HOWMNY is 'A'.
            real(real64), intent(out) :: d(*)
                !! An NEV-element array containing the eigenvalues.
            real(real64), intent(out) :: z(ldz,*)
                !! An LDZ-by-NEV matrix containing the eigenvectors.
            integer(int32), intent(in) :: ldz
                !! The leading dimension of Z.
            real(real64), intent(in) :: sigma
                !! The shift used in shift-invert modes.
            character(len = 1), intent(in) :: bmat
                !! 'I' for a standard problem, or 'G' for a generalized problem.
            integer(int32), intent(in) :: n
                !! The dimension of the problem.
            character(len = 2), intent(in) :: which
                !! Identifies which eigenvalues were computed.
            integer(int32), intent(inout) :: nev
                !! The number of eigenvalues requested.
            real(real64), intent(inout) :: tol
                !! The convergence tolerance.
            real(real64), intent(inout) :: resid(*)
                !! The N-element residual vector from DSAUPD.
            integer(int32), intent(in) :: ncv
                !! The number of Lanczos basis vectors used.
            real(real64), intent(inout) :: v(ldv,*)
                !! The LDV-by-NCV Lanczos basis from DSAUPD.
            integer(int32), intent(in) :: ldv
                !! The leading dimension of V.
            integer(int32), intent(inout) :: iparam(*)
                !! The 11-element iteration parameter array from DSAUPD.
            integer(int32), intent(inout) :: ipntr(*)
                !! The 11-element pointer array from DSAUPD.
            real(real64), intent(inout) :: workd(*)
                !! The reverse communication workspace array from DSAUPD.
            real(real64), intent(inout) :: workl(*)
                !! The LWORKL-element workspace array from DSAUPD.
            integer(int32), intent(in) :: lworkl
                !! The length of WORKL.
            integer(int32), intent(inout) :: info
                !! The status flag.
        end subroutine

        subroutine dnaupd(ido, bmat, n, which, nev, tol, resid, ncv, v, ldv, &
            iparam, ipntr, workd, workl, lworkl, info)
            !! Implements the implicitly restarted Arnoldi iteration for
            !! non-symmetric problems by means of a reverse communication
            !! interface.
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(inout) :: ido
                !! The reverse communication flag.
            character(len = 1), intent(in) :: bmat
                !! 'I' for a standard problem, or 'G' for a generalized problem.
            integer(int32), intent(in) :: n
                !! The dimension of the problem.
            character(len = 2), intent(in) :: which
                !! Identifies which eigenvalues to compute ('LM', 'SM', 'LR',
                !! 'SR', 'LI', or 'SI').
            integer(int32), intent(inout) :: nev
                !! The number of eigenvalues to compute.
            real(real64), intent(inout) :: tol
                !! The convergence tolerance.  If zero, machine precision is
                !! used.
            real(real64), intent(inout) :: resid(*)
                !! An N-element residual vector.
            integer(int32), intent(in) :: ncv
                !! The number of Arnoldi basis vectors to use.
            real(real64), intent(inout) :: v(ldv,*)
                !! An LDV-by-NCV matrix containing the Arnoldi basis vectors.
            integer(int32), intent(in) :: ldv
                !! The leading dimension of V.
            integer(int32), intent(inout) :: iparam(*)
                !! An 11-element array containing the iteration parameters.
            integer(int32), intent(inout) :: ipntr(*)
                !! A 14-element array containing pointers into WORKD and WORKL.
            real(real64), intent(inout) :: workd(*)
                !! A 3*N-element reverse communication workspace array.
            real(real64), intent(inout) :: workl(*)
                !! An LWORKL-element workspace array.
            integer(int32), intent(in) :: lworkl
                !! The length of WORKL.  Must be at least 3 * NCV**2 + 6 * NCV.
            integer(int32), intent(inout) :: info
                !! The status flag.
        end subroutine

        subroutine dneupd(rvec, howmny, sel, dr, di, z, ldz, sigmar, sigmai, &
            workev, bmat, n, which, nev, tol, resid, ncv, v, ldv, iparam, &
            ipntr, workd, workl, lworkl, info)
            !! Extracts the eigenvalues and, optionally, the eigenvectors, from
            !! the results of DNAUPD.
            use iso_fortran_env, only : int32, real64
            logical, intent(in) :: rvec
                !! Set to true to compute the eigenvectors; else, false.
            character(len = 1), intent(in) :: howmny
                !! 'A' to compute all NEV eigenvectors, 'P' for the Schur
                !! basis, or 'S' to compute those flagged in SEL.
            logical, intent(inout) :: sel(*)
                !! An NCV-element array used when HOWMNY is 'S'.  The array is
                !! also used as workspace when HOWMNY is 'A'.
            real(real64), intent(out) :: dr(*)
                !! An (NEV+1)-element array containing the real components of
                !! the eigenvalues.
            real(real64), intent(out) :: di(*)
                !! An (NEV+1)-element array containing the imaginary components
                !! of the eigenvalues.
            real(real64), intent(out) :: z(ldz,*)
                !! An LDZ-by-(NEV+1) matrix containing the eigenvectors stored
                !! in the compact form used by LAPACK's DGEEV.
            integer(int32), intent(in) :: ldz
                !! The leading dimension of Z.
            real(real64), intent(in) :: sigmar
                !! The real component of the shift used in shift-invert modes.
            real(real64), intent(in) :: sigmai
                !! The imaginary component of the shift used in shift-invert
                !! modes.
            real(real64), intent(out) :: workev(*)
                !! A 3*NCV-element workspace array.
            character(len = 1), intent(in) :: bmat
                !! 'I' for a standard problem, or 'G' for a generalized problem.
            integer(int32), intent(in) :: n
                !! The dimension of the problem.
            character(len = 2), intent(in) :: which
                !! Identifies which eigenvalues were computed.
            integer(int32), intent(inout) :: nev
                !! The number of eigenvalues requested.
            real(real64), intent(inout) :: tol
                !! The convergence tolerance.
            real(real64), intent(inout) :: resid(*)
                !! The N-element residual vector from DNAUPD.
            integer(int32), intent(in) :: ncv
                !! The number of Arnoldi basis vectors used.
            real(real64), intent(inout) :: v(ldv,*)
                !! The LDV-by-NCV Arnoldi basis from DNAUPD.
            integer(int32), intent(in) :: ldv
                !! The leading dimension of V.
            integer(int32), intent(inout) :: iparam(*)
                !! The 11-element iteration parameter array from DNAUPD.
            integer(int32), intent(inout) :: ipntr(*)
                !! The 14-element pointer array from DNAUPD.
            real(real64), intent(inout) :: workd(*)
                !! The reverse communication workspace array from DNAUPD.
            real(real64), intent(inout) :: workl(*)
                !! The LWORKL-element workspace array from DNAUPD.
            integer(int32), intent(in) :: lworkl
                !! The length of WORKL.
            integer(int32), intent(inout) :: info
                !! The status flag.
        end subroutine
    end interface
end module
