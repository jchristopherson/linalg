c Copyright (C) 2008, 2009  VZLU Prague, a.s., Czech Republic
c
c Author: Jaroslav Hajek <highegg@gmail.com>
c
c This file is part of qrupdate.
c
c qrupdate is free software; you can redistribute it and/or modify
c it under the terms of the GNU General Public License as published by
c the Free Software Foundation; either version 3 of the License, or
c (at your option) any later version.
c
c This program is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c GNU General Public License for more details.
c
c You should have received a copy of the GNU General Public License
c along with this software; see the file COPYING.  If not, see
c <http://www.gnu.org/licenses/>.
c
      subroutine schinx(n,R,ldr,j,u,w,info)
c purpose:      given an upper triangular matrix R that is a Cholesky
c               factor of a symmetric positive definite matrix A, i.e.
c               A = R'*R, this subroutine updates R -> R1 so that
c               R1'*R1 = A1, A1(jj,jj) = A, A(j,:) = u', A(:,j) = u,
c               jj = [1:j-1,j+1:n+1].
c               (real version)
c arguments:
c n (in)        the order of matrix R.
c R (io)        on entry, the original upper trapezoidal matrix R.
c               on exit, the updated matrix R1.
c ldr (in)      leading dimension of R. ldr >= n+1.
c j (in)        the position of the inserted row/column
c u (io)        on entry, the inserted row/column.
c               on exit, u is destroyed.
c w (out)       workspace vector of size n+1.
c info (out)    on exit, error code:
c                info = 1: update violates positive-definiteness.
c                info = 2: R is singular.
c
      integer n,j,ldr,info
      real R(ldr,*),u(*),w(*)
      external xerbla,scopy,snrm2,strsv,sqrtv1,sqrqh
      real snrm2,t,rho
      integer i

c check arguments
      info = 0
      if (n < 0) then
        info = -1
      else if (j < 1 .or. j > n+1) then
        info = -4
      end if
      if (info /= 0) then
        call xerbla('SCHINX',-info)
        return
      end if

c shift vector.
      t = u(j)
      do i = j,n
        u(i) = u(i+1)
      end do

c check for singularity of R.
      do i = 1,n
        if (R(i,i) == 0e0) goto 20
      end do
c form R' \ u
      call strsv('U','T','N',n,R,ldr,u,1)
      rho = snrm2(n,u,1)
c check positive definiteness.
      rho = t - rho**2
      if (rho <= 0e0) goto 10
c shift columns
      do i = n,j,-1
        call scopy(i,R(1,i),1,R(1,i+1),1)
        R(i+1,i+1) = 0e0
      end do
      call scopy(n,u,1,R(1,j),1)
      R(n+1,j) = sqrt(rho)
c retriangularize
      if (j < n+1) then
c eliminate the introduced spike.
        call sqrtv1(n+2-j,R(j,j),w)
c apply rotations to R
        call sqrqh(n+2-j,n+1-j,R(j,j+1),ldr,w,R(j+1,j))
c zero spike.
        do i = j+1,n+1
          R(i,j) = 0e0
        end do
      end if
c normal return.
      return
c error returns.
 10   info = 1
      return
 20   info = 2
      return
      end subroutine
