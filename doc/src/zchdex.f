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
      subroutine zchdex(n,R,ldr,j,rw)
c purpose:      given an upper triangular matrix R that is a Cholesky
c               factor of a hermitian positive definite matrix A, i.e.
c               A = R'*R, this subroutine updates R -> R1 so that
c               R1'*R1 = A(jj,jj), where jj = [1:j-1,j+1:n+1].
c               (complex version)
c arguments:
c n (in)        the order of matrix R.
c R (io)        on entry, the original upper trapezoidal matrix R.
c               on exit, the updated matrix R1.
c ldr (in)      leading dimension of R. ldr >= n.
c j (in)        the position of the deleted row/column.
c rw (out)      a real workspace vector of size n.
c
      integer n,ldr,j
      double complex R(ldr,*)
      double precision rw(*)
      integer info,i
      external xerbla,zcopy,zqhqr

c quick return if possible
      if (n == 1) return

c check arguments
      info = 0
      if (n < 0) then
        info = 1
      else if (j < 1 .or. j > n) then
        info = 4
      end if
      if (info /= 0) then
        call xerbla('ZCHDEX',info)
        return
      end if

c delete the j-th column.
      do i = j,n-1
        call zcopy(n,R(1,i+1),1,R(1,i),1)
      end do
c retriangularize.
      if (j < n) then
        call zqhqr(n+1-j,n-j,R(j,j),ldr,rw,R(1,n))
      end if
      end subroutine
