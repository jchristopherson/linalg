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

      subroutine zqrder(m,n,Q,ldq,R,ldr,j,w,rw)
c purpose:      updates a QR factorization after deleting a row.
c               i.e., given an m-by-m unitary matrix Q, an m-by-n
c               upper trapezoidal matrix R and index j in the range
c               1:m, this subroutine updates Q ->Q1 and an R -> R1
c               so that Q1 is again unitary, R1 upper trapezoidal,
c               and Q1*R1 = [A(1:j-1,:); A(j+1:m,:)], where A = Q*R.
c               (complex version)
c
c arguments:
c m (in)        number of rows of the matrix Q.
c n (in)        number of columns of the matrix R.
c Q (io)        on entry, the unitary matrix Q.
c               on exit, the updated matrix Q1.
c ldq (in)      leading dimension of Q. ldq >= m.
c R (io)        on entry, the original matrix R.
c               on exit, the updated matrix R1.
c ldr (in)      leading dimension of R. ldr >= m.
c j (in)        the position of the deleted row.
c w (out)       a workspace vector of size m.
c rw (out)      a real workspace vector of size m.
c
      integer m,n,j,ldq,ldr
      double complex Q(ldq,*),R(ldr,*),w(*)
      double precision rw(*)
      external xerbla,zcopy,zqrtv1,zqrot,zqrqh
      integer info,i,k
c quick return if possible
      if (m == 1) return
c check arguments
      info = 0
      if (m < 1) then
        info = 1
      else if (j < 1 .or. j > m) then
        info = 7
      end if
      if (info /= 0) then
        call xerbla('ZQRDER',info)
        return
      end if
c eliminate Q(j,2:m).
      do k = 1,m
        w(k) = conjg(Q(j,k))
      end do
      call zqrtv1(m,w,rw)
c apply rotations to Q.
      call zqrot('B',m,m,Q,ldq,rw,w(2))
c form Q1.
      do k = 1,m-1
        if (j > 1) call zcopy(j-1,Q(1,k+1),1,Q(1,k),1)
        if (j < m) call zcopy(m-j,Q(j+1,k+1),1,Q(j,k),1)
      end do
c apply rotations to R.
      call zqrqh(m,n,R,ldr,rw,w(2))
c form R1.
      do k = 1,n
        do i = 1,m-1
          R(i,k) = R(i+1,k)
        end do
      end do
      end subroutine
