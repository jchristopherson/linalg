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
      subroutine dqrinr(m,n,Q,ldq,R,ldr,j,x,w)
c purpose:      updates a QR factorization after inserting a new
c               row.
c               i.e., given an m-by-m unitary matrix Q, an m-by-n
c               upper trapezoidal matrix R and index j in the range
c               1:m+1, this subroutine updates Q -> Q1  and R -> R1
c               so that Q1 is again unitary, R1 upper trapezoidal,
c               and Q1*R1 = [A(1:j-1,:); x; A(j:m,:)], where A = Q*R.
c               (real version)
c arguments:
c m (in)        number of rows of the matrix Q.
c n (in)        number of columns of the matrix R.
c Q (io)        on entry, the unitary matrix Q.
c               on exit, the updated matrix Q1.
c ldq (in)      leading dimension of Q. ldq >= m+1.
c R (io)        on entry, the original matrix R.
c               on exit, the updated matrix R1.
c ldr (in)      leading dimension of R. ldr >= m+1.
c j (in)        the position of the new row in R1
c x (io)        on entry, the row being added
c               on exit, x is destroyed.
c w (out)       a workspace vector of size min(m,n).
c
      integer m,n,j,ldq,ldr
      double precision Q(ldq,*),R(ldr,*),x(*),w(*)
      external xerbla,dcopy,dqhqr,dqrot
      integer info,i,k
c check arguments
      info = 0
      if (n < 0) then
        info = 2
      else if (j < 1 .or. j > m+1) then
        info = 7
      end if
      if (info /= 0) then
        call xerbla('DQRINR',info)
        return
      end if
c permute the columns of Q1 and rows of R1 so that c the new row ends
c up being the topmost row of R1.
      do i = m,1,-1
        if (j > 1) then
          call dcopy(j-1,Q(1,i),1,Q(1,i+1),1)
        end if
        Q(j,i+1) = 0d0
        if (j <= m) then
          call dcopy(m+1-j,Q(j,i),1,Q(j+1,i+1),1)
        end if
      end do
c set up the 1st column
      do i = 1,j-1
        Q(i,1) = 0d0
      end do
      Q(j,1) = 1d0
      do i = j+1,m+1
        Q(i,1) = 0d0
      end do
c set up the new matrix R1
      do k = 1,n
        if (k < m) R(m+1,k) = 0d0
        do i = min(m,k),1,-1
          R(i+1,k) = R(i,k)
        end do
        R(1,k) = x(k)
      end do
c retriangularize R
      call dqhqr(m+1,n,R,ldr,w,x)
c apply rotations to Q
      call dqrot('F',m+1,min(m,n)+1,Q,ldq,w,x)
      end subroutine
