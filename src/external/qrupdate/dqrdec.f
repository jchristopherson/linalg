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
      subroutine dqrdec(m,n,k,Q,ldq,R,ldr,j,w)
c purpose:      updates a QR factorization after deleting
c               a column.
c               i.e., given an m-by-k orthogonal matrix Q, an k-by-n
c               upper trapezoidal matrix R and index j in the range
c               1:n+1, this subroutine updates the matrix Q -> Q1 and
c               R -> R1 so that Q1 remains orthogonal, R1 is upper
c               trapezoidal, and Q1*R1 = [A(:,1:j-1) A(:,j+1:n)],
c               where A = Q*R.
c               (real version)
c arguments:
c m (in)        number of rows of the matrix Q.
c n (in)        number of columns of the matrix R.
c k (in)        number of columns of Q, and rows of R. Must be
c               either k = m (full Q) or k = n < m (economical form,
c               basis dimension will decrease).
c Q (io)        on entry, the unitary m-by-k matrix Q.
c               on exit, the updated matrix Q1.
c ldq (in)      leading dimension of Q. ldq >= m.
c R (io)        on entry, the original matrix R.
c               on exit, the updated matrix R1.
c ldr (in)      leading dimension of R. ldr >= k.
c j (in)        the position of the deleted column in R.
c               1 <= j <= n.
c w (o)         a workspace vector of size k-j.
c
      integer m,n,k,ldq,ldr,j
      double precision Q(ldq,*),R(ldr,*),w(*)
      external xerbla,dcopy,dqhqr,dqrot
      integer info,i
c quick return if possible.
      if (m == 0 .or. n == 0 .or. j == n) return
c check arguments.
      info = 0
      if (m < 0) then
        info = 1
      else if (n < 0) then
        info = 2
      else if (k /= m .and. (k /= n .or. n >= m)) then
        info = 3
      else if (ldq < m) then
        info = 5
      else if (ldr < k) then
        info = 7
      else if (j < 1 .or. j > n+1) then
        info = 8
      end if
      if (info /= 0) then
        call xerbla('DQRDEC',info)
        return
      end if

c delete the j-th column.
      do i = j,n-1
        call dcopy(k,R(1,i+1),1,R(1,i),1)
      end do
c retriangularize.
      if (j < k) then
        call dqhqr(k+1-j,n-j,R(j,j),ldr,w,R(1,n))
c apply rotations to Q.
        call dqrot('F',m,min(k,n)+1-j,Q(1,j),ldq,w,R(1,n))
      end if
      end subroutine
