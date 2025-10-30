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
      subroutine cqrinc(m,n,k,Q,ldq,R,ldr,j,x,rw)
c purpose:      updates a QR factorization after inserting a new
c               column.
c               i.e., given an m-by-k unitary matrix Q, an m-by-n upper
c               trapezoidal matrix R and index j in the range 1:n+1,
c               this subroutine updates the matrix Q -> Q1 and R -> R1
c               so that Q1 is again unitary, R1 upper trapezoidal, and
c               Q1*R1 = [A(:,1:j-1); x; A(:,j:n)], where A = Q*R.
c               (complex version)
c arguments:
c m (in)        number of rows of the matrix Q.
c n (in)        number of columns of the matrix R.
c k (in)        number of columns of Q, and rows of R. Must be
c               either k = m (full Q) or k = n <= m (economical form,
c               basis dimension will increase).
c Q (io)        on entry, the unitary m-by-k matrix Q.
c               on exit, the updated matrix Q1.
c ldq (in)      leading dimension of Q. ldq >= m.
c R (io)        on entry, the original matrix R.
c               on exit, the updated matrix R1.
c ldr (in)      leading dimension of R. ldr >= min(m,n+1).
c j (in)        the position of the new column in R1
c x (in)        the column being inserted
c rw (out)      a real workspace vector of size k.
c
      integer m,n,k,ldq,ldr,j
      complex Q(ldq,*),R(ldr,*),x(*)
      real rw(*)
      external cqrtv1,cqrqh,cqrot
      external xerbla,ccopy,cdotc,caxpy,csscal,scnrm2
      complex cdotc
      real scnrm2,rx
      integer info,i,k1
      logical full
c quick return if possible.
      if (m == 0) return
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
      else if (ldr < min(m,k+1)) then
        info = 7
      else if (j < 1 .or. j > n+1) then
        info = 8
      end if
      if (info /= 0) then
        call xerbla('CQRINC',info)
        return
      end if

      full = k == m
c insert empty column at j-th position
      do i = n,j,-1
        call ccopy(k,R(1,i),1,R(1,i+1),1)
      end do
c insert Q'*u into R. In the nonfull case, form also u-Q*Q'*u.
      if (full) then
        k1 = k
        do i = 1,k
          R(i,j) = cdotc(m,Q(1,i),1,x,1)
        end do
      else
        k1 = k + 1
c zero last row of R
        do i = 1,n+1
          R(k1,i) = 0e0
        end do
        call ccopy(m,x,1,Q(1,k1),1)
        do i = 1,k
          R(i,j) = cdotc(m,Q(1,i),1,Q(1,k1),1)
          call caxpy(m,-R(i,j),Q(1,i),1,Q(1,k1),1)
        end do
c get norm of the inserted column
        rx = scnrm2(m,Q(1,k1),1)
        R(k1,j) = rx
        if (rx == 0e0) then
c in the rare case when rx is exact zero, we still need to provide
c a valid orthogonal unit vector. The details are boring, so handle
c that elsewhere.
          call cgqvec(m,k,Q,ldq,Q(1,k1))
        else
c otherwise, just normalize the added column.
          call csscal(m,1e0/rx,Q(1,k1),1)
        end if
      end if
c maybe we're finished.
      if (j > k) return
c eliminate the spike.
      call cqrtv1(k1+1-j,R(j,j),rw)
c apply rotations to R(j:k,j:n).
      if (j <= n) call cqrqh(k1+1-j,n+1-j,R(j,j+1),ldr,rw,R(j+1,j))
c apply rotations to Q(:,j:k).
      call cqrot('B',m,k1+1-j,Q(1,j),ldq,rw,R(j+1,j))
c zero spike.
      do i = j+1,k1
        R(i,j) = 0e0
      end do
      end subroutine
