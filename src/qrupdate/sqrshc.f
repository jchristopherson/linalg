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
      subroutine sqrshc(m,n,k,Q,ldq,R,ldr,i,j,w)
c purpose:      updates a QR factorization after circular shift of
c               columns.
c               i.e., given an m-by-k orthogonal matrix Q, an k-by-n
c               upper trapezoidal matrix R and index j in the range
c               1:n+1, this subroutine updates the matrix Q -> Q1 and
c               R -> R1 so that Q1 is again orthogonal, R1 upper
c               trapezoidal, and
c               Q1*R1 = A(:,p), where A = Q*R and p is the permutation
c               [1:i-1,shift(i:j,-1),j+1:n] if i < j  or
c               [1:j-1,shift(j:i,+1),i+1:n] if j < i.
c               (real version)
c arguments:
c m (in)        number of rows of the matrix Q.
c n (in)        number of columns of the matrix R.
c k (in)        number of columns of Q1, and rows of R1. Must be
c               either k = m (full Q) or k = n <= m (economical form).
c Q (io)        on entry, the unitary m-by-k matrix Q.
c               on exit, the updated matrix Q1.
c ldq (in)      leading dimension of Q. ldq >= m.
c R (io)        on entry, the original matrix R.
c               on exit, the updated matrix R1.
c ldr (in)      leading dimension of R. ldr >= k.
c i (in)        the first index determining the range (see above)
c j (in)        the second index determining the range (see above)
c w (o)         a workspace vector of size 2*k.
c
      integer m,n,k,ldq,ldr,i,j
      real Q(ldq,*),R(ldr,*),w(*)
      external xerbla,scopy,sqrtv1,sqrqh,sqhqr
      integer info,jj,kk,l
c quick return if possible.
      if (m == 0 .or. n == 1) return
      info = 0
c check arguments.
      if (m < 0) then
        info = 1
      else if (n < 0) then
        info = 2
      else if (k /= m .and. (k /= n .or. n > m)) then
        info = 3
      else if (i < 1 .or. i > n) then
        info = 6
      else if (j < 1 .or. j > n) then
        info = 7
      end if
      if (info /= 0) then
        call xerbla('SQRSHC',info)
        return
      end if

      if (i < j) then
c shift columns
        call scopy(k,R(1,i),1,w,1)
        do l = i,j-1
          call scopy(k,R(1,l+1),1,R(1,l),1)
        end do
        call scopy(k,w,1,R(1,j),1)
c retriangularize
        if (i < k) then
          kk = min(k,j)
          call sqhqr(kk+1-i,n+1-i,R(i,i),ldr,w(k+1),w)
c apply rotations to Q.
          call sqrot('F',m,kk+1-i,Q(1,i),ldq,w(k+1),w)
        end if
      else if (j < i) then
c shift columns
        call scopy(k,R(1,i),1,w,1)
        do l = i,j+1,-1
          call scopy(k,R(1,l-1),1,R(1,l),1)
        end do
        call scopy(k,w,1,R(1,j),1)
c retriangularize
        if (j < k) then
          jj = min(j+1,n)
          kk = min(k,i)
c eliminate the introduced spike.
          call sqrtv1(kk+1-j,R(j,j),w(k+1))
c apply rotations to R
          call sqrqh(kk+1-j,n-j,R(j,jj),ldr,w(k+1),R(j+1,j))
c apply rotations to Q
          call sqrot('B',m,kk+1-j,Q(1,j),ldq,w(k+1),R(j+1,j))
c zero spike.
          do l = j+1,kk
            R(l,j) = 0e0
          end do
        end if
      end if
      end subroutine
