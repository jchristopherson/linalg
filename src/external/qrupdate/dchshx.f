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
      subroutine dchshx(n,R,ldr,i,j,w)
c purpose:      given an upper triangular matrix R that is a Cholesky
c               factor of a symmetric positive definite matrix A, i.e.
c               A = R'*R, this subroutine updates R -> R1 so that
c               R1'*R1 = A(p,p), where p is the permutation
c               [1:i-1,shift(i:j,-1),j+1:n] if i < j  or
c               [1:j-1,shift(j:i,+1),i+1:n] if j < i.
c               (real version)
c arguments:
c n (in)        the order of matrix R
c R (io)        on entry, the upper triangular matrix R
c               on exit, the updated matrix R1
c ldr (in)      leading dimension of R. ldr >= n.
c i (in)        the first index determining the range (see above).
c j (in)        the second index determining the range (see above).
c w (o)         a workspace vector of size 2*n.
c
      integer n,ldr,i,j
      double precision R(ldr,*),w(*)
      external xerbla,dcopy,dqrtv1,dqrqh,dqhqr
      integer info,l
c quick return if possible.
      if (n == 0 .or. n == 1) return
      info = 0
c check arguments.
      if (n < 0) then
        info = 1
      else if (i < 1 .or. i > n) then
        info = 4
      else if (j < 1 .or. j > n) then
        info = 5
      end if
      if (info /= 0) then
        call xerbla('DCHSHX',info)
        return
      end if

      if (i < j) then
c shift columns
        call dcopy(n,R(1,i),1,w,1)
        do l = i,j-1
          call dcopy(n,R(1,l+1),1,R(1,l),1)
        end do
        call dcopy(n,w,1,R(1,j),1)
c retriangularize
        call dqhqr(n+1-i,n+1-i,R(i,i),ldr,w(n+1),w)
      else if (j < i) then
c shift columns
        call dcopy(n,R(1,i),1,w,1)
        do l = i,j+1,-1
          call dcopy(n,R(1,l-1),1,R(1,l),1)
        end do
        call dcopy(n,w,1,R(1,j),1)
c eliminate the introduced spike.
        call dqrtv1(n+1-j,R(j,j),w(n+1))
c apply rotations to R
        call dqrqh(n+1-j,n-j,R(j,j+1),ldr,w(n+1),R(j+1,j))
c zero spike.
        do l = j+1,n
          R(l,j) = 0d0
        end do
      end if
      end subroutine
