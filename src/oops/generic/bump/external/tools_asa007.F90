!----------------------------------------------------------------------
! Module: tools_asa007
! Purpose: inverse of symmetric positive definite matrix routines
! Source: https://people.sc.fsu.edu/~jburkardt/f_src/asa007/asa007.html
! Author: Michael Healy
! Original licensing: none
! Modified by Alan Miller
! Fortran 90 version by John Burkardt
! Modified by Benjamin Menetrier for BUMP
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
!----------------------------------------------------------------------
module tools_asa007

use tools_kinds, only: kind_real
use tools_repro, only: inf,infeq
use type_mpl, only: mpl_type

implicit none

real(kind_real),parameter :: eta = 1.0e-9_kind_real   ! Small parameter for the Cholesky decomposition

private
public :: asa007_cholesky,asa007_syminv

contains

!----------------------------------------------------------------------
! Subroutine: asa007_cholesky
! Purpose: compute cholesky decomposition
!----------------------------------------------------------------------
subroutine asa007_cholesky(mpl,n,nn,a,u)

implicit none

! Passed variables
type(mpl_type),intent(inout) :: mpl  ! MPI data
integer,intent(in) :: n              ! Matrix rank
integer,intent(in) :: nn             ! Half-matrix size (n*(n-1)/2)
real(kind_real),intent(in) :: a(nn)  ! Matrix
real(kind_real),intent(out) :: u(nn) ! Matrix square-root

! Local variables
integer :: i,icol,ii,irow,j,k,kk,l,m
real(kind_real) :: w,x
character(len=1024),parameter :: subr = 'asa007_cholesky'

! Initialization
ii = 0
j = 1
k = 0
if (nn/=(n*(n+1))/2) then
   call mpl%abort(subr,'wrong size in Cholesky decomposition')
end if
w = 0.0

! Factorize column by column, ICOL = column number
do icol=1,n
   ii = ii+icol
   x = eta**2*a(ii)
   l = 0
   kk = 0

   ! IROW = row number within column ICOL
   do irow=1,icol
      kk = kk+irow
      k = k+1
      w = a(k)
      m = j
      do i=1,irow-1
        l = l+1
        w = w-u(l)*u(m)
        m = m+1
      end do
      l = l+1
      if (irow==icol) exit
      if (abs(u(l))>0.0) then
         u(k) = w/u(l)
      else
         u(k) = 0.0
         if (inf(abs(x*a(k)),w**2)) call mpl%abort(subr,'A is not positive semi-definite')
      end if
   end do

   ! End of row, estimate relative accuracy of diagonal element
   if (infeq(abs(w),abs(eta*a(k)))) then
      u(k) = 0.0
   else
      if (w<0.0) call mpl%abort(subr,'A is not positive semi-definite')
      u(k) = sqrt(w)
   end if
   j = j+icol
end do

end subroutine asa007_cholesky

!----------------------------------------------------------------------
! Subroutine: asa007_syminv
! Purpose: compute inverse of a symmetric matrix
!----------------------------------------------------------------------
subroutine asa007_syminv(mpl,n,nn,a,c)

implicit none

! Passed variables
type(mpl_type),intent(inout) :: mpl  ! MPI data
integer,intent(in) :: n              ! Matrix rank
integer,intent(in) :: nn             ! Half-matrix size (n*(n-1)/2)
real(kind_real),intent(in) :: a(nn)  ! Matrix
real(kind_real),intent(out) :: c(nn) ! Matrix inverse

! Local variables
integer :: i,icol,irow,j,jcol,k,l,mdiag,ndiag,nrow
real(kind_real) :: w(n),x
character(len=1024),parameter :: subr = 'asa007_syminv'

! Initialization
nrow = n
if (nn/=(n*(n+1))/2) then
   call mpl%abort(subr,'wrong size in Cholesky decomposition')
end if
w = 0.0

! Compute the Cholesky factorization of A
call asa007_cholesky(mpl,n,nn,a,c)

! Invert C and form the product (Cinv)' * Cinv, where Cinv is the inverse of C, row by row starting with the last row
irow = nrow
ndiag = nn

do
   if (abs(c(ndiag))>0.0) then
      ! General case
      l = ndiag
      do i=irow,nrow
         w(i) = c(l)
         l = l+i
      end do
      icol = nrow
      jcol = nn
      mdiag = nn
      do
         l = jcol
         if (icol==irow) then
            x = 1.0/w(irow)
         else
            x = 0.0
         end if
         k = nrow
         do while (irow<k)
            x = x-w(k)*c(l)
            k = k-1
            l = l-1
            if (mdiag<l) l = l-k+1
         end do
         c(l) = x/w(irow)
         if (icol<=irow) exit
         mdiag = mdiag-icol
         icol = icol-1
         jcol = jcol-1
      end do
   else
      ! Special case, zero diagonal element
      l = ndiag
      do j=irow,nrow
         c(l) = 0.0
         l = l+j
      end do
   end if
   ndiag = ndiag-irow
   irow = irow-1
   if (irow<=0) exit
end do

end subroutine asa007_syminv

end module tools_asa007
