!----------------------------------------------------------------------
! Module: tools_dirac
!> Purpose: Dirac routines
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module tools_dirac

use tools_const, only: deg2rad
use tools_kinds, only: kind_real
use type_ctree, only: ctreetype,create_ctree,delete_ctree,find_nearest_neighbors

implicit none

private
public :: setup_dirac_points

contains

!----------------------------------------------------------------------
! Subroutine: setup_dirac_points
!> Purpose: setup Dirac points with a cover tree
!----------------------------------------------------------------------
subroutine setup_dirac_points(n,lon,lat,lmask,ndir,londir,latdir,ic0dir)

implicit none

! Passed variables
integer,intent(in) :: n
real(kind_real),intent(in) :: lon(n)
real(kind_real),intent(in) :: lat(n)
logical,intent(in) :: lmask(n)
integer,intent(in) :: ndir
real(kind_real),intent(in) :: londir(ndir)
real(kind_real),intent(in) :: latdir(ndir)
integer,intent(out) :: ic0dir(ndir)

! Local variables
integer :: i,idir
integer :: mask_ctree(n)
real(kind_real) :: dum(1)
type(ctreetype) :: ctree

! Compute cover tree
do i=1,n
   if (lmask(i)) then
      mask_ctree(i) = 1
   else
      mask_ctree(i) = 0
   end if
end do
ctree = create_ctree(n,dble(lon),dble(lat),mask_ctree)

! Compute nearest neighbors
do idir=1,ndir
   call find_nearest_neighbors(ctree,dble(londir(idir)*deg2rad), &
 & dble(latdir(idir)*deg2rad),1,ic0dir(idir:idir),dum)
end do

! Release memory
call delete_ctree(ctree)

end subroutine setup_dirac_points

end module tools_dirac
