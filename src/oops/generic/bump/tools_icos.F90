!----------------------------------------------------------------------
! Module: tools_icos
!> Purpose: icosahedron routines
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2015-... UCAR, CERFACS and METEO-FRANCE
!----------------------------------------------------------------------
module tools_icos

use tools_kinds, only: kind_real
use tools_stripack, only: scoord

implicit none

private
public :: closest_icos,build_icos

contains

!----------------------------------------------------------------------
! Function: closest_icos
!> Purpose: compute the closest (larger) number of nodes for a regular icosahedron
!----------------------------------------------------------------------
subroutine closest_icos(n,fac,np)

implicit none

! Passed variable
integer,intent(in) :: n    !< Target number of nodes
integer,intent(out) :: fac !< Division factor
integer,intent(out) :: np  !< Number of nodes

! Initialization
fac = 1
np = 12

do while (np<n)
   ! Update fac and compute number of nodes
   fac = fac+1
   np = 12+10*3*(fac-1)+10*(fac-2)*(fac-1)
end do

end subroutine closest_icos

!----------------------------------------------------------------------
! Subroutine: build_icos
!> Purpose: build regular icosahedron
!----------------------------------------------------------------------
subroutine build_icos(fac,np,lon,lat)

implicit none

! Passed variable
integer,intent(in) :: fac                !< Division factor
integer,intent(in) :: np                 !< Number of nodes
real(kind_real),intent(inout) :: lon(np) !< Nodes longitudes
real(kind_real),intent(inout) :: lat(np) !< Nodes latitudes

! Local variables
integer :: a,b,c,e,f,f1,f2,i
integer :: edge_point(2,30),face_order(20),face_point(3,20)
real(kind_real) :: aa,bb,dum,phi,zz
real(kind_real) :: point_coord(3,12),node_xyz(3,np)

! Set the point coordinates
phi = 0.5*(sqrt(5.0)+1.0)
aa = phi/sqrt(1.0+phi**2)
bb = 1.0/sqrt(1.0+phi**2)
zz = 0.0
point_coord = reshape((/aa,bb,zz,aa,-bb,zz,bb,zz,aa,bb,zz,-aa,zz,aa,bb,zz,aa,-bb,zz,-aa,bb,zz,-aa,-bb,-bb,zz,aa,-bb,zz,-aa,-aa, &
            & bb,zz,-aa,-bb,zz/),(/3,12/))

! Set the edges
edge_point = reshape((/1,2,1,3,1,4,1,5,1,6,2,3,2,4,2,7,2,8,3,5,3,7,3,9,4,6,4,8,4,10,5,6,5,9,5,11,6,10,6,11,7,8,7,9,7,12,8,10,8, &
           & 12,9,11,9,12,10,11,10,12,11,12/),(/2,30/))

! Set the face orders
face_order = (/3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3/)

!  Set the faces
face_point = reshape((/1,2,4,1,3,2,1,4,6,1,5,3,1,6,5,2,3,7,2,7,8,2,8,4,3,5,9,3,9,7,4,8,10,4,10,6,5,6,11,5,11,9,6,10,11,7,9,12,7, &
           & 12,8,8,12,10,9,11,12,10,12,11/),(/3,20/))

! Generate the point coordinates
i = 0
node_xyz(:,1:12) = point_coord
i = 12
do e=1,30
   a = edge_point(1,e)
   b = edge_point(2,e)
   do f=1,fac-1
      i = i+1
      node_xyz(:,i) = (real(fac-f,kind_real)*point_coord(:,a)+real(f,kind_real)*point_coord(:,b))/real(fac,kind_real)
      node_xyz(:,i) = node_xyz(:,i)/sqrt(sum(node_xyz(:,i)**2))
   end do
end do
do f=1,20
   a = face_point(1,f)
   b = face_point(2,f)
   c = face_point(3,f)
   do f1=1,fac-1
      do f2=1,fac-f1-1
         i = i+1
         node_xyz(:,i) = (real(fac-f1-f2,kind_real)*point_coord(:,a)+real(f1,kind_real)*point_coord(:,b) &
                       & +real(f2,kind_real)*point_coord(:,c))/real(fac,kind_real)
         node_xyz(:,i) = node_xyz(:,i)/sqrt(sum(node_xyz(:,i)**2))
      end do
   end do
end do

! Transform to spherical coordinates
do i=1,np
   call scoord(node_xyz(1,i),node_xyz(2,i),node_xyz(3,i),lat(i),lon(i),dum)
end do

end subroutine build_icos

end module tools_icos
