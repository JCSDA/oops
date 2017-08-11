!----------------------------------------------------------------------
! Module: tools_const
!> Purpose: usual constants
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module tools_const

use tools_kinds,only: kind_real
use tools_missing, only: msr,isnotmsr
implicit none

! Constants
real(kind_real),parameter :: pi=acos(-1.0)    !< Pi
real(kind_real),parameter :: deg2rad=pi/180.0 !< Degree to radian
real(kind_real),parameter :: rad2deg=180.0/pi !< Radian to degree
real(kind_real),parameter :: req=6.371e6      !< Earth radius
real(kind_real),parameter :: epsilon=1.0e-12  !< Small parameter
real(kind_real),parameter :: ps=101325.0      !< Reference surface pressure

private
public :: pi,deg2rad,rad2deg,req,ps,epsilon
public :: lonmod,sphere_dist,vector_product,vector_triple_product

contains

!----------------------------------------------------------------------
! Function: lonmod
!> Purpose: set longitude between -pi and pi
!----------------------------------------------------------------------
real(kind_real) function lonmod(lon)

implicit none

! Passed variables
real(kind_real),intent(in) :: lon !< Longitude

! Check bounds
lonmod = lon
if (lonmod>pi) then
   lonmod = lonmod-2.0*pi
elseif (lonmod<-pi) then
   lonmod = lonmod+2.0*pi
end if

end function lonmod

!----------------------------------------------------------------------
! Function: sphere_dist
!> Purpose: compute the great-circle distance between two points
!----------------------------------------------------------------------
subroutine sphere_dist(lon_i,lat_i,lon_f,lat_f,dist)

implicit none

! Passed variable
real(kind_real),intent(in) :: lon_i !< Initial point longitude (radian)
real(kind_real),intent(in) :: lat_i !< Initial point latitude (radian)
real(kind_real),intent(in) :: lon_f !< Final point longitude (radian)
real(kind_real),intent(in) :: lat_f !< Final point longilatitudetude (radian)
real(kind_real),intent(out) :: dist !< Great-circle distance

! Check that there is no missing value
if (isnotmsr(lon_i).and.isnotmsr(lat_i).and.isnotmsr(lon_f).and.isnotmsr(lat_f)) then
   ! Great-circle distance using Vincenty formula on the sphere
    dist = req*atan2(sqrt((cos(lat_f)*sin(lon_f-lon_i))**2 &
         & +(cos(lat_i)*sin(lat_f)-sin(lat_i)*cos(lat_f)*cos(lon_f-lon_i))**2), & 
         & sin(lat_i)*sin(lat_f)+cos(lat_i)*cos(lat_f)*cos(lon_f-lon_i))
else
   call msr(dist)
end if

end subroutine sphere_dist

!----------------------------------------------------------------------
! Subroutine: vector_product
!> Purpose: compute normalized vector product
!----------------------------------------------------------------------
subroutine vector_product(v1,v2,vp)

implicit none

! Passed variables
real(kind_real),intent(in) :: v1(3)  !< First vector
real(kind_real),intent(in) :: v2(3)  !< Second vector
real(kind_real),intent(out) :: vp(3) !< Vector product

! Local variable
real(kind_real) :: r

! Vector product
vp(1) = v1(2)*v2(3)-v1(3)*v2(2)
vp(2) = v1(3)*v2(1)-v1(1)*v2(3)
vp(3) = v1(1)*v2(2)-v1(2)*v2(1)

! Normalization
r = sqrt(sum(vp**2))
if (r>0.0) vp = vp/r

end subroutine vector_product

!----------------------------------------------------------------------
! Subroutine: vector_triple_product
!> Purpose: compute vector triple product
!----------------------------------------------------------------------
subroutine vector_triple_product(v1,v2,v3,p)

implicit none

! Passed variables
real(kind_real),intent(in) :: v1(3) !< First vector
real(kind_real),intent(in) :: v2(3) !< Second vector
real(kind_real),intent(in) :: v3(3) !< Third vector
real(kind_real),intent(out) :: p    !< Triple product

! Local variable
real(kind_real) :: vp(3)

! Vector product
vp(1) = v1(2)*v2(3)-v1(3)*v2(2)
vp(2) = v1(3)*v2(1)-v1(1)*v2(3)
vp(3) = v1(1)*v2(2)-v1(2)*v2(1)

! Scalar product
p = sum(vp*v3)

end subroutine vector_triple_product

end module tools_const
