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

use tools_display, only: msgerror
use tools_kinds, only: kind_real
use tools_missing, only: msr,isnotmsr
use tools_qsort, only: qsort

implicit none

! Constants
real(kind_real),parameter :: pi=acos(-1.0)    !< Pi
real(kind_real),parameter :: deg2rad=pi/180.0 !< Degree to radian
real(kind_real),parameter :: rad2deg=180.0/pi !< Radian to degree
real(kind_real),parameter :: req=6371229.0    !< Earth radius (m)
real(kind_real),parameter :: reqkm=6371.229   !< Earth radius (km)
real(kind_real),parameter :: ps=101325.0      !< Reference surface pressure

! Internal parameters
real(kind_real),parameter :: qtrim = 0.05     !< Fraction for which upper and lower quantiles are removed in trimmed averages
integer,parameter :: ntrim = 1                !< Minimum number of remaining points for the trimmed average

private
public :: pi,deg2rad,rad2deg,req,reqkm,ps
public :: lonmod,sphere_dist,reduce_arc,vector_product,vector_triple_product,gc99,median,taverage,add,divide,fac

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
   ! Great-circle distance using Vincenty formula on the unit sphere
   dist = atan2(sqrt((cos(lat_f)*sin(lon_f-lon_i))**2 &
        & +(cos(lat_i)*sin(lat_f)-sin(lat_i)*cos(lat_f)*cos(lon_f-lon_i))**2), &
        & sin(lat_i)*sin(lat_f)+cos(lat_i)*cos(lat_f)*cos(lon_f-lon_i))
else
   call msr(dist)
end if

end subroutine sphere_dist

!----------------------------------------------------------------------
! Subroutine: reduce_arc
!> Purpose: reduce arc to a given distance
!----------------------------------------------------------------------
subroutine reduce_arc(lon_i,lat_i,lon_f,lat_f,maxdist,dist)

implicit none

! Passed variables
real(kind_real),intent(in) :: lon_i    !< Initial point longitude
real(kind_real),intent(in) :: lat_i    !< Initial point latitude
real(kind_real),intent(inout) :: lon_f !< Final point longitude
real(kind_real),intent(inout) :: lat_f !< Final point latitude
real(kind_real),intent(in) :: maxdist  !< Maximum distance
real(kind_real),intent(out) :: dist    !< Effective distance

! Local variable
real(kind_real) :: theta

! Compute distance
call sphere_dist(lon_i,lat_i,lon_f,lat_f,dist)

! Check with the maximum distance
if (dist>maxdist) then
   ! Compute bearing
   theta = atan2(sin(lon_f-lon_i)*cos(lat_f),cos(lat_i)*sin(lat_f)-sin(lat_i)*cos(lat_f)*cos(lon_f-lon_i))

   ! Reduce distance
   dist = maxdist

   ! Compute new point
   lat_f = asin(sin(lat_i)*cos(dist)+cos(lat_i)*sin(dist)*cos(theta))
   lon_f = lon_i+atan2(sin(theta)*sin(dist)*cos(lat_i),cos(dist)-sin(lat_i)*sin(lat_f))
end if

end subroutine reduce_arc

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

!----------------------------------------------------------------------
! Function: gc99
!> Purpose: Gaspari and Cohn (1999) function, with the support radius as a parameter
!----------------------------------------------------------------------
function gc99(distnorm)

! Passed variables
real(kind_real),intent(in) :: distnorm !< Normalized distance

! Returned variable
real(kind_real) :: gc99

! Distance check bound
if (distnorm<0.0) call msgerror('negative normalized distance')

if (.true.) then
   ! Gaspari and Cohn (1999) function
   if (distnorm<0.5) then
      gc99 = 1.0-8.0*distnorm**5+8.0*distnorm**4+5.0*distnorm**3-20.0/3.0*distnorm**2
   else if (distnorm<1.0) then
      gc99 = 8.0/3.0*distnorm**5-8.0*distnorm**4+5.0*distnorm**3+20.0/3.0*distnorm**2-10.0*distnorm+4.0-1.0/(3.0*distnorm)
   else
      gc99 = 0.0
   end if
else
   ! Gaussian equivalent
   if (distnorm<1.0) then
      gc99 = exp(-0.5*(3.53*distnorm)**2)
   else
      gc99 = 0.0
   end if
end if

return

end function gc99

!----------------------------------------------------------------------
! Function: median
!> Purpose: compute median of a list
!----------------------------------------------------------------------
function median(n,list)

implicit none

! Passed variables
integer,intent(in) :: n               !< Size of the list
real(kind_real),intent(in) :: list(n) !< List

! Returned variable
real(kind_real) :: median

! Local variables
integer :: order(n)
real(kind_real) :: list_copy(n)

! Copy list
list_copy = list

! Order array
call qsort(n,list_copy,order)

! Get median
call msr(median)
if (mod(n,2)==0) then
   ! Even number of values
   median = 0.5*(list_copy(n/2)+list_copy(n/2+1))
else
   ! Odd number of values
   median = list_copy((n+1)/2)
end if

return

end function median

!----------------------------------------------------------------------
! Function: taverage
!> Purpose: compute the trimmed average
!----------------------------------------------------------------------
real(kind_real) function taverage(n,list)

implicit none

! Passed variables
integer,intent(in) :: n        !< Number of values
real(kind_real),intent(in) :: list(n)     !< List values

! Local variable
integer :: nrm,nvalid
integer :: order(n)
real(kind_real) :: list_copy(n)

! Copy list
list_copy = list

! Compute the number of values to remove
nrm = floor(n*qtrim)

if (n-2*nrm>=ntrim) then
   ! Order array
   call qsort(n,list_copy,order)

   ! Compute trimmed average
   nvalid = count(isnotmsr(list_copy(1+nrm:n-nrm)))
   if (nvalid>0) then
      taverage = sum(list_copy(1+nrm:n-nrm),mask=isnotmsr(list_copy(1+nrm:n-nrm)))/float(nvalid)
   else
      call msr(taverage)
   end if
else
   ! Missing value
   call msr(taverage)
end if

end function taverage

!----------------------------------------------------------------------
! Subroutine: add
!> Purpose: check if missing and add
!----------------------------------------------------------------------
subroutine add(value,cumul,num,wgt)

implicit none

! Passed variables
real(kind_real),intent(in) :: value        !< Value to add
real(kind_real),intent(inout) :: cumul     !< Cumul
real(kind_real),intent(inout) :: num       !< Number of values
real(kind_real),intent(in),optional :: wgt !< Weight

! Local variables
real(kind_real) :: lwgt

! Initialize weight
lwgt = 1.0
if (present(wgt)) lwgt = wgt

! Add value to cumul
if (isnotmsr(value)) then
   cumul = cumul+lwgt*value
   num = num+1.0
end if

end subroutine add

!----------------------------------------------------------------------
! Subroutine: divide
!> Purpose: check if missing and divide
!----------------------------------------------------------------------
subroutine divide(value,num)

implicit none

! Passed variables
real(kind_real),intent(inout) :: value !< Value to divide
real(kind_real),intent(in) :: num      !< Divider

! Divide cumul by num
if (abs(num)>0.0) then
   value = value/num
else
   call msr(value)
end if

end subroutine divide

!----------------------------------------------------------------------
! Function: fac
!> Purpose: factorial
!----------------------------------------------------------------------
function fac(n)

implicit none

! Result
integer :: fac

! Passed variables
integer,intent(in) :: n

! Local variables
integer :: j

if (n<0) call msgerror('factorial requires a non-negative argument')

fac = 1
do j=2,n
   fac = fac*j
end do

end function fac

end module tools_const
