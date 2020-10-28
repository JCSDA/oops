! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

module qg_projection_mod

use kinds
use qg_constants_mod

implicit none

private
public :: xy_to_lonlat,lonlat_to_xy
! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Convert x/y to lon/lat
subroutine xy_to_lonlat(x,y,lon,lat,mapfac)

! Passed variables
real(kind_real), intent(in) :: x
real(kind_real), intent(in) :: y
real(kind_real), intent(out) :: lon
real(kind_real), intent(out) :: lat
real(kind_real), intent(out), optional :: mapfac

! Local variables
real(kind_real) :: rlon,rlat,rlon_c,rlat_c,x_c,y_c

! Test x/y
if ((x<0.0).or.(x>domain_zonal)) call abor1_ftn('xy_to_lonlat: x point out of the domain')
if ((y<0.0).or.(y>domain_meridional)) call abor1_ftn('xy_to_lonlat: y point out of the domain')

! Define central point
rlon_c = 0.0
rlat_c = asin(f0/(2.0*omega))
x_c = 0.5*domain_zonal
y_c = 0.5*domain_meridional

! Define longitude/latitude/map factor
rlon = (x-x_c)/req+rlon_c
rlat = 2.0*atan2(exp(y/req),exp(y_c/req))-0.5*pi+rlat_c
if (present(mapfac)) mapfac = 1.0/cos(rlat)

! Convert lon to radians
lon = rlon*rad_to_deg
lat = rlat*rad_to_deg

end subroutine xy_to_lonlat
! ------------------------------------------------------------------------------
!> Convert lon/lat to x/y
subroutine lonlat_to_xy(lon,lat,x,y)

! Passed variables
real(kind_real), intent(in) :: lon
real(kind_real), intent(in) :: lat
real(kind_real), intent(out) :: x
real(kind_real), intent(out) :: y

! Local variables
real(kind_real) :: rlon,rlat,rlon_c,rlat_c,x_c,y_c

! Convert to radians
rlon = lon*deg_to_rad
rlat = lat*deg_to_rad

! Define central point
rlon_c = 0.0
rlat_c = asin(f0/(2.0*omega))
x_c = 0.5*domain_zonal
y_c = 0.5*domain_meridional

! Define x/y
x = x_c+(rlon-rlon_c)*req
y = y_c+log(tan(0.25*pi+0.5*(rlat-rlat_c)))*req

! Test x/y
if ((x<0.0).or.(x>domain_zonal)) call abor1_ftn('lonlat_to_xy: x point out of the domain')
if ((y<0.0).or.(y>domain_meridional)) call abor1_ftn('lonlat_to_xy: y point out of the domain')

end subroutine lonlat_to_xy
! ------------------------------------------------------------------------------
end module qg_projection_mod
