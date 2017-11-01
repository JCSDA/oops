! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

!> Calculate meridional wind component - Tangent Linear

!> Since the operator is linear, we just call meridional_wind.

subroutine meridional_wind_tl (v,x,nx,ny,deltax)

!--- calculate meridional wind component

use kinds

implicit none
integer, intent(in) :: nx       !< Zonal grid dimension
integer, intent(in) :: ny       !< Meridional grid dimension
real(kind=kind_real), intent(out) :: v(nx,ny,2) !< Meridional wind increment
real(kind=kind_real), intent(in)  :: x(nx,ny,2) !< Streamfunction increment
real(kind=kind_real), intent(in) :: deltax      !< Zonal grid spacing (non-dimensional)

!--- meridional_wind is a linear operator

call meridional_wind (v,x,nx,ny,deltax)

end subroutine meridional_wind_tl
