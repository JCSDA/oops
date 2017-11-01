! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

!> Calculate meridional wind component from the streamfunction

!> Nothing fancy.
!! It's just a standard centred finite difference.

subroutine meridional_wind (v,x,nx,ny,deltax)

!--- calculate meridional wind component

use kinds

implicit none
integer, intent(in) :: nx       !< Zonal grid dimension
integer, intent(in) :: ny       !< Meridional grid dimension
real(kind=kind_real), intent(out) :: v(nx,ny,2) !< Meridional wind
real(kind=kind_real), intent(in)  :: x(nx,ny,2) !< Streamfunction
real(kind=kind_real), intent(in) :: deltax      !< Zonal grid spacing (non-dimensional)

v(1:nx-1,:,:) = (0.5_kind_real/deltax)*x(2:nx,:,:)
v(nx    ,:,:) = (0.5_kind_real/deltax)*x(1   ,:,:)
v(2:nx  ,:,:) = v(2:nx,:,:) - (0.5_kind_real/deltax)*x(1:nx-1,:,:)
v(1     ,:,:) = v(1   ,:,:) - (0.5_kind_real/deltax)*x(nx    ,:,:)

end subroutine meridional_wind
