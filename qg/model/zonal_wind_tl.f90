! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

!> Calculate zonal wind - Tangent Linear

subroutine zonal_wind_tl (u,x,nx,ny,deltay)

!--- calculate zonal wind component

use kinds
implicit none

integer, intent(in) :: nx       !< Zonal grid dimension
integer, intent(in) :: ny       !< Meridional grid dimension
real(kind=kind_real), intent(out) :: u(nx,ny,2) !< Zonal wind increment
real(kind=kind_real), intent(in)  :: x(nx,ny,2) !< Streamfunction increment
real(kind=kind_real), intent(in) :: deltay      !< Meridional grid spacing (non-dimensional)

u(:,2:ny  ,:) =                 (0.5_kind_real/deltay)*x(:,1:ny-1,:)
u(:,1     ,1) = 0.0_kind_real
u(:,1     ,2) = 0.0_kind_real
u(:,1:ny-1,:) = u(:,1:ny-1,:) - (0.5_kind_real/deltay)*x(:,2:ny  ,:)

end subroutine zonal_wind_tl
