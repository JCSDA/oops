! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

!> Calculate zonal wind - Adjoint

subroutine zonal_wind_ad (u,x,nx,ny,deltay)

use kinds
implicit none

integer, intent(in) :: nx         !< Zonal grid dimension
integer, intent(in) :: ny         !< Meridional grid dimension
real(kind=kind_real), intent(inout) :: u(nx,ny,2) !< Zonal wind adjoint variable
real(kind=kind_real), intent(inout) :: x(nx,ny,2) !< Streamfunction adjoint variable
real(kind=kind_real), intent(in) :: deltay        !< Meridional grid spacing (non-dimensional)

x(:,2:ny  ,:) = x(:,2:ny  ,:) - (0.5_kind_real/deltay)*u(:,1:ny-1,:)
x(:,1:ny-1,:) = x(:,1:ny-1,:) + (0.5_kind_real/deltay)*u(:,2:ny  ,:)
u(:,:,:) = 0.0_kind_real

end subroutine zonal_wind_ad
