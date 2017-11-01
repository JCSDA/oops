! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

!> Potential vorticity operator - Adjoint

subroutine pv_operator_ad (x,pv,nx,ny,F1,F2,deltax,deltay)

use kinds

implicit none
integer, intent(in) :: nx          !< Zonal grid dimension
integer, intent(in) :: ny          !< Meridional grid dimension
real(kind=kind_real), intent(inout) :: x(nx,ny,2)  !< Streamfunction adjoint variable
real(kind=kind_real), intent(inout) :: pv(nx,ny,2) !< PV-operator adjoint variable
real(kind=kind_real), intent(in) :: F1             !< Parameter in the PV operator
real(kind=kind_real), intent(in) :: F2             !< Parameter in the PV operator
real(kind=kind_real), intent(in) :: deltax         !< Zonal grid spacing (non-dimensional)
real(kind=kind_real), intent(in) :: deltay         !< Meridional grid spacing (non-dimensional)

!--- vertical differences:

x(:,:,1) = x(:,:,1) - F1*pv(:,:,1)
x(:,:,2) = x(:,:,2) + F1*pv(:,:,1)
x(:,:,1) = x(:,:,1) + F2*pv(:,:,2)
x(:,:,2) = x(:,:,2) - F2*pv(:,:,2)

!--- del-squared of the streamfunction (5-point laplacian)

x(:,1:ny-1,:) = x(:,1:ny-1,:) + (1.0_kind_real/(deltay*deltay))*pv(:,2:ny  ,:)
x(:,2:ny  ,:) = x(:,2:ny  ,:) + (1.0_kind_real/(deltay*deltay))*pv(:,1:ny-1,:)

x(nx    ,:,:) = x(nx    ,:,:) + (1.0_kind_real/(deltax*deltax))*pv(1     ,:,:)

x(1:nx-1,:,:) = x(1:nx-1,:,:)  + (1.0_kind_real/(deltax*deltax))*pv(2:nx  ,:,:) 

x(1     ,:,:) = x(1     ,:,:)  + (1.0_kind_real/(deltax*deltax))*pv(nx    ,:,:)
x(2:nx  ,:,:) = x(2:nx  ,:,:)  + (1.0_kind_real/(deltax*deltax))*pv(1:nx-1,:,:)


x(:,:,:) = x(:,:,:) -2.0_kind_real*( 1.0_kind_real/(deltax*deltax) &
                                  & +1.0_kind_real/(deltay*deltay))*pv(:,:,:)

pv(:,:,:) = 0.0_kind_real

end subroutine pv_operator_ad
