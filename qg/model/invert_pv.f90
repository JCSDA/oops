! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

!> Invert potential vorticity, returning streamfunction

!> Streamfunction is determined from potential vorticity by subtracting
!! the beta and orography terms and then solving the resulting elliptic
!! equation.

subroutine invert_pv (x,pv,x_north,x_south,rs,nx,ny,deltax,deltay,F1,F2,bet)

!--- invert potential vorticity to get streamfunction

use kinds

implicit none
integer, intent(in) :: nx          !< Zonal grid dimension
integer, intent(in) :: ny          !< Meridional grid dimension
real(kind=kind_real), intent(out) :: x(nx,ny,2)    !< Streamfunction
real(kind=kind_real), intent(in)  :: pv(nx,ny,2)   !< Potential vorticity
real(kind=kind_real), intent(in)  :: x_north(2)    !< Streamfunction on northern wall
real(kind=kind_real), intent(in)  :: x_south(2)    !< Streamfunction on southern wall
real(kind=kind_real), intent(in)  :: rs(nx,ny)     !< Orography
real(kind=kind_real), intent(in) :: deltax         !< Zonal grid spacing (non-dimensional)
real(kind=kind_real), intent(in) :: deltay         !< Meridional grid spacing (non-dimensional)
real(kind=kind_real), intent(in) :: F1             !< Coefficient in PV operator
real(kind=kind_real), intent(in) :: F2             !< Coefficient in PV operator
real(kind=kind_real), intent(in) :: bet            !< NS Gradient of Coriolis parameter

real(kind=kind_real) :: y
real(kind=kind_real) :: pv_nobc(nx,ny,2)
integer :: jj

!--- subtract the beta term and the orography/heating term

do jj=1,ny
  y = real(jj-(ny+1)/2,kind_real)*deltay;
  pv_nobc(:,jj,1) = pv(:,jj,1) - bet*y
  pv_nobc(:,jj,2) = pv(:,jj,2) - bet*y -rs(:,jj)
enddo

!--- subtract the contribution from the boundaries

pv_nobc(:,1 ,1) = pv_nobc(:,1 ,1) - (1.0_kind_real/(deltay*deltay))*x_south(1)
pv_nobc(:,1 ,2) = pv_nobc(:,1 ,2) - (1.0_kind_real/(deltay*deltay))*x_south(2)
pv_nobc(:,ny,1) = pv_nobc(:,ny,1) - (1.0_kind_real/(deltay*deltay))*x_north(1)
pv_nobc(:,ny,2) = pv_nobc(:,ny,2) - (1.0_kind_real/(deltay*deltay))*x_north(2)

!--- Solve the elliptic system

call solve_elliptic_system (x,pv_nobc,nx,ny,deltax,deltay,F1,F2)

end subroutine invert_pv
