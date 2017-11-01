! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

!> Invert potential vorticity - Adjoint

subroutine invert_pv_ad (x,pv,nx,ny,deltax,deltay,F1,F2)

use kinds

implicit none
integer, intent(in) :: nx           !< Zonal grid dimension
integer, intent(in) :: ny           !< Meridional grid dimension
real(kind=kind_real), intent(inout) :: x(nx,ny,2)   !< Streamfunction adjoint variable
real(kind=kind_real), intent(inout)  :: pv(nx,ny,2) !< PV adjoint variable
real(kind=kind_real), intent(in) :: deltax          !< Zonal grid spacing (non-dimensional)
real(kind=kind_real), intent(in) :: deltay          !< Meridional grid spacing (non-dimensional)
real(kind=kind_real), intent(in) :: F1              !< Coefficient in PV operator
real(kind=kind_real), intent(in) :: F2              !< Coefficient in PV operator

real(kind=kind_real) :: pv_nobc(nx,ny,2)

!--- Solve the elliptic system. Here, we use the symmetry of the operator,
!--- rather than coding the adjoint of the solver.

x(:,:,1) = -F1*x(:,:,1)
x(:,:,2) = -F2*x(:,:,2)

call solve_elliptic_system (pv_nobc,x,nx,ny,deltax,deltay,F1,F2)

x(:,:,:) = 0.0_kind_real
pv(:,:,1) = pv(:,:,1) + pv_nobc(:,:,1)*(-1.0_kind_real/F1)
pv(:,:,2) = pv(:,:,2) + pv_nobc(:,:,2)*(-1.0_kind_real/F2)

end subroutine invert_pv_ad
