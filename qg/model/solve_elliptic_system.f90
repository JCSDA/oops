! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

!> Solves the elliptic system that arises when inverting potential vorticity.

!> Here, we are solving the coupled system:
!! \f{eqnarray}{
!! \nabla^2 \psi_1 - F_1 (\psi_1 -\psi_2 ) &=& p_1 \\\\
!! \nabla^2 \psi_2 - F_2 (\psi_2 -\psi_1 ) &=& p_2 \nonumber
!! \f}
!! where the subscript refers to model level, and
!! \f$ \nabla^2 \f$ is the 5-point finite-difference two-dimensional Laplacian.
!!
!! We reduce this to a 2D problem, which we solve for
!! \f$ \nabla^2 \psi_1 \f$:
!! \f[
!!  \nabla^2 \left(\nabla^2 \psi_1 \right) - (F_1 + F_2 ) \nabla^2 \psi_1 =
!!   \nabla^2 p_1 - F_2 p_1 - F_1 p_2
!! \f]
!!
!! Having found \f$ \nabla^2 \psi_1 \f$, we invert the Laplacian
!! to determine \f$ \psi_1 \f$, and then get \f$ \psi_2 \f$ by substitution
!! in equation 1 above.

subroutine solve_elliptic_system (x,pv,nx,ny,deltax,deltay,F1,F2)

!--- invert potential vorticity to get streamfunction

use kinds
implicit none

integer, intent(in) :: nx        !< Zonal grid dimension
integer, intent(in) :: ny        !< Meridional grid dimension
real(kind=kind_real), intent(out) :: x(nx,ny,2)  !< Streamfunction
real(kind=kind_real), intent(in)  :: pv(nx,ny,2) !< Right hand side
real(kind=kind_real), intent(in) :: deltax       !< Zonal grid spacing (non-dimensional)
real(kind=kind_real), intent(in) :: deltay       !< Meridional grid spacing (non-dimensional)
real(kind=kind_real), intent(in) :: F1           !< Coefficient in the PV operator
real(kind=kind_real), intent(in) :: F2           !< Coefficient in the PV operator

real(kind=kind_real) :: rhs(nx,ny), del2x1(nx,ny)

!--- Solve the 2D problem

call laplacian_2d (pv(:,:,1),rhs,nx,ny,deltax,deltay)
rhs(:,:) = rhs(:,:) - F2*pv(:,:,1) - F1*pv(:,:,2)

call solve_helmholz (del2x1,rhs,-(F1+F2),nx,ny,deltax,deltay)

!--- Invert Laplacian to get x on level 1

call solve_helmholz (x(:,:,1),del2x1,0.0_kind_real,nx,ny,deltax,deltay)

!--- Calculate x on level 2

x(:,:,2) = x(:,:,1) + (pv(:,:,1)-del2x1(:,:))*(1.0_kind_real/F1)

end subroutine solve_elliptic_system
