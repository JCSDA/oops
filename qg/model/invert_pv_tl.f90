! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

!> Invert potential vorticity - Tangent Linear

subroutine invert_pv_tl (x,pv,nx,ny,deltax,deltay,F1,F2)

!--- invert potential vorticity to get streamfunction

use kinds

implicit none
integer, intent(in) :: nx         !< Zonal grid dimension
integer, intent(in) :: ny         !< Meridional grid dimension
real(kind=kind_real), intent(out) :: x(nx,ny,2)   !< Streamfunction increment
real(kind=kind_real), intent(in)  :: pv(nx,ny,2)  !< Potential vorticity increment
real(kind=kind_real), intent(in) :: deltax        !< Zonal grid spacing (non-dimensional)
real(kind=kind_real), intent(in) :: deltay        !< Meridional grid spacing (non-dimensional)
real(kind=kind_real), intent(in) :: F1            !< Coefficient in PV operator
real(kind=kind_real), intent(in) :: F2            !< Coefficient in PV operator

!--- Solve the elliptic system

call solve_elliptic_system (x,pv,nx,ny,deltax,deltay,F1,F2)

end subroutine invert_pv_tl
