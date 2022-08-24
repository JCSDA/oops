! (C) Copyright 2009-2016 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

module qg_convert_x_to_q_mod

use differential_solver_mod
use kinds
!$ use omp_lib
use qg_geom_mod

implicit none

private
public :: convert_x_to_q,convert_x_to_q_tl,convert_x_to_q_ad
! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Convert streamfunction to potential vorticity
subroutine convert_x_to_q(geom,x,x_north,x_south,q)

implicit none

! Passed variables
type(qg_geom),intent(in) :: geom                          !< Geometry
real(kind_real),intent(in) :: x(geom%nx,geom%ny,geom%nz)  !< Streamfunction
real(kind_real),intent(in) :: x_north(geom%nz)            !< Streamfunction on northern wall
real(kind_real),intent(in) :: x_south(geom%nz)            !< Streamfunction on southern wall
real(kind_real),intent(inout) :: q(geom%nx,geom%ny,geom%nz) !< Potential vorticity

! Local variables
integer :: ix,iy,iz
real(kind_real) :: del2x(geom%nx,geom%ny,geom%nz)
real(kind_real) :: zz

! Laplacian of the streamfunction
call laplacian_2d(geom,x,del2x)

! Vertical differences
!$omp parallel do schedule(static) private(iz,iy,ix)
do iz=1,geom%nz
  do iy=1,geom%ny
    do ix=1,geom%nx
      q(ix,iy,iz) = del2x(ix,iy,iz)+sum(geom%f(iz,:)*x(ix,iy,:))
    end do
  end do
end do
!$omp end parallel do

! Add the contribution from the boundaries
zz = 1.0 / (geom%deltay * geom%deltay)
!$omp parallel do schedule(static) private(iz)
do iz=1,geom%nz
  q(:,1,iz) = q(:,1,iz)+x_south(iz)*zz
  q(:,geom%ny,iz) = q(:,geom%ny,iz)+x_north(iz)*zz
enddo
!$omp end parallel do

! Add the beta term and the heating term
!$omp parallel do schedule(static) private(iy)
do iy=1,geom%ny
  q(:,iy,:) = q(:,iy,:)+geom%bet(iy)
  if(geom%ph_coeff/=0) then
    q(:,iy,1) = q(:,iy,1)+geom%heat(:,iy)*tanh(geom%ph_coeff*q(:,iy,1))
  else
    q(:,iy,1) = q(:,iy,1)+geom%heat(:,iy)
  end if 
enddo
!$omp end parallel do

end subroutine convert_x_to_q
! ------------------------------------------------------------------------------
!> Convert streamfunction to potential vorticity - tangent linear
subroutine convert_x_to_q_tl(geom,x,q)

implicit none

! Passed variables
type(qg_geom),intent(in) :: geom                          !< Geometry
real(kind_real),intent(in) :: x(geom%nx,geom%ny,geom%nz)  !< Streamfunction
real(kind_real),intent(inout) :: q(geom%nx,geom%ny,geom%nz) !< Potential vorticity

! Local variables
integer :: ix,iy,iz
real(kind_real) :: del2x(geom%nx,geom%ny,geom%nz)

! Laplacian of the streamfunction
call laplacian_2d(geom,x,del2x)

! Vertical differences
!$omp parallel do schedule(static) private(iz,iy,ix)
do iz=1,geom%nz
  do iy=1,geom%ny
    do ix=1,geom%nx
      q(ix,iy,iz) = del2x(ix,iy,iz)+sum(geom%f(iz,:)*x(ix,iy,:))
    end do
  end do
end do
!$omp end parallel do

! Add the contribution from the boundaries (=> TL of this is identity)

! Add the beta term and the heating term (=> TL of this is identity)

end subroutine convert_x_to_q_tl
! ------------------------------------------------------------------------------
!> Convert streamfunction to potential vorticity - adjoint
subroutine convert_x_to_q_ad(geom,q,x)

implicit none

! Passed variables
type(qg_geom),intent(in) :: geom                            !< Geometry
real(kind_real),intent(in) :: q(geom%nx,geom%ny,geom%nz)    !< Potential vorticity
real(kind_real),intent(inout) :: x(geom%nx,geom%ny,geom%nz) !< Streamfunction

! Local variables
integer :: ix,iy,iz
real(kind_real) :: del2x(geom%nx,geom%ny,geom%nz)

! Add the beta term and the heating term (=> AD of this is identity)

! Add the contribution from the boundaries (=> AD of this is identity)

! Initialization
del2x = 0.0

! Vertical differences
do iz=1,geom%nz
  do iy=1,geom%ny
    do ix=1,geom%nx
      del2x(ix,iy,iz) = del2x(ix,iy,iz)+q(ix,iy,iz)
      x(ix,iy,iz) = x(ix,iy,iz)+sum(geom%f(:,iz)*q(ix,iy,:))
    end do
  end do
end do

! Laplacian of the streamfunction
call laplacian_2d_ad(geom,del2x,x)

end subroutine convert_x_to_q_ad
! ------------------------------------------------------------------------------
end module qg_convert_x_to_q_mod
