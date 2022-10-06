! (C) Copyright 2009-2016 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

module qg_convert_q_to_x_mod

use differential_solver_mod
use kinds
!$ use omp_lib
use qg_geom_mod

implicit none

private
public :: convert_q_to_x,convert_q_to_x_tl,convert_q_to_x_ad
! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Convert potential vorticity to streamfunction
subroutine convert_q_to_x(geom,q,x_north,x_south,x)

implicit none

! Passed variables
type(qg_geom),intent(in) :: geom                          !< Geometry
real(kind_real),intent(in) :: q(geom%nx,geom%ny,geom%nz)  !< Potential vorticity
real(kind_real),intent(in) :: x_north(geom%nz)            !< Streamfunction on northern wall
real(kind_real),intent(in) :: x_south(geom%nz)            !< Streamfunction on southern wall
real(kind_real),intent(inout) :: x(geom%nx,geom%ny,geom%nz) !< Streamfunction

! Local variables
integer :: ix,iy,iz
real(kind_real) :: q_nobc(geom%nx,geom%ny,geom%nz),pinv_q(geom%nx,geom%ny,geom%nz),pinv_x(geom%nx,geom%ny,geom%nz)
real(kind_real) :: zz

! Subtract the beta term and the heating term
!$omp parallel do schedule(static) private(iy)
do iy=1,geom%ny
  q_nobc(:,iy,:) = q(:,iy,:)-geom%bet(iy)
  if(geom%ph_coeff/=0) then
    q_nobc(:,iy,1) = q_nobc(:,iy,1)-geom%heat(:,iy)*tanh(geom%ph_coeff*q_nobc(:,iy,1))
  else
    q_nobc(:,iy,1) = q_nobc(:,iy,1)-geom%heat(:,iy)
  endif
enddo
!$omp end parallel do

! Subtract the contribution from the boundaries
zz = 1.0 / (geom%deltay * geom%deltay)
!$omp parallel do schedule(static) private(iz)
do iz=1,geom%nz
  q_nobc(:,1,iz) = q_nobc(:,1,iz)-x_south(iz)*zz
  q_nobc(:,geom%ny,iz) = q_nobc(:,geom%ny,iz)-x_north(iz)*zz
enddo
!$omp end parallel do

! Apply ff_pinv
!$omp parallel do schedule(static) private(iz,iy,ix)
do iz=1,geom%nz
  do iy=1,geom%ny
    do ix=1,geom%nx
      pinv_q(ix,iy,iz) = sum(geom%f_pinv(iz,:)*q_nobc(ix,iy,:))
    end do
  end do
end do
!$omp end parallel do

! Solve Helmholz equation for each layer
do iz=1,geom%nz
  call solve_helmholz(geom,geom%f_d(iz),pinv_q(:,:,iz),pinv_x(:,:,iz))
end do

! Apply ff_p
!$omp parallel do schedule(static) private(iz,iy,ix)
do iz=1,geom%nz
  do iy=1,geom%ny
    do ix=1,geom%nx
      x(ix,iy,iz) = sum(geom%f_p(iz,:)*pinv_x(ix,iy,:))
    end do
  end do
end do
!$omp end parallel do

end subroutine convert_q_to_x
! ------------------------------------------------------------------------------
!> Convert potential vorticity to streamfunction - tangent linear
subroutine convert_q_to_x_tl(geom,q,x)

implicit none

! Passed variables
type(qg_geom),intent(in) :: geom                          !< Geometry
real(kind_real),intent(in) :: q(geom%nx,geom%ny,geom%nz)  !< Potential vorticity
real(kind_real),intent(inout) :: x(geom%nx,geom%ny,geom%nz) !< Streamfunction

! Local variables
integer :: ix,iy,iz
real(kind_real) :: pinv_q(geom%nx,geom%ny,geom%nz),pinv_x(geom%nx,geom%ny,geom%nz)

! Subtract the beta term and the heating term (=> TL of this is identity)

! Subtract the contribution from the boundaries (=> TL of this is identity)

! Apply ff_pinv
!$omp parallel do schedule(static) private(iz,iy,ix)
do iz=1,geom%nz
  do iy=1,geom%ny
    do ix=1,geom%nx
      pinv_q(ix,iy,iz) = sum(geom%f_pinv(iz,:)*q(ix,iy,:))
    end do
  end do
end do
!$omp end parallel do

! Solve Helmholz equation for each layer
do iz=1,geom%nz
  call solve_helmholz(geom,geom%f_d(iz),pinv_q(:,:,iz),pinv_x(:,:,iz))
end do

! Apply ff_p
!$omp parallel do schedule(static) private(iz,iy,ix)
do iz=1,geom%nz
  do iy=1,geom%ny
    do ix=1,geom%nx
      x(ix,iy,iz) = sum(geom%f_p(iz,:)*pinv_x(ix,iy,:))
    end do
  end do
end do
!$omp end parallel do

end subroutine convert_q_to_x_tl
! ------------------------------------------------------------------------------
!> Convert potential vorticity to streamfunction - adjoint
subroutine convert_q_to_x_ad(geom,x,q)

implicit none

! Passed variables
type(qg_geom),intent(in) :: geom                            !< Geometry
real(kind_real),intent(in) :: x(geom%nx,geom%ny,geom%nz)    !< Streamfunction
real(kind_real),intent(inout) :: q(geom%nx,geom%ny,geom%nz) !< Potential vorticity

! Local variables
integer :: ix,iy,iz
real(kind_real) :: pinv_q(geom%nx,geom%ny,geom%nz),pinv_x(geom%nx,geom%ny,geom%nz),qtmp(geom%nx,geom%ny,geom%nz)

! Apply ff_p
!$omp parallel do schedule(static) private(iz,iy,ix)
do iz=1,geom%nz
  do iy=1,geom%ny
    do ix=1,geom%nx
      pinv_x(ix,iy,iz) = sum(geom%f_p(:,iz)*x(ix,iy,:))
    end do
  end do
end do
!$omp end parallel do

! Solve Helmholz equation for each layer
pinv_q = 0.0
do iz=1,geom%nz
  call solve_helmholz_ad(geom,geom%f_d(iz),pinv_x(:,:,iz),pinv_q(:,:,iz))
end do

! Apply ff_pinv
!$omp parallel do schedule(static) private(iz,iy,ix)
do iz=1,geom%nz
  do iy=1,geom%ny
    do ix=1,geom%nx
      qtmp(ix,iy,iz) = sum(geom%f_pinv(:,iz)*pinv_q(ix,iy,:))
    end do
  end do
end do
!$omp end parallel do
q = q+qtmp

! Subtract the beta term and the heating term (=> AD of this is identity)

! Subtract the contribution from the boundaries (=> AD of this is identity)

end subroutine convert_q_to_x_ad
! ------------------------------------------------------------------------------
end module qg_convert_q_to_x_mod
