! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

module qg_advect_q_mod

use kinds
!$ use omp_lib
use qg_constants_mod
use qg_geom_mod
use qg_interp_mod

implicit none

private
public :: advect_q,advect_q_tl,advect_q_ad
! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Advect potential vorticity
subroutine advect_q(geom,dt,u,v,q,q_north,q_south,qnew)

implicit none

! Passed variables
type(qg_geom),intent(in) :: geom                             !< Geometry
real(kind_real),intent(in)  :: dt                            !< Time step
real(kind_real),intent(in)  :: u(geom%nx,geom%ny,geom%nz)    !< Zonal wind
real(kind_real),intent(in)  :: v(geom%nx,geom%ny,geom%nz)    !< Meridional wind
real(kind_real),intent(in)  :: q(geom%nx,geom%ny,geom%nz)    !< Input potential vorticity
real(kind_real),intent(in)  :: q_north(geom%nx,geom%nz)      !< Potential vorticity on northern wall
real(kind_real),intent(in)  :: q_south(geom%nx,geom%nz)      !< Potential vorticity on southern wall
real(kind_real),intent(out) :: qnew(geom%nx,geom%ny,geom%nz) !< Output potential vorticity

! Local variables
integer :: ix,iy,iz
real(kind_real) :: x,y
real(kind_real) :: qext(geom%nx,-1:geom%ny+2,geom%nz)

! Extend q field
qext(:,1:geom%ny,:) = q
qext(:,-1,:) = q_south
qext(:,0,:) = q_south
qext(:,geom%ny+1,:) = q_north
qext(:,geom%ny+2,:) = q_north

! Advect q
!$omp parallel do schedule(static) private(iz,iy,ix,x,y)
do iz=1,geom%nz
  do iy=1,geom%ny
    do ix=1,geom%nx
      ! Find the interpolation point
      x = geom%x(ix)-u(ix,iy,iz)*dt
      y = geom%y(iy)-v(ix,iy,iz)*dt

      ! Interpolate
      call qg_interp_bicubic(geom,x,y,qext(:,:,iz),qnew(ix,iy,iz))
    enddo
  enddo
enddo
!$omp end parallel do

end subroutine advect_q
! ------------------------------------------------------------------------------
!> Advect potential vorticity - tangent linear
subroutine advect_q_tl(geom,dt,u_traj,v_traj,q_traj,q_traj_north,q_traj_south,u,v,q,qnew)

implicit none

! Passed variables
type(qg_geom),intent(in) :: geom                               !< Geometry
real(kind_real),intent(in) :: dt                               !< Time step
real(kind_real),intent(in)  :: u_traj(geom%nx,geom%ny,geom%nz) !< Zonal wind (trajectory)
real(kind_real),intent(in)  :: v_traj(geom%nx,geom%ny,geom%nz) !< Meridional wind (trajectory)
real(kind_real),intent(in)  :: q_traj(geom%nx,geom%ny,geom%nz) !< Potential vorticity (trajectory)
real(kind_real),intent(in)  :: q_traj_north(geom%nx,geom%nz)   !< Potential vorticity on northern wall (trajectory)
real(kind_real),intent(in)  :: q_traj_south(geom%nx,geom%nz)   !< Potential vorticity on southern wall (trajectory)
real(kind_real),intent(in)  :: u(geom%nx,geom%ny,geom%nz)      !< Zonal wind (perturbation)
real(kind_real),intent(in)  :: v(geom%nx,geom%ny,geom%nz)      !< Meridional wind (perturbation)
real(kind_real),intent(in)  :: q(geom%nx,geom%ny,geom%nz)      !< Input potential vorticity (perturbation)
real(kind_real),intent(out) :: qnew(geom%nx,geom%ny,geom%nz)   !< Output potential vorticity (perturbation)

! Local variables
integer :: ix,iy,iz
real(kind_real) :: x,y,dx,dy
real(kind_real) :: qext_traj(geom%nx,-1:geom%ny+2,geom%nz),qext(geom%nx,-1:geom%ny+2,geom%nz)

! Extend q field (trajectory)
qext_traj(:,1:geom%ny,:) = q_traj
qext_traj(:,-1,:) = q_traj_south
qext_traj(:,0,:) = q_traj_south
qext_traj(:,geom%ny+1,:) = q_traj_north
qext_traj(:,geom%ny+2,:) = q_traj_north

! Extend q field (perturbation)
qext(:,1:geom%ny,:) = q
qext(:,-1,:) = 0.0
qext(:,0,:) = 0.0
qext(:,geom%ny+1,:) = 0.0
qext(:,geom%ny+2,:) = 0.0

! Advect q
!$omp parallel do schedule(static) private(iz,iy,ix,x,y,dx,dy)
do iz=1,geom%nz
  do iy=1,geom%ny
    do ix=1,geom%nx
      ! Find the interpolation point
      x = geom%x(ix)-u_traj(ix,iy,iz)*dt
      y = geom%y(iy)-v_traj(ix,iy,iz)*dt

      ! Find interpolation point perturbation
      dx = -u(ix,iy,iz)*dt
      dy = -v(ix,iy,iz)*dt

      ! Interpolate
      call qg_interp_bicubic_tl(geom,x,y,qext_traj(:,:,iz),dx,dy,qext(:,:,iz),qnew(ix,iy,iz)) 
    enddo
  enddo
enddo
!$omp end parallel do

end subroutine advect_q_tl
! ------------------------------------------------------------------------------
!> Advect potential vorticity - adjoint
subroutine advect_q_ad(geom,dt,u_traj,v_traj,q_traj,q_traj_north,q_traj_south,qnew,u,v,q)

implicit none

! Passed variables
type(qg_geom),intent(in) :: geom                              !< Geometry
real(kind_real),intent(in) :: dt                              !< Time step
real(kind_real),intent(in) :: u_traj(geom%nx,geom%ny,geom%nz) !< Zonal wind (trajectory)
real(kind_real),intent(in) :: v_traj(geom%nx,geom%ny,geom%nz) !< Meridional wind (trajectory)
real(kind_real),intent(in) :: q_traj(geom%nx,geom%ny,geom%nz) !< Potential vorticity (trajectory)
real(kind_real),intent(in) :: q_traj_north(geom%nx,geom%nz)   !< Potential vorticity on northern wall (trajectory)
real(kind_real),intent(in) :: q_traj_south(geom%nx,geom%nz)   !< Potential vorticity on southern wall (trajectory)
real(kind_real),intent(in) :: qnew(geom%nx,geom%ny,geom%nz)   !< Output potential vorticity (perturbation)
real(kind_real),intent(inout) :: u(geom%nx,geom%ny,geom%nz)   !< Zonal wind (perturbation)
real(kind_real),intent(inout) :: v(geom%nx,geom%ny,geom%nz)   !< Meridional wind (perturbation)
real(kind_real),intent(inout) :: q(geom%nx,geom%ny,geom%nz)   !< Input potential vorticity (perturbation)

! Local variables
integer :: ix,iy,iz
real(kind_real) :: x,y,dx,dy
real(kind_real) :: qext_traj(geom%nx,-1:geom%ny+2,geom%nz),qext(geom%nx,-1:geom%ny+2,geom%nz)

! Extend q field (trajectory)
qext_traj(:,1:geom%ny,:) = q_traj
qext_traj(:,-1,:) = q_traj_south
qext_traj(:,0,:) = q_traj_south
qext_traj(:,geom%ny+1,:) = q_traj_north
qext_traj(:,geom%ny+2,:) = q_traj_north

! Initialization
qext = 0.0

! Advect q
do iz=geom%nz,1,-1
  do iy=geom%ny,1,-1
    do ix=geom%nx,1,-1
      ! Find the interpolation point
      x = geom%x(ix)-u_traj(ix,iy,iz)*dt
      y = geom%y(iy)-v_traj(ix,iy,iz)*dt

      ! Initialization
      dx = 0.0
      dy = 0.0

      ! Interpolate,adjoint
      call qg_interp_bicubic_ad(geom,x,y,qext_traj(:,:,iz),qnew(ix,iy,iz),dx,dy,qext(:,:,iz))

      ! Find interpolation point perturbation,adjoint
      u(ix,iy,iz) = u(ix,iy,iz)-dx*dt
      v(ix,iy,iz) = v(ix,iy,iz)-dy*dt
    enddo
  enddo
enddo

! Extend q field (perturbation)
q = q+qext(:,1:geom%ny,:)

end subroutine advect_q_ad
! ------------------------------------------------------------------------------
end module qg_advect_q_mod
