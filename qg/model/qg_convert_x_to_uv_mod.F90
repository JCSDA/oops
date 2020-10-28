! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to ageom%ny jurisdiction.

module qg_convert_x_to_uv_mod

use kinds
use qg_geom_mod

implicit none

private
public :: convert_x_to_uv,convert_x_to_uv_tl,convert_x_to_uv_ad
! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Convert streafunction to wind components
subroutine convert_x_to_uv(geom,x,x_north,x_south,u,v)

implicit none

! Passed variables
type(qg_geom),intent(in) :: geom                          !< Geometry
real(kind_real),intent(in) :: x(geom%nx,geom%ny,geom%nz)  !< Streamfunction
real(kind_real),intent(in) :: x_north(geom%nz)            !< Streamfunction on northern wall
real(kind_real),intent(in) :: x_south(geom%nz)            !< Streamfunction on southern wall
real(kind_real),intent(out) :: u(geom%nx,geom%ny,geom%nz) !< Zonal wind
real(kind_real),intent(out) :: v(geom%nx,geom%ny,geom%nz) !< Meridional wind

! Local variables
integer :: iz

! Zonal wind
!$omp parallel do schedule(static) private(iz)
do iz=1,geom%nz
  u(:,2:geom%ny,iz) = 0.5*x(:,1:geom%ny-1,iz)/geom%deltay
  u(:,1,iz) = 0.5*x_south(iz)/geom%deltay
  u(:,1:geom%ny-1,iz) = u(:,1:geom%ny-1,iz)-0.5*x(:,2:geom%ny,iz)/geom%deltay
  u(:,geom%ny,iz) = u(:,geom%ny,iz)-0.5*x_north(iz)/geom%deltay
enddo
!$omp end parallel do

! Meridional wind
!$omp parallel do schedule(static) private(iz)
do iz=1,geom%nz
   v(1:geom%nx-1,:,iz) = 0.5*x(2:geom%nx,:,iz)/geom%deltax
   v(geom%nx,:,iz) = 0.5*x(1,:,iz)/geom%deltax
   v(2:geom%nx,:,iz) = v(2:geom%nx,:,iz)-0.5*x(1:geom%nx-1,:,iz)/geom%deltax
   v(1,:,iz) = v(1,:,iz)-0.5*x(geom%nx,:,iz)/geom%deltax
enddo
!$omp end parallel do

end subroutine convert_x_to_uv
! ------------------------------------------------------------------------------
!> Convert streafunction to wind components - tangent Linear
subroutine convert_x_to_uv_tl(geom,x,u,v)

implicit none

! Passed variables
type(qg_geom),intent(in) :: geom                          !< Geometry
real(kind_real),intent(in) :: x(geom%nx,geom%ny,geom%nz)  !< Streamfunction
real(kind_real),intent(out) :: u(geom%nx,geom%ny,geom%nz) !< Zonal wind
real(kind_real),intent(out) :: v(geom%nx,geom%ny,geom%nz) !< Meridional wind

! Local variables
integer :: iz

! Zonal wind
!$omp parallel do schedule(static) private(iz)
do iz=1,geom%nz
  u(:,2:geom%ny,iz) = 0.5*x(:,1:geom%ny-1,iz)/geom%deltay
  u(:,1,iz) = 0.0
  u(:,1:geom%ny-1,iz) = u(:,1:geom%ny-1,iz)-0.5*x(:,2:geom%ny,iz)/geom%deltay
enddo
!$omp end parallel do

! Meridional wind
!$omp parallel do schedule(static) private(iz)
do iz=1,geom%nz
   v(1:geom%nx-1,:,iz) = 0.5*x(2:geom%nx,:,iz)/geom%deltax
   v(geom%nx,:,iz) = 0.5*x(1,:,iz)/geom%deltax
   v(2:geom%nx,:,iz) = v(2:geom%nx,:,iz)-0.5*x(1:geom%nx-1,:,iz)/geom%deltax
   v(1,:,iz) = v(1,:,iz)-0.5*x(geom%nx,:,iz)/geom%deltax
enddo
!$omp end parallel do

end subroutine convert_x_to_uv_tl
! ------------------------------------------------------------------------------
!> Convert streafunction to wind components - adjoint
subroutine convert_x_to_uv_ad(geom,u,v,x)

implicit none

! Passed variables
type(qg_geom),intent(in) :: geom                             !< Geometry
real(kind_real),intent(in) :: u(geom%nx,geom%ny,geom%nz)     !< Zonal wind
real(kind_real),intent(in) :: v(geom%nx,geom%ny,geom%nz)     !< Meridional wind
real(kind_real),intent(inout)  :: x(geom%nx,geom%ny,geom%nz) !< Streamfunction

! Local variables
integer :: iz

! Zonal wind
do iz=1,geom%nz
   x(:,2:geom%ny,iz) = x(:,2:geom%ny,iz)-0.5/geom%deltay*u(:,1:geom%ny-1,iz)
   x(:,1:geom%ny-1,iz) = x(:,1:geom%ny-1,iz)+0.5/geom%deltay*u(:,2:geom%ny,iz)
enddo

! Meridional wind
do iz=1,geom%nz
   x(geom%nx,:,iz) = x(geom%nx,:,iz)-0.5/geom%deltax*v(1,:,iz)
   x(1:geom%nx-1,:,iz) = x(1:geom%nx-1,:,iz)-0.5/geom%deltax*v(2:geom%nx,:,iz)
   x(1,:,iz) = x(1,:,iz)+0.5/geom%deltax*v(geom%nx,:,iz)
   x(2:geom%nx,:,iz) = x(2:geom%nx,:,iz)+0.5/geom%deltax*v(1:geom%nx-1,:,iz)
enddo

end subroutine convert_x_to_uv_ad
! ------------------------------------------------------------------------------
end module qg_convert_x_to_uv_mod
