! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

module differential_solver_mod

use iso_c_binding
use fckit_log_module, only: fckit_log
use fft_mod
use kinds
!$ use omp_lib
use qg_constants_mod
use qg_geom_mod

implicit none

private
public :: solve_helmholz,solve_helmholz_ad,laplacian_2d,laplacian_2d_ad
! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Solve a Helmholz equation
subroutine solve_helmholz(geom,c,b,x)

implicit none

! Passed variables
type(qg_geom),intent(in) :: geom                  !< Geometry
real(kind_real),intent(in) :: c                   !< Coefficient in the linear operator
real(kind_real),intent(in) :: b(geom%nx,geom%ny)  !< Right hand side
real(kind_real),intent(out) :: x(geom%nx,geom%ny) !< Solution

! Local variables
integer :: kx,iri,ix,iy
real(kind_real) :: v(geom%ny),w(geom%ny),z(geom%ny)
real(kind_real) :: bext(geom%nx+2,geom%ny),xext(geom%nx+2,geom%ny)
real(kind_real) :: am,bm
real(kind_real) :: tmp_in(geom%nx,geom%ny),tmp_out(geom%nx,geom%ny)
character(len=2014) :: record

! Transform
do iy=1,geom%ny
  call fft_fwd(geom%nx,b(:,iy),bext(:,iy))
enddo

! Solve the tri-diagonal systems
do kx=0,geom%nx/2
  ! am and bm parameters
  am = c+2.0*(cos(2.0*real(kx,kind_real)*pi/real(geom%nx,kind_real))-1.0)/geom%deltax**2-2.0/geom%deltay**2
  bm = 1.0/geom%deltay**2

  do iri=1,2
    ! v and w parameters
    v(1) = bm/am
    do iy=2,geom%ny
      w(iy) = am-bm*v(iy-1)
      v(iy) = bm/w(iy)
    enddo

    ! bext to z
    z(1) = bext(2*kx+iri,1)/am
    do iy=2,geom%ny
      z(iy) = (bext(2*kx+iri,iy)-bm*z(iy-1))/w(iy)
    enddo

    ! z to xext
    xext(2*kx+iri,geom%ny) = z(geom%ny)
    do iy=geom%ny-1,1,-1
       xext(2*kx+iri,iy) = z(iy)-v(iy)*xext(2*kx+iri,iy+1) 
    enddo
  enddo
enddo

! Transform back
do iy=1,geom%ny
  call fft_inv(geom%nx,xext(:,iy),x(:,iy))
enddo

end subroutine solve_helmholz
! ------------------------------------------------------------------------------
!> Solve a Helmholz equation - adjoint
subroutine solve_helmholz_ad(geom,c,x,b)

implicit none

! Passed variables
type(qg_geom),intent(in) :: geom                    !< Geometry
real(kind_real),intent(in) :: c                     !< Coefficient in the linear operator
real(kind_real),intent(in) :: x(geom%nx,geom%ny)    !< Solution
real(kind_real),intent(inout) :: b(geom%nx,geom%ny) !< Right hand side

! Local variables
real(kind_real) :: v(geom%ny),w(geom%ny),z(geom%ny)
real(kind_real) :: bext(geom%nx+2,geom%ny),xext(geom%nx+2,geom%ny)
real(kind_real) :: btmp(geom%nx,geom%ny)
real(kind_real) :: am,bm
integer :: kx,iri,iy

! Initialization
bext = 0.0

! Transform back
do iy=geom%ny,1,-1
 call fft_fwd(geom%nx,x(:,iy),xext(:,iy))
enddo

! Solve the tri-diagonal systems
do kx=geom%nx/2,0,-1
  ! am and bm parameters
  am = c+2.0*(cos(2.0*real(kx,kind_real)*pi/real(geom%nx,kind_real))-1.0)/geom%deltax**2-2.0/geom%deltay**2
  bm = 1.0/geom%deltay**2
  
  do iri=2,1,-1
    ! v and w parameters
    v(1) = bm/am
    do iy=2,geom%ny
      w(iy) = am-bm*v(iy-1)
      v(iy) = bm/w(iy)
    enddo
  
    ! Initialization
    z = 0.0
  
    ! z to xext
    do iy=1,geom%ny-1
      xext(2*kx+iri,iy+1) = xext(2*kx+iri,iy+1)-v(iy)*xext(2*kx+iri,iy)
      z(iy) = z(iy)+xext(2*kx+iri,iy)
    enddo
    z(geom%ny) = z(geom%ny)+xext(2*kx+iri,geom%ny)
  
    ! bext to z
    do iy=geom%ny,2,-1
      bext(2*kx+iri,iy) = bext(2*kx+iri,iy)+z(iy)/w(iy)
      z(iy-1) = z(iy-1)-bm*z(iy)/w(iy)
    enddo
    bext(2*kx+iri,1) = bext(2*kx+iri,1)+z(1)/am
  enddo
enddo

! Transform
do iy=geom%ny,1,-1
  call fft_inv(geom%nx,bext(:,iy),btmp(:,iy))
enddo
b = b+btmp

end subroutine solve_helmholz_ad
! ------------------------------------------------------------------------------
!> Horizontal Laplacian operator
subroutine laplacian_2d(geom,x,del2x)

implicit none

! Passed variables
type(qg_geom),intent(in) :: geom                              !< Geometry
real(kind_real),intent(in) :: x(geom%nx,geom%ny,geom%nz)      !< Streamfunction
real(kind_real),intent(out) :: del2x(geom%nx,geom%ny,geom%nz) !< Result of applying Laplacian to x

! Local variables
integer :: iz

!$omp parallel do schedule(static) private(iz)
do iz=1,geom%nz
  ! 5-point laplacian
  del2x(:,:,iz) = -2.0*x(:,:,iz)*(1.0/geom%deltax**2+1.0/geom%deltay**2)
  del2x(1:geom%nx-1,:,iz) = del2x(1:geom%nx-1,:,iz)+x(2:geom%nx,:,iz)/geom%deltax**2
  del2x(geom%nx,:,iz) = del2x(geom%nx,:,iz)+x(1,:,iz)/geom%deltax**2
  del2x(2:geom%nx,:,iz) = del2x(2:geom%nx,:,iz)+x(1:geom%nx-1,:,iz)/geom%deltax**2
  del2x(1,:,iz) = del2x(1,:,iz)+x(geom%nx,:,iz)/geom%deltax**2
  del2x(:,1:geom%ny-1,iz) = del2x(:,1:geom%ny-1,iz)+x(:,2:geom%ny,iz)/geom%deltay**2
  del2x(:,2:geom%ny,iz) = del2x(:,2:geom%ny,iz)+x(:,1:geom%ny-1,iz)/geom%deltay**2
enddo
!$omp end parallel do

end subroutine laplacian_2d
! ------------------------------------------------------------------------------
!> Horizontal Laplacian operator - adjoint
subroutine laplacian_2d_ad(geom,del2x,x)

implicit none

! Passed variables
type(qg_geom),intent(in) :: geom                             !< Geometry
real(kind_real),intent(in) :: del2x(geom%nx,geom%ny,geom%nz) !< Result of applying Laplacian to x
real(kind_real),intent(inout) :: x(geom%nx,geom%ny,geom%nz)  !< Streamfunction

! Local variables
integer :: iz

do iz=1,geom%nz
  ! 5-point laplacian
  x(:,1:geom%ny-1,iz) = x(:,1:geom%ny-1,iz)+del2x(:,2:geom%ny,iz)/geom%deltay**2
  x(:,2:geom%ny,iz) = x(:,2:geom%ny,iz)+del2x(:,1:geom%ny-1,iz)/geom%deltay**2
  x(geom%nx,:,iz) = x(geom%nx,:,iz)+del2x(1,:,iz)/geom%deltax**2
  x(1:geom%nx-1,:,iz) = x(1:geom%nx-1,:,iz)+del2x(2:geom%nx,:,iz)/geom%deltax**2
  x(1,:,iz) = x(1,:,iz)+del2x(geom%nx,:,iz)/geom%deltax**2
  x(2:geom%nx,:,iz) = x(2:geom%nx,:,iz)+del2x(1:geom%nx-1,:,iz)/geom%deltax**2
  x(:,:,iz) = x(:,:,iz)-2.0*del2x(:,:,iz)*(1.0/geom%deltax**2+1.0/geom%deltay**2)
enddo

end subroutine laplacian_2d_ad
! ------------------------------------------------------------------------------
end module differential_solver_mod
