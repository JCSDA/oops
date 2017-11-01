! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

!> Solve a Helmholz equation.

!> This routine solves the Helmholz equation:
!! \f$ \nabla^2 \psi + c \psi = b \f$,
!! where \f$ c \f$ is a constant, and where \f$ \nabla^2 \f$ is the
!! two-dimensional, 5-point finite-difference Laplacian.
!!
!! The solution method is to apply an FFT in the zonal direction. This
!! reduces the problem to a set of un-coupled tri-diagonal systems that
!! are solved with the standard back-substitution algorithm.

subroutine solve_helmholz (x,b,c,nx,ny,deltax,deltay)
use fft_mod
use kinds
implicit none

integer, intent(in) :: nx     !< Zonal grid dimension
integer, intent(in) :: ny     !< Meridional grid dimension
real(kind=kind_real), intent(out) :: x(nx,ny) !< Solution
real(kind=kind_real), intent(in) :: b(nx,ny)  !< Right hand side
real(kind=kind_real), intent(in) :: c         !< Coefficient in the linear operator
real(kind=kind_real), intent(in) :: deltax    !< Zonal grid spacing (non-dimensional)
real(kind=kind_real), intent(in) :: deltay    !< Meridional grid spacing (non-dimensional)

real(kind=kind_real) :: v(ny), z(ny), bext(nx+2,ny), xext(nx+2,ny)
real(kind=kind_real) :: pi, am, bm, w
integer :: k, i, iri, jj

x(:,:) = 0.0_kind_real

if (maxval(abs(b))>0.0_kind_real) then

  !--- transform
  do jj=1,ny
    call fft_fwd(nx,b(:,jj),bext(:,jj))
  enddo

  pi = 4.0_kind_real*atan(1.0_kind_real)

  !--- loop over wavenumber

  do k=0,nx/2

    !--- solve the tri-diagonal systems

    am = c + 2.0_kind_real*( cos( 2.0_kind_real*real(k,kind_real)*pi &
                                 /real(nx,kind_real)) &
                            -1.0_kind_real)/(deltax*deltax) &
           - 2.0_kind_real/(deltay*deltay)

    bm = 1.0_kind_real/(deltay*deltay)

    do iri=1,2
      v(1) = bm/am
      z(1) = bext(2*k+iri,1)/am

      do i=2,ny
        w    = am-bm*v(i-1)
        v(i) = bm/w
        z(i) = (bext(2*k+iri,i)-bm*z(i-1))/w
      enddo

      xext(2*k+iri,ny) = z(ny)
      do i=ny-1,1,-1
         xext(2*k+iri,i) = z(i) - v(i)*xext(2*k+iri,i+1) 
      enddo
    enddo
  enddo

  !--- transform back
  do jj=1,ny
    call fft_inv(nx,xext(:,jj),x(:,jj))
  enddo

endif

end subroutine solve_helmholz
