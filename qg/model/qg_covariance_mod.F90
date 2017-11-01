! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

!> Structure holding configuration variables for the 3d error
!! covariance matrices of the QG analysis.

module qg_covariance_mod

use kinds
implicit none

!> Fortran derived type to hold configuration data for the QG background/model covariance
type :: qg_3d_covar_config
  integer :: nx !< Zonal grid dimension
  integer :: ny !< Meridional grid dimension
  real(kind=kind_real)    :: sigma        !< Standard deviation
  real(kind=kind_real)    :: vert_corr    !< Vertical correlation between levels
  real(kind=kind_real), allocatable :: sqrt_merid(:,:) !< sqrt(meridional correlation matrix)
  real(kind=kind_real), allocatable :: sqrt_inv_merid(:,:) !< Its inverse
  real(kind=kind_real), allocatable :: sqrt_zonal(:) !< Spectral weights for sqrt(zonal corr)

end type qg_3d_covar_config

#define LISTED_TYPE qg_3d_covar_config

!> Linked list interface - defines registry_t type
#include "linkedList_i.f"

!> Global registry
type(registry_t) :: qg_3d_cov_registry

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Linked list implementation
#include "linkedList_c.f"
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

!> Setup for the QG model's 3d error covariance matrices (B and Q_i)

!> This routine queries the configuration for the parameters that define the
!! covariance matrix, and stores the relevant values in the
!! error covariance structure.

subroutine qg_3d_covar_setup(c_model, geom, config)

use qg_constants
use qg_geom_mod
use iso_c_binding
use config_mod
use fft_mod
use kinds
use fckit_log_module, only : fckit_log

implicit none
type(c_ptr), intent(in)   :: c_model  !< The configuration
type(qg_geom), intent(in) :: geom     !< Geometry
type(qg_3d_covar_config), intent(inout) :: config !< The covariance structure

integer :: ifax(13)
real(kind=kind_real), allocatable :: trigs(:), struct_fn(:), evals(:), &
                      work(:), evects(:,:), revals(:)
character(len=160) :: record

real(kind=kind_real) :: corr_length_scale, dy, z, dx, condition_number, threshold
integer :: i, j, k, info

config%nx         = geom%nx
config%ny         = geom%ny
config%sigma      = config_get_real(c_model,"standard_deviation")
config%vert_corr  = config_get_real(c_model,"vertical_correlation")
corr_length_scale = config_get_real(c_model,"horizontal_length_scale")
condition_number  = config_get_real(c_model,"maximum_condition_number")

dx = domain_zonal/real(config%nx,kind_real)
dy = domain_meridional/real(1+config%ny,kind_real)

if (mod(config%nx,2)/=0) then
  write(record,*) "c_qg_3d_covar_setup: number of zonal gridpoints nx=",config%nx
  call fckit_log%error(record)
  call abor1_ftn("c_qg_3d_covar_setup: odd number of zonal grid points")
endif

!--- Generate the lower triangle of the meridional correlation matrix

allocate(evects(config%ny,config%ny))

do i=1,config%ny
  z = (sqrt(0.5_kind_real)*dy/corr_length_scale)*real(i-1,kind_real)
  evects(i,1) = exp(-z*z)
enddo

do j=2,config%ny
  do i=j,config%ny
    evects(i,j) = evects(i-j+1,1)
  enddo
enddo

!--- Calculate the eigen-decomposition (NB: assuming 8-byte reals)

allocate(evals(config%ny))
allocate(work((config%ny+3)*config%ny))

call DSYEV ('V','L',config%ny,evects,config%ny,evals,work,size(work),info)

if (info/=0) then
  write(record,*) "c_qg_3d_covar_setup: DSYEV returns info=",info
  call fckit_log%error(record)
  call abor1_ftn ("c_qg_3d_covar_setup: DSYEV failed to sqrt matrix")
endif

deallocate(work)

!--- Calculate the lower triangle of the symmetric square-root and its inverse
!--- (NB: assuming 8-byte reals)

allocate(revals(config%ny))

threshold = maxval(evals(:))/condition_number

do j=1,config%ny
  evals(j) = sqrt(max(threshold,evals(j)))
  revals(j) = 1.0_kind_real/evals(j)
enddo

allocate(config%sqrt_merid(config%ny,config%ny))
allocate(config%sqrt_inv_merid(config%ny,config%ny))

do j=1,config%ny
  do i=j,config%ny
    config%sqrt_merid(i,j) = 0.0_kind_real
    config%sqrt_inv_merid(i,j) = 0.0_kind_real

    do k=1,config%ny
      config%sqrt_merid(i,j) = config%sqrt_merid(i,j) &
                 & + evects(i,k)*evects(j,k)*evals(k)
      config%sqrt_inv_merid(i,j) = config%sqrt_inv_merid(i,j) &
                 & + evects(i,k)*evects(j,k)*revals(k)
    enddo
  enddo
enddo

deallocate(revals)
deallocate(evects)
deallocate(evals)

!--- Calculate spectral weights for zonal correlations:
!--- First we construct the structure function in grid space and FFT it

allocate(struct_fn(config%nx))
do j=1,config%nx
  z = (sqrt(0.5_kind_real)*dx/corr_length_scale) &
   & *real(min(j-1,1+config%nx-j),kind_real)
  struct_fn(j) = exp(-z*z)
enddo

allocate(work(config%nx+2))

call fft_fwd(config%nx,struct_fn,work)
work(:) = work(:) / real(config%nx, kind_real)

!--- We are after sqrt(B) or sqrt(Q), so square-root the coefficents.
!--- The factor N ensures we get the struture function back if we apply
!--- B or Q to a unit delta function.

allocate(config%sqrt_zonal(0:config%nx/2))

threshold = 0.0_kind_real
do i=0,config%nx/2
  threshold = max(threshold,work(1+2*i))
enddo
threshold = threshold/condition_number

do i=0,config%nx/2
  config%sqrt_zonal(i) = sqrt(max(threshold,real(config%nx,kind_real)*work(1+2*i)))
enddo

deallocate(struct_fn)
deallocate(work)

return
end subroutine qg_3d_covar_setup

! ------------------------------------------------------------------------------

!> Delete for the QG model's 3d error covariance matrices

subroutine qg_3d_covar_delete(c_key_conf)

use iso_c_binding

implicit none
integer(c_int), intent(inout) :: c_key_conf !< The model covariance structure

type(qg_3d_covar_config), pointer :: conf !< covar structure

call qg_3d_cov_registry%get(c_key_conf, conf)

deallocate(conf%sqrt_zonal)
deallocate(conf%sqrt_merid)
deallocate(conf%sqrt_inv_merid)
call qg_3d_cov_registry%remove(c_key_conf)

end subroutine qg_3d_covar_delete

! ------------------------------------------------------------------------------

!> Multiply streamfunction by inverse(sqrt(C)), where C is 3d covariance matrix

subroutine qg_3d_covar_sqrt_inv_mult(kx,ky,xctl,xincr,config)
use iso_c_binding
use fft_mod
use kinds
use qg_fields

implicit none
integer(c_int), intent(in)    :: kx            !< Zonal grid dimension
integer(c_int), intent(in)    :: ky            !< Meridional grid dimension
real(c_double), intent(inout) :: xctl(kx,ky,2) !< inv(sqrt(C)) times psi
type(qg_field), intent(in)    :: xincr         !< Streamfunction: psi
type(qg_3d_covar_config), intent(in) :: config !< covar config structure

real(kind=kind_real) :: zfour(kx+2), work(ky)
integer :: i, j, k, iri, m
real(kind=kind_real) :: zc, zero, one

!--- multiply by standard deviation

zc = 1.0_kind_real/config%sigma
do k=1,2
  do j=1,ky
    do i=1,kx
      xctl(i,j,k) = zc * xincr%x(i,j,k)
    enddo
  enddo
enddo

!--- multiply by inverse square-root of zonal correlation matrix

do k=1,2
  do j=1,ky
    call fft_fwd(kx,xctl(:,j,k),zfour)
    do m=0,kx/2
      do iri=1,2
        zfour(2*m+iri) = zfour(2*m+iri) / config%sqrt_zonal(m)
      enddo
    enddo
    call fft_inv(kx,zfour,xctl(:,j,k))
  enddo
enddo

!--- multiply by inverse square-root of meridional correlation matrix

zero = 0.0_kind_real
one = 1.0_kind_real
do k=1,2
  do i=1,kx
    call DSYMV('L',ky,one,config%sqrt_inv_merid,ky,xctl(i,1,k),kx,zero,work,1)
    do j=1,ky
      xctl(i,j,k) = work(j)
    enddo
  enddo
enddo

!--- multiply by inverse symmetric square-root of vertical correlation matrix

zc = 1.0_kind_real / sqrt(1.0_kind_real-config%vert_corr*config%vert_corr)
do j=1,ky
  do i=1,kx
    xctl(i,j,2) = zc * (xctl(i,j,2) - config%vert_corr * xctl(i,j,1))
  enddo
enddo

end subroutine qg_3d_covar_sqrt_inv_mult

! ------------------------------------------------------------------------------

!> Multiply streamfunction by inverse(sqrt(C)) - Adjoint

subroutine qg_3d_covar_sqrt_inv_mult_ad(kx,ky,xctl,xincr,config)
use iso_c_binding
use fft_mod
use kinds
use qg_fields

implicit none
integer(c_int), intent(in) :: kx               !< Zonal grid dimension
integer(c_int), intent(in) :: ky               !< Meridional grid dimension
type(qg_field), intent(inout) :: xincr         !< sqrt(C) times streamfunction
real(c_double), intent(in) :: xctl(kx,ky,2)    !< Streamfunction
type(qg_3d_covar_config), intent(in) :: config !< covariance config structure

real(kind=kind_real), allocatable :: xout(:,:,:)
real(kind=kind_real) :: zfour(kx+2), work(ky)
integer :: i, j, k, iri, m
real(kind=kind_real) :: zc, zero, one

!--- adoint of multiplication by inverse symmetric square-root of vertical
!--- correlation matrix

allocate(xout(kx,ky,2))

zc = sqrt(1.0_kind_real-config%vert_corr*config%vert_corr)
do j=1,ky
  do i=1,kx
    xout(i,j,1) = xctl(i,j,1) -config%vert_corr * xctl(i,j,2) &
                    & *(1.0_kind_real/zc)
    xout(i,j,2) = xctl(i,j,2)*(1.0_kind_real/zc)
  enddo
enddo

!--- adjoint multiplication by inverse sqrt of meridional correlation matrix

zero = 0.0_kind_real
one = 1.0_kind_real
do k=1,2
  do i=1,kx
    call DSYMV('L',ky,one,config%sqrt_inv_merid,ky,xout(i,1,k), &
        &      kx,zero,work,1)
    do j=1,ky
      xout(i,j,k) = work(j)
    enddo
  enddo
enddo

!--- multiply by inverse square-root of zonal correlation matrix

do k=1,2
  do j=1,ky
    call fft_fwd(kx,xout(:,j,k),zfour)
    do m=0,kx/2
      do iri=1,2
        zfour(2*m+iri) = zfour(2*m+iri) / config%sqrt_zonal(m)
      enddo
    enddo
    call fft_inv(kx,zfour,xout(:,j,k))
  enddo
enddo

!--- adjoint of multiplication by standard deviation

do k=1,2
  do j=1,ky
    do i=1,kx
      xincr%x(i,j,k) = xincr%x(i,j,k) + xout(i,j,k) &
                    & *(1.0_kind_real/config%sigma)
    enddo
  enddo
enddo

deallocate(xout)

end subroutine qg_3d_covar_sqrt_inv_mult_ad

! ------------------------------------------------------------------------------

!> Multiply streamfunction by sqrt(C), where C is a 3d covariance matrix

subroutine qg_3d_covar_sqrt_mult(kx,ky,xincr,xctl,config)
use iso_c_binding
use fft_mod
use kinds
use qg_fields

implicit none
integer(c_int), intent(in) :: kx               !< Zonal grid dimension
integer(c_int), intent(in) :: ky               !< Meridional grid dimension
type(qg_field), intent(inout) :: xincr         !< sqrt(C) times streamfunction
real(c_double), intent(in) :: xctl(kx,ky,2)    !< Streamfunction
type(qg_3d_covar_config), intent(in) :: config !< covariance config structure

integer :: i, j, k, iri, m
real(kind=kind_real) :: zc, zero, one
real(kind=kind_real) :: zfour(kx+2), work(ky)

!--- multiply by symmetric square-root of vertical correlation matrix

zc = sqrt(1.0_kind_real-config%vert_corr*config%vert_corr)
do j=1,ky
  do i=1,kx
    xincr%x(i,j,1) = xctl(i,j,1)
    xincr%x(i,j,2) = config%vert_corr * xctl(i,j,1) + zc * xctl(i,j,2)
  enddo
enddo

!--- multiply by square-root of meridional correlation matrix

zero = 0.0_kind_real
one = 1.0_kind_real
do k=1,2
  do i=1,kx
    call DSYMV('L',ky,one,config%sqrt_merid,ky,xincr%x(i,1,k),kx,zero,work,1)
    do j=1,ky
      xincr%x(i,j,k) = work(j)
    enddo
  enddo
enddo

!--- multiply by square-root of zonal correlation matrix

do k=1,2
  do j=1,ky
    call fft_fwd(kx,xincr%x(:,j,k),zfour)
    do m=0,kx/2
      do iri=1,2
        zfour(2*m+iri) = zfour(2*m+iri) * config%sqrt_zonal(m)
      enddo
    enddo
    call fft_inv(kx,zfour,xincr%x(:,j,k))
  enddo
enddo

!--- multiply by standard deviation

do k=1,2
  do j=1,ky
    do i=1,kx
      xincr%x(i,j,k) = xincr%x(i,j,k) * config%sigma
    enddo
  enddo
enddo

end subroutine qg_3d_covar_sqrt_mult

! ------------------------------------------------------------------------------

!> Multiply streamfunction by sqrt(C) - Adjoint

subroutine qg_3d_covar_sqrt_mult_ad(kx,ky,xincr,xctl,config)
use iso_c_binding
use fft_mod
use kinds
use qg_fields

implicit none
integer(c_int), intent(in)    :: kx            !< Zonal grid spacing
integer(c_int), intent(in)    :: ky            !< Meridional grid spacing
real(c_double), intent(inout) :: xctl(kx,ky,2) !< Result
type(qg_field), intent(in)    :: xincr         !< Streamfunction: psi
type(qg_3d_covar_config), intent(in) :: config !< covar config structure

real(kind=kind_real), allocatable :: xout(:,:,:)
integer :: i, j, k, iri, m
real(kind=kind_real) :: zc, zero, one
real(kind=kind_real) :: zgrid(kx), zfour(kx+2), work(ky)

!--- adjoint of multiplication by standard deviation

allocate(xout(kx,ky,2))

do k=1,2
  do j=1,ky
    zgrid(:) = xincr%x(:,j,k) * config%sigma
    call fft_fwd(kx,zgrid,zfour)
    do m=0,kx/2
      do iri=1,2
        zfour(2*m+iri) = zfour(2*m+iri) * config%sqrt_zonal(m)
      enddo
    enddo
    call fft_inv(kx,zfour,xout(:,j,k))
  enddo
enddo

!--- adjoint of multiplication by square-root of meridional correlation matrix

zero = 0.0_kind_real
one = 1.0_kind_real
do k=1,2
  do i=1,kx
    call DSYMV('L',ky,ONE,CONFIg%sqrt_merid,ky,xout(i,1,k),kx,zero,work,1)
    do j=1,ky
      xout(i,j,k) = work(j)
    enddo
  enddo
enddo

!--- adjoint multiplication by symmetric square-root of vert correlation matrix

zc = sqrt(1.0_kind_real-config%vert_corr*config%vert_corr)
do j=1,ky
  do i=1,kx
    xctl(i,j,1) = xctl(i,j,1) + xout(i,j,1) &
              & +config%vert_corr * xout(i,j,2)
    xctl(i,j,2) = xctl(i,j,2) + xout(i,j,2) * zc
  enddo
enddo

deallocate(xout)

end subroutine qg_3d_covar_sqrt_mult_ad

! ------------------------------------------------------------------------------

end module qg_covariance_mod
