! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

module qg_error_covariance_mod

use fckit_configuration_module, only: fckit_configuration
use fckit_log_module, only: fckit_log
use fft_mod
use iso_c_binding
use kinds
!$ use omp_lib
use oops_variables_mod
use qg_constants_mod
use qg_fields_mod
use qg_geom_mod
use random_mod

implicit none

private
public :: qg_error_covariance_config
public :: qg_error_covariance_registry
public :: qg_error_covariance_setup,qg_error_covariance_delete,qg_error_covariance_mult, &
        & qg_error_covariance_randomize
! ------------------------------------------------------------------------------
type :: qg_error_covariance_config
  integer :: nx                                      !< Number of points in the zonal direction
  integer :: ny                                      !< Number of points in the meridional direction
  integer :: nz                                      !< Number of vertical levels
  real(kind_real) :: sigma                           !< Standard deviation
  real(kind_real),allocatable :: sqrt_zonal(:)       !< Spectral weights for the spectral of the zonal correlation matrix
  real(kind_real),allocatable :: sqrt_merid(:,:)     !< Square-root of the meridional correlation matrix
  real(kind_real),allocatable :: sqrt_vert(:,:)      !< Square-root of the meridional correlation matrix
  real(kind_real),allocatable :: norm(:,:)           !< Normalization factor
  integer :: seed                                    !< Randomization seed
end type qg_error_covariance_config

real(kind_real),parameter :: eps_ad = 1.0e-10 !< Epsilon value for adjoint tests

#define LISTED_TYPE qg_error_covariance_config

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry
type(registry_t) :: qg_error_covariance_registry
! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
! Public
! ------------------------------------------------------------------------------
!> Linked list implementation
#include "oops/util/linkedList_c.f"
! ------------------------------------------------------------------------------
!> Setup error covariance matrix
subroutine qg_error_covariance_setup(self,f_conf,geom)

implicit none

! Passed variables
type(qg_error_covariance_config),intent(inout) :: self !< Error covariance configuration
type(fckit_configuration),intent(in) :: f_conf         !< FCKIT configuration
type(qg_geom),intent(in) :: geom                       !< Geometry

! Local variables
integer :: ix,iy,jy,ky,iz,jz,kz,info
real(kind_real) :: horizontal_length_scale,vertical_length_scale
real(kind_real) :: distx,condition_number,threshold
real(kind_real),allocatable :: struct_fn(:),workx(:)
real(kind_real),allocatable :: evalsy(:),worky(:),evectsy(:,:),revalsy(:)
real(kind_real),allocatable :: evalsz(:),workz(:),evectsz(:,:),revalsz(:)
real(kind_real),allocatable :: norm(:,:)
character(len=160) :: record
type(qg_fields) :: fld_in,fld_out
type(oops_variables) :: vars

! Get parameters
call f_conf%get_or_die("standard_deviation",self%sigma)
call f_conf%get_or_die("horizontal_length_scale",horizontal_length_scale)
call f_conf%get_or_die("vertical_length_scale",vertical_length_scale)
call f_conf%get_or_die("maximum_condition_number",condition_number)
if (f_conf%has("randomization_seed")) then
   call f_conf%get_or_die("randomization_seed",self%seed)
else
   self%seed = rseed
end if

! Check nx
if (mod(geom%nx,2)/=0) then
  write(record,*) 'qg_error_covariance_setup: number of zonal gridpoints nx=',geom%nx
  call fckit_log%error(record)
  call abor1_ftn('qg_error_covariance_setup: odd number of zonal grid points')
endif

! Copy grid size
self%nx = geom%nx
self%ny = geom%ny
self%nz = geom%nz

! Allocation
allocate(self%sqrt_merid(geom%ny,geom%ny))
allocate(self%sqrt_vert(geom%nz,geom%nz))
allocate(self%norm(geom%ny,geom%nz))
allocate(struct_fn(geom%nx))
allocate(workx(geom%nx+2))
allocate(self%sqrt_zonal(0:geom%nx/2))
allocate(evectsy(geom%ny,geom%ny))
allocate(evalsy(geom%ny))
allocate(worky((geom%ny+3)*geom%ny))
allocate(revalsy(geom%ny))
allocate(evectsz(geom%nz,geom%nz))
allocate(evalsz(geom%nz))
allocate(workz((geom%nz+3)*geom%nz))
allocate(revalsz(geom%nz))
allocate(norm(geom%ny,geom%nz))

! Calculate spectral weights for zonal correlations:
! First we construct the structure function in grid space and FFT it
do ix=1,geom%nx
  distx = geom%x(ix)-geom%x(1)
  if (distx>0.5*domain_zonal) distx = domain_zonal-distx
  struct_fn(ix) = exp(-0.5*(distx/horizontal_length_scale)**2)
enddo
call fft_fwd(geom%nx,struct_fn,workx)
workx = workx/real(geom%nx,kind_real)

! We are after sqrt(B) or sqrt(Q),so square-root the coefficents.
! The factor N ensures we get the struture function back if we apply
! B or Q to a unit delta function.
threshold = 0.0
do ix=0,geom%nx/2
  threshold = max(threshold,workx(1+2*ix))
enddo
threshold = threshold/condition_number
do ix=0,geom%nx/2
  self%sqrt_zonal(ix) = sqrt(max(threshold,real(geom%nx,kind_real)*workx(1+2*ix)))
enddo

! Generate the lower triangle of the meridional correlation matrix
do jy=1,geom%ny
  do iy=jy,geom%ny
    evectsy(iy,jy) = exp(-0.5*((geom%y(iy)-geom%y(jy))/horizontal_length_scale)**2)
    evectsy(jy,iy) = evectsy(iy,jy)
  enddo
enddo

! Calculate the eigen-decomposition
call dsyev('V','L',geom%ny,evectsy,geom%ny,evalsy,worky,size(worky),info)
if (info/=0) then
  write(record,*) 'qg_error_covariance_setup: dsyev returns info=',info
  call fckit_log%error(record)
  call abor1_ftn ('qg_error_covariance_setup: dsyev failed to sqrt matrix')
endif

! Calculate the lower triangle of the symmetric square-root
threshold = maxval(evalsy)/condition_number
do iy=1,geom%ny
  evalsy(iy) = sqrt(max(threshold,evalsy(iy)))
  revalsy(iy) = 1.0/evalsy(iy)
enddo
do jy=1,geom%ny
  do iy=jy,geom%ny
    self%sqrt_merid(iy,jy) = 0.0
    do ky=1,geom%ny
      self%sqrt_merid(iy,jy) = self%sqrt_merid(iy,jy)+evectsy(iy,ky)*evectsy(jy,ky)*evalsy(ky)
    enddo
  enddo
enddo

! Generate the lower triangle of the vertical correlation matrix
do jz=1,geom%nz
  do iz=jz,geom%nz
    evectsz(iz,jz) = exp(-0.5*((geom%z(iz)-geom%z(jz))/vertical_length_scale)**2)
    evectsz(jz,iz) = evectsz(iz,jz)
  enddo
enddo

! Calculate the eigen-decomposition
call dsyev('V','L',geom%nz,evectsz,geom%nz,evalsz,workz,size(workz),info)
if (info/=0) then
  write(record,*) 'qg_error_covariance_setup: dsyev returns info=',info
  call fckit_log%error(record)
  call abor1_ftn ('qg_error_covariance_setup: dsyev failed to sqrt matrix')
endif

! Calculate the lower triangle of the symmetric square-root
threshold = maxval(evalsz)/condition_number
do iz=1,geom%nz
  evalsz(iz) = sqrt(max(threshold,evalsz(iz)))
  revalsz(iz) = 1.0/evalsz(iz)
enddo
do jz=1,geom%nz
  do iz=jz,geom%nz
    self%sqrt_vert(iz,jz) = 0.0
    do kz=1,geom%nz
      self%sqrt_vert(iz,jz) = self%sqrt_vert(iz,jz)+evectsz(iz,kz)*evectsz(jz,kz)*evalsz(kz)
    enddo
  enddo
enddo

! Compute normalization factor
vars = oops_variables()
call vars%push_back('x')
call qg_fields_create(fld_in,geom,vars,.false.)
call qg_fields_create(fld_out,geom,vars,.false.)
self%norm = 1.0
do iz=1,geom%nz
  do iy=1,geom%ny
    call qg_fields_zero(fld_in)
    fld_in%x(1,iy,iz) = 1.0
    call qg_error_covariance_mult(self,fld_in,fld_out)
    norm(iy,iz) = 1.0/sqrt(fld_out%x(1,iy,iz))
  end do
end do
self%norm = norm*self%sigma

! Release memory
deallocate(struct_fn)
deallocate(workx)
deallocate(worky)
deallocate(revalsy)
deallocate(evectsy)
deallocate(evalsy)
deallocate(workz)
deallocate(revalsz)
deallocate(evectsz)
deallocate(evalsz)
call vars%destruct()
call qg_fields_delete(fld_in)
call qg_fields_delete(fld_out)

end subroutine qg_error_covariance_setup
! ------------------------------------------------------------------------------
!> Delete error covariance matrix
subroutine qg_error_covariance_delete(self)

implicit none

! Passed variables
type(qg_error_covariance_config),intent(inout) :: self !< Error covariance configuration

! Release memory
deallocate(self%sqrt_zonal)
deallocate(self%sqrt_merid)
deallocate(self%sqrt_vert)
deallocate(self%norm)

end subroutine qg_error_covariance_delete
! ------------------------------------------------------------------------------
!> Multiply by error covariance matrix
subroutine qg_error_covariance_mult(self,fld_in,fld_out)

implicit none

! Passed variables
type(qg_error_covariance_config),intent(in) :: self !< Error covariance configuration
type(qg_fields),intent(in) :: fld_in                !< Input field
type(qg_fields),intent(inout) :: fld_out            !< Output field

! Local variables
type(qg_fields) :: fld_tmp

! Initialization
call qg_fields_create_from_other(fld_tmp,fld_in,fld_in%geom)
call qg_fields_copy(fld_tmp,fld_in)

! Apply covariance matrix
call qg_error_covariance_sqrt_mult_ad(self,fld_tmp,fld_out)
call qg_fields_copy(fld_tmp,fld_out)
call qg_error_covariance_sqrt_mult(self,fld_tmp,fld_out)

! Release memory
call qg_fields_delete(fld_tmp)

end subroutine qg_error_covariance_mult
! ------------------------------------------------------------------------------
!> Randomize error covariance
subroutine qg_error_covariance_randomize(self,fld_out)

implicit none

! Passed variables
type(qg_error_covariance_config),intent(in) :: self !< Error covariance configuration
type(qg_fields),intent(inout) :: fld_out            !< Output field

! Local variables
type(qg_fields) :: fld_tmp

! Initialize temporary field
call qg_fields_create_from_other(fld_tmp,fld_out,fld_out%geom)
call qg_fields_random(fld_tmp,self%seed)

! Apply square-root of the covariance matrix
call qg_error_covariance_sqrt_mult(self,fld_tmp,fld_out)

! Release memory
call qg_fields_delete(fld_tmp)

end subroutine qg_error_covariance_randomize
! ------------------------------------------------------------------------------
! Private
! ------------------------------------------------------------------------------
!> Multiply by error covariance matrix square-root, zonal part
subroutine qg_error_covariance_sqrt_mult_zonal(self,fld)

implicit none

! Passed variables
type(qg_error_covariance_config),intent(in) :: self           !< Error covariance configuration
real(kind_real),intent(inout) :: fld(self%nx,self%ny,self%nz) !< Field

! Local variables
integer :: iy,iz,m,iri
real(kind_real) :: zfour(self%nx+2)

do iz=1,self%nz
  do iy=1,self%ny
    call fft_fwd(self%nx,fld(:,iy,iz),zfour)
    !$omp parallel do schedule(static) private(m,iri)
    do m=0,self%nx/2
      do iri=1,2
        zfour(2*m+iri) = zfour(2*m+iri)*self%sqrt_zonal(m)
      enddo
    enddo
    !$omp end parallel do
    call fft_inv(self%nx,zfour,fld(:,iy,iz))
  enddo
enddo

end subroutine qg_error_covariance_sqrt_mult_zonal
! ------------------------------------------------------------------------------
!> Multiply by error covariance matrix square-root - meridional part
subroutine qg_error_covariance_sqrt_mult_meridional(self,fld)

implicit none

! Passed variables
type(qg_error_covariance_config),intent(in) :: self           !< Error covariance configuration
real(kind_real),intent(inout) :: fld(self%nx,self%ny,self%nz) !< Field

! Local variables
integer :: ix,iz
real(kind_real),allocatable :: arr_in(:),arr_out(:)

!$omp parallel do schedule(static) private(iz,ix) firstprivate(arr_in,arr_out)
do iz=1,self%nz
  do ix=1,self%nx
    ! Allocation
    allocate(arr_in(self%ny))
    allocate(arr_out(self%ny))

    ! Initialize
    arr_in = fld(ix,:,iz)

    ! Apply transform
    call dsymv('L',self%ny,1.0_kind_real,self%sqrt_merid,self%ny,arr_in,1,0.0_kind_real,arr_out,1)

    ! Copy
    fld(ix,:,iz) = arr_out

    ! Release memory
    deallocate(arr_in)
    deallocate(arr_out)
  enddo
enddo
!$omp end parallel do

end subroutine qg_error_covariance_sqrt_mult_meridional
! ------------------------------------------------------------------------------
!> Multiply by error covariance matrix square-root - vertical part
subroutine qg_error_covariance_sqrt_mult_vertical(self,fld)

implicit none

! Passed variables
type(qg_error_covariance_config),intent(in) :: self           !< Error covariance configuration
real(kind_real),intent(inout) :: fld(self%nx,self%ny,self%nz) !< Field

! Local variables
integer :: ix,iy
real(kind_real),allocatable :: arr_in(:),arr_out(:)

!$omp parallel do schedule(static) private(iy,ix) firstprivate(arr_in,arr_out)
do iy=1,self%ny
  do ix=1,self%nx
    ! Allocation
    allocate(arr_in(self%nz))
    allocate(arr_out(self%nz))

    ! Initialize
    arr_in = fld(ix,iy,:)

    ! Apply transform
    call dsymv('L',self%nz,1.0_kind_real,self%sqrt_vert,self%nz,arr_in,1,0.0_kind_real,arr_out,1)

    ! Copy
    fld(ix,iy,:) = arr_out

    ! Release memory
    deallocate(arr_in)
    deallocate(arr_out)
  enddo
enddo
!$omp end parallel do

end subroutine qg_error_covariance_sqrt_mult_vertical
! ------------------------------------------------------------------------------
!> Multiply by error covariance matrix square-root
subroutine qg_error_covariance_sqrt_mult(self,fld_in,fld_out)

implicit none

! Passed variables
type(qg_error_covariance_config),intent(in) :: self !< Error covariance configuration
type(qg_fields),intent(in) :: fld_in                !< Input field
type(qg_fields),intent(inout) :: fld_out            !< Output field

! Local variables
integer :: ix

! Check input/output
if (.not.allocated(fld_in%x)) call abor1_ftn("qg_error_covariance_sqrt_mult: x required as input")
if (.not.allocated(fld_out%x)) call abor1_ftn("qg_error_covariance_sqrt_mult: x required as output")

! Copy field
call qg_fields_copy(fld_out,fld_in)

! Multiply by symmetric square-root of vertical correlation matrix
call qg_error_covariance_sqrt_mult_vertical(self,fld_out%x)

! Multiply by square-root of meridional correlation matrix
call qg_error_covariance_sqrt_mult_meridional(self,fld_out%x)

! Multiply by square-root of zonal correlation matrix
call qg_error_covariance_sqrt_mult_zonal(self,fld_out%x)

! Multiply by normalization factor
!$omp parallel do schedule(static) private(ix)
do ix=1,fld_out%geom%nx
  fld_out%x(ix,:,:) = fld_out%x(ix,:,:)*self%norm
end do
!$omp end parallel do

! Multiply by standard deviation
fld_out%x = fld_out%x*self%sigma

end subroutine qg_error_covariance_sqrt_mult
! ------------------------------------------------------------------------------
!> Multiply by error covariance matrix square-root - adjoint
subroutine qg_error_covariance_sqrt_mult_ad(self,fld_in,fld_out)

implicit none

! Passed variables
type(qg_error_covariance_config),intent(in) :: self !< Error covariance configuration
type(qg_fields),intent(in) :: fld_in                !< Input field
type(qg_fields),intent(inout) :: fld_out            !< Output field

! Local variables
integer :: ix

! Check input/output
if (.not.allocated(fld_in%x)) call abor1_ftn("qg_error_covariance_sqrt_mult: x required as input")
if (.not.allocated(fld_out%x)) call abor1_ftn("qg_error_covariance_sqrt_mult: x required as output")

! Copy field
call qg_fields_copy(fld_out,fld_in)

! Multiply by standard deviation
fld_out%x = fld_out%x*self%sigma

! Multiply by normalization factor
!$omp parallel do schedule(static) private(ix)
do ix=1,fld_out%geom%nx
  fld_out%x(ix,:,:) = fld_out%x(ix,:,:)*self%norm
end do
!$omp end parallel do

! Multiply by square-root of zonal correlation matrix
call qg_error_covariance_sqrt_mult_zonal(self,fld_out%x)

! Multiply by square-root of meridional correlation matrix
call qg_error_covariance_sqrt_mult_meridional(self,fld_out%x)

! Multiply by symmetric square-root of vertical correlation matrix
call qg_error_covariance_sqrt_mult_vertical(self,fld_out%x)

end subroutine qg_error_covariance_sqrt_mult_ad
! ------------------------------------------------------------------------------
end module qg_error_covariance_mod
