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
use qg_constants_mod
use qg_fields_mod
use qg_geom_mod
use oops_variables_mod
use random_mod

implicit none

private
public :: qg_error_covariance_config
public :: qg_error_covariance_registry
public :: qg_error_covariance_setup,qg_error_covariance_delete,qg_error_covariance_mult, &
        & qg_error_covariance_randomize,qg_error_covariance_test
! ------------------------------------------------------------------------------
type :: qg_error_covariance_config
  real(kind_real) :: sigma                           !< Standard deviation
  real(kind_real),allocatable :: sqrt_zonal(:)       !< Spectral weights for the spectral of the zonal correlation matrix
  real(kind_real),allocatable :: sqrt_merid(:,:)     !< Square-root of the meridional correlation matrix
  real(kind_real),allocatable :: sqrt_vert(:,:)      !< Square-root of the meridional correlation matrix
  real(kind_real),allocatable :: norm(:,:)           !< Normalization factor
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

! Get parameters
call f_conf%get_or_die("standard_deviation",self%sigma)
call f_conf%get_or_die("horizontal_length_scale",horizontal_length_scale)
call f_conf%get_or_die("vertical_length_scale",vertical_length_scale)
call f_conf%get_or_die("maximum_condition_number",condition_number)

! Check nx
if (mod(geom%nx,2)/=0) then
  write(record,*) 'qg_error_covariance_setup: number of zonal gridpoints nx=',geom%nx
  call fckit_log%error(record)
  call abor1_ftn('qg_error_covariance_setup: odd number of zonal grid points')
endif

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
call qg_fields_create_default(fld_in,geom,.false.)
call qg_fields_create_default(fld_out,geom,.false.)
self%norm = 1.0
do iz=1,geom%nz
  do iy=1,geom%ny
    call qg_fields_zero(fld_in)
    fld_in%gfld3d(1,iy,iz) = 1.0
    call qg_error_covariance_mult(self,fld_in,fld_out)
    norm(iy,iz) = 1.0/sqrt(fld_out%gfld3d(1,iy,iz))
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
subroutine qg_error_covariance_mult(conf,fld_in,fld_out)

implicit none

! Passed variables
type(qg_error_covariance_config),intent(in) :: conf !< Error covariance configuration
type(qg_fields),intent(in) :: fld_in                !< Input field
type(qg_fields),intent(inout) :: fld_out            !< Output field

! Local variables
type(qg_fields) :: fld_tmp

! Initialization
call qg_fields_create_from_other(fld_tmp,fld_in)
call qg_fields_copy(fld_tmp,fld_in)

! Apply covariance matrix
call qg_error_covariance_sqrt_mult_ad(conf,fld_tmp,fld_out)
call qg_fields_copy(fld_tmp,fld_out)
call qg_error_covariance_sqrt_mult(conf,fld_tmp,fld_out)

end subroutine qg_error_covariance_mult
! ------------------------------------------------------------------------------
!> Randomize error covariance
subroutine qg_error_covariance_randomize(conf,fld_out)

implicit none

! Passed variables
type(qg_error_covariance_config),intent(in) :: conf !< Error covariance configuration
type(qg_fields),intent(inout) :: fld_out            !< Output field

! Local variables
type(qg_fields) :: fld_tmp

! Initialize temporary field
call qg_fields_create_from_other(fld_tmp,fld_out)
call qg_fields_random(fld_tmp)

! Apply square-root of the covariance matrix
call qg_error_covariance_sqrt_mult(conf,fld_tmp,fld_out)

end subroutine qg_error_covariance_randomize
! ------------------------------------------------------------------------------
!> Test error covariance adjoint
subroutine qg_error_covariance_test(conf,fld_in)

implicit none

! Passed variables
type(qg_error_covariance_config),intent(in) :: conf !< Error covariance configuration
type(qg_fields),intent(in) :: fld_in                !< Input field

! Local variables
type(qg_fields) :: fld1,fld2,fld3,fld1_save,fld2_save
real(kind_real) :: dot1,dot2
character(len=1024) :: record

! Allocation
call qg_fields_create_from_other(fld1,fld_in)
call qg_fields_create_from_other(fld2,fld_in)
call qg_fields_create_from_other(fld3,fld_in)
call qg_fields_create_from_other(fld1_save,fld_in)
call qg_fields_create_from_other(fld2_save,fld_in)

! Initialization
call qg_fields_random(fld1_save)
call qg_fields_random(fld2_save)

! Zonal correlation square-root auto-adjointness
call qg_error_covariance_sqrt_mult_zonal(conf,fld1_save,fld2)
call qg_error_covariance_sqrt_mult_zonal(conf,fld2_save,fld1)
call qg_fields_dot_prod(fld1,fld1_save,dot1)
call qg_fields_dot_prod(fld2,fld2_save,dot2)
if (abs(dot1-dot2)/abs(dot1)>eps_ad) then
  write(record,*) 'qg_error_covariance_test: sqrt_mult_zonal/sqrt_mult_zonal failed: ',abs(dot1-dot2)/abs(dot1)
  call abor1_ftn(record)
endif

! Meridional correlation square-root auto-adjointness
call qg_error_covariance_sqrt_mult_meridional(conf,fld1_save,fld2)
call qg_error_covariance_sqrt_mult_meridional(conf,fld2_save,fld1)
call qg_fields_dot_prod(fld1,fld1_save,dot1)
call qg_fields_dot_prod(fld2,fld2_save,dot2)
if (abs(dot1-dot2)/abs(dot1)>eps_ad) then
  write(record,*) 'qg_error_covariance_test: sqrt_mult_meridional/sqrt_mult_meridional failed: ',abs(dot1-dot2)/abs(dot1)
  call abor1_ftn(record)
endif

! Vertical correlation square-root auto-adjointness
call qg_error_covariance_sqrt_mult_vertical(conf,fld1_save,fld2)
call qg_error_covariance_sqrt_mult_vertical(conf,fld2_save,fld1)
call qg_fields_dot_prod(fld1,fld1_save,dot1)
call qg_fields_dot_prod(fld2,fld2_save,dot2)
if (abs(dot1-dot2)/abs(dot1)>eps_ad) then
  write(record,*) 'qg_error_covariance_test: sqrt_mult_vertical/sqrt_mult_vertical failed: ',abs(dot1-dot2)/abs(dot1)
  call abor1_ftn(record)
endif

! Correlation square-root adjointness
call qg_error_covariance_sqrt_mult(conf,fld1_save,fld2)
call qg_error_covariance_sqrt_mult_ad(conf,fld2_save,fld1)
call qg_fields_dot_prod(fld1,fld1_save,dot1)
call qg_fields_dot_prod(fld2,fld2_save,dot2)
if (abs(dot1-dot2)/abs(dot1)>eps_ad) then
  write(record,*) 'qg_error_covariance_test: sqrt_mult/sqrt_mult_ad failed: ',abs(dot1-dot2)/abs(dot1)
  call abor1_ftn(record)
endif

! Correlation auto-adjointness
call qg_error_covariance_mult(conf,fld1_save,fld2)
call qg_error_covariance_mult(conf,fld2_save,fld1)
call qg_fields_dot_prod(fld1,fld1_save,dot1)
call qg_fields_dot_prod(fld2,fld2_save,dot2)
if (abs(dot1-dot2)/abs(dot1)>eps_ad) then
  write(record,*) 'qg_error_covariance_test: mult/mult failed: ',abs(dot1-dot2)/abs(dot1)
  call abor1_ftn(record)
endif

end subroutine qg_error_covariance_test
! ------------------------------------------------------------------------------
! Private
! ------------------------------------------------------------------------------
!> Multiply by error covariance matrix square-root, zonal part
subroutine qg_error_covariance_sqrt_mult_zonal(conf,fld_in,fld_out)

implicit none

! Passed variables
type(qg_error_covariance_config),intent(in) :: conf !< Error covariance configuration
type(qg_fields),intent(in) :: fld_in                !< Input field
type(qg_fields),intent(inout) :: fld_out            !< Output field

! Local variables
integer :: iy,iz,m,iri
real(kind_real) :: zfour(fld_in%geom%nx+2)

do iz=1,fld_in%geom%nz
  do iy=1,fld_in%geom%ny
    call fft_fwd(fld_in%geom%nx,fld_in%gfld3d(:,iy,iz),zfour)
    !$omp parallel do schedule(static) private(m,iri)
    do m=0,fld_in%geom%nx/2
      do iri=1,2
        zfour(2*m+iri) = zfour(2*m+iri)*conf%sqrt_zonal(m)
      enddo
    enddo
    !$omp end parallel do
    call fft_inv(fld_in%geom%nx,zfour,fld_out%gfld3d(:,iy,iz))
  enddo
enddo

end subroutine qg_error_covariance_sqrt_mult_zonal
! ------------------------------------------------------------------------------
!> Multiply by error covariance matrix square-root - meridional part
subroutine qg_error_covariance_sqrt_mult_meridional(conf,fld_in,fld_out)

implicit none

! Passed variables
type(qg_error_covariance_config),intent(in) :: conf !< Error covariance configuration
type(qg_fields),intent(in) :: fld_in                !< Input field
type(qg_fields),intent(inout) :: fld_out            !< Output field

! Local variables
integer :: ix,iz
real(kind_real),allocatable :: arr_in(:),arr_out(:)

!$omp parallel do schedule(static) private(iz,ix) firstprivate(arr_in,arr_out)
do iz=1,fld_in%geom%nz
  do ix=1,fld_in%geom%nx
    ! Allocation
    allocate(arr_in(fld_in%geom%ny))
    allocate(arr_out(fld_in%geom%ny))

    ! Initialize
    arr_in = fld_in%gfld3d(ix,:,iz)

    ! Apply transform
    call dsymv('L',fld_in%geom%ny,1.0_kind_real,conf%sqrt_merid,fld_in%geom%ny,arr_in,1,0.0_kind_real,arr_out,1)

    ! Copy
    fld_out%gfld3d(ix,:,iz) = arr_out

    ! Release memory
    deallocate(arr_in)
    deallocate(arr_out)
  enddo
enddo
!$omp end parallel do

end subroutine qg_error_covariance_sqrt_mult_meridional
! ------------------------------------------------------------------------------
!> Multiply by error covariance matrix square-root - vertical part
subroutine qg_error_covariance_sqrt_mult_vertical(conf,fld_in,fld_out)

implicit none

! Passed variables
type(qg_error_covariance_config),intent(in) :: conf !< Error covariance configuration
type(qg_fields),intent(in) :: fld_in                !< Input field
type(qg_fields),intent(inout) :: fld_out            !< Output field

! Local variables
integer :: ix,iy
real(kind_real),allocatable :: arr_in(:),arr_out(:)

!$omp parallel do schedule(static) private(iy,ix) firstprivate(arr_in,arr_out)
do iy=1,fld_in%geom%ny
  do ix=1,fld_in%geom%nx
    ! Allocation
    allocate(arr_in(fld_in%geom%nz))
    allocate(arr_out(fld_in%geom%nz))

    ! Initialize
    arr_in = fld_in%gfld3d(ix,iy,:)

    ! Apply transform
    call dsymv('L',fld_in%geom%nz,1.0_kind_real,conf%sqrt_vert,fld_in%geom%nz,arr_in,1,0.0_kind_real,arr_out,1)

    ! Copy
    fld_out%gfld3d(ix,iy,:) = arr_out

    ! Release memory
    deallocate(arr_in)
    deallocate(arr_out)
  enddo
enddo
!$omp end parallel do

end subroutine qg_error_covariance_sqrt_mult_vertical
! ------------------------------------------------------------------------------
!> Multiply by error covariance matrix square-root
subroutine qg_error_covariance_sqrt_mult(conf,fld_in,fld_out)

implicit none

! Passed variables
type(qg_error_covariance_config),intent(in) :: conf !< Error covariance configuration
type(qg_fields),intent(in) :: fld_in                !< Input field
type(qg_fields),intent(inout) :: fld_out            !< Output field

! Local variables
integer :: ix
type(qg_fields) :: fld_tmp

! Initialization
call qg_fields_create_from_other(fld_tmp,fld_in)
call qg_fields_copy(fld_tmp,fld_in)

! Multiply by symmetric square-root of vertical correlation matrix
call qg_error_covariance_sqrt_mult_vertical(conf,fld_tmp,fld_out)
call qg_fields_copy(fld_tmp,fld_out)

! Multiply by square-root of meridional correlation matrix
call qg_error_covariance_sqrt_mult_meridional(conf,fld_tmp,fld_out)
call qg_fields_copy(fld_tmp,fld_out)

! Multiply by square-root of zonal correlation matrix
call qg_error_covariance_sqrt_mult_zonal(conf,fld_tmp,fld_out)
call qg_fields_copy(fld_tmp,fld_out)

! Multiply by normalization factor
!$omp parallel do schedule(static) private(ix)
do ix=1,fld_in%geom%nx
  fld_out%gfld3d(ix,:,:) = fld_tmp%gfld3d(ix,:,:)*conf%norm
end do
!$omp end parallel do
call qg_fields_copy(fld_tmp,fld_out)

! Multiply by standard deviation
fld_out%gfld3d = fld_tmp%gfld3d*conf%sigma

end subroutine qg_error_covariance_sqrt_mult
! ------------------------------------------------------------------------------
!> Multiply by error covariance matrix square-root - adjoint
subroutine qg_error_covariance_sqrt_mult_ad(conf,fld_in,fld_out)

implicit none

! Passed variables
type(qg_error_covariance_config),intent(in) :: conf !< Error covariance configuration
type(qg_fields),intent(in) :: fld_in                !< Input field
type(qg_fields),intent(inout) :: fld_out            !< Output field

! Local variables
integer :: ix
type(qg_fields) :: fld_tmp

! Initialization
call qg_fields_create_from_other(fld_tmp,fld_in)
call qg_fields_copy(fld_tmp,fld_in)

! Multiply by standard deviation
fld_out%gfld3d = fld_tmp%gfld3d*conf%sigma
call qg_fields_copy(fld_tmp,fld_out)

! Multiply by normalization factor
!$omp parallel do schedule(static) private(ix)
do ix=1,fld_in%geom%nx
  fld_out%gfld3d(ix,:,:) = fld_tmp%gfld3d(ix,:,:)*conf%norm
end do
!$omp end parallel do
call qg_fields_copy(fld_tmp,fld_out)

! Multiply by square-root of zonal correlation matrix
call qg_error_covariance_sqrt_mult_zonal(conf,fld_tmp,fld_out)
call qg_fields_copy(fld_tmp,fld_out)

! Multiply by square-root of meridional correlation matrix
call qg_error_covariance_sqrt_mult_meridional(conf,fld_tmp,fld_out)
call qg_fields_copy(fld_tmp,fld_out)

! Multiply by symmetric square-root of vertical correlation matrix
call qg_error_covariance_sqrt_mult_vertical(conf,fld_tmp,fld_out)

end subroutine qg_error_covariance_sqrt_mult_ad
! ------------------------------------------------------------------------------
end module qg_error_covariance_mod
