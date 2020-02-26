! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

module qg_geom_mod

use atlas_module
use fckit_configuration_module, only: fckit_configuration
use fckit_log_module,only: fckit_log
use kinds
use iso_c_binding
use qg_constants_mod
use qg_projection_mod

implicit none

private
public :: qg_geom
public :: qg_geom_registry
public :: qg_geom_setup,qg_geom_create_atlas_grid_conf,qg_geom_fill_atlas_fieldset,qg_geom_clone,qg_geom_delete,qg_geom_info
! ------------------------------------------------------------------------------
type :: qg_geom
  integer :: nx                                                 !< Number of points in the zonal direction
  integer :: ny                                                 !< Number of points in the meridional direction
  integer :: nz                                                 !< Number of vertical levels
  real(kind_real) :: deltax                                     !< Zonal cell size
  real(kind_real) :: deltay                                     !< Meridional cell size
  real(kind_real),allocatable :: x(:)                           !< Zonal coordinate
  real(kind_real),allocatable :: y(:)                           !< Meridional coordinate
  real(kind_real),allocatable :: z(:)                           !< Altitude
  real(kind_real),allocatable :: lat(:,:)                       !< Latitude
  real(kind_real),allocatable :: lon(:,:)                       !< Longitude
  real(kind_real),allocatable :: area(:,:)                      !< Area
  real(kind_real),allocatable :: f(:,:)                         !< Coefficients of PV operator
  real(kind_real),allocatable :: f_p(:,:)                       !< Coefficients of PV operator, right eigenvectors
  real(kind_real),allocatable :: f_pinv(:,:)                    !< Coefficients of PV operator, right eigenvectors inverse
  real(kind_real),allocatable :: f_d(:)                         !< Coefficients of PV operator, eigenvalues
  real(kind_real),allocatable :: bet(:)                         !< Beta coefficient
  real(kind_real),allocatable :: heat(:,:)                      !< Heating term
  type(atlas_functionspace_structuredcolumns) :: afunctionspace !< ATLAS function space
  type(atlas_fieldset) :: afieldset                             !< ATLAS fieldset
end type qg_geom

#define LISTED_TYPE qg_geom

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry
type(registry_t) :: qg_geom_registry
! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
! Public
! ------------------------------------------------------------------------------
!> Linked list implementation
#include "oops/util/linkedList_c.f"
! ------------------------------------------------------------------------------
!> Setup geometry
subroutine qg_geom_setup(self,f_conf)

! Passed variables
type(qg_geom),intent(inout) :: self            !< Geometry
type(fckit_configuration),intent(in) :: f_conf !< FCKIT configuration

! Local variables
integer :: ix,iy,iz,ix_c,iy_c,lwork,info
integer,allocatable :: ipiv(:),ipivsave(:)
real(kind_real) :: mapfac,distx,disty,f
real(kind_real),allocatable :: real_array(:),depths(:),wi(:),vl(:,:),work(:)
real(kind_real),allocatable :: fsave(:,:),vrlu(:,:),vrlusave(:,:)
character(len=1024) :: htype,record
character(len=:),allocatable :: str

! Get horizontal resolution data
call f_conf%get_or_die("nx",self%nx)
call f_conf%get_or_die("ny",self%ny)
self%nz = f_conf%get_size("depths")

! Allocation
allocate(depths(self%nz))
allocate(self%x(self%nx))
allocate(self%y(self%ny))
allocate(self%z(self%nz))
allocate(self%lon(self%nx,self%ny))
allocate(self%lat(self%nx,self%ny))
allocate(self%area(self%nx,self%ny))
allocate(self%f(self%nz,self%nz))
allocate(self%f_p(self%nz,self%nz))
allocate(self%f_pinv(self%nz,self%nz))
allocate(self%f_d(self%nz))
allocate(self%bet(self%ny))
allocate(self%heat(self%nx,self%ny))
allocate(wi(self%nz))
allocate(vl(self%nz,self%nz))
allocate(fsave(self%nz,self%nz))
allocate(vrlu(self%nz,self%nz))
allocate(vrlusave(self%nz,self%nz))
allocate(ipiv(self%nz))
allocate(ipivsave(self%nz))

! Get depths
call f_conf%get_or_die("depths",real_array)
depths = real_array

! Define dx/dy
self%deltax = domain_zonal/real(self%nx,kind_real)
self%deltay = domain_meridional/real(self%ny+1,kind_real)

! Print sizes and dx/dy
write(record,'(a,i3,a,i3,a,i3)') 'qg_geom_create: nx/ny/nz = ',self%nx,' /',self%ny,' /',self%nz
call fckit_log%info(record)
write(record,'(a,f7.2,a,f7.2,a)') '                deltax/deltay = ',self%deltax*1.0e-3,' km / ',self%deltay*1.0e-3,' km'
call fckit_log%info(record)

! Define x/y
do ix=1,self%nx
  self%x(ix) = (real(ix,kind_real)-0.5)*self%deltax
enddo
do iy=1,self%ny
  self%y(iy) = real(iy,kind_real)*self%deltay
enddo

! Define lon/lat/area
do iy=1,self%ny
  do ix=1,self%nx
    call xy_to_lonlat(self%x(ix),self%y(iy),self%lon(ix,iy),self%lat(ix,iy),mapfac)
    self%area(ix,iy) = self%deltax*self%deltay/mapfac
  end do
end do

! Set heights
self%z(1) = 0.5*depths(1)
do iz=2,self%nz
  self%z(iz) = sum(depths(1:iz-1))+0.5*depths(iz)
end do

! Coefficients of PV operator
self%f = 0.0
do iz=1,self%nz
  f = f0**2*real(self%nz-1,kind_real)/(g*dlogtheta*depths(iz))
  if (iz>1) then
    self%f(iz,iz-1) = f
    self%f(iz,iz) = self%f(iz,iz)-f
  end if
  if (iz<self%nz) then
    self%f(iz,iz+1) = f
    self%f(iz,iz) = self%f(iz,iz)-f
  end if
enddo

! Compute eigendecomposition of ff
fsave = self%f
allocate(work(1))
call dgeev('V','V',self%nz,fsave,self%nz,self%f_d,wi,vl,self%nz,self%f_p,self%nz,work,-1,info)
if (info/=0) call abor1_ftn('error in dgeev, first pass')
lwork = int(work(1))
deallocate(work)
allocate(work(lwork))
fsave = self%f
call dgeev('V','V',self%nz,fsave,self%nz,self%f_d,wi,vl,self%nz,self%f_p,self%nz,work,lwork,info)
if (info/=0) call abor1_ftn('error in dgeev, second pass')
deallocate(work)

! Compute inverse of right eigenvectors of ff
vrlu = self%f_p
call dgetrf(self%nz,self%nz,vrlu,self%nz,ipiv,info)
if (info/=0) call abor1_ftn('error in dgetrf')
allocate(work(1))
vrlusave = vrlu
ipivsave = ipiv
call dgetri(self%nz,vrlusave,self%nz,ipivsave,work,-1,info)
if (info/=0) call abor1_ftn('error in dgetri, first pass')
lwork = int(work(1))
deallocate(work)
allocate(work(lwork))
self%f_pinv = vrlu
ipivsave = ipiv
call dgetri(self%nz,self%f_pinv,self%nz,ipivsave,work,lwork,info)
if (info/=0) call abor1_ftn('error in dgetri, second pass')
deallocate(work)

! Beta coefficient
do iy=1,self%ny
  self%bet(iy) = real(iy-(self%ny+1)/2,kind_real)*self%deltay*bet0
enddo

! Set heating term 
if (f_conf%has("heating")) then
  call f_conf%get_or_die("heating",str)
  htype = str
else
  ! This default value is for backward compatibility
  htype = 'gaussian'
endif
call fckit_log%info('qg_geom_setup: heating type = '//trim(htype))
select case (trim(htype))
case ('none')
  ! No heating term
  self%heat = 0.0
case ('gaussian')
  ! Gaussian source
  ix_c = self%nx/4
  iy_c = 3*self%ny/4
  do iy=1,self%ny
    do ix=1,self%nx
      distx = abs(self%x(ix)-self%x(ix_c))
      if (distx>0.5*domain_zonal) distx = domain_zonal-distx
      disty = abs(self%y(iy)-self%y(iy_c))
      self%heat(ix,iy) = heating_amplitude*exp(-(distx**2+disty**2)/heating_scale**2)
    enddo
  enddo
case default
  call abor1_ftn('qg_geom_setup: wrong heating type')
endselect

end subroutine qg_geom_setup
! ------------------------------------------------------------------------------
!> Create ATLAS grid configuration
subroutine qg_geom_create_atlas_grid_conf(self,fconf)

! Passed variables
type(qg_geom),intent(inout) :: self              !< Geometry
type(fckit_configuration),intent(inout) :: fconf !< ATLAS grid configuration

! Local variables
type(fckit_configuration) :: xspace,yspace

! Create xspace configuration
xspace = fckit_configuration()
call xspace%set("type","linear")
call xspace%set("N",self%nx)
call xspace%set("start",self%lon(1,1))
call xspace%set("end",self%lon(self%nx,1))
call xspace%set("endpoint",.true.)

! Create yspace configuration
yspace = fckit_configuration()
call yspace%set("type","custom")
call yspace%set("N",self%ny)
call yspace%set("values",self%lat(1,:))

! Create ATLAS grid configuration
call fconf%set("type","structured")
call fconf%set("xspace",xspace)
call fconf%set("yspace",yspace)

end subroutine qg_geom_create_atlas_grid_conf
! ------------------------------------------------------------------------------
!> Fill ATLAS fieldset
subroutine qg_geom_fill_atlas_fieldset(self,afieldset)

! Passed variables
type(qg_geom),intent(inout) :: self                                      !< Geometry
type(atlas_fieldset),intent(inout) :: afieldset                          !< ATLAS fieldset

! Local variables
integer :: i,j,k,node
real(kind_real),pointer :: real_ptr_1(:),real_ptr_2(:,:)
type(atlas_field) :: afield

! Add area
afield = self%afunctionspace%create_field(name='area',kind=atlas_real(kind_real),levels=0)
call afield%data(real_ptr_1)
node = 0
do j=self%afunctionspace%j_begin(),self%afunctionspace%j_end()
  do i=self%afunctionspace%i_begin(j),self%afunctionspace%i_end(j)
    node = node+1
    real_ptr_1(node) = self%area(i,j)
  enddo
enddo
call afieldset%add(afield)
call afield%final()

! Add vertical unit
afield = self%afunctionspace%create_field(name='vunit',kind=atlas_real(kind_real),levels=self%nz)
call afield%data(real_ptr_2)
do k=1,self%nz
  real_ptr_2(k,1:self%afunctionspace%size_owned()) = self%z(k)
end do
call afieldset%add(afield)
call afield%final()

end subroutine qg_geom_fill_atlas_fieldset
! ------------------------------------------------------------------------------
!> Clone geometry
subroutine qg_geom_clone(self,other)

! Passed variables
type(qg_geom),intent(inout) :: self !< Geometry
type(qg_geom),intent(in) :: other   !< Other geometry

! Copy dimensions
self%nx = other%nx
self%ny = other%ny
self%nz = other%nz

! Allocation
allocate(self%x(self%nx))
allocate(self%y(self%ny))
allocate(self%z(self%nz))
allocate(self%lon(self%nx,self%ny))
allocate(self%lat(self%nx,self%ny))
allocate(self%area(self%nx,self%ny))
allocate(self%f(self%nz,self%nz))
allocate(self%f_p(self%nz,self%nz))
allocate(self%f_pinv(self%nz,self%nz))
allocate(self%f_d(self%nz))
allocate(self%bet(self%ny))
allocate(self%heat(self%nx,self%ny))

! Copy data
self%deltax = other%deltax
self%deltay = other%deltay
self%x = other%x
self%y = other%y
self%z = other%z
self%lon = other%lon
self%lat = other%lat
self%area = other%area
self%f = other%f
self%f_p = other%f_p
self%f_pinv = other%f_pinv
self%f_d = other%f_d
self%bet = other%bet
self%heat = other%heat
self%afunctionspace = atlas_functionspace_nodecolumns(other%afunctionspace%c_ptr())
self%afieldset = atlas_functionspace_nodecolumns(other%afieldset%c_ptr())

end subroutine qg_geom_clone
! ------------------------------------------------------------------------------
!> Delete geometry
subroutine qg_geom_delete(self)

! Passed variables
type(qg_geom),intent(inout) :: self !< Geometry

! Release memory
deallocate(self%x)
deallocate(self%y)
deallocate(self%z)
deallocate(self%lon)
deallocate(self%lat)
deallocate(self%area)
deallocate(self%f)
deallocate(self%f_p)
deallocate(self%f_pinv)
deallocate(self%f_d)
deallocate(self%bet)
deallocate(self%heat)
call self%afunctionspace%final()
call self%afieldset%final()

end subroutine qg_geom_delete
! ------------------------------------------------------------------------------
!> Get geometry info
subroutine qg_geom_info(self,nx,ny,nz,deltax,deltay)

! Passed variables
type(qg_geom),intent(in) :: self      !< Geometry
integer,intent(out) :: nx             !< Number of points in the zonal direction
integer,intent(out) :: ny             !< Number of points in the meridional direction
integer,intent(out) :: nz             !< Number of vertical levels
real(kind_real),intent(out) :: deltax !< Zonal cell size
real(kind_real),intent(out) :: deltay !< Meridional cell size

! Copy data
nx = self%nx
ny = self%ny
nz = self%nz
deltax = self%deltax
deltay = self%deltay

end subroutine qg_geom_info
! ------------------------------------------------------------------------------
end module qg_geom_mod
