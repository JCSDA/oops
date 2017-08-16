! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

!> Handle fields for the QG model

module qg_fields

use config_mod
use qg_geom_mod
use qg_locs_mod
use qg_vars_mod
use qg_goms_mod
use kinds

implicit none
private

public :: qg_field, &
        & create, delete, zeros, dirac, random, copy, &
        & self_add, self_schur, self_sub, self_mul, axpy, &
        & dot_prod, add_incr, diff_incr, &
        & read_file, write_file, gpnorm, fldrms, &
        & change_resol, interp_tl, interp_ad, convert_to_ug, convert_from_ug

! ------------------------------------------------------------------------------

public :: qg_field_registry

! ------------------------------------------------------------------------------

!> Fortran derived type to hold QG fields
type :: qg_field
  integer :: nx                     !< Zonal grid dimension
  integer :: ny                     !< Meridional grid dimension
  real(kind=kind_real), allocatable :: lons(:)  !< Longitudes
  real(kind=kind_real), allocatable :: lats(:)  !< Latitudes
  integer :: nl                     !< Number of levels
  integer :: nf                     !< Number of fields
  logical :: lbc                    !< North-South boundary is present
  real(kind=kind_real), pointer :: gfld3d(:,:,:)    !< 3D fields
  real(kind=kind_real), pointer :: x(:,:,:)         !< Stream function
  real(kind=kind_real), pointer :: q(:,:,:)         !< Potential vorticity
  real(kind=kind_real), pointer :: u(:,:,:)         !< Zonal wind
  real(kind=kind_real), pointer :: v(:,:,:)         !< Meridional wind
  real(kind=kind_real), pointer :: xbound(:)        !< Streamfunction on walls
  real(kind=kind_real), pointer :: x_north(:)       !< Streamfunction on northern wall
  real(kind=kind_real), pointer :: x_south(:)       !< Streamfunction on southern wall
  real(kind=kind_real), pointer :: qbound(:,:)      !< PV on walls
  real(kind=kind_real), pointer :: q_north(:,:)     !< PV on northern wall
  real(kind=kind_real), pointer :: q_south(:,:)     !< PV on southern wall
  character(len=1), allocatable :: fldnames(:)      !< Variable identifiers
end type qg_field

#define LISTED_TYPE qg_field

!> Linked list interface - defines registry_t type
#include "util/linkedList_i.f"

!> Global registry
type(registry_t) :: qg_field_registry

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Linked list implementation
#include "util/linkedList_c.f"

! ------------------------------------------------------------------------------

subroutine create(self, geom, vars)
implicit none
type(qg_field), intent(inout) :: self
type(qg_geom),  intent(in)    :: geom
type(qg_vars),  intent(in)    :: vars
integer :: ioff

self%nx = geom%nx
self%ny = geom%ny
allocate(self%lons(self%nx))
allocate(self%lats(self%ny))
self%lons = geom%lons
self%lats = geom%lats
self%nl = 2
self%nf = vars%nv
self%lbc = vars%lbc

allocate(self%gfld3d(self%nx,self%ny,self%nl*self%nf))
self%gfld3d(:,:,:)=0.0_kind_real

ioff=0
self%x => self%gfld3d(:,:,ioff+1:ioff+self%nl)
ioff =ioff + self%nl

if (self%nf>1) then
  self%q => self%gfld3d(:,:,ioff+1:ioff+self%nl)
  ioff =ioff + self%nl
  self%u => self%gfld3d(:,:,ioff+1:ioff+self%nl)
  ioff =ioff + self%nl
  self%v => self%gfld3d(:,:,ioff+1:ioff+self%nl)
  ioff =ioff + self%nl
else
  self%q => null()
  self%u => null()
  self%v => null()
endif
if (ioff/=self%nf*self%nl)  call abor1_ftn ("qg_fields:create error number of fields")

allocate(self%fldnames(self%nf))
self%fldnames(:)=vars%fldnames(:)

if (self%lbc) then
  allocate(self%xbound(4))
  self%x_north => self%xbound(1:2)
  self%x_south => self%xbound(3:4)
  allocate(self%qbound(self%nx,4))
  self%q_north => self%qbound(:,1:2)
  self%q_south => self%qbound(:,3:4)
else
  self%xbound => null()
  self%x_north => null()
  self%x_south => null()
  self%qbound => null()
  self%q_north => null()
  self%q_south => null()
endif

call check(self)

end subroutine create

! ------------------------------------------------------------------------------

subroutine delete(self)
implicit none
type(qg_field), intent(inout) :: self

call check(self)

if (associated(self%gfld3d)) deallocate(self%gfld3d)
if (associated(self%xbound)) deallocate(self%xbound)
if (associated(self%qbound)) deallocate(self%qbound)
if (allocated(self%fldnames)) deallocate(self%fldnames)

end subroutine delete

! ------------------------------------------------------------------------------

subroutine zeros(self)
implicit none
type(qg_field), intent(inout) :: self

call check(self)

self%gfld3d(:,:,:) = 0.0_kind_real
if (self%lbc) then
  self%xbound(:)   = 0.0_kind_real
  self%qbound(:,:) = 0.0_kind_real
endif

end subroutine zeros

! ------------------------------------------------------------------------------

subroutine dirac(self, c_conf)
use iso_c_binding
implicit none
type(qg_field), intent(inout) :: self
type(c_ptr), intent(in)       :: c_conf   !< Configuration
integer :: ndir,idir,ildir,ifdir,ioff
integer,allocatable :: ixdir(:),iydir(:)
character(len=3) :: idirchar

call check(self)

! Get Diracs positions
ndir = config_get_int(c_conf,"ndir")
allocate(ixdir(ndir))
allocate(iydir(ndir))
do idir=1,ndir
   write(idirchar,'(i3)') idir
   ixdir(idir) = config_get_int(c_conf,"ixdir("//trim(adjustl(idirchar))//")")
   iydir(idir) = config_get_int(c_conf,"iydir("//trim(adjustl(idirchar))//")")
end do
ildir = config_get_int(c_conf,"ildir")
ifdir = config_get_int(c_conf,"ifdir")

! Check
if (ndir<1) call abor1_ftn("qg_fields:dirac non-positive ndir")
if (any(ixdir<1).or.any(ixdir>self%nx)) call abor1_ftn("qg_fields:dirac invalid ixdir")
if (any(iydir<1).or.any(iydir>self%ny)) call abor1_ftn("qg_fields:dirac invalid iydir")
if ((ildir<1).or.(ildir>self%nl)) call abor1_ftn("qg_fields:dirac invalid ildir")
if ((ifdir<1).or.(ifdir>self%nf)) call abor1_ftn("qg_fields:dirac invalid ifdir")

! Setup Diracs
call zeros(self)
ioff = (ifdir-1)*self%nl
do idir=1,ndir
   self%gfld3d(ixdir(idir),iydir(idir),ioff+ildir) = 1.0
end do

end subroutine dirac

! ------------------------------------------------------------------------------

subroutine random(self)
use random_vectors_mod
implicit none
type(qg_field), intent(inout) :: self

call check(self)

call random_vector(self%gfld3d(:,:,:))
if (self%lbc) then
  call random_vector(self%xbound(:))
  call random_vector(self%qbound(:,:))
endif

end subroutine random

! ------------------------------------------------------------------------------

subroutine copy(self,rhs)
implicit none
type(qg_field), intent(inout) :: self
type(qg_field), intent(in)    :: rhs
integer :: nf

call check_resolution(self, rhs)

nf = common_vars(self, rhs)
self%gfld3d(:,:,1:nf) = rhs%gfld3d(:,:,1:nf)
if (self%nf>nf) self%gfld3d(:,:,nf+1:self%nf) = 0.0_kind_real

if (self%lbc) then
  if (rhs%lbc) then
    self%xbound(:) = rhs%xbound(:)
    self%qbound(:,:) = rhs%qbound(:,:)
  else
    self%xbound(:) = 0.0_kind_real
    self%qbound(:,:) = 0.0_kind_real
  endif
endif

return
end subroutine copy

! ------------------------------------------------------------------------------

subroutine self_add(self,rhs)
implicit none
type(qg_field), intent(inout) :: self
type(qg_field), intent(in)    :: rhs
integer :: nf

call check_resolution(self, rhs)

nf = common_vars(self, rhs)
self%gfld3d(:,:,1:nf) = self%gfld3d(:,:,1:nf) + rhs%gfld3d(:,:,1:nf)

if (self%lbc .and. rhs%lbc) then
  self%xbound(:)   = self%xbound(:)   + rhs%xbound(:)
  self%qbound(:,:) = self%qbound(:,:) + rhs%qbound(:,:)
endif

return
end subroutine self_add

! ------------------------------------------------------------------------------

subroutine self_schur(self,rhs)
implicit none
type(qg_field), intent(inout) :: self
type(qg_field), intent(in)    :: rhs
integer :: nf

call check_resolution(self, rhs)

nf = common_vars(self, rhs)
self%gfld3d(:,:,1:nf) = self%gfld3d(:,:,1:nf) * rhs%gfld3d(:,:,1:nf)

if (self%lbc .and. rhs%lbc) then
  self%xbound(:)   = self%xbound(:)   * rhs%xbound(:)
  self%qbound(:,:) = self%qbound(:,:) * rhs%qbound(:,:)
endif

return
end subroutine self_schur

! ------------------------------------------------------------------------------

subroutine self_sub(self,rhs)
implicit none
type(qg_field), intent(inout) :: self
type(qg_field), intent(in)    :: rhs
integer :: nf

call check_resolution(self, rhs)

nf = common_vars(self, rhs)
self%gfld3d(:,:,1:nf) = self%gfld3d(:,:,1:nf) - rhs%gfld3d(:,:,1:nf)

if (self%lbc .and. rhs%lbc) then
  self%xbound(:)   = self%xbound(:)   - rhs%xbound(:)
  self%qbound(:,:) = self%qbound(:,:) - rhs%qbound(:,:)
endif

return
end subroutine self_sub

! ------------------------------------------------------------------------------

subroutine self_mul(self,zz)
implicit none
type(qg_field), intent(inout) :: self
real(kind=kind_real), intent(in) :: zz

call check(self)

self%gfld3d(:,:,:) = zz * self%gfld3d(:,:,:)
if (self%lbc) then
  self%xbound(:)   = zz * self%xbound(:)
  self%qbound(:,:) = zz * self%qbound(:,:)
endif

return
end subroutine self_mul

! ------------------------------------------------------------------------------

subroutine axpy(self,zz,rhs)
implicit none
type(qg_field), intent(inout) :: self
real(kind=kind_real), intent(in) :: zz
type(qg_field), intent(in)    :: rhs
integer :: nf

call check_resolution(self, rhs)

nf = common_vars(self, rhs)
self%gfld3d(:,:,1:nf) = self%gfld3d(:,:,1:nf) + zz * rhs%gfld3d(:,:,1:nf)

if (self%lbc .and. rhs%lbc) then
  self%xbound(:)   = self%xbound(:)   + zz * rhs%xbound(:)
  self%qbound(:,:) = self%qbound(:,:) + zz * rhs%qbound(:,:)
endif

return
end subroutine axpy

! ------------------------------------------------------------------------------

subroutine dot_prod(fld1,fld2,zprod)
implicit none
type(qg_field), intent(in) :: fld1, fld2
real(kind=kind_real), intent(inout) :: zprod

integer :: jx,jy,jz,jj
real(kind=kind_real), allocatable :: zz(:)

call check_resolution(fld1, fld2)
if (fld1%nf /= fld2%nf .or. fld1%nl /= fld2%nl) then
  call abor1_ftn("qg_fields:field_prod error number of fields")
endif
!if (fld1%lbc .or. fld2%lbc) then
!  call abor1_ftn("qg_fields:field_prod should never dot_product full state")
!endif

allocate(zz(fld1%nl*fld1%nf))
zz(:)=0.0_kind_real
do jy=1,fld1%ny
do jx=1,fld1%nx
do jz=1,fld1%nl*fld1%nf
  zz(jz) = zz(jz) + fld1%gfld3d(jx,jy,jz) * fld2%gfld3d(jx,jy,jz)
enddo
enddo
enddo

zprod=0.0_kind_real
do jj=1,fld1%nl*fld1%nf
  zprod=zprod+zz(jj)
enddo
deallocate(zz)

return
end subroutine dot_prod

! ------------------------------------------------------------------------------

subroutine add_incr(self,rhs)
implicit none
type(qg_field), intent(inout) :: self
type(qg_field), intent(in)    :: rhs

call check(self)
call check(rhs)

if (self%nx==rhs%nx .and. self%ny==rhs%ny) then
  self%x(:,:,:) = self%x(:,:,:) + rhs%x(:,:,:)
else
  call abor1_ftn("qg_fields:add_incr: not coded for low res increment yet")
endif

if (self%nf>1) self%x(:,:,self%nl+1:) = 0.0_kind_real

return
end subroutine add_incr

! ------------------------------------------------------------------------------

subroutine diff_incr(lhs,x1,x2)
implicit none
type(qg_field), intent(inout) :: lhs
type(qg_field), intent(in)    :: x1
type(qg_field), intent(in)    :: x2

call check(lhs)
call check(x1)
call check(x2)

call zeros(lhs)
if (x1%nx==x2%nx .and. x1%ny==x2%ny) then
  if (lhs%nx==x1%nx .and. lhs%ny==x1%ny) then
    lhs%x(:,:,:) = x1%x(:,:,:) - x2%x(:,:,:)
  else
    call abor1_ftn("qg_fields:diff_incr: not coded for low res increment yet")
  endif
else
  call abor1_ftn("qg_fields:diff_incr: states not at same resolution")
endif

return
end subroutine diff_incr

! ------------------------------------------------------------------------------

subroutine change_resol(fld,rhs)
implicit none
type(qg_field), intent(inout) :: fld
type(qg_field), intent(in)    :: rhs
real(kind=kind_real), allocatable :: ztmp(:,:)
real(kind=kind_real) :: dy1, dy2, ya, yb, dx1, dx2, xa, xb
integer :: jx, jy, jf, iy, ia, ib

call check(fld)
call check(rhs)

if (fld%nx==rhs%nx .and. fld%ny==rhs%ny) then
  call copy(fld, rhs)
else
  call abor1_ftn("qg_fields:field_resol: untested code")

allocate(ztmp(rhs%nx,rhs%nl*rhs%nf))
dy1=1.0_kind_real/real(fld%ny,kind_real)
dy2=1.0_kind_real/real(rhs%ny,kind_real)
dx1=1.0_kind_real/real(fld%nx,kind_real)
dx2=1.0_kind_real/real(rhs%nx,kind_real)
do jy=1,fld%ny
! North-south interpolation
  if (jy==1 .or. jy==fld%ny) then
    if (jy==1) then
      iy=1
    else
      iy=rhs%ny
    endif
    do jx=1,rhs%nx
      do jf=1,rhs%nl*rhs%nf
        ztmp(jx,jf)=rhs%gfld3d(jx,iy,jf)
      enddo
    enddo
  else
    call lin_weights(jy,dy1,dy2,ia,ib,ya,yb)
    if (ib>rhs%ny) call abor1_ftn("qg_fields:field_resol: ib too large")
    do jx=1,rhs%nx
      do jf=1,rhs%nl*rhs%nf
        ztmp(jx,jf)=ya*rhs%gfld3d(jx,ia,jf)+yb*rhs%gfld3d(jx,ib,jf)
      enddo
    enddo
  endif
! East-west interpolation (periodic)
  do jx=1,fld%nx
    call lin_weights(jx,dx1,dx2,ia,ib,xa,xb)
    if (ib>rhs%nx) ib=ib-rhs%nx
    do jf=1,fld%nl*rhs%nf
      fld%gfld3d(jx,jy,jf)=xa*ztmp(jx,jf)+xb*ztmp(jx+1,jf)
    enddo
  enddo
enddo
deallocate(ztmp)

endif

return
end subroutine change_resol

! ------------------------------------------------------------------------------

subroutine read_file(fld, c_conf, vdate)
! Needs more interface clean-up here...
use iso_c_binding
use datetime_mod
use fckit_log_module, only : fckit_log

implicit none
type(qg_field), intent(inout) :: fld      !< Fields
type(c_ptr), intent(in)       :: c_conf   !< Configuration
type(datetime), intent(inout) :: vdate    !< DateTime

integer, parameter :: iunit=10
integer, parameter :: max_string_length=800 ! Yuk!
character(len=max_string_length+50) :: record
character(len=max_string_length) :: filename
character(len=20) :: sdate, fmtn
character(len=4)  :: cnx
character(len=11) :: fmt1='(X,ES24.16)'
character(len=1024)  :: buf
integer :: ic, iy, il, ix, is, jx, jy, jf, iread, nf
real(kind=kind_real), allocatable :: zz(:)

iread = 1
if (config_element_exists(c_conf,"read_from_file")) then
  iread = config_get_int(c_conf,"read_from_file")
endif
if (iread==0) then
  call fckit_log%warning("qg_fields:read_file: Inventing State")
  call invent_state(fld,c_conf)
  sdate = config_get_string(c_conf,len(sdate),"date")
  WRITE(buf,*) 'validity date is: '//sdate
  call fckit_log%info(buf)
  call datetime_set(sdate, vdate)
else
  call zeros(fld)
  filename = config_get_string(c_conf,len(filename),"filename")
  WRITE(buf,*) 'qg_field:read_file: opening '//filename
  call fckit_log%info(buf)
  open(unit=iunit, file=trim(filename), form='formatted', action='read')

  read(iunit,*) ix, iy, il, ic, is
  if (ix /= fld%nx .or. iy /= fld%ny .or. il /= fld%nl) then
    write (record,*) "qg_fields:read_file: ", &
                   & "input fields have wrong dimensions: ",ix,iy,il
    call fckit_log%error(record)
    write (record,*) "qg_fields:read_file: expected: ",fld%nx,fld%ny,fld%nl
    call fckit_log%error(record)
    call abor1_ftn("qg_fields:read_file: input fields have wrong dimensions")
  endif

  read(iunit,*) sdate
  WRITE(buf,*) 'validity date is: '//sdate
  call fckit_log%info(buf)
  call datetime_set(sdate, vdate)

  if (fld%nx>9999)  call abor1_ftn("Format too small")
  write(cnx,'(I4)')fld%nx
  fmtn='('//trim(cnx)//fmt1//')'

  nf = min(fld%nf, ic)
  do jf=1,il*nf
    do jy=1,fld%ny
      read(iunit,fmtn) (fld%gfld3d(jx,jy,jf), jx=1,fld%nx)
    enddo
  enddo
! Skip un-necessary data from file if any
  allocate(zz(fld%nx))
  do jf=nf*il+1, ic*il
    do jy=1,fld%ny
      read(iunit,fmtn) (zz(jx), jx=1,fld%nx)
    enddo
  enddo
  deallocate(zz)

  if (fld%lbc) then
    do jf=1,4
      read(iunit,fmt1) fld%xbound(jf)
    enddo
    do jf=1,4
      read(iunit,fmtn) (fld%qbound(jx,jf), jx=1,fld%nx)
    enddo
  endif

  close(iunit)
endif

call check(fld)

return
end subroutine read_file

! ------------------------------------------------------------------------------

subroutine write_file(fld, c_conf, vdate)
use iso_c_binding
use datetime_mod
use fckit_log_module, only : fckit_log

implicit none
type(qg_field), intent(in) :: fld    !< Fields
type(c_ptr), intent(in)    :: c_conf !< Configuration
type(datetime), intent(in) :: vdate  !< DateTime

integer, parameter :: iunit=11
integer, parameter :: max_string_length=800 ! Yuk!
character(len=max_string_length+50) :: record
character(len=max_string_length) :: filename
character(len=20) :: sdate, fmtn
character(len=4)  :: cnx
character(len=11) :: fmt1='(X,ES24.16)'
character(len=1024):: buf
integer :: jf, jy, jx, is

call check(fld)

filename = genfilename(c_conf,max_string_length,vdate)
WRITE(buf,*) 'qg_field:write_file: writing '//filename
call fckit_log%info(buf)
open(unit=iunit, file=trim(filename), form='formatted', action='write')

is=0
if (fld%lbc) is=1

write(iunit,*) fld%nx, fld%ny, fld%nl, fld%nf, is

call datetime_to_string(vdate, sdate)
write(iunit,*) sdate

if (fld%nx>9999)  call abor1_ftn("Format too small")
write(cnx,'(I4)')fld%nx
fmtn='('//trim(cnx)//fmt1//')'

do jf=1,fld%nl*fld%nf
  do jy=1,fld%ny
    write(iunit,fmtn) (fld%gfld3d(jx,jy,jf), jx=1,fld%nx)
  enddo
enddo

if (fld%lbc) then
  do jf=1,4
    write(iunit,fmt1) fld%xbound(jf)
  enddo
  do jf=1,4
    write(iunit,fmtn) (fld%qbound(jx,jf), jx=1,fld%nx)
  enddo
endif

close(iunit)

return
end subroutine write_file

! ------------------------------------------------------------------------------

subroutine gpnorm(fld, nf, pstat)
implicit none
type(qg_field), intent(in) :: fld
integer, intent(in) :: nf
real(kind=kind_real), intent(inout) :: pstat(3, nf)
integer :: jj,joff

call check(fld)

do jj=1,fld%nf
  joff=(jj-1)*fld%nl
  pstat(1,jj)=minval(fld%gfld3d(:,:,joff+1:joff+fld%nl))
  pstat(2,jj)=maxval(fld%gfld3d(:,:,joff+1:joff+fld%nl))
  pstat(3,jj)=sqrt(sum(fld%gfld3d(:,:,joff+1:joff+fld%nl)**2) &
               & /real(fld%nl*fld%nx*fld%ny,kind_real))
enddo
jj=jj-1

if (fld%lbc) then
  jj=jj+1
  pstat(1,jj)=minval(fld%xbound(:))
  pstat(2,jj)=maxval(fld%xbound(:))
  pstat(3,jj)=sqrt(sum(fld%xbound(:)**2)/real(4,kind_real))

  jj=jj+1
  pstat(1,jj)=minval(fld%qbound(:,:))
  pstat(2,jj)=maxval(fld%qbound(:,:))
  pstat(3,jj)=sqrt(sum(fld%qbound(:,:)**2)/real(4*fld%nx,kind_real))
endif

if (jj /= nf) call abor1_ftn("qg_fields_gpnorm: error number of fields")

return
end subroutine gpnorm

! ------------------------------------------------------------------------------

subroutine fldrms(fld, prms)
implicit none
type(qg_field), intent(in) :: fld
real(kind=kind_real), intent(out) :: prms
integer :: jf,jy,jx,ii
real(kind=kind_real) :: zz

call check(fld)

zz = 0.0

do jf=1,fld%nl*fld%nf
  do jy=1,fld%ny
    do jx=1,fld%nx
      zz = zz + fld%gfld3d(jx,jy,jf)*fld%gfld3d(jx,jy,jf)
    enddo
  enddo
enddo

ii = fld%nl*fld%nf*fld%ny*fld%nx

if (fld%lbc) then
  do jf=1,4
    zz = zz + fld%xbound(jf)*fld%xbound(jf)
  enddo
  ii = ii + 4
  do jf=1,4
    do jx=1,fld%nx
     zz = zz + fld%qbound(jx,jf)*fld%qbound(jx,jf)
    enddo
  enddo
  ii = ii + 4*fld%nx
endif

prms = sqrt(zz/real(ii,kind_real))

end subroutine fldrms

! ------------------------------------------------------------------------------

subroutine interp_tl(fld, locs, gom)
implicit none
type(qg_field), intent(in)   :: fld
type(qg_locs), intent(in)    :: locs
type(qg_goms), intent(inout) :: gom
real(kind=kind_real) :: di, dj, ai, aj, valb, valt
integer :: jloc,joff,jvar,ii,jj,kk,iright,jsearch
character(len=1) :: cvar

call check(fld)

do jloc=1,locs%nloc
! Convert horizontal coordinate stored in ODB to grid locations
  di=locs%xyz(1,jloc)*real(fld%nx,kind_real)
  ii=int(di)
  ai=di-real(ii,kind_real);
  ii=ii+1
  if (ii<1.or.ii>fld%nx) call abor1_ftn("qg_fields_interp_tl: error ii")
  iright=ii+1
  if (iright==fld%nx+1) iright=1

  dj=locs%xyz(2,jloc)*real(fld%ny,kind_real)
  jj=int(dj)
  aj=dj-real(jj,kind_real);
  jj=jj+1
  if (jj<1.or.jj>fld%ny) call abor1_ftn("qg_fields_interp_tl: error jj")

  kk=nint(locs%xyz(3,jloc))

! Loop on required variables
  gom%used=gom%used+1
  do jvar=1,gom%nvar
    cvar=gom%variables(jvar)
    joff=-1
    do jsearch=1,fld%nf
      if (fld%fldnames(jsearch)==cvar) joff=jsearch-1
    enddo
    if (joff<0) call abor1_ftn("qg_fields_interp_tl: unknown variable")
    joff=fld%nl*joff

!   Interpolate east-west
    if (jj>1) then
      valb=(1.0_kind_real-ai)*fld%gfld3d(ii    ,jj,joff+kk) &
                       & +ai *fld%gfld3d(iright,jj,joff+kk)
    else
      select case (cvar)
      case ("x")
        valb=0.0_kind_real
      case ("u")
        valb=2.0_kind_real*( (1.0_kind_real-ai)*fld%gfld3d(ii,1,joff+kk) &
             &                             +ai *fld%gfld3d(iright,1,joff+kk) ) &
             & -           ( (1.0_kind_real-ai)*fld%gfld3d(ii,2,joff+kk) &
             &                             +ai *fld%gfld3d(iright,2,joff+kk) )
      case ("v")
        valb=0.0_kind_real
      end select
    endif

    if (jj<fld%ny) then
      valt=(1.0_kind_real-ai)*fld%gfld3d(ii    ,jj+1,joff+kk) &
                       & +ai *fld%gfld3d(iright,jj+1,joff+kk)
    else
      select case (cvar)
      case ("x")
        valt = 0.0_kind_real
      case ("u")
        valt=2.0_kind_real*( (1.0_kind_real-ai)*fld%gfld3d(ii,fld%ny,joff+kk) &
             &                             +ai *fld%gfld3d(iright,fld%ny,joff+kk) ) &
             & -           ( (1.0_kind_real-ai)*fld%gfld3d(ii,fld%ny-1,joff+kk) &
             &                             +ai *fld%gfld3d(iright,fld%ny-1,joff+kk) )
      case ("v")
        valt = 0.0_kind_real
      end select
    endif

!   Interpolate north-south
    gom%values(jvar,gom%used) = (1.0_kind_real-aj)*valb + aj*valt
  enddo
enddo

return
end subroutine interp_tl

! ------------------------------------------------------------------------------

subroutine interp_ad(fld, locs, gom)
implicit none
type(qg_field), intent(inout) :: fld
type(qg_locs), intent(in)     :: locs
type(qg_goms), intent(inout)  :: gom
real(kind=kind_real) :: di, dj, ai, aj, valb, valt
integer :: jloc,joff,jvar,ii,jj,kk,iright,jsearch
character(len=1) :: cvar

call check(fld)

do jloc=locs%nloc,1,-1
! Convert horizontal coordinate stored in ODB to grid locations
  di=locs%xyz(1,jloc)*real(fld%nx,kind_real)
  ii=int(di)
  ai=di-real(ii,kind_real);
  ii=ii+1
  if (ii<1.or.ii>fld%nx) call abor1_ftn("qg_fields_interp_ad: error ii")
  iright=ii+1
  if (iright==fld%nx+1) iright=1

  dj=locs%xyz(2,jloc)*real(fld%ny,kind_real)
  jj=int(dj)
  aj=dj-real(jj,kind_real);
  jj=jj+1
  if (jj<1.or.jj>fld%ny) call abor1_ftn("qg_fields_interp_ad: error jj")

  kk=nint(locs%xyz(3,jloc))

! Loop on required variables
  gom%used=gom%used+1
  do jvar=gom%nvar,1,-1
    cvar=gom%variables(jvar)
    joff=-1
    do jsearch=1,fld%nf
      if (fld%fldnames(jsearch)==cvar) joff=jsearch-1
    enddo
    if (joff<0) call abor1_ftn("qg_fields_interp_ad: unknown variable")
    joff=fld%nl*joff

!   Interpolate north-south
    valb = (1.0_kind_real-aj)*gom%values(jvar,gom%nobs-gom%used+1)
    valt = aj*gom%values(jvar,gom%nobs-gom%used+1)

!   Interpolate east-west
    if (jj>1) then
      fld%gfld3d(ii    ,jj,joff+kk) = &
             & fld%gfld3d(ii    ,jj,joff+kk)+(1.0_kind_real-ai)*valb
      fld%gfld3d(iright,jj,joff+kk) = &
             & fld%gfld3d(iright,jj,joff+kk)+ai*valb
    else
      select case (cvar)
      case ("u")
        fld%gfld3d(ii,1,joff+kk) = &
             & fld%gfld3d(ii,1,joff+kk)+2.0_kind_real*(1.0_kind_real-ai)*valb
        fld%gfld3d(iright,1,joff+kk) = &
             & fld%gfld3d(iright,1,joff+kk)+2.0_kind_real*ai*valb
        fld%gfld3d(ii,2,joff+kk) = &
             & fld%gfld3d(ii,2,joff+kk)-(1.0_kind_real-ai)*valb
        fld%gfld3d(iright,2,joff+kk) = &
             & fld%gfld3d(iright,2,joff+kk)-ai*valb
      end select
    endif

    if (jj<fld%ny) then
      fld%gfld3d(ii    ,jj+1,joff+kk) = &
             & fld%gfld3d(ii    ,jj+1,joff+kk)+(1.0_kind_real-ai)*valt
      fld%gfld3d(iright,jj+1,joff+kk) = &
             & fld%gfld3d(iright,jj+1,joff+kk)+ai*valt
    else
      select case (cvar)
      case ("u")
        fld%gfld3d(ii,fld%ny,joff+kk) = &
             & fld%gfld3d(ii,fld%ny,joff+kk)+2.0_kind_real*(1.0_kind_real-ai)*valt
        fld%gfld3d(iright,fld%ny,joff+kk) = &
             & fld%gfld3d(iright,fld%ny,joff+kk)+2.0_kind_real*ai*valt
        fld%gfld3d(ii,fld%ny-1,joff+kk) = &
             & fld%gfld3d(ii,fld%ny-1,joff+kk)-(1.0_kind_real-ai)*valt
        fld%gfld3d(iright,fld%ny-1,joff+kk)= &
             & fld%gfld3d(iright,fld%ny-1,joff+kk)-ai*valt
      end select
    endif
  enddo
enddo

return
end subroutine interp_ad

! ------------------------------------------------------------------------------

subroutine lin_weights(kk,delta1,delta2,k1,k2,w1,w2)
implicit none
integer, intent(in)  :: kk
real(kind=kind_real), intent(in)     :: delta1,delta2
integer, intent(out) :: k1,k2
real(kind=kind_real), intent(out)    :: w1,w2

integer :: ii
real(kind=kind_real) :: zz

zz=real(kk-1,kind_real)*delta1
zz=zz/delta2
ii=int(zz)
w1=zz-real(ii,kind_real)
w2=1.0_kind_real-w1
k1=ii+1
k2=ii+2

return
end subroutine lin_weights

! ------------------------------------------------------------------------------

subroutine convert_to_ug(self, ug)
use unstructured_grid_mod
implicit none
type(qg_field), intent(in) :: self
type(unstructured_grid), intent(inout) :: ug
real(kind=kind_real) :: zz(2)
integer :: jx,jy,cmask(2)

zz(1) = 0.0
zz(2) = 1.0
call create_unstructured_grid(ug, 2, zz)

cmask(1) = 1
cmask(2) = 1
do jy=1,self%ny
  do jx=1,self%nx
    call add_column(ug, self%lats(jy), self%lons(jx), 2, 1, 0, cmask, 1)
    ug%last%column%cols(1) = self%x(jx,jy,1)
    ug%last%column%cols(2) = self%x(jx,jy,2)
  enddo
enddo

end subroutine convert_to_ug

! ------------------------------------------------------------------------------

subroutine convert_from_ug(self, ug)
use unstructured_grid_mod
implicit none
type(qg_field), intent(inout) :: self
type(unstructured_grid), intent(in) :: ug
type(column_element), pointer :: current
integer :: jx,jy

current => ug%head
do jy=1,self%ny
  do jx=1,self%nx
    self%x(jx,jy,1) = current%column%cols(1)
    self%x(jx,jy,2) = current%column%cols(2)
    current => current%next
  enddo
enddo

end subroutine convert_from_ug

! ------------------------------------------------------------------------------

function genfilename (c_conf,length,vdate)
use iso_c_binding
use datetime_mod
use duration_mod
type(c_ptr), intent(in)    :: c_conf  !< Configuration
integer, intent(in) :: length
character(len=length) :: genfilename
type(datetime), intent(in) :: vdate

character(len=length) :: fdbdir, expver, typ, validitydate, referencedate, sstep, &
                       & prefix, mmb
type(datetime) :: rdate
type(duration) :: step
integer lenfn

! here we should query the length and then allocate "string".
! But Fortran 90 does not allow variable-length allocatable strings.
! config_get_string checks the string length and aborts if too short.
fdbdir = config_get_string(c_conf,len(fdbdir),"datadir")
expver = config_get_string(c_conf,len(expver),"exp")
typ    = config_get_string(c_conf,len(typ)   ,"type")

if (typ=="ens") then
  mmb = config_get_string(c_conf, len(mmb), "member")
  lenfn = LEN_TRIM(fdbdir) + 1 + LEN_TRIM(expver) + 1 + LEN_TRIM(typ) + 1 + LEN_TRIM(mmb)
  prefix = TRIM(fdbdir) // "/" // TRIM(expver) // "." // TRIM(typ) // "." // TRIM(mmb)
else
  lenfn = LEN_TRIM(fdbdir) + 1 + LEN_TRIM(expver) + 1 + LEN_TRIM(typ)
  prefix = TRIM(fdbdir) // "/" // TRIM(expver) // "." // TRIM(typ)
endif

if (typ=="fc" .or. typ=="ens") then
  referencedate = config_get_string(c_conf,len(referencedate),"date")
  call datetime_to_string(vdate, validitydate)
  call datetime_create(TRIM(referencedate), rdate)
  call datetime_diff(vdate, rdate, step)
  call duration_to_string(step, sstep)
  lenfn = lenfn + 1 + LEN_TRIM(referencedate) + 1 + LEN_TRIM(sstep)
  genfilename = TRIM(prefix) // "." // TRIM(referencedate) // "." // TRIM(sstep)
endif

if (typ=="an") then
  call datetime_to_string(vdate, validitydate)
  lenfn = lenfn + 1 + LEN_TRIM(validitydate)
  genfilename = TRIM(prefix) // "." // TRIM(validitydate)
endif

if (lenfn>length) &
  & call abor1_ftn("qg_fields:genfilename: filename too long")

end function genfilename

! ------------------------------------------------------------------------------

function common_vars(x1, x2)

implicit none
type(qg_field), intent(in) :: x1, x2
integer :: common_vars
integer :: jf

! We assume here that one set of fields is a subset of the other,
! that fields are always in the same order starting with x,
! and that the common fields are the first ones.

common_vars = min(x1%nf, x2%nf)
do jf = 1, common_vars
  if (x1%fldnames(jf)/=x2%fldnames(jf)) &
    & call abor1_ftn("common_vars: fields do not match")
enddo
if (x1%nl /= x2%nl) call abor1_ftn("common_vars: error number of levels")
common_vars = x1%nl * common_vars

end function common_vars

! ------------------------------------------------------------------------------

subroutine check_resolution(x1, x2)

implicit none
type(qg_field), intent(in) :: x1, x2

if (x1%nx /= x2%nx .or.  x1%ny /= x2%ny .or.  x1%nl /= x2%nl) then
  call abor1_ftn ("qg_fields: resolution error")
endif
call check(x1)
call check(x2)

end subroutine check_resolution

! ------------------------------------------------------------------------------

subroutine check(self)
implicit none
type(qg_field), intent(in) :: self
logical :: bad

bad = .not.associated(self%gfld3d)

bad = bad .or. (size(self%gfld3d, 1) /= self%nx)
bad = bad .or. (size(self%gfld3d, 2) /= self%ny)
bad = bad .or. (size(self%gfld3d, 3) /= self%nl*self%nf)

bad = bad .or. .not.associated(self%x)

if (self%nf>1) then
  bad = bad .or. .not.associated(self%q)
  bad = bad .or. .not.associated(self%u)
  bad = bad .or. .not.associated(self%v)
else
  bad = bad .or. associated(self%q)
  bad = bad .or. associated(self%u)
  bad = bad .or. associated(self%v)
endif

!allocate(self%fldnames(self%nf))

if (self%lbc) then
  bad = bad .or. .not.associated(self%xbound)
  bad = bad .or. (size(self%xbound) /= 4)

  bad = bad .or. .not.associated(self%qbound)
  bad = bad .or. (size(self%qbound, 1) /= self%nx)
  bad = bad .or. (size(self%qbound, 2) /= 4)
else
  bad = bad .or. associated(self%xbound)
  bad = bad .or. associated(self%qbound)
endif

if (bad) then
  write(0,*)'nx, ny, nf, nl, lbc = ',self%nx,self%ny,self%nf,self%nl,self%lbc
  if (associated(self%gfld3d)) write(0,*)'shape(gfld3d) = ',shape(self%gfld3d)
  if (associated(self%xbound)) write(0,*)'shape(xbound) = ',shape(self%xbound)
  if (associated(self%qbound)) write(0,*)'shape(qbound) = ',shape(self%qbound)
  call abor1_ftn ("qg_fields: field not consistent")
endif

end subroutine check

! ------------------------------------------------------------------------------

end module qg_fields
