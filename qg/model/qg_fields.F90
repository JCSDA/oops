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
use calculate_pv, only : calc_pv
use netcdf
use tools_nc, only: ncerr
use kinds

implicit none
private

public :: qg_field, &
        & create, delete, zeros, random, copy, &
        & self_add, self_schur, self_sub, self_mul, axpy, &
        & dot_prod, add_incr, diff_incr, &
        & read_file, write_file, gpnorm, fldrms, &
        & change_resol, interp_tl, interp_ad, &
        & analytic_init
public :: qg_field_registry

! ------------------------------------------------------------------------------

!> Fortran derived type to hold QG fields
type :: qg_field
  type(qg_geom), pointer :: geom                    !< Geometry
  integer :: nl                                     !< Number of levels
  integer :: nf                                     !< Number of fields
  logical :: lbc                                    !< North-South boundary is present
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

logical :: netcdfio = .true.

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
type(qg_geom), target, intent(in) :: geom
type(qg_vars),  intent(in)    :: vars
integer :: ioff

self%geom => geom
self%nl = 2
self%nf = vars%nv
self%lbc = vars%lbc

allocate(self%gfld3d(self%geom%nx,self%geom%ny,self%nl*self%nf))
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
  allocate(self%qbound(self%geom%nx,4))
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
do jy=1,fld1%geom%ny
do jx=1,fld1%geom%nx
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

if (self%geom%nx==rhs%geom%nx .and. self%geom%ny==rhs%geom%ny) then
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
if (x1%geom%nx==x2%geom%nx .and. x1%geom%ny==x2%geom%ny) then
  if (lhs%geom%nx==x1%geom%nx .and. lhs%geom%ny==x1%geom%ny) then
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

if (fld%geom%nx==rhs%geom%nx .and. fld%geom%ny==rhs%geom%ny) then
  call copy(fld, rhs)
else
  call abor1_ftn("qg_fields:field_resol: untested code")

allocate(ztmp(rhs%geom%nx,rhs%nl*rhs%nf))
dy1=1.0_kind_real/real(fld%geom%ny,kind_real)
dy2=1.0_kind_real/real(rhs%geom%ny,kind_real)
dx1=1.0_kind_real/real(fld%geom%nx,kind_real)
dx2=1.0_kind_real/real(rhs%geom%nx,kind_real)
do jy=1,fld%geom%ny
! North-south interpolation
  if (jy==1 .or. jy==fld%geom%ny) then
    if (jy==1) then
      iy=1
    else
      iy=rhs%geom%ny
    endif
    do jx=1,rhs%geom%nx
      do jf=1,rhs%nl*rhs%nf
        ztmp(jx,jf)=rhs%gfld3d(jx,iy,jf)
      enddo
    enddo
  else
    call lin_weights(jy,dy1,dy2,ia,ib,ya,yb)
    if (ib>rhs%geom%ny) call abor1_ftn("qg_fields:field_resol: ib too large")
    do jx=1,rhs%geom%nx
      do jf=1,rhs%nl*rhs%nf
        ztmp(jx,jf)=ya*rhs%gfld3d(jx,ia,jf)+yb*rhs%gfld3d(jx,ib,jf)
      enddo
    enddo
  endif
! East-west interpolation (periodic)
  do jx=1,fld%geom%nx
    call lin_weights(jx,dx1,dx2,ia,ib,xa,xb)
    if (ib>rhs%geom%nx) ib=ib-rhs%geom%nx
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
integer :: ncid, nx_id, ny_id, nl_id, nc_id, gfld3d_id, xbound_id, qbound_id
real(kind=kind_real), allocatable :: zz(:)
character(len=1024) :: subr='read_file'

iread = 1
if (config_element_exists(c_conf,"read_from_file")) then
  iread = config_get_int(c_conf,"read_from_file")
endif
if (iread==0) then
  call fckit_log%warning("qg_fields:read_file: Inventing State")
  call invent_state(fld,c_conf)
  sdate = config_get_string(c_conf,len(sdate),"date")
  write(buf,*) 'validity date is: '//sdate
  call fckit_log%info(buf)
  call datetime_set(sdate, vdate)
else
  call zeros(fld)
  filename = config_get_string(c_conf,len(filename),"filename")
  write(buf,*) 'qg_field:read_file: opening '//filename
  call fckit_log%info(buf)

  ! Read data
  if (netcdfio) then

    call ncerr(subr,nf90_open(trim(filename)//'.nc',nf90_nowrite,ncid))
    call ncerr(subr,nf90_inq_dimid(ncid,'nx',nx_id))
    call ncerr(subr,nf90_inq_dimid(ncid,'ny',ny_id))
    call ncerr(subr,nf90_inq_dimid(ncid,'nl',nl_id))
    call ncerr(subr,nf90_inq_dimid(ncid,'nc',nc_id))
    call ncerr(subr,nf90_inquire_dimension(ncid,nx_id,len=ix))
    call ncerr(subr,nf90_inquire_dimension(ncid,ny_id,len=iy))
    call ncerr(subr,nf90_inquire_dimension(ncid,nl_id,len=il))
    call ncerr(subr,nf90_inquire_dimension(ncid,nl_id,len=ic))
    if (ix /= fld%geom%nx .or. iy /= fld%geom%ny .or. il /= fld%nl) then
      write (record,*) "qg_fields:read_file: ", &
                     & "input fields have wrong dimensions: ",ix,iy,il
      call fckit_log%error(record)
      write (record,*) "qg_fields:read_file: expected: ",fld%geom%nx,fld%geom%ny,fld%nl
      call fckit_log%error(record)
      call abor1_ftn("qg_fields:read_file: input fields have wrong dimensions")
    endif
    call ncerr(subr,nf90_get_att(ncid,nf90_global,'lbc',is))
  
    call ncerr(subr,nf90_get_att(ncid,nf90_global,'sdate',sdate))
    write(buf,*) 'validity date is: '//sdate
    call fckit_log%info(buf)
  
    nf = min(fld%nf, ic)
    call ncerr(subr,nf90_inq_varid(ncid,'gfld3d',gfld3d_id))
    do jf=1,il*nf
      call ncerr(subr,nf90_get_var(ncid,gfld3d_id,fld%gfld3d(:,:,jf),(/1,1,jf/),(/ix,iy,1/)))
    enddo
  
    if (fld%lbc) then
      call ncerr(subr,nf90_inq_varid(ncid,'xbound',xbound_id))
      call ncerr(subr,nf90_get_var(ncid,xbound_id,fld%xbound))
      call ncerr(subr,nf90_inq_varid(ncid,'qbound',qbound_id))
      call ncerr(subr,nf90_get_var(ncid,qbound_id,fld%qbound))
    endif
  
    call ncerr(subr,nf90_close(ncid))

  else

    open(unit=iunit, file=trim(filename), form='formatted', action='read')
  
    read(iunit,*) ix, iy, il, ic, is
    if (ix /= fld%geom%nx .or. iy /= fld%geom%ny .or. il /= fld%nl) then
      write (record,*) "qg_fields:read_file: ", &
                     & "input fields have wrong dimensions: ",ix,iy,il
      call fckit_log%error(record)
      write (record,*) "qg_fields:read_file: expected: ",fld%geom%nx,fld%geom%ny,fld%nl
      call fckit_log%error(record)
      call abor1_ftn("qg_fields:read_file: input fields have wrong dimensions")
    endif
  
    read(iunit,*) sdate
    write(buf,*) 'validity date is: '//sdate
    call fckit_log%info(buf)
    call datetime_set(sdate, vdate)
  
    if (fld%geom%nx>9999)  call abor1_ftn("Format too small")
    write(cnx,'(I4)')fld%geom%nx
    fmtn='('//trim(cnx)//fmt1//')'
  
    nf = min(fld%nf, ic)
    do jf=1,il*nf
      do jy=1,fld%geom%ny
        read(iunit,fmtn) (fld%gfld3d(jx,jy,jf), jx=1,fld%geom%nx)
      enddo
    enddo
! Skip un-necessary data from file if any
    allocate(zz(fld%geom%nx))
    do jf=nf*il+1, ic*il
      do jy=1,fld%geom%ny
        read(iunit,fmtn) (zz(jx), jx=1,fld%geom%nx)
      enddo
    enddo
    deallocate(zz)
  
    if (fld%lbc) then
      do jf=1,4
        read(iunit,fmt1) fld%xbound(jf)
      enddo
      do jf=1,4
        read(iunit,fmtn) (fld%qbound(jx,jf), jx=1,fld%geom%nx)
      enddo
    endif
  
    close(iunit)

  endif

  ! Set date
  call datetime_set(sdate, vdate)
endif

call check(fld)

return
end subroutine read_file

! ------------------------------------------------------------------------------
!> Analytic Initialization for the QG model
!!
!! \details **analytic_init()** initializes the qg field using one of several alternative idealized analytic models.
!! This is intended to faciltitate testing by eliminating the need to read in the initial state
!! from a file and by providing exact expressions to test interpolations.
!! This function is activated by setting the initial.analytic_init field in the configuration file.
!!
!! Initialization options that begin with "dcmip" refer to tests defined by the multi-institutional
!! 2012 [Dynamical Core Intercomparison Project](https://earthsystealcmcog.org/projects/dcmip-2012)
!! and the associated Summer School, sponsored by NOAA, NSF, DOE, NCAR, and the University of Michigan.
!!
!! Currently implemented options for analytic_init include:
!! * baroclinic-instability: A two-layer flow that is baroclinically unstable in the presence of non-zero orography
!! * dcmip-test-1-1: 3D deformational flow
!! * dcmip-test-1-2: 3D Hadley-like meridional circulation
!!
!! Relevant parameters specified in the **State.StateGenerate** section of the config file include:
!! * **analytic_init** the type of initialization as described in analytic_init())
!! * **variables** currently must be set to ["x","q","u","v"]
!! * **date** a reference date and time for the state
!! * **top_layer_depth**, **bottom_layer_depth** are needed in order to convert the normalized height
!!     coordinate into meters so it can be used in the analytic formulae
!! 
!! \author M. Miesch (JCSDA, adapted from a pre-existing call to invent_state)
!! \date March 7, 2018: Created
!!
!! \warning If you use the dcmip routines, be sure to include u and v in the variables list.
!!
!! \warning For dcmip initial conditions, this routine computes the streamfunction x and the potential
!! vorticity q from the horizontal winds u and v.  Since q (or alternatively x) is the state variable
!! for the evolution of the QG model, only the non-divergent (vortical) component of the horizontal
!! flow is captured.  In other words, u and v are projected onto a flow field with zero horizontal
!! divergence.  Furthermore, this projection neglects geometric curvature factors.  This projection
!! is used to calculate the streamfunction X and the potential vorticity q.  However, on exiting this
!! function, the u and v components of the qg_field structure are still set to their original
!! (non-projected) values.  This is for use with the State interpolation test.  The winds u and v
!! are re-calculated from X before doing a forecast by the routine c_qg_prepare_integration() in
!! c_qg_model.F90.
!!
!! \warning For all but the baroclinic-instability option, the PV here is computed assuming zero
!! orography.  For particular model configurations, the orography is defined in c_qg_setup()
!! [c_qg_model.F90] and the PV is recalculated in c_qg_prepare_integration() [c_qg_model.F90]
!! before the integration of the model commences.
!!

subroutine analytic_init(fld, geom, config, vdate)

  use iso_c_binding
  use datetime_mod
  use fckit_log_module, only : fckit_log
  use dcmip_initial_conditions_test_1_2_3, only : test1_advection_deformation, test1_advection_hadley
  use qg_constants

  implicit none
  type(qg_field), intent(inout) :: fld      !< Fields
  type(qg_geom), intent(inout)  :: geom     !< Grid information
  type(c_ptr), intent(in)       :: config   !< Configuration
  type(datetime), intent(inout) :: vdate    !< DateTime

  character(len=30) :: ic
  character(len=20) :: sdate
  character(len=1024)  :: buf
  real(kind=kind_real) :: d1, d2, height(2), z, f1, f2, xdum(2)
  real(kind=kind_real) :: deltax,deltay
  real(kind=kind_real) :: p0,u0,v0,w0,t0,phis0,ps0,rho0,hum0,q1,q2,q3,q4
  real(kind=kind_real), allocatable :: rs(:,:)
  integer :: ix, iy, il

  if (config_element_exists(config,"analytic_init")) then
     ic = trim(config_get_string(config,len(ic),"analytic_init"))
  else
     ! this is mainly for backward compatibility
     ic = "baroclinic-instability"
  endif

  call fckit_log%warning("qg_fields:analytic_init: "//IC)
  sdate = config_get_string(config,len(sdate),"date")
  write(buf,*) 'validity date is: '//sdate
  call fckit_log%info(buf)
  call datetime_set(sdate,vdate)
    
  ! To use the DCMIP routines we need a vertical coordinate
  ! assume 2 levels for now but easily to generalize as needed
  d1  = config_get_real(config,"top_layer_depth")
  d2  = config_get_real(config,"bottom_layer_depth")

  height(2) = 0.5_kind_real*d2
  height(1) = d2 + 0.5_kind_real*d1
  xdum(:) = 0.0_kind_real

  deltax = domain_zonal/(scale_length*real(geom%nx,kind_real))
  deltay = domain_meridional/(scale_length*real(geom%ny,kind_real))

  f1 = f0*f0*scale_length*scale_length/(g*dlogtheta*d1)
  f2 = f0*f0*scale_length*scale_length/(g*dlogtheta*d2)

  allocate(rs(geom%nx,geom%ny))
  rs = 0.0_kind_real
  
  !--------------------------------------
  init_option: select case (ic)

     case ("baroclinic-instability")
        call invent_state(fld,config)   

     case ("dcmip-test-1-1")

        do il = 1, fld%nl
           z = height(il)
           do iy = 1,fld%geom%ny ; do ix = 1, fld%geom%nx 

              call test1_advection_deformation(geom%lon(ix),geom%lat(iy),p0,z,1,&
                   u0,v0,w0,t0,phis0,ps0,rho0,hum0,q1,q2,q3,q4)
              
              fld%u(ix,iy,il) = u0 / ubar
              fld%v(ix,iy,il) = v0 / ubar
              
           enddo; enddo
        enddo

        call calc_streamfunction(geom%nx,geom%ny,fld%x,fld%u,fld%v,deltax,deltay)

        ! Some compilers don't tolerate passing a null pointer
        if (fld%lbc) then
           call calc_pv(geom%nx,geom%ny,fld%q,fld%x,fld%x_north,fld%x_south,&
                        f1,f2,deltax,deltay,bet,rs,fld%lbc)
        else
           call calc_pv(geom%nx,geom%ny,fld%q,fld%x,xdum,xdum,&
                        f1,f2,deltax,deltay,bet,rs,fld%lbc)
        endif
             
     case ("dcmip-test-1-2")

        do il = 1, fld%nl
           z = height(il)
           do iy = 1,fld%geom%ny ; do ix = 1, fld%geom%nx 

              call test1_advection_hadley(geom%lon(ix),geom%lat(iy),p0,z,1,&
                   u0,v0,w0,t0,phis0,ps0,rho0,hum0,q1)
           
              fld%u(ix,iy,il) = u0 / ubar
              fld%v(ix,iy,il) = v0 / ubar
              
           enddo; enddo
        enddo

        call calc_streamfunction(geom%nx,geom%ny,fld%x,fld%u,fld%v,deltax,deltay)

        ! Some compilers don't tolerate passing a null pointer
        if (fld%lbc) then
           call calc_pv(geom%nx,geom%ny,fld%q,fld%x,fld%x_north,fld%x_south,&
                        f1,f2,deltax,deltay,bet,rs,fld%lbc)
        else
           call calc_pv(geom%nx,geom%ny,fld%q,fld%x,xdum,xdum,&
                        f1,f2,deltax,deltay,bet,rs,fld%lbc)
        endif

     case default

        call abor1_ftn ("qg_fields:analytic_init not properly defined")
        
  end select init_option
  
  !--------------------------------------

  call check(fld)

  return
  
end subroutine analytic_init

! ------------------------------------------------------------------------------
!> Calculate Streamfunction from u and v
!!
!! \details **calc_streamfunction()** computes the streamfunction x from the velocities u and v.
!! It was originally developed to support the dcmip initial conditions in analytic_init()
!! but is potentially of broader applicability.
!! The algorithm first computes the vertical vorticity and then inverts the
!! resulting Laplacian for x.
!!
!! \author M. Miesch 
!! \date March 7, 2018: Created
!!
!! \warning This routine currently does not take into account latitudinal boundary conditions on x.  

subroutine calc_streamfunction(kx,ky,x,u,v,dx,dy)

  integer, intent(in) :: kx           !< Zonal grid dimension
  integer, intent(in) :: ky           !< Meridional grid dimension
  real(kind=kind_real), intent(out)   :: x(kx,ky,2)   !< Streamfunction
  real(kind=kind_real), intent(in)    :: u(kx,ky,2)   !< zonal velocity (non-dimensional)
  real(kind=kind_real), intent(in)    :: v(kx,ky,2)   !< meridional velocity (non-dimensional)
  real(kind=kind_real), intent(in)    :: dx           !< Zonal grid spacing (non-dimensional)
  real(kind=kind_real), intent(in)    :: dy           !< Meridional grid spacing (non-dimensional)
  
  real(kind=kind_real), allocatable :: vort(:,:), psi(:,:)
  real(kind=kind_real), allocatable :: dudy(:), dvdx(:)
  integer :: ix, iy, il

  allocate(vort(kx,ky),psi(kx,ky))
  allocate(dudy(ky),dvdx(kx))

  ! assume two levels - easy to generalize if needed
  do il = 1, 2

     ! compute vertical vorticity from u and v
     do iy=1,ky
        call deriv_1D(dvdx,v(:,iy,il),kx,dx,periodic=.true.)
        vort(:,iy) = dvdx
     enddo

     do ix=1,kx
        call deriv_1D(dudy,u(ix,:,il),ky,dy)
        vort(ix,:) = vort(ix,:) - dudy
     enddo
        
     ! invert the Laplacian to get the streamfunction
     call solve_helmholz(psi,vort,0.0_kind_real,kx,ky,dx,dy)
     x(:,:,il) = psi
        
  enddo

  deallocate(vort,psi,dudy,dvdx)

end subroutine calc_streamfunction

! ------------------------------------------------------------------------------
!> Finite Difference Derivative
!!
!! \details **deriv()** computes the first derivative of a function using a simple
!! fourth-order, central-difference scheme
!!
!! \author M. Miesch 
!! \date March 7, 2018: Created

subroutine deriv_1d(df,f,n,delta,periodic)
  
  real(kind=kind_real), intent(In) :: f(:)     !< 1D function
  real(kind=kind_real), intent(In) :: delta    !< grid spacing (assumed uniform)
  Integer, intent(In) :: n                     !< Number of grid points
  real(kind=kind_real), intent(InOut) :: df(:) !< result
  logical, optional, intent(In) :: periodic    !< Set to true if the domain is periodic
  
  integer :: i
  real(kind=kind_real) :: ho2, ho12
  logical :: wrap

  if (present(periodic)) then
     wrap=Periodic
  else
     wrap=.false.
  endIf
  
  ho2 = 1.0_kind_real/(2.0_kind_real*delta)
  ho12 = 1.0_kind_real/(12.0_kind_real*delta)

  ! interior points
  do i=3,n-2 
     df(i) = ho12*(f(i-2)-8.0_kind_real*f(i-1)+8.0_kind_real*f(i+1)-f(i+2))
  enddo

  if (wrap) then

     df(1) = ho12*(f(n-1)-8.0_kind_real*f(n)+8.0_kind_real*f(2)-f(3))
     df(2) = ho12*(f(n)-8.0_kind_real*f(1)+8.0_kind_real*f(3)-f(4))
     df(N-1) = ho12*(f(n-3)-8.0_kind_real*f(n-2)+8.0_kind_real*f(n)-f(1))
     df(N) = ho12*(f(n-2)-8.0_kind_real*f(n-1)+8.0_kind_real*f(1)-f(2))

  else

     ! left edge (lower order)
     df(1) = ho2*(-3.0_kind_real*f(1)+4.0_kind_real*f(2)-f(3))
     df(2) = ho12*(-3.0_kind_real*f(1)-10.0_kind_real*f(2)+18.0_kind_real*f(3)-6.0_kind_real*f(4)+f(5))                 

     ! right edge (lower order)
     df(n-1) = ho12*(3.0_kind_real*f(n)+10.0_kind_real*f(n-1)-18.0_kind_real*f(n-2)+6.0_kind_real*f(n-3)-f(n-4))                 
     df(n) = ho2*(3.0_kind_real*f(n)-4.0_kind_real*f(n-1)+f(n-2))

  endIf
     
end subroutine deriv_1d

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
integer :: ncid, nx_id, ny_id, nl_id, nc_id, ntot_id, four_id, gfld3d_id, xbound_id, qbound_id
character(len=max_string_length+50) :: record
character(len=max_string_length) :: filename
character(len=20) :: sdate, fmtn
character(len=4)  :: cnx
character(len=11) :: fmt1='(X,ES24.16)'
character(len=1024):: buf
integer :: jf, jy, jx, is
character(len=1024) :: subr='write_file'

call check(fld)

filename = genfilename(c_conf,max_string_length,vdate)
WRITE(buf,*) 'qg_field:write_file: writing '//filename
call fckit_log%info(buf)

is=0
if (fld%lbc) is=1
call datetime_to_string(vdate, sdate)

if (netcdfio) then
  
  call ncerr(subr,nf90_create(trim(filename)//'.nc',or(nf90_clobber,nf90_64bit_offset),ncid))
  call ncerr(subr,nf90_def_dim(ncid,'nx',fld%geom%nx,nx_id))
  call ncerr(subr,nf90_def_dim(ncid,'ny',fld%geom%ny,ny_id))
  call ncerr(subr,nf90_def_dim(ncid,'nl',fld%nl,nl_id))
  call ncerr(subr,nf90_def_dim(ncid,'nc',fld%nf,nc_id))
  call ncerr(subr,nf90_def_dim(ncid,'ntot',fld%nl*fld%nf,ntot_id))
  call ncerr(subr,nf90_def_dim(ncid,'four',4,four_id))
  call ncerr(subr,nf90_put_att(ncid,nf90_global,'lbc',is))
  
  call ncerr(subr,nf90_put_att(ncid,nf90_global,'sdate',sdate))
  call ncerr(subr,nf90_def_var(ncid,'gfld3d',nf90_double,(/nx_id,ny_id,ntot_id/),gfld3d_id))
  
  if (fld%lbc) then
    call ncerr(subr,nf90_def_var(ncid,'xbound',nf90_double,(/four_id/),xbound_id))
    call ncerr(subr,nf90_def_var(ncid,'qbound',nf90_double,(/nx_id,four_id/),qbound_id))
  end if
  
  call ncerr(subr,nf90_enddef(ncid))
  
  call ncerr(subr,nf90_put_var(ncid,gfld3d_id,fld%gfld3d))
  
  if (fld%lbc) then
    call ncerr(subr,nf90_put_var(ncid,xbound_id,fld%xbound))
    call ncerr(subr,nf90_put_var(ncid,qbound_id,fld%qbound))
  end if
  
  call ncerr(subr,nf90_close(ncid))

else

  open(unit=iunit, file=trim(filename), form='formatted', action='write')
  
  write(iunit,*) fld%geom%nx, fld%geom%ny, fld%nl, fld%nf, is
  write(iunit,*) sdate
  
  if (fld%geom%nx>9999)  call abor1_ftn("Format too small")
  write(cnx,'(I4)')fld%geom%nx
  fmtn='('//trim(cnx)//fmt1//')'
  
  do jf=1,fld%nl*fld%nf
    do jy=1,fld%geom%ny
      write(iunit,fmtn) (fld%gfld3d(jx,jy,jf), jx=1,fld%geom%nx)
    enddo
  enddo
  
  if (fld%lbc) then
    do jf=1,4
      write(iunit,fmt1) fld%xbound(jf)
    enddo
    do jf=1,4
      write(iunit,fmtn) (fld%qbound(jx,jf), jx=1,fld%geom%nx)
    enddo
  endif
  
  close(iunit) 
endif

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
               & /real(fld%nl*fld%geom%nx*fld%geom%ny,kind_real))
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
  pstat(3,jj)=sqrt(sum(fld%qbound(:,:)**2)/real(4*fld%geom%nx,kind_real))
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
  do jy=1,fld%geom%ny
    do jx=1,fld%geom%nx
      zz = zz + fld%gfld3d(jx,jy,jf)*fld%gfld3d(jx,jy,jf)
    enddo
  enddo
enddo

ii = fld%nl*fld%nf*fld%geom%ny*fld%geom%nx

if (fld%lbc) then
  do jf=1,4
    zz = zz + fld%xbound(jf)*fld%xbound(jf)
  enddo
  ii = ii + 4
  do jf=1,4
    do jx=1,fld%geom%nx
     zz = zz + fld%qbound(jx,jf)*fld%qbound(jx,jf)
    enddo
  enddo
  ii = ii + 4*fld%geom%nx
endif

prms = sqrt(zz/real(ii,kind_real))

end subroutine fldrms

! ------------------------------------------------------------------------------

subroutine interp_tl(fld, locs, vars, gom)
implicit none
type(qg_field), intent(in)   :: fld
type(qg_locs), intent(in)    :: locs
type(qg_vars), intent(in)    :: vars
type(qg_goms), intent(inout) :: gom
real(kind=kind_real) :: di, dj, ai, aj, valb, valt
integer :: jloc,joff,jvar,ii,jj,kk,iright,jsearch
character(len=1) :: cvar

call check(fld)

do jloc=1,locs%nloc
! Convert horizontal coordinate stored in ODB to grid locations
  di=locs%xyz(1,jloc)*real(fld%geom%nx,kind_real) 
  ii=int(di)
  ai=di-real(ii,kind_real);
  ii=ii+1
  if (ii<1.or.ii>fld%geom%nx) call abor1_ftn("qg_fields_interp_tl: error ii")
  iright=ii+1
  if (iright==fld%geom%nx+1) iright=1

  dj=locs%xyz(2,jloc)*real(fld%geom%ny,kind_real) 
  jj=int(dj)
  aj=dj-real(jj,kind_real);
  jj=jj+1
  if (jj<1.or.jj>fld%geom%ny) call abor1_ftn("qg_fields_interp_tl: error jj")

  kk=nint(locs%xyz(3,jloc))

! Loop on required variables
  gom%used=gom%used+1
  if (vars%nv/=gom%nvar) call abor1_ftn ("qg_fields_interp_tl: inconsistent var and gom")
  do jvar=1,vars%nv
    cvar=vars%fldnames(jvar)
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

    if (jj<fld%geom%ny) then
      valt=(1.0_kind_real-ai)*fld%gfld3d(ii    ,jj+1,joff+kk) &
                       & +ai *fld%gfld3d(iright,jj+1,joff+kk)
    else
      select case (cvar)
      case ("x")
        valt = 0.0_kind_real
      case ("u")
        valt=2.0_kind_real*( (1.0_kind_real-ai)*fld%gfld3d(ii,fld%geom%ny,joff+kk) &
             &                             +ai *fld%gfld3d(iright,fld%geom%ny,joff+kk) ) &
             & -           ( (1.0_kind_real-ai)*fld%gfld3d(ii,fld%geom%ny-1,joff+kk) &
             &                             +ai *fld%gfld3d(iright,fld%geom%ny-1,joff+kk) )
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

subroutine interp_ad(fld, locs, vars, gom)
implicit none
type(qg_field), intent(inout) :: fld
type(qg_locs), intent(in)     :: locs
type(qg_vars), intent(in)     :: vars
type(qg_goms), intent(inout)  :: gom
real(kind=kind_real) :: di, dj, ai, aj, valb, valt
integer :: jloc,joff,jvar,ii,jj,kk,iright,jsearch
character(len=1) :: cvar

call check(fld)

do jloc=locs%nloc,1,-1
! Convert horizontal coordinate stored in ODB to grid locations
  di=locs%xyz(1,jloc)*real(fld%geom%nx,kind_real)
  ii=int(di)
  ai=di-real(ii,kind_real);
  ii=ii+1
  if (ii<1.or.ii>fld%geom%nx) call abor1_ftn("qg_fields_interp_ad: error ii")
  iright=ii+1
  if (iright==fld%geom%nx+1) iright=1

  dj=locs%xyz(2,jloc)*real(fld%geom%ny,kind_real)
  jj=int(dj)
  aj=dj-real(jj,kind_real);
  jj=jj+1
  if (jj<1.or.jj>fld%geom%ny) call abor1_ftn("qg_fields_interp_ad: error jj")

  kk=nint(locs%xyz(3,jloc))

! Loop on required variables
  gom%used=gom%used+1
  if (vars%nv/=gom%nvar) call abor1_ftn ("qg_fields_interp_tl: inconsistent var and gom")
  do jvar=vars%nv,1,-1
    cvar=vars%fldnames(jvar)
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

    if (jj<fld%geom%ny) then
      fld%gfld3d(ii    ,jj+1,joff+kk) = &
             & fld%gfld3d(ii    ,jj+1,joff+kk)+(1.0_kind_real-ai)*valt
      fld%gfld3d(iright,jj+1,joff+kk) = &
             & fld%gfld3d(iright,jj+1,joff+kk)+ai*valt
    else
      select case (cvar)
      case ("u")
        fld%gfld3d(ii,fld%geom%ny,joff+kk) = &
             & fld%gfld3d(ii,fld%geom%ny,joff+kk)+2.0_kind_real*(1.0_kind_real-ai)*valt
        fld%gfld3d(iright,fld%geom%ny,joff+kk) = &
             & fld%gfld3d(iright,fld%geom%ny,joff+kk)+2.0_kind_real*ai*valt
        fld%gfld3d(ii,fld%geom%ny-1,joff+kk) = &
             & fld%gfld3d(ii,fld%geom%ny-1,joff+kk)-(1.0_kind_real-ai)*valt
        fld%gfld3d(iright,fld%geom%ny-1,joff+kk)= &
             & fld%gfld3d(iright,fld%geom%ny-1,joff+kk)-ai*valt
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

if (x1%geom%nx /= x2%geom%nx .or.  x1%geom%ny /= x2%geom%ny .or.  x1%nl /= x2%nl) then
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

bad = bad .or. (size(self%gfld3d, 1) /= self%geom%nx)
bad = bad .or. (size(self%gfld3d, 2) /= self%geom%ny)
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
  bad = bad .or. (size(self%qbound, 1) /= self%geom%nx)
  bad = bad .or. (size(self%qbound, 2) /= 4)
else
  bad = bad .or. associated(self%xbound)
  bad = bad .or. associated(self%qbound)
endif

if (bad) then
  write(0,*)'nx, ny, nf, nl, lbc = ',self%geom%nx,self%geom%ny,self%nf,self%nl,self%lbc
  if (associated(self%gfld3d)) write(0,*)'shape(gfld3d) = ',shape(self%gfld3d)
  if (associated(self%xbound)) write(0,*)'shape(xbound) = ',shape(self%xbound)
  if (associated(self%qbound)) write(0,*)'shape(qbound) = ',shape(self%qbound)
  call abor1_ftn ("qg_fields: field not consistent")
endif

end subroutine check

! ------------------------------------------------------------------------------

end module qg_fields
