! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

module qg_gom_mod

use config_mod
use fckit_log_module, only: fckit_log
use iso_c_binding
use kinds
use netcdf
use qg_constants_mod
use qg_geom_mod
use qg_locs_mod
use qg_projection_mod
use qg_tools_mod
use qg_vars_mod
use random_mod

implicit none
private
public :: qg_gom
public :: qg_gom_registry
public :: qg_gom_setup,qg_gom_create,qg_gom_delete,qg_gom_copy,qg_gom_zero,qg_gom_abs,qg_gom_random,qg_gom_mult, &
        & qg_gom_add,qg_gom_diff,qg_gom_schurmult,qg_gom_divide,qg_gom_rms,qg_gom_dotprod,qg_gom_stats,qg_gom_maxloc, &
        & qg_gom_read_file, qg_gom_write_file,qg_gom_analytic_init
! ------------------------------------------------------------------------------
type :: qg_gom
  integer :: nobs                             !< Number of observations
  integer :: used                             !< Index of used observation
  integer :: ix                               !< Streamfunction index
  integer :: iq                               !< Potential vorticity index
  integer :: iu                               !< Zonal wind index
  integer :: iv                               !< Meridian wind index
  integer :: nv                               !< Number of variables
  integer,allocatable :: indx(:)              !< Observations index
  real(kind_real), allocatable :: values(:,:) !< Observations values
  logical :: lalloc                           !< Allocation flag
end type qg_gom

#define LISTED_TYPE qg_gom

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry
type(registry_t) :: qg_gom_registry
! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
! Public
! ------------------------------------------------------------------------------
!> Linked list implementation
#include "oops/util/linkedList_c.f"
! ------------------------------------------------------------------------------
!> Setup GOM
subroutine qg_gom_setup(self,kobs,vars)

implicit none

! Passed variables
type(qg_gom),intent(inout) :: self !< GOM
integer,intent(in) :: kobs(:)      !< Observations index
type(qg_vars),intent(in) :: vars   !< Variables

! Set attributes
self%nobs = size(kobs)
self%used = 0
self%nv = 0
if (vars%lx) then
  self%nv = self%nv+1
  self%ix = self%nv
else
  self%ix = 0
endif
if (vars%lq) then
  self%nv = self%nv+1
  self%iq = self%nv
else
  self%iq = 0
endif
if (vars%lu) then
  self%nv = self%nv+1
  self%iu = self%nv
else
  self%iu = 0
endif
if (vars%lv) then
  self%nv = self%nv+1
  self%iv = self%nv
else
  self%iv = 0
endif

! Allocation
allocate(self%indx(self%nobs))
allocate(self%values(self%nv,self%nobs))
self%lalloc = .true.

! Initialization
self%indx = kobs

end subroutine qg_gom_setup
! ------------------------------------------------------------------------------
!> Create GOM
subroutine qg_gom_create(self)

implicit none

! Passed variables
type(qg_gom),intent(inout) :: self !< GOM

! Set allocation flag
self%lalloc = .false.

end subroutine qg_gom_create
! ------------------------------------------------------------------------------
!> Delete GOM
subroutine qg_gom_delete(self)

implicit none

! Passed variables
type(qg_gom),intent(inout) :: self !< GOM

! Release memory
if (self%lalloc) then
  deallocate(self%values)
  deallocate(self%indx)
endif

end subroutine qg_gom_delete
! ------------------------------------------------------------------------------
!> Copy GOM
subroutine qg_gom_copy(self,other)

implicit none

! Passed variables
type(qg_gom),intent(inout) :: self !< GOM
type(qg_gom),intent(in) :: other   !< Other GOM

! Copy attribues
self%nobs = other%nobs
self%ix = other%ix
self%iq = other%iq
self%iu = other%iu
self%iv = other%iv
self%nv = other%nv
self%used = other%used

! Allocation
if (.not.self%lalloc) then
   allocate(self%values(self%nv,self%nobs))
   allocate(self%indx(self%nobs))
   self%lalloc = .true.
endif

! Copy
self%values = other%values
self%indx = other%indx

end subroutine qg_gom_copy
! ------------------------------------------------------------------------------
!> Set GOM to zero
subroutine qg_gom_zero(self)

implicit none

! Passed variables
type(qg_gom),intent(inout) :: self !< GOM

! Set to zero
self%values = 0.0

end subroutine qg_gom_zero
! ------------------------------------------------------------------------------
!> Get GOM absolute value
subroutine qg_gom_abs(self)

implicit none

! Passed variables
type(qg_gom),intent(inout) :: self !< GOM

! Get absolute value
self%values = abs(self%values)

end subroutine qg_gom_abs
! ------------------------------------------------------------------------------
!> Generate random GOM values
subroutine qg_gom_random(self)

implicit none

! Passed variables
type(qg_gom),intent(inout) :: self !< GOM

! Generate random GOM values
call normal_distribution(self%values,0.0_kind_real,1.0_kind_real)

end subroutine qg_gom_random
! ------------------------------------------------------------------------------
!> Multiply GOM with a scalar
subroutine qg_gom_mult(self,zz)

implicit none

! Passed variables
type(qg_gom),intent(inout) :: self !< GOM
real(kind_real),intent(in) :: zz   !< Multiplier

! Local variables
integer :: jo,jv

! Multiply GOM with a scalar
do jo=1,self%nobs
  do jv=1,self%nv
    self%values(jv,jo) = zz*self%values(jv,jo)
  enddo
enddo

end subroutine qg_gom_mult
! ------------------------------------------------------------------------------
!> Add GOM
subroutine qg_gom_add(self,other)

implicit none

! Passed variables
type(qg_gom),intent(inout) :: self !< GOM
type(qg_gom),intent(in) :: other   !< Other GOM

! Local variables
integer :: jo,jv

! Add GOM
do jo=1,self%nobs
  do jv=1,self%nv
    self%values(jv,jo) = self%values(jv,jo)+other%values(jv,jo)
  enddo
enddo

end subroutine qg_gom_add
! ------------------------------------------------------------------------------
!> Subtract GOM
subroutine qg_gom_diff(self,other)

implicit none

! Passed variables
type(qg_gom),intent(inout) :: self !< GOM
type(qg_gom),intent(in) :: other   !< Other GOM

! Local variables
integer :: jo,jv

! Subtract GOM
do jo=1,self%nobs
  do jv=1,self%nv
    self%values(jv,jo) = self%values(jv,jo)-other%values(jv,jo)
  enddo
enddo

end subroutine qg_gom_diff
! ------------------------------------------------------------------------------
!> Schur product for GOM
subroutine qg_gom_schurmult(self,other)

implicit none

! Passed variables
type(qg_gom),intent(inout) :: self !< GOM
type(qg_gom),intent(in) :: other   !< Other GOM

! Local variables
integer :: jo,jv

! Add GOM
do jo=1,self%nobs
  do jv=1,self%nv
    self%values(jv,jo) = self%values(jv,jo)*other%values(jv,jo)
  enddo
enddo

end subroutine qg_gom_schurmult
! ------------------------------------------------------------------------------
!> Schur division for GOM
subroutine qg_gom_divide(self,other)

implicit none

! Passed variables
type(qg_gom),intent(inout) :: self !< GOM
type(qg_gom),intent(in) :: other   !< Other GOM

! Local variables
real(kind_real) :: tol
integer :: jloc,jvar

! Set tolerance
tol = epsilon(tol)

! Conditional division
do jvar=1,self%nv
  do jloc=1,self%nobs
    if (abs(other%values(jvar,jloc))>tol) then
      self%values(jvar,jloc) = self%values(jvar,jloc)/other%values(jvar,jloc)
    else
      self%values(jvar,jloc) = 0.0
    endif
  enddo
enddo

end subroutine qg_gom_divide
! ------------------------------------------------------------------------------
!> Compute GOM RMS
subroutine qg_gom_rms(self,rms)

implicit none

! Passed variables
type(qg_gom),intent(inout) :: self   !< GOM
real(kind_real),intent(inout) :: rms !< RMS

! Local variables
integer :: jo,jv

! Initialization
rms = 0.0

! Loop over values
do jo=1,self%nobs
  do jv=1,self%nv
    rms = rms+self%values(jv,jo)**2
  enddo
enddo

! Normalize and take square-root
rms = sqrt(rms/(self%nobs*self%nv))

end subroutine qg_gom_rms
! ------------------------------------------------------------------------------
!> GOM dot product
subroutine qg_gom_dotprod(gom1,gom2,prod)

implicit none

! Passed variables
type(qg_gom),intent(in) :: gom1       !< GOM 1
type(qg_gom),intent(in) :: gom2       !< GOM 2
real(kind_real),intent(inout) :: prod !< Dot product

! Local variables
integer :: jo,jv

! Check
if ((gom1%nv/=gom2%nv).or.(gom1%nobs/=gom2%nobs)) call abor1_ftn('qg_gom_dotprod: inconsistent GOM sizes')

! Initialization
prod = 0.0

! Loop over values
do jo=1,gom1%nobs
  do jv=1,gom1%nv
    prod = prod+gom1%values(jv,jo)*gom2%values(jv,jo)
  enddo
enddo

end subroutine qg_gom_dotprod
! ------------------------------------------------------------------------------
!> Compute GOM stats
subroutine qg_gom_stats(self,kobs,scaling,pmin,pmax,prms)

implicit none

! Passed variables
type(qg_gom),intent(inout) :: self       !< GOM
integer,intent(inout) :: kobs            !< Number of observations
real(kind_real),intent(inout) :: scaling !< Scaling value
real(kind_real),intent(inout) :: pmin    !< Minimum value
real(kind_real),intent(inout) :: pmax    !< Maximum value
real(kind_real),intent(inout) :: prms    !< RMS

! Local variables
real(kind_real) :: expo

! Compute GOM stats
kobs = self%nobs
pmin = minval(self%values)
pmax = maxval(self%values)
prms = sqrt(sum(self%values**2)/real(self%nobs*self%nv,kind_real))

! Scaling
if (abs(prms)>0.0) then
  expo = aint(log(abs(prms))/log(10.0_kind_real))
  scaling = 10.0**expo
else
  scaling = 1.0
endif
pmin = pmin/scaling
pmax = pmax/scaling
prms = prms/scaling


end subroutine qg_gom_stats
! ------------------------------------------------------------------------------
!> Find and locate GOM max. value
subroutine qg_gom_maxloc(self,mxval,iloc,ivar)

implicit none

! Passed variables
type(qg_gom),intent(inout) :: self     !< GOM
real(kind_real),intent(inout) :: mxval !< Maximum value
integer,intent(inout) :: iloc          !< Location of maximum value
integer,intent(inout) :: ivar          !< Variable with maximum value

! Local variables
integer :: mxloc(2)

! Find GOM max. value
mxval = maxval(self%values)
mxloc = maxloc(self%values)

! Locate GOM max. value
ivar = mxloc(1)
iloc = mxloc(2)

end subroutine qg_gom_maxloc
! ------------------------------------------------------------------------------
!> Read GOM from file
subroutine qg_gom_read_file(self,conf)

implicit none

! Passed variables
type(qg_gom),intent(inout) :: self !< GOM
type(c_ptr),intent(in) :: conf     !< Configuration

! Local variables
integer :: ncid,nobs_id,nv_id,indx_id,values_id
character(len=1024) :: filename

! Check allocation
if (self%lalloc) call abor1_ftn('qg_gom_read_file: gom alredy allocated')

! Get filename
filename = trim(config_get_string(conf,len(filename),'filename'))
call fckit_log%info('qg_gom_read_file: reading '//trim(filename))

! Open NetCDF file
call ncerr(nf90_open(trim(filename)//'.nc',nf90_nowrite,ncid))

! Get dimensions ids
call ncerr(nf90_inq_dimid(ncid,'nobs',nobs_id))
call ncerr(nf90_inq_dimid(ncid,'nv',nv_id))

! Get dimensions
call ncerr(nf90_inquire_dimension(ncid,nobs_id,len=self%nobs))
call ncerr(nf90_inquire_dimension(ncid,nv_id,len=self%nv))

! Get attributes
call ncerr(nf90_get_att(ncid,nf90_global,'used',self%used))
call ncerr(nf90_get_att(ncid,nf90_global,'ix',self%ix))
call ncerr(nf90_get_att(ncid,nf90_global,'iq',self%iq))
call ncerr(nf90_get_att(ncid,nf90_global,'iu',self%iu))
call ncerr(nf90_get_att(ncid,nf90_global,'iv',self%iv))
call ncerr(nf90_get_att(ncid,nf90_global,'nv',self%nv))

! Allocation
allocate(self%indx(self%nobs))
allocate(self%values(self%nv,self%nobs))
self%lalloc = .true.

! Initialization
call qg_gom_zero(self)

! Get variables ids
call ncerr(nf90_inq_varid(ncid,'indx',indx_id))
call ncerr(nf90_inq_varid(ncid,'values',values_id))

! Get variables
call ncerr(nf90_get_var(ncid,indx_id,self%indx))
call ncerr(nf90_get_var(ncid,values_id,self%values))

! Close NetCDF file
call ncerr(nf90_close(ncid))

end subroutine qg_gom_read_file
! ------------------------------------------------------------------------------
!> Write GOM to file
subroutine qg_gom_write_file(self,conf)

implicit none

! Passed variables
type(qg_gom),intent(inout) :: self !< GOM
type(c_ptr),intent(in) :: conf     !< Configuration

! Local variables
integer :: ncid,nobs_id,nv_id,indx_id,values_id
character(len=1024) :: filename

! Check allocation
if (.not.self%lalloc) call abor1_ftn('qg_gom_write_file: gom not allocated')

! Set filename
filename = trim(config_get_string(conf,len(filename),'filename'))
call fckit_log%info('qg_gom_write_file: writing '//trim(filename))

! Create NetCDF file
call ncerr(nf90_create(trim(filename)//'.nc',or(nf90_clobber,nf90_64bit_offset),ncid))

! Define dimensions
call ncerr(nf90_def_dim(ncid,'nobs',self%nobs,nobs_id))
call ncerr(nf90_def_dim(ncid,'nv',self%nv,nv_id))

! Define attributes
call ncerr(nf90_put_att(ncid,nf90_global,'used',self%used))
call ncerr(nf90_put_att(ncid,nf90_global,'ix',self%ix))
call ncerr(nf90_put_att(ncid,nf90_global,'iq',self%iq))
call ncerr(nf90_put_att(ncid,nf90_global,'iu',self%iu))
call ncerr(nf90_put_att(ncid,nf90_global,'iv',self%iv))
call ncerr(nf90_put_att(ncid,nf90_global,'nv',self%nv))

! Define variables
call ncerr(nf90_def_var(ncid,'indx',nf90_int,(/nobs_id/),indx_id))
call ncerr(nf90_def_var(ncid,'values',nf90_double,(/nv_id,nobs_id/),values_id))

! End definitions
call ncerr(nf90_enddef(ncid))

! Put variables
call ncerr(nf90_put_var(ncid,indx_id,self%indx))
call ncerr(nf90_put_var(ncid,values_id,self%values))

! Close NetCDF file
call ncerr(nf90_close(ncid))

end subroutine qg_gom_write_file
! ------------------------------------------------------------------------------
!> GOM analytic initialization
subroutine qg_gom_analytic_init(self,locs,conf)

implicit none

! Passed variables
type(qg_gom),intent(inout) :: self  !< GOM
type(qg_locs),intent(inout) :: locs !< Locations
type(c_ptr),intent(in)    :: conf   !< Configuration

! Local variables
integer :: iloc
real(kind_real) :: x,y
character(len=30) :: ic

! Check allocation
if (.not. self%lalloc) call abor1_ftn('qg_gom_analytic init: gom not allocated')

! Get analytic configuration
ic = config_get_string(conf,len(ic),'analytic_init')
call fckit_log%info('qg_gom_analytic_init: ic = '//trim(ic))
do iloc=1,locs%nlocs
  select case (trim(ic))
  case ('baroclinic-instability')
    ! Go to cartesian coordinates
    call lonlat_to_xy(locs%lon(iloc),locs%lat(iloc),x,y)

    ! Compute values for baroclinic instability
    if (self%ix>0) call baroclinic_instability(x,y,locs%z(iloc),'x',self%values(self%ix,iloc))
    if (self%iq>0) call baroclinic_instability(x,y,locs%z(iloc),'q',self%values(self%iq,iloc))
    if (self%iu>0) call baroclinic_instability(x,y,locs%z(iloc),'u',self%values(self%iu,iloc))
    if (self%iv>0) call baroclinic_instability(x,y,locs%z(iloc),'v',self%values(self%iv,iloc))
  case ('large-vortices')
    ! Go to cartesian coordinates
    call lonlat_to_xy(locs%lon(iloc),locs%lat(iloc),x,y)

    ! Compute values for large vortices
    if (self%ix>0) call large_vortices(x,y,locs%z(iloc),'x',self%values(self%ix,iloc))
    if (self%iq>0) call large_vortices(x,y,locs%z(iloc),'q',self%values(self%iq,iloc))
    if (self%iu>0) call large_vortices(x,y,locs%z(iloc),'u',self%values(self%iu,iloc))
    if (self%iv>0) call large_vortices(x,y,locs%z(iloc),'v',self%values(self%iv,iloc))
  case default
    call abor1_ftn('qg_gom_analytic_init: unknown initialization')
  endselect
enddo

end subroutine qg_gom_analytic_init
! ------------------------------------------------------------------------------
end module qg_gom_mod
