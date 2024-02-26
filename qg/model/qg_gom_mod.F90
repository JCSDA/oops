! (C) Copyright 2009-2016 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

module qg_gom_mod

use atlas_module, only: atlas_field
use fckit_configuration_module, only: fckit_configuration
use fckit_log_module, only: fckit_log
use iso_c_binding
use datetime_mod
use kinds
use netcdf
use oops_variables_mod
use qg_constants_mod
use qg_geom_mod
use qg_locs_mod
use qg_projection_mod
use qg_tools_mod
use random_mod
use string_f_c_mod

implicit none
private
public :: qg_gom
public :: qg_gom_registry
public :: qg_gom_setup,qg_gom_alloc,qg_gom_delete,qg_gom_copy,qg_gom_zero,qg_gom_abs,qg_gom_random,qg_gom_mult, &
        & qg_gom_add,qg_gom_diff,qg_gom_schurmult,qg_gom_divide,qg_gom_rms,qg_gom_dotprod,qg_gom_stats,qg_gom_maxloc, &
        & qg_gom_fill, qg_gom_fillad, qg_gom_read_file, qg_gom_write_file,qg_gom_analytic_init
! ------------------------------------------------------------------------------
type :: qg_gom
  integer :: nobs                               !< Number of observations
  real(kind_real), pointer :: x(:,:) => null()  !< Streamfunction observations values
  real(kind_real), pointer :: q(:,:) => null()  !< Potential vorticity observations values
  real(kind_real), pointer :: u(:,:) => null()  !< Zonal wind observations values
  real(kind_real), pointer :: v(:,:) => null()  !< Meridian wind observations values
  real(kind_real), pointer :: z(:,:) => null()  !< Height values
  logical :: lalloc = .false.                   !< Allocation flag
  type(oops_variables) :: vars                  !< Variables
  integer :: levs                               !< Number of levels
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
subroutine qg_gom_setup(self,npaths,vars,levs)

implicit none

! Passed variables
type(qg_gom),intent(inout) :: self      !< GOM
integer(c_int), intent(in) :: npaths    !< Number of interpolation paths
type(oops_variables), intent(in) :: vars
integer(c_int),intent(in)  :: levs

! Set attributes
self%vars = vars

call qg_gom_alloc(self, npaths, levs)

end subroutine qg_gom_setup
! ------------------------------------------------------------------------------
!> Setup GOM
subroutine qg_gom_alloc(self, nobs, levs)

implicit none
type(qg_gom),intent(inout) :: self      !< GOM
integer, intent(in) :: nobs
integer, intent(in) :: levs

self%levs=levs
self%nobs=nobs
if (self%vars%has('x')) allocate(self%x(self%levs, self%nobs))
if (self%vars%has('q')) allocate(self%q(self%levs, self%nobs))
if (self%vars%has('u')) allocate(self%u(self%levs, self%nobs))
if (self%vars%has('v')) allocate(self%v(self%levs, self%nobs))
if (self%vars%has('z')) allocate(self%z(self%levs, self%nobs))
self%lalloc = .true.

end subroutine qg_gom_alloc
! ------------------------------------------------------------------------------
!> Delete GOM
subroutine qg_gom_delete(self)

implicit none

! Passed variables
type(qg_gom),intent(inout) :: self !< GOM

! Release memory
if (associated(self%x)) deallocate(self%x)
if (associated(self%q)) deallocate(self%q)
if (associated(self%u)) deallocate(self%u)
if (associated(self%v)) deallocate(self%v)
if (associated(self%z)) deallocate(self%z)
self%lalloc = .false.

end subroutine qg_gom_delete
! ------------------------------------------------------------------------------
!> Copy GOM
subroutine qg_gom_copy(self,other)

implicit none

! Passed variables
type(qg_gom),intent(inout) :: self            !< GOM
type(qg_gom),intent(in) :: other              !< Other GOM

! Copy attributes
self%nobs = other%nobs
self%levs = other%levs
self%vars = other%vars

! Allocation
if (.not.self%lalloc) then
  if (self%vars%has('x')) allocate(self%x(self%levs, self%nobs))
  if (self%vars%has('q')) allocate(self%q(self%levs, self%nobs))
  if (self%vars%has('u')) allocate(self%u(self%levs, self%nobs))
  if (self%vars%has('v')) allocate(self%v(self%levs, self%nobs))
  if (self%vars%has('z')) allocate(self%z(self%levs, self%nobs))
  self%lalloc = .true.
endif

! Copy
if (self%vars%has('x')) self%x = other%x
if (self%vars%has('q')) self%q = other%q
if (self%vars%has('u')) self%u = other%u
if (self%vars%has('v')) self%v = other%v
if (self%vars%has('z')) self%z = other%z

end subroutine qg_gom_copy
! ------------------------------------------------------------------------------
subroutine qg_gom_fill(self, lvar, c_var, c_nloc, c_indx, c_nlev, c_vals)
implicit none
type(qg_gom), intent(inout) :: self
integer(c_int), intent(in) :: lvar
character(kind=c_char, len=1), intent(in) :: c_var(lvar+1)
integer(c_int), intent(in) :: c_nloc
integer(c_int), intent(in) :: c_indx(c_nloc)
integer(c_int), intent(in) :: c_nlev
real(c_double), intent(in) :: c_vals(c_nloc, c_nlev)

character(len=1024) :: fieldname
real(kind_real),pointer :: gval(:,:)
integer :: jlev, jloc, iloc, ii

if (.not.self%lalloc) call abor1_ftn('qg_gom_fill: gom not allocated')
if (self%levs /= c_nlev) call abor1_ftn('qg_gom_fillad: incorrect number of levels')

call c_f_string(c_var, fieldname)

select case (trim(fieldname))
case ('x')
  gval => self%x(:,:)
case ('q')
  gval => self%q(:,:)
case ('u')
  gval => self%u(:,:)
case ('v')
  gval => self%v(:,:)
case ('z')
  gval => self%z(:,:)
case default
  call abor1_ftn('qg_gom_fill: wrong variable')
endselect

do jlev = 1, self%levs
  do jloc=1,c_nloc
    iloc = c_indx(jloc)
    gval(jlev,iloc) = c_vals(jloc, jlev)
  enddo
enddo

end subroutine qg_gom_fill
! ------------------------------------------------------------------------------
subroutine qg_gom_fillad(self, lvar, c_var, c_nloc, c_indx, c_nlev, c_vals)
implicit none
type(qg_gom), intent(in) :: self
integer(c_int), intent(in) :: lvar
character(kind=c_char, len=1), intent(in) :: c_var(lvar+1)
integer(c_int), intent(in) :: c_nloc
integer(c_int), intent(in) :: c_indx(c_nloc)
integer(c_int), intent(in) :: c_nlev
real(c_double), intent(inout) :: c_vals(c_nloc, c_nlev)

character(len=1024) :: fieldname
real(kind_real),pointer :: gval(:,:)
integer :: jlev, jloc, iloc

if (.not.self%lalloc) call abor1_ftn('qg_gom_fillad: gom not allocated')
if (self%levs /= c_nlev) call abor1_ftn('qg_gom_fillad: incorrect number of levels')

call c_f_string(c_var, fieldname)
select case (trim(fieldname))
case ('x')
  gval => self%x(:,:)
case ('q')
  gval => self%q(:,:)
case ('u')
  gval => self%u(:,:)
case ('v')
  gval => self%v(:,:)
case ('z')
  gval => self%z(:,:)
case default
  call abor1_ftn('qg_gom_fillad: wrong variable')
endselect

do jlev = 1, self%levs
  do jloc=1,c_nloc
    iloc = c_indx(jloc)
    c_vals(jloc,jlev) = gval(jlev,iloc)
  enddo
enddo

end subroutine qg_gom_fillad
! ------------------------------------------------------------------------------
!> Set GOM to zero
subroutine qg_gom_zero(self)

implicit none

! Passed variables
type(qg_gom),intent(inout) :: self !< GOM

! Set to zero
if (self%vars%has('x')) self%x = 0.0
if (self%vars%has('q')) self%q = 0.0
if (self%vars%has('u')) self%u = 0.0
if (self%vars%has('v')) self%v = 0.0

end subroutine qg_gom_zero
! ------------------------------------------------------------------------------
!> Get GOM absolute value
subroutine qg_gom_abs(self)

implicit none

! Passed variables
type(qg_gom),intent(inout) :: self !< GOM

! Get absolute value
if (self%vars%has('x')) self%x = abs(self%x)
if (self%vars%has('q')) self%q = abs(self%q)
if (self%vars%has('u')) self%u = abs(self%u)
if (self%vars%has('v')) self%v = abs(self%v)

end subroutine qg_gom_abs
! ------------------------------------------------------------------------------
!> Generate random GOM values
subroutine qg_gom_random(self)

implicit none

! Passed variables
type(qg_gom),intent(inout) :: self !< GOM

! Local variables
integer :: nv
real(kind_real),allocatable :: values(:,:,:)

! TODO(Benjamin): change that in a following PR
nv = 0
if (self%vars%has('x')) nv = nv+1
if (self%vars%has('q')) nv = nv+1
if (self%vars%has('u')) nv = nv+1
if (self%vars%has('v')) nv = nv+1
allocate(values(nv,self%levs,self%nobs))

! Generate random GOM values
call normal_distribution(values,0.0_kind_real,1.0_kind_real)

! Split random values
nv = 0
if (self%vars%has('x')) then
   nv = nv+1
   self%x = values(nv,:,:)
endif
if (self%vars%has('q')) then
   nv = nv+1
   self%q = values(nv,:,:)
endif
if (self%vars%has('u')) then
   nv = nv+1
   self%u = values(nv,:,:)
endif
if (self%vars%has('v')) then
   nv = nv+1
   self%v = values(nv,:,:)
endif

end subroutine qg_gom_random
! ------------------------------------------------------------------------------
!> Multiply GOM with a scalar
subroutine qg_gom_mult(self,zz)

implicit none

! Passed variables
type(qg_gom),intent(inout) :: self !< GOM
real(kind_real),intent(in) :: zz   !< Multiplier

! Multiply GOM with a scalar
if (self%vars%has('x')) self%x = zz*self%x
if (self%vars%has('q')) self%q = zz*self%q
if (self%vars%has('u')) self%u = zz*self%u
if (self%vars%has('v')) self%v = zz*self%v

end subroutine qg_gom_mult
! ------------------------------------------------------------------------------
!> Add GOM
subroutine qg_gom_add(self,other)

implicit none

! Passed variables
type(qg_gom),intent(inout) :: self !< GOM
type(qg_gom),intent(in) :: other   !< Other GOM

! Add GOM
if (self%vars%has('x')) self%x = self%x+other%x
if (self%vars%has('q')) self%q = self%q+other%q
if (self%vars%has('u')) self%u = self%u+other%u
if (self%vars%has('v')) self%v = self%v+other%v

end subroutine qg_gom_add
! ------------------------------------------------------------------------------
!> Subtract GOM
subroutine qg_gom_diff(self,other)

implicit none

! Passed variables
type(qg_gom),intent(inout) :: self !< GOM
type(qg_gom),intent(in) :: other   !< Other GOM

! Subtract GOM
if (self%vars%has('x')) self%x = self%x-other%x
if (self%vars%has('q')) self%q = self%q-other%q
if (self%vars%has('u')) self%u = self%u-other%u
if (self%vars%has('v')) self%v = self%v-other%v

end subroutine qg_gom_diff
! ------------------------------------------------------------------------------
!> Schur product for GOM
subroutine qg_gom_schurmult(self,other)

implicit none

! Passed variables
type(qg_gom),intent(inout) :: self !< GOM
type(qg_gom),intent(in) :: other   !< Other GOM

! Multiply GOM
if (self%vars%has('x')) self%x = self%x*other%x
if (self%vars%has('q')) self%q = self%q*other%q
if (self%vars%has('u')) self%u = self%u*other%u
if (self%vars%has('v')) self%v = self%v*other%v

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
integer :: jloc, jlev

! Set tolerance
tol = epsilon(tol)

! Conditional division
do jloc=1,self%nobs
 do jlev=1,self%levs
  if (self%vars%has('x')) then
    if (abs(other%x(jlev,jloc))>tol) then
      self%x(jlev,jloc) = self%x(jlev,jloc)/other%x(jlev,jloc)
    else
      self%x(jlev,jloc) = 0.0
    endif
  endif
  if (self%vars%has('q')) then
    if (abs(other%q(jlev,jloc))>tol) then
      self%q(jlev,jloc) = self%q(jlev,jloc)/other%q(jlev,jloc)
    else
      self%q(jlev,jloc) = 0.0
    endif
  endif
  if (self%vars%has('u')) then
    if (abs(other%u(jlev,jloc))>tol) then
      self%u(jlev,jloc) = self%u(jlev,jloc)/other%u(jlev,jloc)
    else
      self%u(jlev,jloc) = 0.0
    endif
  endif
  if (self%vars%has('v')) then
    if (abs(other%v(jlev,jloc))>tol) then
      self%v(jlev,jloc) = self%v(jlev,jloc)/other%v(jlev,jloc)
    else
      self%v(jlev,jloc) = 0.0
    endif
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
integer :: nv

! Initialization
rms = 0.0
nv = 0

! Loop over values
if (self%vars%has('x')) then
  rms = rms+sum(self%x**2)
  nv = nv+1
endif
if (self%vars%has('q')) then
  rms = rms+sum(self%q**2)
  nv = nv+1
endif
if (self%vars%has('u')) then
  rms = rms+sum(self%u**2)
  nv = nv+1
endif
if (self%vars%has('v')) then
  rms = rms+sum(self%v**2)
  nv = nv+1
endif

! Normalize and take square-root
rms = sqrt(rms/real(self%nobs*nv,kind_real))

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
if (gom1%nobs/=gom2%nobs) call abor1_ftn('qg_gom_dotprod: inconsistent GOM sizes')

! Initialization
prod = 0.0

! Dot product
if (gom1%vars%has('x').and.gom2%vars%has('x')) prod = prod+sum(gom1%x*gom2%x)
if (gom1%vars%has('q').and.gom2%vars%has('q')) prod = prod+sum(gom1%q*gom2%q)
if (gom1%vars%has('u').and.gom2%vars%has('u')) prod = prod+sum(gom1%u*gom2%u)
if (gom1%vars%has('v').and.gom2%vars%has('v')) prod = prod+sum(gom1%v*gom2%v)

end subroutine qg_gom_dotprod
! ------------------------------------------------------------------------------
!> Compute GOM stats
subroutine qg_gom_stats(self,kobs,pmin,pmax,prms)

implicit none

! Passed variables
type(qg_gom),intent(inout) :: self       !< GOM
integer,intent(inout) :: kobs            !< Number of observations
real(kind_real),intent(inout) :: pmin    !< Minimum value
real(kind_real),intent(inout) :: pmax    !< Maximum value
real(kind_real),intent(inout) :: prms    !< RMS

! Local variables
integer :: nv

! Compute GOM stats
kobs = self%nobs
if (self%nobs>0) then
  pmin = huge(1.0)
  pmax = -huge(1.0)
  prms = 0.0
  nv = 0
  if (self%vars%has('x')) then
    pmin = min(pmin,minval(self%x))
    pmax = max(pmax,maxval(self%x))
    prms = prms+sum(self%x**2)
    nv = nv+1
  endif
  if (self%vars%has('q')) then
    pmin = min(pmin,minval(self%q))
    pmax = max(pmax,maxval(self%q))
    prms = prms+sum(self%q**2)
    nv = nv+1
  endif
  if (self%vars%has('u')) then
    pmin = min(pmin,minval(self%u))
    pmax = max(pmax,maxval(self%u))
    prms = prms+sum(self%u**2)
    nv = nv+1
  endif
  if (self%vars%has('v')) then
    pmin = min(pmin,minval(self%v))
    pmax = max(pmax,maxval(self%v))
    prms = prms+sum(self%v**2)
    nv = nv+1
  endif
  prms = sqrt(prms/real(self%nobs*nv,kind_real))
else
  pmin = 0.0
  pmax = 0.0
  prms = 0.0
end if

end subroutine qg_gom_stats
! ------------------------------------------------------------------------------
!> Find and locate GOM max. value
subroutine qg_gom_maxloc(self,mxval,mxloc,mxvar)

implicit none

! Passed variables
type(qg_gom),intent(inout) :: self          !< GOM
real(kind_real),intent(inout) :: mxval      !< Maximum value
integer,intent(inout) :: mxloc              !< Location of maximum value
type(oops_variables),intent(inout) :: mxvar !< Variable of maximum value

! Local variables
integer :: mxloc_arr(2),mxval_tmp
character(len=1) :: var

! Initialization
mxval = -huge(1.0)

! Find GOM max. value
if (self%vars%has('x')) then
  mxval_tmp = maxval(self%x)
  if (mxval_tmp>mxval) then
    mxval = mxval
    mxloc_arr = maxloc(self%x)
    var = 'x'
  endif
endif
if (self%vars%has('q')) then
  mxval_tmp = maxval(self%q)
  if (mxval_tmp>mxval) then
    mxval = mxval
    mxloc_arr = maxloc(self%q)
    var = 'q'
  endif
endif
if (self%vars%has('u')) then
  mxval_tmp = maxval(self%u)
  if (mxval_tmp>mxval) then
    mxval = mxval
    mxloc_arr = maxloc(self%u)
    var = 'u'
  endif
endif
if (self%vars%has('v')) then
  mxval_tmp = maxval(self%v)
  if (mxval_tmp>mxval) then
    mxval = mxval
    mxloc_arr = maxloc(self%v)
    var = 'v'
  endif
endif

! Locate GOM max. value
mxloc = mxloc_arr(2)

! Set GOM max. variable
call mxvar%push_back(var)

end subroutine qg_gom_maxloc
! ------------------------------------------------------------------------------
!> Read GOM from file
subroutine qg_gom_read_file(self,vars,f_conf)

implicit none

! Passed variables
type(qg_gom),intent(inout) :: self             !< GOM
type(oops_variables), intent(in) :: vars
type(fckit_configuration),intent(in) :: f_conf !< FCKIT configuration

! Local variables
integer :: ncid,nobs_id,nlev_id,nobs,levs,x_id,q_id,u_id,v_id,z_id
character(len=1024) :: filename
character(len=:),allocatable :: str

! Get filename
call f_conf%get_or_die("filename",str)
filename = str
call fckit_log%info('qg_gom_read_file: reading '//trim(filename))

! Open NetCDF file
call ncerr(nf90_open(trim(filename)//'.nc',nf90_nowrite,ncid))

! Get dimension id
call ncerr(nf90_inq_dimid(ncid,'nobs',nobs_id))
call ncerr(nf90_inq_dimid(ncid,'levs',nlev_id))

! Get dimension
call ncerr(nf90_inquire_dimension(ncid,nobs_id,len=nobs))
call ncerr(nf90_inquire_dimension(ncid,nlev_id,len=levs))

! GOM setup
self%vars = vars
call qg_gom_alloc(self, nobs, levs)

! Get variables ids
if (self%vars%has('x')) call ncerr(nf90_inq_varid(ncid,'x',x_id))
if (self%vars%has('q')) call ncerr(nf90_inq_varid(ncid,'q',q_id))
if (self%vars%has('u')) call ncerr(nf90_inq_varid(ncid,'u',u_id))
if (self%vars%has('v')) call ncerr(nf90_inq_varid(ncid,'v',v_id))
if (self%vars%has('z')) call ncerr(nf90_inq_varid(ncid,'z',z_id))

! Get variables
if (self%vars%has('x')) call ncerr(nf90_get_var(ncid,x_id,self%x))
if (self%vars%has('q')) call ncerr(nf90_get_var(ncid,q_id,self%q))
if (self%vars%has('u')) call ncerr(nf90_get_var(ncid,u_id,self%u))
if (self%vars%has('v')) call ncerr(nf90_get_var(ncid,v_id,self%v))
if (self%vars%has('z')) call ncerr(nf90_get_var(ncid,z_id,self%z))

! Close NetCDF file
call ncerr(nf90_close(ncid))

end subroutine qg_gom_read_file
! ------------------------------------------------------------------------------
!> Write GOM to file
subroutine qg_gom_write_file(self,f_conf)

implicit none

! Passed variables
type(qg_gom),intent(inout) :: self !< GOM
type(fckit_configuration),intent(in) :: f_conf !< FCKIT configuration

! Local variables
integer :: ncid,nobs_id,nlev_id,x_id,q_id,u_id,v_id,z_id
character(len=1024) :: filename
character(len=:),allocatable :: str

! Check allocation
if (.not.self%lalloc) call abor1_ftn('qg_gom_write_file: gom not allocated')

! Set filename
call f_conf%get_or_die("filename",str)
filename = str
call fckit_log%info('qg_gom_write_file: writing '//trim(filename))

! Create NetCDF file
call ncerr(nf90_create(trim(filename)//'.nc',ior(nf90_clobber,nf90_64bit_offset),ncid))

! Define dimensions
call ncerr(nf90_def_dim(ncid,'nobs',self%nobs,nobs_id))
call ncerr(nf90_def_dim(ncid,'levs',self%levs,nlev_id))

! Define variables
if (self%vars%has('x')) call ncerr(nf90_def_var(ncid,'x',nf90_double,(/nlev_id,nobs_id/),x_id))
if (self%vars%has('q')) call ncerr(nf90_def_var(ncid,'q',nf90_double,(/nlev_id,nobs_id/),q_id))
if (self%vars%has('u')) call ncerr(nf90_def_var(ncid,'u',nf90_double,(/nlev_id,nobs_id/),u_id))
if (self%vars%has('v')) call ncerr(nf90_def_var(ncid,'v',nf90_double,(/nlev_id,nobs_id/),v_id))
if (self%vars%has('z')) call ncerr(nf90_def_var(ncid,'z',nf90_double,(/nlev_id,nobs_id/),z_id))

! End definitions
call ncerr(nf90_enddef(ncid))

! Put variables
if (self%vars%has('x')) call ncerr(nf90_put_var(ncid,x_id,self%x))
if (self%vars%has('q')) call ncerr(nf90_put_var(ncid,q_id,self%q))
if (self%vars%has('u')) call ncerr(nf90_put_var(ncid,u_id,self%u))
if (self%vars%has('v')) call ncerr(nf90_put_var(ncid,v_id,self%v))
if (self%vars%has('z')) call ncerr(nf90_put_var(ncid,z_id,self%z))

! Close NetCDF file
call ncerr(nf90_close(ncid))

end subroutine qg_gom_write_file
! ------------------------------------------------------------------------------
!> GOM analytic initialization
subroutine qg_gom_analytic_init(self,locs,f_conf)

implicit none

! Passed variables
type(qg_gom),intent(inout) :: self             !< GOM
type(qg_locs),intent(inout) :: locs            !< Locations
type(fckit_configuration),intent(in) :: f_conf !< FCKIT configuration

! Local variables
integer :: iloc,ilev,iz,nz
real(kind_real) :: x,y
character(len=30) :: ic
character(len=:),allocatable :: str
real(kind_real), pointer :: lonlat(:,:)
real(kind_real), allocatable :: depths(:)
type(atlas_field) :: lonlat_field
real(kind_real) :: z(self%levs)

! Get depths
!nz = f_conf%get_size("depths")
!if (nz /= self%levs) call abor1_ftn('qg_gom_analytic init: wrong number of levels')
!allocate(depths(nz))
!call f_conf%get_or_die("depths",depths)

! Should get depths from yaml as above, however, since the GetValues is going to be
! re-writen it is not worth spending the time getting it through the Parameters now.
! Values are from the geometry section of getvalues.yaml
nz = 2
allocate(depths(nz))
depths(1) = 4500.0
depths(2) = 5500.0

! Set heights
z(1) = 0.5*depths(1)
do iz=2,nz
  z(iz) = sum(depths(1:iz-1))+0.5*depths(iz)
end do

! get locations
lonlat_field = locs%lonlat()
call lonlat_field%data(lonlat)

! Check allocation
if (.not.self%lalloc) call abor1_ftn('qg_gom_analytic init: gom not allocated')
if (.not.associated(self%z)) allocate(self%z(self%levs,locs%nlocs()))

! Get analytic configuration
call f_conf%get_or_die("method",str)
ic = str
call fckit_log%info('qg_gom_analytic_init: ic = '//trim(ic))
do iloc=1,locs%nlocs()
 do ilev=1,self%levs
  select case (trim(ic))
  case ('baroclinic-instability')
    ! Go to cartesian coordinates
    call lonlat_to_xy(lonlat(1,iloc),lonlat(2,iloc),x,y)

    ! Compute values for baroclinic instability
    if (self%vars%has('x')) call baroclinic_instability(x,y,z(ilev),'x',self%x(ilev,iloc))
    if (self%vars%has('q')) call baroclinic_instability(x,y,z(ilev),'q',self%q(ilev,iloc))
    if (self%vars%has('u')) call baroclinic_instability(x,y,z(ilev),'u',self%u(ilev,iloc))
    if (self%vars%has('v')) call baroclinic_instability(x,y,z(ilev),'v',self%v(ilev,iloc))
  case ('large-vortices')
    ! Go to cartesian coordinates
    call lonlat_to_xy(lonlat(1,iloc),lonlat(2,iloc),x,y)

    ! Compute values for large vortices
    if (self%vars%has('x')) call large_vortices(x,y,z(ilev),'x',self%x(ilev,iloc))
    if (self%vars%has('q')) call large_vortices(x,y,z(ilev),'q',self%q(ilev,iloc))
    if (self%vars%has('u')) call large_vortices(x,y,z(ilev),'u',self%u(ilev,iloc))
    if (self%vars%has('v')) call large_vortices(x,y,z(ilev),'v',self%v(ilev,iloc))
  case default
    call abor1_ftn('qg_gom_analytic_init: unknown initialization')
  endselect
 enddo
 self%z(:,iloc) = z(:)
enddo

call lonlat_field%final()

end subroutine qg_gom_analytic_init
! ------------------------------------------------------------------------------
end module qg_gom_mod
