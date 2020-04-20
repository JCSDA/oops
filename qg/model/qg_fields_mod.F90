! (C) Copyright 2009-2016 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

module qg_fields_mod

use atlas_module
use fckit_configuration_module, only: fckit_configuration
use datetime_mod
use duration_mod
use fckit_log_module,only: fckit_log
use iso_c_binding
use kinds
use missing_values_mod
use netcdf
!$ use omp_lib
use oops_variables_mod
use qg_constants_mod
use qg_convert_q_to_x_mod
use qg_convert_x_to_q_mod
use qg_convert_x_to_uv_mod
use qg_geom_mod
use qg_geom_iter_mod
use qg_gom_mod
use qg_interp_mod
use qg_locs_mod
use qg_tools_mod
use oops_variables_mod
use random_mod

implicit none

private
public :: qg_fields
public :: qg_fields_registry
public :: qg_fields_create,qg_fields_create_default,qg_fields_create_from_other,qg_fields_delete, &
        & qg_fields_zero,qg_fields_dirac,qg_fields_random, &
        & qg_fields_copy,qg_fields_self_add,qg_fields_self_sub,qg_fields_self_mul,qg_fields_axpy,qg_fields_self_schur, &
        & qg_fields_dot_prod,qg_fields_add_incr,qg_fields_diff_incr,qg_fields_change_resol,qg_fields_read_file, &
        & qg_fields_write_file,qg_fields_analytic_init,qg_fields_gpnorm,qg_fields_rms,qg_fields_sizes,qg_fields_vars, &
        & qg_fields_set_atlas,qg_fields_to_atlas, qg_fields_from_atlas, &
        & qg_fields_getpoint,qg_fields_setpoint,qg_fields_serialize,qg_fields_deserialize, qg_fields_check, &
        & qg_fields_check_resolution,qg_fields_check_variables
! ------------------------------------------------------------------------------
integer,parameter :: rseed = 7 !< Random seed (for reproducibility)

type :: qg_fields
  type(qg_geom),pointer :: geom                !< Geometry
  logical :: lq                                !< PV as main variable (streamfunction if false)
  logical :: lbc                               !< Boundaries are present
  real(kind_real),allocatable :: gfld3d(:,:,:) !< 3d field
  real(kind_real),allocatable :: x_north(:)    !< Streamfunction on northern wall
  real(kind_real),allocatable :: x_south(:)    !< Streamfunction on southern wall
  real(kind_real),allocatable :: q_north(:,:)  !< PV on northern wall
  real(kind_real),allocatable :: q_south(:,:)  !< PV on southern wall
end type qg_fields

#define LISTED_TYPE qg_fields

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry
type(registry_t) :: qg_fields_registry
! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
! Public
! ------------------------------------------------------------------------------
!> Linked list implementation
#include "oops/util/linkedList_c.f"
! ------------------------------------------------------------------------------
!> Create fields from geometry and variables
subroutine qg_fields_create(self,geom,vars,lbc)

implicit none

! Passed variables
type(qg_fields),intent(inout) :: self   !< Fields
type(qg_geom),target,intent(in) :: geom !< Geometry
type(oops_variables),intent(in) :: vars !< Variables
logical,intent(in) :: lbc               !< Boundaries flag

! Local variables
character(len=1024) :: record

! Associate geometry
self%geom => geom

! Set variables
if (vars%has('x') .and. vars%has('q')) then
  call abor1_ftn('qg_fields_create: x and q cannot be set as fields together')
elseif (vars%has('u') .or. vars%has('v')) then
  call abor1_ftn('qg_fieldsçcreate: u and v cannot be set as fields')
elseif (vars%has('x')) then
  self%lq = .false.
elseif (vars%has('q')) then
  self%lq = .true.
else
  call abor1_ftn('qg_fields_create: x or q should be set as fields')
endif

! Set boundaries
self%lbc = lbc

! Allocate 3d field
allocate(self%gfld3d(self%geom%nx,self%geom%ny,self%geom%nz))

! Allocate boundaries
if (self%lbc) then
  ! Allocation
  allocate(self%x_north(self%geom%nz))
  allocate(self%x_south(self%geom%nz))
  allocate(self%q_north(self%geom%nx,self%geom%nz))
  allocate(self%q_south(self%geom%nx,self%geom%nz))
endif

! Initialize
call qg_fields_zero(self)

end subroutine qg_fields_create

!> Create fields from geometry (x)
subroutine qg_fields_create_default(self,geom,lbc)

implicit none

! Passed variables
type(qg_fields),intent(inout) :: self   !< Fields
type(qg_geom),target,intent(in) :: geom !< Geometry
logical,intent(in) :: lbc               !< Boundaries flag

! Local variables
character(len=1024) :: record

! Associate geometry
self%geom => geom

! Set variables
self%lq = .false.

! Set boundaries
self%lbc = lbc

! Allocate 3d field
allocate(self%gfld3d(self%geom%nx,self%geom%ny,self%geom%nz))

! Allocate boundaries
if (self%lbc) then
  ! Allocation
  allocate(self%x_north(self%geom%nz))
  allocate(self%x_south(self%geom%nz))
  allocate(self%q_north(self%geom%nx,self%geom%nz))
  allocate(self%q_south(self%geom%nx,self%geom%nz))
endif

! Initialize
call qg_fields_zero(self)

end subroutine qg_fields_create_default
! ------------------------------------------------------------------------------
!> Create fields from another one
subroutine qg_fields_create_from_other(self,other)

implicit none

! Passed variables
type(qg_fields),intent(inout) :: self !< Fields
type(qg_fields),intent(in) :: other   !< Other fields

! Associate geometry
self%geom => other%geom

! Copy attributes
self%lq = other%lq
self%lbc = other%lbc

! Allocate 3d field
allocate(self%gfld3d(self%geom%nx,self%geom%ny,self%geom%nz))

! Allocate boundaries
if (self%lbc) then
  ! Allocation
  allocate(self%x_north(self%geom%nz))
  allocate(self%x_south(self%geom%nz))
  allocate(self%q_north(self%geom%nx,self%geom%nz))
  allocate(self%q_south(self%geom%nx,self%geom%nz))
endif

! Initialize
call qg_fields_zero(self)

end subroutine qg_fields_create_from_other
! ------------------------------------------------------------------------------
!> Delete fields
subroutine qg_fields_delete(self)

implicit none

! Passed variables
type(qg_fields),intent(inout) :: self !< Fields

! Release memory
if (allocated(self%gfld3d)) deallocate(self%gfld3d)
if (allocated(self%x_north)) deallocate(self%x_north)
if (allocated(self%x_south)) deallocate(self%x_south)
if (allocated(self%q_north)) deallocate(self%q_north)
if (allocated(self%q_south)) deallocate(self%q_south)

end subroutine qg_fields_delete
! ------------------------------------------------------------------------------
!> Set fields to zero
subroutine qg_fields_zero(self)

implicit none

! Passed variables
type(qg_fields),intent(inout) :: self

! Check field
call qg_fields_check(self)

! Set fields to zero
self%gfld3d = 0.0
if (self%lbc) then
  self%x_north = 0.0
  self%x_south = 0.0
  self%q_north = 0.0
  self%q_south = 0.0
endif

end subroutine qg_fields_zero
! ------------------------------------------------------------------------------
!> Set fields to Diracs
subroutine qg_fields_dirac(self,f_conf)

implicit none

! Passed variables
type(qg_fields),intent(inout) :: self          !< Fields
type(fckit_configuration),intent(in) :: f_conf !< FCKIT configuration

! Local variables
integer :: ndir,idir
integer,allocatable :: ixdir(:),iydir(:),izdir(:)

! Check field
call qg_fields_check(self)

! Get Diracs size
ndir = f_conf%get_size('ixdir')
if ((f_conf%get_size('iydir')/=ndir).or.(f_conf%get_size('izdir')/=ndir)) &
 & call abor1_ftn('qg_fields_dirac: inconsistent sizes for ixdir, iydir and izdir')

! Allocation
allocate(ixdir(ndir))
allocate(iydir(ndir))
allocate(izdir(ndir))

! Get Diracs positions
call f_conf%get_or_die("ixdir",ixdir)
call f_conf%get_or_die("iydir",iydir)
call f_conf%get_or_die("izdir",izdir)

! Check Diracs positions
if (any(ixdir<1).or.any(ixdir>self%geom%nx)) call abor1_ftn('qg_fields_dirac: invalid ixdir')
if (any(iydir<1).or.any(iydir>self%geom%ny)) call abor1_ftn('qg_fields_dirac: invalid iydir')
if (any(izdir<1).or.any(izdir>self%geom%nz)) call abor1_ftn('qg_fields_dirac: invalid izdir')

! Setup Diracs
call qg_fields_zero(self)
do idir=1,ndir
   self%gfld3d(ixdir(idir),iydir(idir),izdir(idir)) = 1.0
end do

end subroutine qg_fields_dirac
! ------------------------------------------------------------------------------
!> Generate random fields
subroutine qg_fields_random(self)

implicit none

! Passed variables
type(qg_fields),intent(inout) :: self !< Fields

! Check field
call qg_fields_check(self)

! Set at random value
call normal_distribution(self%gfld3d,0.0_kind_real,1.0_kind_real,rseed)

end subroutine qg_fields_random
! ------------------------------------------------------------------------------
!> Copy fields
subroutine qg_fields_copy(self,other)

implicit none

! Passed variables
type(qg_fields),intent(inout) :: self  !< Fields
type(qg_fields),intent(in)    :: other !< Other fields

! Check resolution
call qg_fields_check_resolution(self,other)
call qg_fields_check_variables(self,other)

! Copy field
self%gfld3d = other%gfld3d
if (self%lbc) then
  if (other%lbc) then
    self%x_north = other%x_north
    self%x_south = other%x_south
    self%q_north = other%q_north
    self%q_south = other%q_south
  else
    self%x_north = 0.0
    self%x_south = 0.0
    self%q_north = 0.0
    self%q_south = 0.0
  endif
endif

end subroutine qg_fields_copy
! ------------------------------------------------------------------------------
!> Add fields
subroutine qg_fields_self_add(self,rhs)

implicit none

! Passed variables
type(qg_fields),intent(inout) :: self !< Fields
type(qg_fields),intent(in)    :: rhs  !< Right-hand side

! Check resolution
call qg_fields_check_resolution(self,rhs)
call qg_fields_check_variables(self,rhs)

! Add field
self%gfld3d = self%gfld3d+rhs%gfld3d
if (self%lbc.and.rhs%lbc) then
  self%x_north = self%x_north+rhs%x_north
  self%x_south = self%x_south+rhs%x_south
  self%q_north = self%q_north+rhs%q_north
  self%q_south = self%q_south+rhs%q_south
endif

end subroutine qg_fields_self_add
! ------------------------------------------------------------------------------
!> Subtract fields
subroutine qg_fields_self_sub(self,rhs)

implicit none

! Passed variables
type(qg_fields),intent(inout) :: self !< Fields
type(qg_fields),intent(in)    :: rhs  !< Right-hand side

! Check resolution
call qg_fields_check_resolution(self,rhs)
call qg_fields_check_variables(self,rhs)

! Subtract field
self%gfld3d = self%gfld3d-rhs%gfld3d
if (self%lbc.and.rhs%lbc) then
  self%x_north = self%x_north-rhs%x_north
  self%x_south = self%x_south-rhs%x_south
  self%q_north = self%q_north-rhs%q_north
  self%q_south = self%q_south-rhs%q_south
endif

end subroutine qg_fields_self_sub
! ------------------------------------------------------------------------------
!> Multiply fields by a scalar
subroutine qg_fields_self_mul(self,zz)

implicit none

! Passed variables
type(qg_fields),intent(inout) :: self !< Fields
real(kind_real),intent(in) :: zz      !< Multiplier

! Check field
call qg_fields_check(self)

! Multiply with a scalar
self%gfld3d = zz*self%gfld3d
if (self%lbc) then
  self%x_north = zz*self%x_north
  self%x_south = zz*self%x_south
  self%q_north = zz*self%q_north
  self%q_south = zz*self%q_south
endif

end subroutine qg_fields_self_mul
! ------------------------------------------------------------------------------
!> Apply axpy operator to fields
subroutine qg_fields_axpy(self,zz,rhs)

implicit none

! Passed variables
type(qg_fields),intent(inout) :: self !< Fields
real(kind_real),intent(in) :: zz      !< Multiplier
type(qg_fields),intent(in) :: rhs     !< Right-hand side

! Check resolution
call qg_fields_check_resolution(self,rhs)
call qg_fields_check_variables(self,rhs)

! Apply apxy
self%gfld3d = self%gfld3d+zz*rhs%gfld3d
if (self%lbc.and.rhs%lbc) then
  self%x_north = self%x_north+zz*rhs%x_north
  self%x_south = self%x_south+zz*rhs%x_south
  self%q_north = self%q_north+zz*rhs%q_north
  self%q_south = self%q_south+zz*rhs%q_south
endif

end subroutine qg_fields_axpy
! ------------------------------------------------------------------------------
!> Schur product of fields
subroutine qg_fields_self_schur(self,rhs)

implicit none

! Passed variables
type(qg_fields),intent(inout) :: self !< Fields
type(qg_fields),intent(in)    :: rhs  !< Right-hand side

! Check resolution
call qg_fields_check_resolution(self,rhs)
call qg_fields_check_variables(self,rhs)

! Schur product
self%gfld3d = self%gfld3d*rhs%gfld3d
if (self%lbc.and.rhs%lbc) then
  self%x_north = self%x_north*rhs%x_north
  self%x_south = self%x_south*rhs%x_south
  self%q_north = self%q_north*rhs%q_north
  self%q_south = self%q_south*rhs%q_south
endif

end subroutine qg_fields_self_schur
! ------------------------------------------------------------------------------
!> Compute dot product for fields
subroutine qg_fields_dot_prod(fld1,fld2,zprod)

implicit none

! Passed variables
type(qg_fields),intent(in) :: fld1   !< First fields
type(qg_fields),intent(in) :: fld2   !< Second fields
real(kind_real),intent(out) :: zprod !< Dot product

! Check resolution
call qg_fields_check_resolution(fld1,fld2)
call qg_fields_check_variables(fld1,fld2)

! Compute dot product
zprod = sum(fld1%gfld3d*fld2%gfld3d)

end subroutine qg_fields_dot_prod
! ------------------------------------------------------------------------------
!> Add increment to fields
subroutine qg_fields_add_incr(self,rhs)

implicit none

! Passed variables
type(qg_fields),intent(inout) :: self !< Fields
type(qg_fields),intent(in)    :: rhs  !< Right-hand side

! Check fields
call qg_fields_check(self)
call qg_fields_check(rhs)

if (self%lq.eqv.rhs%lq) then
  if ((self%geom%nx==rhs%geom%nx).and.(self%geom%ny==rhs%geom%ny).and.(self%geom%nz==rhs%geom%nz)) then
    ! Same resolution
    self%gfld3d = self%gfld3d+rhs%gfld3d
  else
    ! Different resolutions
    call abor1_ftn('qg_fields_add_incr: not coded for low res increment yet')
  endif
else
  call abor1_ftn('qg_fields_add_incr: different variables')
endif

end subroutine qg_fields_add_incr
! ------------------------------------------------------------------------------
!> Compute increment from the difference of two fields
subroutine qg_fields_diff_incr(lhs,fld1,fld2)

implicit none

! Passed variables
type(qg_fields),intent(inout) :: lhs !< Left-hand side
type(qg_fields),intent(in) :: fld1   !< First fields
type(qg_fields),intent(in) :: fld2   !< Second fields

! Check resolution
call qg_fields_check_resolution(fld1,fld2)
call qg_fields_check_variables(fld1,fld2)
call qg_fields_check(lhs)

! Initialization
call qg_fields_zero(lhs)

if (lhs%lq.eqv.fld1%lq) then
  if ((fld1%geom%nx==lhs%geom%nx).and.(fld1%geom%ny==lhs%geom%ny).and.(fld1%geom%nz==lhs%geom%nz)) then
    ! Same resolution
    lhs%gfld3d = fld1%gfld3d-fld2%gfld3d
  else
    ! Different resolutions
    call abor1_ftn('qg_fields_diff_incr: not coded for low res increment yet')
  endif
else
  call abor1_ftn('qg_fields_diff_incr: different variables')
endif

end subroutine qg_fields_diff_incr
! ------------------------------------------------------------------------------
!> Change fields resolution
subroutine qg_fields_change_resol(fld,rhs)

implicit none

! Passed variables
type(qg_fields),intent(inout) :: fld !< Fields
type(qg_fields),intent(in)    :: rhs !< Right-hand side

integer :: ix,iy,iz
real(kind_real), allocatable, dimension(:,:,:) :: q1, q2

! Check fields
call qg_fields_check(fld)
call qg_fields_check(rhs)

if (fld%lq.eqv.rhs%lq) then
  if ((fld%geom%nx==rhs%geom%nx).and.(fld%geom%ny==rhs%geom%ny).and.(fld%geom%nz==rhs%geom%nz)) then
    ! Same resolution
    call qg_fields_copy(fld,rhs)
  else
    do ix = 1,fld%geom%nx
      do iy = 1,fld%geom%ny
        do iz = 1,fld%geom%nz
          call qg_interp_trilinear( rhs%geom,fld%geom%lon(ix,iy),fld%geom%lat(ix,iy),fld%geom%z(iz), &
                                    rhs%gfld3d,fld%gfld3d(ix,iy,iz) )
        enddo
      enddo
    enddo
    if (fld%lbc) then
      if (rhs%lbc) then
        allocate(q1(rhs%geom%nx,rhs%geom%ny,rhs%geom%nz))
        allocate(q2(fld%geom%nx,fld%geom%ny,fld%geom%nz))
        do iy = 1,rhs%geom%ny
          q1(:,iy,:) = rhs%q_south
        enddo
        do ix = 1,fld%geom%nx
          do iz = 1,fld%geom%nz
            call qg_interp_trilinear( rhs%geom,fld%geom%lon(ix,1),fld%geom%lat(ix,1),fld%geom%z(iz), &
                                      q1,q2(ix,1,iz) )
          enddo
        enddo
        fld%q_south = q2(:,1,:)
        do iy = 1,rhs%geom%ny
          q1(:,iy,:) = rhs%q_north
        enddo
        do ix = 1,fld%geom%nx
          do iz = 1,fld%geom%nz
            call qg_interp_trilinear( rhs%geom,fld%geom%lon(ix,1),fld%geom%lat(ix,1),fld%geom%z(iz), &
                                      q1,q2(ix,1,iz) )
          enddo
        enddo
        fld%q_north = q2(:,1,:)
        deallocate(q1,q2)
        fld%x_north = rhs%x_north
        fld%x_south = rhs%x_south
      else
        fld%x_north = 0.0
        fld%x_south = 0.0
        fld%q_north = 0.0
        fld%q_south = 0.0
      endif
    endif
  endif
else
  call abor1_ftn('qg_fields_change_resol: different variables')
endif

end subroutine qg_fields_change_resol
! ------------------------------------------------------------------------------
!> Read fields from file
subroutine qg_fields_read_file(fld,f_conf,vdate)
use string_utils

implicit none

! Passed variables
type(qg_fields),intent(inout) :: fld           !< Fields
type(fckit_configuration),intent(in) :: f_conf !< FCKIT configuration
type(datetime),intent(inout) :: vdate          !< Date and time

! Local variables
integer :: iread,nx,ny,nz,bc
integer :: ncid,nx_id,ny_id,nz_id,gfld3d_id,x_north_id,x_south_id,q_north_id,q_south_id
logical :: lbc
character(len=20) :: sdate
character(len=1024) :: record,filename
character(len=:),allocatable :: str

! Check field
call qg_fields_check(fld)

! Check whether the field should be invented or read from file
iread = 1
if (f_conf%has("read_from_file")) call f_conf%get_or_die("read_from_file",iread)

if (iread==0) then
  ! Invent field
  call fckit_log%warning('qg_fields_read_file: inventing field')
  call qg_fields_analytic_init(fld,f_conf,vdate)
else
  ! Read field from file

  ! Get filename
  call f_conf%get_or_die("filename",str)
  call swap_name_member(f_conf, str)
  filename = str
  call fckit_log%info('qg_fields_read_file: opening '//trim(filename))

  ! Initialize field
  call qg_fields_zero(fld)

  ! Open NetCDF file
  call ncerr(nf90_open(trim(filename)//'.nc',nf90_nowrite,ncid))

  ! Get dimensions ids
  call ncerr(nf90_inq_dimid(ncid,'nx',nx_id))
  call ncerr(nf90_inq_dimid(ncid,'ny',ny_id))
  call ncerr(nf90_inq_dimid(ncid,'nz',nz_id))

  ! Get dimensions
  call ncerr(nf90_inquire_dimension(ncid,nx_id,len=nx))
  call ncerr(nf90_inquire_dimension(ncid,ny_id,len=ny))
  call ncerr(nf90_inquire_dimension(ncid,nz_id,len=nz))

  ! Test dimensions consistency with the field geometry
  if ((nx/=fld%geom%nx).or.(ny/=fld%geom%ny).or.(nz/=fld%geom%nz)) then
    write (record,*) 'qg_fields_read_file: input fields have wrong dimensions: ',nx,ny,nz
    call fckit_log%error(record)
    write (record,*) 'qg_fields_read_file: expected: ',fld%geom%nx,fld%geom%ny,fld%geom%nz
    call fckit_log%error(record)
    call abor1_ftn('qg_fields_read_file: input fields have wrong dimensions')
  endif

  ! Get attributes
  call ncerr(nf90_get_att(ncid,nf90_global,'bc',bc))
  if (bc==0) then
    lbc = .false.
  elseif (bc==1) then
    lbc = .true.
  else
    call abor1_ftn('qg_fields_read_file: wrong bc value')
  end if
  call ncerr(nf90_get_att(ncid,nf90_global,'sdate',sdate))

  ! Test attributes consistency with the field
  if ((.not.lbc).and.fld%lbc) call abor1_ftn('qg_fields_read_file: LBC are missing in NetCDF file')

  ! Get variables ids
  if (fld%lq) then
    call ncerr(nf90_inq_varid(ncid,'q',gfld3d_id))
  else
    call ncerr(nf90_inq_varid(ncid,'x',gfld3d_id))
  endif
  if (fld%lbc) then
    call ncerr(nf90_inq_varid(ncid,'x_north',x_north_id))
    call ncerr(nf90_inq_varid(ncid,'x_south',x_south_id))
    call ncerr(nf90_inq_varid(ncid,'q_north',q_north_id))
    call ncerr(nf90_inq_varid(ncid,'q_south',q_south_id))
  endif

  ! Get variables
  call ncerr(nf90_get_var(ncid,gfld3d_id,fld%gfld3d))
  if (fld%lbc) then
    call ncerr(nf90_get_var(ncid,x_north_id,fld%x_north))
    call ncerr(nf90_get_var(ncid,x_south_id,fld%x_south))
    call ncerr(nf90_get_var(ncid,q_north_id,fld%q_north))
    call ncerr(nf90_get_var(ncid,q_south_id,fld%q_south))
  endif

  ! Close NetCDF file
  call ncerr(nf90_close(ncid))

  ! Set date
  call fckit_log%info('qg_fields_read_file: validity date is '//sdate)
  call datetime_set(sdate,vdate)
end if

! Check field
call qg_fields_check(fld)

end subroutine qg_fields_read_file
! ------------------------------------------------------------------------------
!> Write fields to file
subroutine qg_fields_write_file(fld,f_conf,vdate)

implicit none

! Passed variables
type(qg_fields),intent(in) :: fld              !< Fields
type(fckit_configuration),intent(in) :: f_conf !< FCKIT configuration
type(datetime),intent(in) :: vdate             !< Date and time

! Local variables
integer :: ncid,nx_id,ny_id,nz_id,lon_id,lat_id,z_id,area_id,heat_id,x_id,q_id,u_id,v_id
integer :: x_north_id,x_south_id,q_north_id,q_south_id
integer :: info
real(kind_real) :: x(fld%geom%nx,fld%geom%ny,fld%geom%nz),q(fld%geom%nx,fld%geom%ny,fld%geom%nz)
real(kind_real) :: u(fld%geom%nx,fld%geom%ny,fld%geom%nz),v(fld%geom%nx,fld%geom%ny,fld%geom%nz)
logical :: lwx,lwq,lwuv,ismpi,mainpe
character(len=20) :: sdate
character(len=1024) :: filename

! Check field
call qg_fields_check(fld)

! Compute streamfunction and potential vorticity
if (fld%lq) then
  q = fld%gfld3d
  lwq = .true.
  if (fld%lbc) then
    call convert_q_to_x(fld%geom,q,fld%x_north,fld%x_south,x)
    lwx = .true.
  else
    lwx = .false.
  endif
else
  x = fld%gfld3d
  lwx = .true.
  if (fld%lbc) then
    call convert_x_to_q(fld%geom,x,fld%x_north,fld%x_south,q)
    lwq = .true.
  else
    lwq = .false.
  endif
endif

! Compute wind
if (fld%lbc) then
  call convert_x_to_uv(fld%geom,x,fld%x_north,fld%x_south,u,v)
  lwuv = .true.
else
  lwuv = .false.
endif

! Set filename
filename = genfilename(f_conf,800,vdate)
call fckit_log%info('qg_fields_write_file: writing '//trim(filename))

! Set date
call datetime_to_string(vdate,sdate)

! Create NetCDF file
call ncerr(nf90_create(trim(filename)//'.nc',or(nf90_clobber,nf90_64bit_offset),ncid))

! Define dimensions
call ncerr(nf90_def_dim(ncid,'nx',fld%geom%nx,nx_id))
call ncerr(nf90_def_dim(ncid,'ny',fld%geom%ny,ny_id))
call ncerr(nf90_def_dim(ncid,'nz',fld%geom%nz,nz_id))

! Define attributes
if (fld%lbc) then
  call ncerr(nf90_put_att(ncid,nf90_global,'bc',1))
else
  call ncerr(nf90_put_att(ncid,nf90_global,'bc',0))
end if
call ncerr(nf90_put_att(ncid,nf90_global,'sdate',sdate))

! Define variables
call ncerr(nf90_def_var(ncid,'lon',nf90_double,(/nx_id,ny_id/),lon_id))
call ncerr(nf90_def_var(ncid,'lat',nf90_double,(/nx_id,ny_id/),lat_id))
call ncerr(nf90_def_var(ncid,'z',nf90_double,(/nz_id/),z_id))
call ncerr(nf90_def_var(ncid,'area',nf90_double,(/nx_id,ny_id/),area_id))
call ncerr(nf90_def_var(ncid,'heat',nf90_double,(/nx_id,ny_id/),heat_id))
if (lwx) then
  call ncerr(nf90_def_var(ncid,'x',nf90_double,(/nx_id,ny_id,nz_id/),x_id))
  call ncerr(nf90_put_att(ncid,x_id,'_FillValue',missing_value(1.0_kind_real)))
endif
if (lwq) then
  call ncerr(nf90_def_var(ncid,'q',nf90_double,(/nx_id,ny_id,nz_id/),q_id))
  call ncerr(nf90_put_att(ncid,q_id,'_FillValue',missing_value(1.0_kind_real)))
endif
if (lwuv) then
  call ncerr(nf90_def_var(ncid,'u',nf90_double,(/nx_id,ny_id,nz_id/),u_id))
  call ncerr(nf90_put_att(ncid,u_id,'_FillValue',missing_value(1.0_kind_real)))
  call ncerr(nf90_def_var(ncid,'v',nf90_double,(/nx_id,ny_id,nz_id/),v_id))
  call ncerr(nf90_put_att(ncid,v_id,'_FillValue',missing_value(1.0_kind_real)))
endif
if (fld%lbc) then
  call ncerr(nf90_def_var(ncid,'x_north',nf90_double,(/nz_id/),x_north_id))
  call ncerr(nf90_put_att(ncid,x_north_id,'_FillValue',missing_value(1.0_kind_real)))
  call ncerr(nf90_def_var(ncid,'x_south',nf90_double,(/nz_id/),x_south_id))
  call ncerr(nf90_put_att(ncid,x_south_id,'_FillValue',missing_value(1.0_kind_real)))
  call ncerr(nf90_def_var(ncid,'q_north',nf90_double,(/nx_id,nz_id/),q_north_id))
  call ncerr(nf90_put_att(ncid,q_north_id,'_FillValue',missing_value(1.0_kind_real)))
  call ncerr(nf90_def_var(ncid,'q_south',nf90_double,(/nx_id,nz_id/),q_south_id))
  call ncerr(nf90_put_att(ncid,q_south_id,'_FillValue',missing_value(1.0_kind_real)))
end if

! End definitions
call ncerr(nf90_enddef(ncid))

! Put variables
call ncerr(nf90_put_var(ncid,lon_id,fld%geom%lon))
call ncerr(nf90_put_var(ncid,lat_id,fld%geom%lat))
call ncerr(nf90_put_var(ncid,z_id,fld%geom%z))
call ncerr(nf90_put_var(ncid,area_id,fld%geom%area))
call ncerr(nf90_put_var(ncid,heat_id,fld%geom%heat))
if (lwx) call ncerr(nf90_put_var(ncid,x_id,x))
if (lwq) call ncerr(nf90_put_var(ncid,q_id,q))
if (lwuv) then
  call ncerr(nf90_put_var(ncid,u_id,u))
  call ncerr(nf90_put_var(ncid,v_id,v))
endif
if (fld%lbc) then
  call ncerr(nf90_put_var(ncid,x_north_id,fld%x_north))
  call ncerr(nf90_put_var(ncid,x_south_id,fld%x_south))
  call ncerr(nf90_put_var(ncid,q_north_id,fld%q_north))
  call ncerr(nf90_put_var(ncid,q_south_id,fld%q_south))
endif

! Close NetCDF file
call ncerr(nf90_close(ncid))

end subroutine qg_fields_write_file
! ------------------------------------------------------------------------------
!> Analytic initialization of fields
subroutine qg_fields_analytic_init(fld,f_conf,vdate)

implicit none

! Passed variables
type(qg_fields),intent(inout) :: fld           !< Fields
type(fckit_configuration),intent(in) :: f_conf !< FCKIT configuration
type(datetime),intent(inout) :: vdate          !< Date and time

! Local variables
integer :: ix,iy,iz
real(kind_real),allocatable :: x(:,:,:),q(:,:,:)
character(len=30) :: ic
character(len=20) :: sdate
character(len=:),allocatable :: str

! Check configuration
if (f_conf%has("analytic_init")) then
  call f_conf%get_or_die("analytic_init",str)
  ic = str
else
  ic = 'baroclinic-instability'
endif
call fckit_log%warning('qg_fields_analytic_init: '//trim(ic))

! Set date
call f_conf%get_or_die("date",str)
sdate = str
call fckit_log%info('qg_fields_analytic_init: validity date is '//sdate)
call datetime_set(sdate,vdate)

! Check boundaries
if (.not.fld%lbc) call abor1_ftn('qg_fields_analytic_init: boundaries required')

! Allocation
allocate(x(fld%geom%nx,fld%geom%ny,fld%geom%nz))
allocate(q(fld%geom%nx,fld%geom%ny,fld%geom%nz))

! Define state
select case (trim(ic))
case ('baroclinic-instability')
  ! Baroclinic instability
  do iz=1,fld%geom%nz
    do iy=1,fld%geom%ny
      do ix=1,fld%geom%nx
        call baroclinic_instability(fld%geom%x(ix),fld%geom%y(iy),fld%geom%z(iz),'x',x(ix,iy,iz))
      enddo
    enddo
    call baroclinic_instability(0.0_kind_real,domain_meridional,fld%geom%z(iz),'x',fld%x_north(iz))
    call baroclinic_instability(0.0_kind_real,0.0_kind_real,fld%geom%z(iz),'x',fld%x_south(iz))
  end do
case ('large-vortices')
  ! Large vortices
  do iz=1,fld%geom%nz
    do iy=1,fld%geom%ny
      do ix=1,fld%geom%nx
        call large_vortices(fld%geom%x(ix),fld%geom%y(iy),fld%geom%z(iz),'x',x(ix,iy,iz))
      enddo
    enddo
    call large_vortices(0.0_kind_real,domain_meridional,fld%geom%z(iz),'x',fld%x_north(iz))
    call large_vortices(0.0_kind_real,0.0_kind_real,fld%geom%z(iz),'x',fld%x_south(iz))
  enddo
case default
  call abor1_ftn ('qg_fields_analytic_init: unknown initialization')
end select

! Compute PV
call convert_x_to_q(fld%geom,x,fld%x_north,fld%x_south,q)
do iz=1,fld%geom%nz
  do ix=1,fld%geom%nx
    fld%q_south(ix,iz) = 2.0*q(ix,1,iz)-q(ix,2,iz)
  enddo
  do ix=1,fld%geom%nx
    fld%q_north(ix,iz) = 2.0*q(ix,fld%geom%ny,iz)-q(ix,fld%geom%ny-1,iz)
  enddo
end do

! Copy 3d field
if (fld%lq) then
  fld%gfld3d = q
else
  fld%gfld3d = x
endif

! Check field
call qg_fields_check(fld)

end subroutine qg_fields_analytic_init
! ------------------------------------------------------------------------------
!> Fields statistics
subroutine qg_fields_gpnorm(fld,nb,pstat)

implicit none

! Passed variables
type(qg_fields),intent(in) :: fld                !< Fields
integer,intent(in) :: nb                         !< Number of boundaries
real(kind_real),intent(inout) :: pstat(4*(1+nb)) !< Statistics

! Local variables
integer :: jj,js,jvb
real(kind_real) :: expo,stat(4,1+nb)

! Check field
call qg_fields_check(fld)

! Check number of stats
if ((fld%lbc.and.(nb/=2)).or.((.not.fld%lbc).and.(nb>0))) call abor1_ftn('qg_fields_gpnorm: error number of fields')

! 3d field
stat(2,1) = minval(fld%gfld3d)
stat(3,1) = maxval(fld%gfld3d)
stat(4,1) = sqrt(sum(fld%gfld3d**2)/real(fld%geom%nx*fld%geom%ny*fld%geom%nz,kind_real))
if (stat(4,1)>0.0) then
  expo = aint(log(stat(4,1))/log(10.0_kind_real))
  stat(1,1) = 10.0**expo
else
  stat(1,1) = 1.0
endif
stat(2:4,1) = stat(2:4,1)/stat(1,1)

! Boundaries
if (nb==2) then
  ! Streamfunction
  stat(2,2) = min(minval(fld%x_north),minval(fld%x_south))
  stat(3,2) = max(maxval(fld%x_north),maxval(fld%x_south))
  stat(4,2) = sqrt(sum(fld%x_north**2+fld%x_south**2)/real(2*fld%geom%nz,kind_real))
  if (stat(4,2)>0.0) then
    expo = aint(log(stat(4,2))/log(10.0_kind_real))
    stat(1,2) = 10.0**expo
  else
    stat(1,2) = 1.0
  endif
  stat(2:4,2) = stat(2:4,2)/stat(1,2)

  ! Potential vorticity
  stat(2,3) = min(minval(fld%q_north),minval(fld%q_south))
  stat(3,3) = max(maxval(fld%q_north),maxval(fld%q_south))
  stat(4,3) = sqrt(sum(fld%q_north**2+fld%q_south**2)/real(2*fld%geom%nx*fld%geom%nz,kind_real))
  if (stat(4,3)>0.0) then
    expo = aint(log(stat(4,3))/log(10.0_kind_real))
    stat(1,3) = 10.0**expo
  else
    stat(1,3) = 1.0
  endif
  stat(2:4,3) = stat(2:4,3)/stat(1,3)
endif

! Pack
jj = 0
do jvb=1,1+nb
  do js=1,4
    jj = jj+1
    pstat(jj) = stat(js,jvb)
  enddo
enddo

end subroutine qg_fields_gpnorm
! ------------------------------------------------------------------------------
!> Fields RMS
subroutine qg_fields_rms(fld,prms)

implicit none

! Passed variables
type(qg_fields),intent(in) :: fld   !< Fields
real(kind_real),intent(out) :: prms !< RMS

! Local variables
integer :: norm
real(kind_real) :: zz

! Check field
call qg_fields_check(fld)

! 3d field
zz = sum(fld%gfld3d**2)
norm = fld%geom%nx*fld%geom%ny*fld%geom%nz

! Boundaries
if (fld%lbc) then
  zz = zz+sum(fld%x_north**2+fld%x_south**2)+sum(fld%q_north**2+fld%q_south**2)
  norm = norm+2*(1+fld%geom%nx)*fld%geom%nz
end if

! Normalize
prms = sqrt(zz/real(norm,kind_real))

end subroutine qg_fields_rms
! ------------------------------------------------------------------------------
!> Get fields geometry
subroutine qg_fields_sizes(fld,nx,ny,nz,nb)

implicit none

! Passed variables
type(qg_fields),intent(in) :: fld !< Fields
integer,intent(out) :: nx         !< X size
integer,intent(out) :: ny         !< Y size
integer,intent(out) :: nz         !< Z size
integer,intent(out) :: nb         !< Number of boundaries

! Copy sizes
nx = fld%geom%nx
ny = fld%geom%ny
nz = fld%geom%nz
if (fld%lbc) then
  ! North and South boundaries
  nb = 2
else
  ! No boundaries
  nb = 0
endif

end subroutine qg_fields_sizes
! ------------------------------------------------------------------------------
!> Get fields variables
subroutine qg_fields_vars(fld,lq,lbc)

implicit none

! Passed variables
type(qg_fields),intent(in) :: fld !< Fields
integer,intent(out) :: lq         !< Potential vorticity flag
integer,intent(out) :: lbc        !< Boundaries flag

! Potential vorticity flag
if (fld%lq) then
  lq = 1
else
  lq = 0
endif

! Boundaries flag
if (fld%lbc) then
  lbc = 1
else
  lbc = 0
endif

end subroutine qg_fields_vars
! ------------------------------------------------------------------------------
!> Set ATLAS field
subroutine qg_fields_set_atlas(self,vars,vdate,afieldset)

implicit none

! Passed variables
type(qg_fields),intent(in) :: self              !< Fields
type(oops_variables),intent(in) :: vars         !< Variables
type(datetime),intent(in) :: vdate              !< Date and time
type(atlas_fieldset),intent(inout) :: afieldset !< ATLAS fieldset

! Local variables
character(len=20) :: sdate
character(len=1024) :: fieldname
type(atlas_field) :: afield

! Set date
call datetime_to_string(vdate,sdate)

! Get or create field
if (vars%has('x')) fieldname = 'x_'//sdate
if (vars%has('q')) fieldname = 'q_'//sdate
if (afieldset%has_field(trim(fieldname))) then
  ! Get afield
  afield = afieldset%field(trim(fieldname))
else
  ! Create field
  afield = self%geom%afunctionspace%create_field(name=trim(fieldname),kind=atlas_real(kind_real),levels=self%geom%nz)

  ! Add field
  call afieldset%add(afield)
end if

! Release pointer
call afield%final()

end subroutine qg_fields_set_atlas
! ------------------------------------------------------------------------------
!> Convert fields to ATLAS
subroutine qg_fields_to_atlas(self,vars,vdate,afieldset)

implicit none

! Passed variables
type(qg_fields),intent(in) :: self              !< Fields
type(oops_variables),intent(in) :: vars         !< Variables
type(datetime),intent(in) :: vdate              !< Date and time
type(atlas_fieldset),intent(inout) :: afieldset !< ATLAS fieldset

! Local variables
integer :: iv,i,j,k,node
integer(kind_int),pointer :: int_ptr_1(:),int_ptr_2(:,:)
real(kind_real) :: gfld3d(self%geom%nx,self%geom%ny,self%geom%nz)
real(kind_real),pointer :: real_ptr_1(:),real_ptr_2(:,:)
character(len=20) :: sdate
character(len=1024) :: fieldname
type(atlas_field) :: afield

! Set date
call datetime_to_string(vdate,sdate)

! Get variable
if (vars%has('x')) then
  if (self%lq) then
    call convert_q_to_x(self%geom,self%gfld3d,self%x_north,self%x_south,gfld3d)
  else
    gfld3d = self%gfld3d
  end if
end if
if (vars%has('q')) then
  if (self%lq) then
    gfld3d = self%gfld3d
  else
    call convert_x_to_q(self%geom,self%gfld3d,self%x_north,self%x_south,gfld3d)
  end if
end if

! Get or create field
if (vars%has('x')) fieldname = 'x_'//sdate
if (vars%has('q')) fieldname = 'q_'//sdate
if (afieldset%has_field(trim(fieldname))) then
  ! Get afield
  afield = afieldset%field(trim(fieldname))
else
  ! Create field
  afield = self%geom%afunctionspace%create_field(name=trim(fieldname),kind=atlas_real(kind_real),levels=self%geom%nz)

  ! Add field
  call afieldset%add(afield)
end if

! Copy field
call afield%data(real_ptr_2)
do k=1,self%geom%nz
  node = 0
  do j=self%geom%afunctionspace%j_begin(),self%geom%afunctionspace%j_end()
    do i=self%geom%afunctionspace%i_begin(j),self%geom%afunctionspace%i_end(j)
      node = node+1
      real_ptr_2(k,node) = gfld3d(i,j,k)
    enddo
  enddo
enddo

! Release pointer
call afield%final()

end subroutine qg_fields_to_atlas
! ------------------------------------------------------------------------------
!> Get fields from ATLAS
subroutine qg_fields_from_atlas(self,vars,vdate,afieldset)

implicit none

! Passed variables
type(qg_fields),intent(inout) :: self           !< Fields
type(oops_variables),intent(in) :: vars         !< Variables
type(datetime),intent(in) :: vdate              !< Date and time
type(atlas_fieldset),intent(inout) :: afieldset !< ATLAS fieldset

! Local variables
integer :: i,j,k,node
real(kind_real) :: gfld3d(self%geom%nx,self%geom%ny,self%geom%nz)
real(kind_real),pointer :: real_ptr_1(:),real_ptr_2(:,:)
character(len=1) :: cgrid
character(len=20) :: sdate
character(len=1024) :: fieldname
type(atlas_field) :: afield

! Set date
call datetime_to_string(vdate,sdate)

! Get field
if (vars%has('x')) fieldname = 'x_'//sdate
if (vars%has('q')) fieldname = 'q_'//sdate
afield = afieldset%field(trim(fieldname))

! Copy field
call afield%data(real_ptr_2)
do k=1,self%geom%nz
  node = 0
  do j=1,self%geom%ny
    do i=1,self%geom%nx
      node = node+1
      gfld3d(i,j,k) = real_ptr_2(k,node)
    enddo
  enddo
enddo

! Get variable
if (vars%has('x')) then
  if (self%lq) then
    call convert_x_to_q(self%geom,gfld3d,self%x_north,self%x_south,self%gfld3d)
  else
    self%gfld3d = gfld3d
  end if
end if
if (vars%has('q')) then
  if (self%lq) then
    self%gfld3d = gfld3d
  else
    call convert_x_to_q(self%geom,gfld3d,self%x_north,self%x_south,self%gfld3d)
  end if
end if

! Release pointer
call afield%final()

end subroutine qg_fields_from_atlas
! ------------------------------------------------------------------------------
!> Get points from fields
subroutine qg_fields_getpoint(fld,iter,nval,vals)

implicit none

! Passed variables
type(qg_fields),intent(in) :: fld           !< Fields
type(qg_geom_iter),intent(in) :: iter       !< Geometry iterator
integer,intent(in) :: nval                  !< Number of values
real(kind_real),intent(inout) :: vals(nval) !< Values

! Local variables
character(len=1024) :: record

! Check
if ((fld%geom%nx/=iter%geom%nx).or.(fld%geom%ny/=iter%geom%ny).or.(fld%geom%nz/=iter%geom%nz)) then
  write(record,*) 'qg_fields_getpoint: resolution inconsistency, ',fld%geom%nx,'/',fld%geom%ny,'/',fld%geom%nz, &
& ' and ',iter%geom%nx,'/',iter%geom%ny,'/',iter%geom%nz
  call abor1_ftn(record)
endif
if (fld%geom%nz/=nval) then
  write(record,*) 'qg_fields_getpoint: array sizes are different: ',fld%geom%nz,'/',nval
  call abor1_ftn(record)
endif

! Get values
vals = fld%gfld3d(iter%ilon,iter%ilat,:)

end subroutine qg_fields_getpoint
! ------------------------------------------------------------------------------
!> Set fields values at a specified gridpoint
subroutine qg_fields_setpoint(fld,iter,nval,vals)

implicit none

! Passed variables
type(qg_fields),intent(inout) :: fld     !< Fields
type(qg_geom_iter),intent(in) :: iter    !< Geometry iterator
integer,intent(in) :: nval               !< Number of values
real(kind_real),intent(in) :: vals(nval) !< Values

! Local variables
character(len=1024) :: record

! Check
if ((fld%geom%nx/=iter%geom%nx).or.(fld%geom%ny/=iter%geom%ny).or.(fld%geom%nz/=iter%geom%nz)) then
  write(record,*) 'qg_fields_getpoint: resolution inconsistency,',fld%geom%nx,'/',fld%geom%ny,'/',fld%geom%nz, &
& ' and ',iter%geom%nx,'/',iter%geom%ny,'/',iter%geom%nz
  call abor1_ftn(record)
endif
if (fld%geom%nz/=nval) then
  write(record,*) 'qg_fields_getpoint: array sizes are different:',fld%geom%nz,'/',nval
  call abor1_ftn(record)
endif

! Set values
fld%gfld3d(iter%ilon,iter%ilat,:) = vals

end subroutine qg_fields_setpoint
! ------------------------------------------------------------------------------
!> Serialize fields
subroutine qg_fields_serialize(fld,vsize,vect_fld)

implicit none

! Passed variables
type(qg_fields),intent(in) :: fld              !< Fields
integer,intent(in) :: vsize                    !< Size
real(kind_real),intent(out) :: vect_fld(vsize) !< Vector

! Local variables
integer :: ix,iy,iz,ind

! Initialize
ind = 0

! Copy
do iz = 1,fld%geom%nz
  do iy = 1,fld%geom%ny
    do ix = 1,fld%geom%nx
      ind = ind + 1
      vect_fld(ind) = fld%gfld3d(ix,iy,iz)
    enddo
  enddo
enddo

if (fld%lbc) then
  do iz=1,fld%geom%nz
    ind = ind + 1
    vect_fld(ind) = fld%x_north(iz)
    ind = ind + 1
    vect_fld(ind) = fld%x_south(iz)
    do ix=1,fld%geom%nx
      ind = ind + 1
      vect_fld(ind) = fld%q_north(ix,iz)
      ind = ind + 1
      vect_fld(ind) = fld%q_south(ix,iz)
    enddo
  enddo
endif

end subroutine qg_fields_serialize
! ------------------------------------------------------------------------------
!> Deserialize fields
subroutine qg_fields_deserialize(self,vsize,vect_fld,index)

implicit none

! Passed variables
type(qg_fields),intent(inout) :: self         !< Fields
integer,intent(in) :: vsize                   !< Size
real(kind_real),intent(in) :: vect_fld(vsize) !< Vector
integer,intent(inout) :: index                !< Index

! Local variables
integer :: ix,iy,iz

! 3d field
index = 1 + index
do iz = 1,self%geom%nz
  do iy = 1,self%geom%ny
    do ix = 1,self%geom%nx
      self%gfld3d(ix,iy,iz) = vect_fld(index)
      index = index+1
    enddo
  enddo
enddo

! Boundaries
if (self%lbc) then
  do iz=1,self%geom%nz
    index = index + 1
    self%x_north(iz) = vect_fld(index)
    index = index + 1
    self%x_south(iz) = vect_fld(index)
    do ix=1,self%geom%nx
      index = index + 1
      self%q_north(ix,iz) = vect_fld(index)
      index = index + 1
      self%q_south(ix,iz) = vect_fld(index)
    enddo
  enddo
endif

index = index - 1

end subroutine qg_fields_deserialize
! ------------------------------------------------------------------------------
!> Check fields
subroutine qg_fields_check(self)

implicit none

! Passed variables
type(qg_fields),intent(in) :: self !< Fields

! Local variables
logical :: bad
character(len=1024) :: record

! Initialization
bad = .false.

! Check 3d field
bad = bad.or.(.not.allocated(self%gfld3d))
bad = bad.or.(size(self%gfld3d,1)/=self%geom%nx)
bad = bad.or.(size(self%gfld3d,2)/=self%geom%ny)
bad = bad.or.(size(self%gfld3d,3)/=self%geom%nz)

! Check boundaries
if (self%lbc) then
  bad = bad.or.(.not.allocated(self%x_north))
  bad = bad.or.(.not.allocated(self%x_south))
  bad = bad.or.(.not.allocated(self%q_north))
  bad = bad.or.(.not.allocated(self%q_south))
  bad = bad.or.(size(self%x_north)/=self%geom%nz)
  bad = bad.or.(size(self%x_south)/=self%geom%nz)
  bad = bad.or.(size(self%q_north,1)/=self%geom%nx)
  bad = bad.or.(size(self%q_north,2)/=self%geom%nz)
  bad = bad.or.(size(self%q_south,1)/=self%geom%nx)
  bad = bad.or.(size(self%q_south,2)/=self%geom%nz)
else
  bad = bad.or.allocated(self%x_north)
  bad = bad.or.allocated(self%x_south)
  bad = bad.or.allocated(self%q_north)
  bad = bad.or.allocated(self%q_south)
endif

if (bad) then
  call fckit_log%info('qg_fields_check: field not consistent')
  write(record,*) '  nx,ny,nz,lbc = ',self%geom%nx,self%geom%ny,self%geom%nz,self%lbc
  call fckit_log%info(record)
  write(record,*) '  shape(gfld3d) = ',shape(self%gfld3d)
  call fckit_log%info(record)
  if (self%lbc) then
    write(record,*) '  shape(x_north) = ',shape(self%x_north)
    call fckit_log%info(record)
    write(record,*) '  shape(x_south) = ',shape(self%x_south)
    call fckit_log%info(record)
    write(record,*) '  shape(q_north) = ',shape(self%q_north)
    call fckit_log%info(record)
    write(record,*) '  shape(q_south) = ',shape(self%q_south)
    call fckit_log%info(record)
    call abor1_ftn ('qg_fields_check: field not consistent')
  endif
endif

end subroutine qg_fields_check
! ------------------------------------------------------------------------------
!> Check fields resolution
subroutine qg_fields_check_resolution(fld1,fld2)

implicit none

! Passed variables
type(qg_fields),intent(in) :: fld1 !< First fields
type(qg_fields),intent(in) :: fld2 !< Second fields

! Local variables
character(len=1024) :: record

! Check fields consistency
if ((fld1%geom%nx/=fld2%geom%nx).or.(fld1%geom%ny/=fld2%geom%ny).or.(fld1%geom%nz/=fld2%geom%nz)) then
  write(record,*) 'qg_fields_check_resolution: resolution inconsistency, ',fld1%geom%nx,'/',fld1%geom%ny,'/',fld1%geom%nz, &
& ' and ',fld2%geom%nx,'/',fld2%geom%ny,'/',fld2%geom%nz
  call abor1_ftn(record)
endif

! Check fields independently
call qg_fields_check(fld1)
call qg_fields_check(fld2)

end subroutine qg_fields_check_resolution
! ------------------------------------------------------------------------------
!> Check fields variables
subroutine qg_fields_check_variables(fld1,fld2)

implicit none

! Passed variables
type(qg_fields),intent(in) :: fld1 !< First fields
type(qg_fields),intent(in) :: fld2 !< Second fields

! Local variables
character(len=1024) :: record

! Check fields consistency
if (fld1%lq.neqv.fld2%lq) then
  write(record,*) 'qg_fields_check_variables: variables inconsistency, ',fld1%lq,' and ',fld2%lq
  call abor1_ftn(record)
endif

! Check fields independently
call qg_fields_check(fld1)
call qg_fields_check(fld2)

end subroutine qg_fields_check_variables
! ------------------------------------------------------------------------------
end module qg_fields_mod
