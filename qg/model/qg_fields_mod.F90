! (C) Copyright 2009-2016 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

module qg_fields_mod

use atlas_module, only: atlas_field, atlas_fieldset, atlas_real
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
use qg_convert_x_to_u_mod
use qg_convert_x_to_v_mod
use qg_geom_mod
use qg_geom_iter_mod
use qg_interp_mod
use qg_locs_mod
use qg_tools_mod
use random_mod

implicit none

private
public :: rseed
public :: qg_fields
public :: qg_fields_registry
public :: qg_fields_create,qg_fields_create_from_other,qg_fields_delete, &
        & qg_fields_zero,qg_fields_ones,qg_fields_dirac,qg_fields_random, &
        & qg_fields_copy,qg_fields_copy_lbc,qg_fields_self_add,qg_fields_self_sub,qg_fields_self_mul,qg_fields_axpy, &
        & qg_fields_self_schur,qg_fields_dot_prod,qg_fields_add_incr,qg_fields_diff_incr,qg_fields_change_resol, &
        & qg_fields_change_resol_ad,qg_fields_read_file,qg_fields_write_file,qg_fields_analytic_init,qg_fields_gpnorm, &
        & qg_fields_rms,qg_fields_sizes,qg_fields_lbc,qg_fields_to_fieldset,qg_fields_to_fieldset_ad,qg_fields_from_fieldset, &
        & qg_fields_getvals, qg_fields_getvalsad, &
        & qg_fields_getpoint,qg_fields_setpoint,qg_fields_serialize,qg_fields_deserialize, &
        & qg_fields_complete,qg_fields_check,qg_fields_check_resolution
! ------------------------------------------------------------------------------
integer,parameter :: rseed = 7 !< Random seed (for reproducibility)

type :: qg_fields
  type(qg_geom),pointer :: geom                !< Geometry
  type(oops_variables) :: vars                 !< List of variables
  logical :: lbc                               !< Boundaries are present
  real(kind_real),allocatable :: x(:,:,:)      !< Streamfunction
  real(kind_real),allocatable :: q(:,:,:)      !< Potential vorticity
  real(kind_real),allocatable :: u(:,:,:)      !< U wind
  real(kind_real),allocatable :: v(:,:,:)      !< V wind
  real(kind_real),allocatable :: x_north(:)    !< Streamfunction on northern wall
  real(kind_real),allocatable :: x_south(:)    !< Streamfunction on southern wall
  real(kind_real),allocatable :: q_north(:,:)  !< q on northern wall
  real(kind_real),allocatable :: q_south(:,:)  !< q on southern wall
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
type(oops_variables),intent(in) :: vars !< List of variables
logical,intent(in) :: lbc               !< Boundaries flag

! Local variables
character(len=1024) :: record

! Associate geometry
self%geom => geom

! Set variables
self%vars = oops_variables(vars)

! Set boundaries
self%lbc = lbc

! Allocate 3d fields
if (self%vars%has('x')) allocate(self%x(self%geom%nx,self%geom%ny,self%geom%nz))
if (self%vars%has('q')) allocate(self%q(self%geom%nx,self%geom%ny,self%geom%nz))
if (self%vars%has('u')) allocate(self%u(self%geom%nx,self%geom%ny,self%geom%nz))
if (self%vars%has('v')) allocate(self%v(self%geom%nx,self%geom%ny,self%geom%nz))

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
! ------------------------------------------------------------------------------
!> Create fields from another one
subroutine qg_fields_create_from_other(self,other,geom)

implicit none

! Passed variables
type(qg_fields),intent(inout) :: self   !< Fields
type(qg_fields),intent(in) :: other     !< Other fields
type(qg_geom),target,intent(in) :: geom !< Geometry

! Associate geometry
self%geom => geom

! Copy variables
self%vars = oops_variables(other%vars)

! Copy attributes
self%lbc = other%lbc

! Allocate 3d fields
if (allocated(other%x)) allocate(self%x(self%geom%nx,self%geom%ny,self%geom%nz))
if (allocated(other%q)) allocate(self%q(self%geom%nx,self%geom%ny,self%geom%nz))
if (allocated(other%u)) allocate(self%u(self%geom%nx,self%geom%ny,self%geom%nz))
if (allocated(other%v)) allocate(self%v(self%geom%nx,self%geom%ny,self%geom%nz))

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
call self%vars%destruct()
if (allocated(self%x)) deallocate(self%x)
if (allocated(self%q)) deallocate(self%q)
if (allocated(self%u)) deallocate(self%u)
if (allocated(self%v)) deallocate(self%v)
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
if (allocated(self%x)) self%x = 0.0_kind_real
if (allocated(self%q)) self%q = 0.0_kind_real
if (allocated(self%u)) self%u = 0.0_kind_real
if (allocated(self%v)) self%v = 0.0_kind_real
if (self%lbc) then
  self%x_north = 0.0_kind_real
  self%x_south = 0.0_kind_real
  self%q_north = 0.0_kind_real
  self%q_south = 0.0_kind_real
endif

end subroutine qg_fields_zero
! ------------------------------------------------------------------------------
!> Set fields to ones
subroutine qg_fields_ones(self)

implicit none

! Passed variables
type(qg_fields),intent(inout) :: self

! Check field
call qg_fields_check(self)

! Set fields to ones
if (allocated(self%x)) self%x = 1.0
if (allocated(self%q)) self%q = 1.0
if (allocated(self%u)) self%u = 1.0
if (allocated(self%v)) self%v = 1.0
if (self%lbc) then
  self%x_north = 1.0
  self%x_south = 1.0
  self%q_north = 1.0
  self%q_south = 1.0
endif

end subroutine qg_fields_ones
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
character(len=1) :: var
character(len=:),allocatable :: str

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
call f_conf%get_or_die("var",str)
var = str

! Check Diracs positions
if (any(ixdir<1).or.any(ixdir>self%geom%nx)) call abor1_ftn('qg_fields_dirac: invalid ixdir')
if (any(iydir<1).or.any(iydir>self%geom%ny)) call abor1_ftn('qg_fields_dirac: invalid iydir')
if (any(izdir<1).or.any(izdir>self%geom%nz)) call abor1_ftn('qg_fields_dirac: invalid izdir')

! Setup Diracs
call qg_fields_zero(self)
do idir=1,ndir
  select case (var)
  case ('x')
     if (.not.allocated(self%x)) call abor1_ftn('qg_fields_dirac: x should be allocated')
     self%x(ixdir(idir),iydir(idir),izdir(idir)) = 1.0
  case ('q')
     if (.not.allocated(self%q)) call abor1_ftn('qg_fields_dirac: q should be allocated')
     self%q(ixdir(idir),iydir(idir),izdir(idir)) = 1.0
  case default
     call abor1_ftn('qg_fields_dirac: wrong variable')
  endselect
enddo

! Complete other fields
call qg_fields_complete(self,var)

end subroutine qg_fields_dirac
! ------------------------------------------------------------------------------
!> Generate random fields
subroutine qg_fields_random(self, seed)

implicit none

! Passed variables
type(qg_fields),intent(inout) :: self !< Fields

! Local variables
logical :: allocate_x
integer,intent(in),optional :: seed   !< Optional seed

! Local variables
integer :: lseed

! Check field
call qg_fields_check(self)

! Local seed
lseed = rseed
if (present(seed)) lseed = seed

! Allocation
allocate_x = .not.allocated(self%x)
if (allocate_x) allocate(self%x(self%geom%nx,self%geom%ny,self%geom%nz))

! Set at random value
call normal_distribution(self%x,0.0_kind_real,1.0_kind_real,lseed)

! Complete other fields
call qg_fields_complete(self,'x')

! Release memory
if (allocate_x) deallocate(self%x)

end subroutine qg_fields_random
! ------------------------------------------------------------------------------
!> Copy fields
subroutine qg_fields_copy(self,other)

implicit none

! Passed variables
type(qg_fields),intent(inout) :: self   !< Fields
type(qg_fields),intent(in)    :: other  !< Other fields

! Check resolution
call qg_fields_check_resolution(self,other)

! Copy 3D field
if (allocated(self%x).and.allocated(other%x)) self%x = other%x
if (allocated(self%q).and.allocated(other%q)) self%q = other%q
if (allocated(self%u).and.allocated(other%u)) self%u = other%u
if (allocated(self%v).and.allocated(other%v)) self%v = other%v

! Copy LBC
call qg_fields_copy_lbc(self,other)

end subroutine qg_fields_copy
! ------------------------------------------------------------------------------
!> Copy fields LBC
subroutine qg_fields_copy_lbc(self,other)

implicit none

! Passed variables
type(qg_fields),intent(inout) :: self  !< Fields
type(qg_fields),intent(in)    :: other !< Other fields

! Check resolution
call qg_fields_check_resolution(self,other)

if (self%lbc) then
  if (other%lbc) then
    self%x_north = other%x_north
    self%x_south = other%x_south
    self%q_north = other%q_north
    self%q_south = other%q_south
  else
    self%x_north = 0.0_kind_real
    self%x_south = 0.0_kind_real
    self%q_north = 0.0_kind_real
    self%q_south = 0.0_kind_real
  endif
endif

end subroutine qg_fields_copy_lbc
! ------------------------------------------------------------------------------
!> Add fields
subroutine qg_fields_self_add(self,rhs)

implicit none

! Passed variables
type(qg_fields),intent(inout) :: self !< Fields
type(qg_fields),intent(in)    :: rhs  !< Right-hand side

! Check resolution
call qg_fields_check_resolution(self,rhs)

! Add field
if (allocated(self%x).and.allocated(rhs%x)) self%x = self%x + rhs%x
if (allocated(self%q).and.allocated(rhs%q)) self%q = self%q + rhs%q
if (allocated(self%u).and.allocated(rhs%u)) self%u = self%u + rhs%u
if (allocated(self%v).and.allocated(rhs%v)) self%v = self%v + rhs%v
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

! Subtract field
if (allocated(self%x).and.allocated(rhs%x)) self%x = self%x - rhs%x
if (allocated(self%q).and.allocated(rhs%q)) self%q = self%q - rhs%q
if (allocated(self%u).and.allocated(rhs%u)) self%u = self%u - rhs%u
if (allocated(self%v).and.allocated(rhs%v)) self%v = self%v - rhs%v
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
if (allocated(self%x)) self%x = zz * self%x
if (allocated(self%q)) self%q = zz * self%q
if (allocated(self%u)) self%u = zz * self%u
if (allocated(self%v)) self%v = zz * self%v
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

! Apply apxy
if (allocated(self%x).and.allocated(rhs%x)) self%x = self%x + zz * rhs%x
if (allocated(self%q).and.allocated(rhs%q)) self%q = self%q + zz * rhs%q
if (allocated(self%u).and.allocated(rhs%u)) self%u = self%u + zz * rhs%u
if (allocated(self%v).and.allocated(rhs%v)) self%v = self%v + zz * rhs%v
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

! Schur product
if (allocated(self%x).and.allocated(rhs%x)) self%x = self%x * rhs%x
if (allocated(self%q).and.allocated(rhs%q)) self%q = self%q * rhs%q
if (allocated(self%u).and.allocated(rhs%u)) self%u = self%u * rhs%u
if (allocated(self%v).and.allocated(rhs%v)) self%v = self%v * rhs%v
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

! Compute dot product
zprod = 0.0_kind_real
if (allocated(fld1%x).and.allocated(fld2%x)) zprod = zprod+sum(fld1%x*fld2%x)
if (allocated(fld1%q).and.allocated(fld2%q)) zprod = zprod+sum(fld1%q*fld2%q)
if (allocated(fld1%u).and.allocated(fld2%u)) zprod = zprod+sum(fld1%u*fld2%u)
if (allocated(fld1%v).and.allocated(fld2%v)) zprod = zprod+sum(fld1%v*fld2%v)

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

if ((self%geom%nx==rhs%geom%nx).and.(self%geom%ny==rhs%geom%ny).and.(self%geom%nz==rhs%geom%nz)) then
  ! Same resolution
  if (allocated(self%x).and.allocated(rhs%x)) self%x = self%x + rhs%x
  if (allocated(self%q).and.allocated(rhs%q)) self%q = self%q + rhs%q
  if (allocated(self%u).and.allocated(rhs%u)) self%u = self%u + rhs%u
  if (allocated(self%v).and.allocated(rhs%v)) self%v = self%v + rhs%v
else
  ! Different resolutions
  call abor1_ftn('qg_fields_add_incr: not coded for low res increment yet')
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
call qg_fields_check(lhs)

! Initialization
call qg_fields_zero(lhs)

if ((fld1%geom%nx==lhs%geom%nx).and.(fld1%geom%ny==lhs%geom%ny).and.(fld1%geom%nz==lhs%geom%nz)) then
  ! Same resolution
  if (allocated(lhs%x).and.allocated(fld1%x)) lhs%x = fld1%x - fld2%x
  if (allocated(lhs%q).and.allocated(fld1%q)) lhs%q = fld1%q - fld2%q
  if (allocated(lhs%u).and.allocated(fld1%u)) lhs%u = fld1%u - fld2%u
  if (allocated(lhs%v).and.allocated(fld1%v)) lhs%v = fld1%v - fld2%v
else
  ! Different resolutions
  call abor1_ftn('qg_fields_diff_incr: not coded for low res increment yet')
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

if ((fld%geom%nx==rhs%geom%nx).and.(fld%geom%ny==rhs%geom%ny).and.(fld%geom%nz==rhs%geom%nz)) then
  ! Same resolution
   call qg_fields_copy(fld,rhs)
else
  ! Trilinear interpolation
  do ix=1,fld%geom%nx
    do iy=1,fld%geom%ny
      do iz=1,fld%geom%nz
        if (allocated(rhs%x).and.allocated(fld%x)) then
          call qg_interp_trilinear(rhs%geom,fld%geom%lon(ix,iy),fld%geom%lat(ix,iy),fld%geom%z(iz), &
                                   rhs%x,fld%x(ix,iy,iz))
        endif
        if (allocated(rhs%q).and.allocated(fld%q)) then
          call qg_interp_trilinear(rhs%geom,fld%geom%lon(ix,iy),fld%geom%lat(ix,iy),fld%geom%z(iz), &
                                   rhs%q,fld%q(ix,iy,iz))
        endif
        if (allocated(rhs%u).and.allocated(fld%u)) then
          call qg_interp_trilinear(rhs%geom,fld%geom%lon(ix,iy),fld%geom%lat(ix,iy),fld%geom%z(iz), &
                                   rhs%u,fld%u(ix,iy,iz))
        endif
        if (allocated(rhs%v).and.allocated(fld%v)) then
          call qg_interp_trilinear(rhs%geom,fld%geom%lon(ix,iy),fld%geom%lat(ix,iy),fld%geom%z(iz), &
                                   rhs%v,fld%v(ix,iy,iz))
        endif
      enddo
    enddo
  enddo

  ! Deal with boundary conditions
  if (fld%lbc) then
    if (rhs%lbc) then
      allocate(q1(rhs%geom%nx,rhs%geom%ny,rhs%geom%nz))
      allocate(q2(fld%geom%nx,fld%geom%ny,fld%geom%nz))
      do iy=1,rhs%geom%ny
        q1(:,iy,:) = rhs%q_south
      enddo
      do ix=1,fld%geom%nx
        do iz=1,fld%geom%nz
          call qg_interp_trilinear(rhs%geom,fld%geom%lon(ix,1),fld%geom%lat(ix,1),fld%geom%z(iz),q1,q2(ix,1,iz) )
        enddo
      enddo
      fld%q_south = q2(:,1,:)
      do iy=1,rhs%geom%ny
        q1(:,iy,:) = rhs%q_north
      enddo
      do ix=1,fld%geom%nx
        do iz=1,fld%geom%nz
          call qg_interp_trilinear(rhs%geom,fld%geom%lon(ix,1),fld%geom%lat(ix,1),fld%geom%z(iz),q1,q2(ix,1,iz))
        enddo
      enddo
      fld%q_north = q2(:,1,:)
      deallocate(q1,q2)
      fld%x_north = rhs%x_north
      fld%x_south = rhs%x_south
    else
      fld%x_north = 0.0_kind_real
      fld%x_south = 0.0_kind_real
      fld%q_north = 0.0_kind_real
      fld%q_south = 0.0_kind_real
    endif
  endif
endif

end subroutine qg_fields_change_resol
! ------------------------------------------------------------------------------
!> Change fields resolution (adjoint)
subroutine qg_fields_change_resol_ad(fld, rhs)

implicit none

! Passed variables
type(qg_fields),intent(inout)    :: fld !< Fields (blank field of target geom)
type(qg_fields),intent(inout)    :: rhs !< Right-hand side (to be mapped to target geom via adjoint interpolation)

integer :: ix,iy,iz
real(kind_real), allocatable, dimension(:,:,:) :: q1, q2

if ((rhs%geom%nx==fld%geom%nx).and.(rhs%geom%ny==fld%geom%ny).and.(rhs%geom%nz==fld%geom%nz)) then
  ! Same resolution
  call qg_fields_copy(fld,rhs)

else
  ! Deal with boundary conditions (adjoint)
  if (rhs%lbc) then
    if (fld%lbc) then
      allocate(q1(fld%geom%nx,fld%geom%ny,fld%geom%nz))
      allocate(q2(rhs%geom%nx,rhs%geom%ny,rhs%geom%nz))
      do iy=1,fld%geom%ny
        q1(:,iy,:) = fld%q_south
      enddo
      do iz=rhs%geom%nz,1,-1
        do ix=rhs%geom%nx,1,-1
          call qg_interp_trilinear_ad(fld%geom,rhs%geom%lon(ix,1),rhs%geom%lat(ix,1),rhs%geom%z(iz),q2(ix,1,iz),q1)
        enddo
      enddo
      rhs%q_south = q2(:,1,:)
      do iy=1,fld%geom%ny
        q1(:,iy,:) = fld%q_north
      enddo
      do iz=rhs%geom%nz,1,-1
        do ix=rhs%geom%nx,1,-1
          call qg_interp_trilinear_ad(fld%geom,rhs%geom%lon(ix,1),rhs%geom%lat(ix,1),rhs%geom%z(iz),q2(ix,1,iz),q1)
        enddo
      enddo
      rhs%q_north = q2(:,1,:)
      deallocate(q1,q2)
      rhs%x_north = fld%x_north
      rhs%x_south = fld%x_south
    else
      rhs%x_north = 0.0_kind_real
      rhs%x_south = 0.0_kind_real
      rhs%q_north = 0.0_kind_real
      rhs%q_south = 0.0_kind_real
    endif
  endif

  ! Trilinear interpolation (adjoint)
  do iz=rhs%geom%nz,1,-1
    do iy=rhs%geom%ny,1,-1
      do ix=rhs%geom%nx,1,-1
        if (allocated(fld%x).and.allocated(rhs%x)) then
          call qg_interp_trilinear_ad(fld%geom,rhs%geom%lon(ix,iy),rhs%geom%lat(ix,iy),rhs%geom%z(iz), &
                                      rhs%x(ix,iy,iz),fld%x)
        endif
        if (allocated(fld%q).and.allocated(rhs%q)) then
          call qg_interp_trilinear_ad(fld%geom,rhs%geom%lon(ix,iy),rhs%geom%lat(ix,iy),rhs%geom%z(iz), &
                                      rhs%q(ix,iy,iz),fld%q)
        endif
        if (allocated(fld%u).and.allocated(rhs%u)) then
          call qg_interp_trilinear_ad(fld%geom,rhs%geom%lon(ix,iy),rhs%geom%lat(ix,iy),rhs%geom%z(iz), &
                                      rhs%u(ix,iy,iz),fld%u)
        endif
        if (allocated(fld%v).and.allocated(rhs%v)) then
          call qg_interp_trilinear_ad(fld%geom,rhs%geom%lon(ix,iy),rhs%geom%lat(ix,iy),rhs%geom%z(iz), &
                                      rhs%v(ix,iy,iz),fld%v)
        endif
      enddo
    enddo
  enddo
endif

end subroutine qg_fields_change_resol_ad
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
integer :: ncid,nx_id,ny_id,nz_id,x_id,q_id,u_id,v_id,x_north_id,x_south_id,q_north_id,q_south_id
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
  call swap_name_member(f_conf, str, 6)
  filename = str
  call fckit_log%info('qg_fields_read_file: opening '//trim(filename))

  ! Initialize field
  call qg_fields_zero(fld)

  ! Open NetCDF file
  call ncerr(nf90_open(trim(filename),nf90_nowrite,ncid))

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
  endif
  call ncerr(nf90_get_att(ncid,nf90_global,'sdate',sdate))

  ! Test attributes consistency with the field
  if ((.not.lbc).and.fld%lbc) call abor1_ftn('qg_fields_read_file: LBC are missing in NetCDF file')

  ! Get variables ids
  if (allocated(fld%x)) call ncerr(nf90_inq_varid(ncid,'x',x_id))
  if (allocated(fld%q)) call ncerr(nf90_inq_varid(ncid,'q',q_id))
  if (allocated(fld%u)) call ncerr(nf90_inq_varid(ncid,'u',u_id))
  if (allocated(fld%v)) call ncerr(nf90_inq_varid(ncid,'v',v_id))
  if (fld%lbc) then
    call ncerr(nf90_inq_varid(ncid,'x_north',x_north_id))
    call ncerr(nf90_inq_varid(ncid,'x_south',x_south_id))
    call ncerr(nf90_inq_varid(ncid,'q_north',q_north_id))
    call ncerr(nf90_inq_varid(ncid,'q_south',q_south_id))
  endif

  ! Get variables
  if (allocated(fld%x)) call ncerr(nf90_get_var(ncid,x_id,fld%x))
  if (allocated(fld%q)) call ncerr(nf90_get_var(ncid,q_id,fld%q))
  if (allocated(fld%u)) call ncerr(nf90_get_var(ncid,u_id,fld%u))
  if (allocated(fld%v)) call ncerr(nf90_get_var(ncid,v_id,fld%v))
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
endif

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
character(len=:),allocatable :: str
character(len=20) :: sdate
character(len=1024) :: typ,filename
type(oops_variables) :: vars
type(qg_fields) :: fld_io
logical :: date_cols

! Get output type
call f_conf%get_or_die("type",str)
typ = str

! Check field
call qg_fields_check(fld)

! Get all variables
vars = oops_variables()
call vars%push_back('x')
call vars%push_back('q')
call vars%push_back('u')
call vars%push_back('v')
call qg_fields_create(fld_io,fld%geom,vars,.true.)
call qg_fields_copy_lbc(fld_io,fld)
if ((trim(typ)=='diag') .or. (trim(typ)=='in')) then
  ! Diagnostic or increment file: don't complete fields
  if (allocated(fld%x)) then
    fld_io%x = fld%x
  else
    fld_io%x = missing_value(1.0_kind_real)
  endif
  if (allocated(fld%q)) then
    fld_io%q = fld%q
  else
    fld_io%q = missing_value(1.0_kind_real)
  endif
  if (allocated(fld%u)) then
    fld_io%u = fld%u
  else
    fld_io%u = missing_value(1.0_kind_real)
  endif
  if (allocated(fld%v)) then
    fld_io%v = fld%v
  else
    fld_io%v = missing_value(1.0_kind_real)
  endif
else
  ! Usual file: complete fields
  if (allocated(fld%x)) then
    fld_io%x = fld%x
    call qg_fields_complete(fld_io,'x')
  elseif (allocated(fld%q)) then
    fld_io%q = fld%q
    call qg_fields_complete(fld_io,'q')
  else
    call abor1_ftn('qg_fields_write_file: x or q required')
  endif
endif

! Get date IO format (colons or not?)
date_cols = .true.
if (f_conf%has("date colons")) then
  call f_conf%get_or_die("date colons", date_cols)
end if

! Set filename
filename = genfilename(f_conf,800,vdate,date_cols)
call fckit_log%info('qg_fields_write_file: writing '//trim(filename))

! Set date
call datetime_to_string(vdate,sdate)

! Create NetCDF file
call ncerr(nf90_create(trim(filename),ior(nf90_clobber,nf90_64bit_offset),ncid))

! Define dimensions
call ncerr(nf90_def_dim(ncid,'nx',fld%geom%nx,nx_id))
call ncerr(nf90_def_dim(ncid,'ny',fld%geom%ny,ny_id))
call ncerr(nf90_def_dim(ncid,'nz',fld%geom%nz,nz_id))

! Define attributes
if (fld%lbc) then
  call ncerr(nf90_put_att(ncid,nf90_global,'bc',1))
else
  call ncerr(nf90_put_att(ncid,nf90_global,'bc',0))
endif
call ncerr(nf90_put_att(ncid,nf90_global,'sdate',sdate))

! Define variables
call ncerr(nf90_def_var(ncid,'lon',nf90_double,(/nx_id,ny_id/),lon_id))
call ncerr(nf90_def_var(ncid,'lat',nf90_double,(/nx_id,ny_id/),lat_id))
call ncerr(nf90_def_var(ncid,'z',nf90_double,(/nz_id/),z_id))
call ncerr(nf90_def_var(ncid,'area',nf90_double,(/nx_id,ny_id/),area_id))
call ncerr(nf90_def_var(ncid,'heat',nf90_double,(/nx_id,ny_id/),heat_id))
call ncerr(nf90_def_var(ncid,'x',nf90_double,(/nx_id,ny_id,nz_id/),x_id))
call ncerr(nf90_put_att(ncid,x_id,'_FillValue',missing_value(1.0_kind_real)))
call ncerr(nf90_def_var(ncid,'q',nf90_double,(/nx_id,ny_id,nz_id/),q_id))
call ncerr(nf90_put_att(ncid,q_id,'_FillValue',missing_value(1.0_kind_real)))
call ncerr(nf90_def_var(ncid,'u',nf90_double,(/nx_id,ny_id,nz_id/),u_id))
call ncerr(nf90_put_att(ncid,u_id,'_FillValue',missing_value(1.0_kind_real)))
call ncerr(nf90_def_var(ncid,'v',nf90_double,(/nx_id,ny_id,nz_id/),v_id))
call ncerr(nf90_put_att(ncid,v_id,'_FillValue',missing_value(1.0_kind_real)))
if (fld%lbc) then
  call ncerr(nf90_def_var(ncid,'x_north',nf90_double,(/nz_id/),x_north_id))
  call ncerr(nf90_put_att(ncid,x_north_id,'_FillValue',missing_value(1.0_kind_real)))
  call ncerr(nf90_def_var(ncid,'x_south',nf90_double,(/nz_id/),x_south_id))
  call ncerr(nf90_put_att(ncid,x_south_id,'_FillValue',missing_value(1.0_kind_real)))
  call ncerr(nf90_def_var(ncid,'q_north',nf90_double,(/nx_id,nz_id/),q_north_id))
  call ncerr(nf90_put_att(ncid,q_north_id,'_FillValue',missing_value(1.0_kind_real)))
  call ncerr(nf90_def_var(ncid,'q_south',nf90_double,(/nx_id,nz_id/),q_south_id))
  call ncerr(nf90_put_att(ncid,q_south_id,'_FillValue',missing_value(1.0_kind_real)))
endif

! End definitions
call ncerr(nf90_enddef(ncid))

! Put variables
call ncerr(nf90_put_var(ncid,lon_id,fld%geom%lon))
call ncerr(nf90_put_var(ncid,lat_id,fld%geom%lat))
call ncerr(nf90_put_var(ncid,z_id,fld%geom%z))
call ncerr(nf90_put_var(ncid,area_id,fld%geom%area))
call ncerr(nf90_put_var(ncid,heat_id,fld%geom%heat))
call ncerr(nf90_put_var(ncid,x_id,fld_io%x))
call ncerr(nf90_put_var(ncid,q_id,fld_io%q))
call ncerr(nf90_put_var(ncid,u_id,fld_io%u))
call ncerr(nf90_put_var(ncid,v_id,fld_io%v))
if (fld%lbc) then
  call ncerr(nf90_put_var(ncid,x_north_id,fld%x_north))
  call ncerr(nf90_put_var(ncid,x_south_id,fld%x_south))
  call ncerr(nf90_put_var(ncid,q_north_id,fld%q_north))
  call ncerr(nf90_put_var(ncid,q_south_id,fld%q_south))
endif

! Close NetCDF file
call ncerr(nf90_close(ncid))

! Release memory
call vars%destruct()
call qg_fields_delete(fld_io)

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
real(kind_real) :: uval
real(kind_real),allocatable :: x(:,:,:),q(:,:,:)
character(len=30) :: ic
character(len=20) :: sdate
character(len=:),allocatable :: str

! Check configuration
if (f_conf%has("analytic init.method")) then
  call f_conf%get_or_die("analytic init.method",str)
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
  enddo
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
case ('uniform_field')
  ! Uniform field
  call f_conf%get_or_die("uval",uval)
  x = uval
case default
  call abor1_ftn('qg_fields_analytic_init: unknown initialization')
endselect

! Compute q
call convert_x_to_q(fld%geom,x,fld%x_north,fld%x_south,q)
do iz=1,fld%geom%nz
  do ix=1,fld%geom%nx
    fld%q_south(ix,iz) = 2.0*q(ix,1,iz)-q(ix,2,iz)
  enddo
  do ix=1,fld%geom%nx
    fld%q_north(ix,iz) = 2.0*q(ix,fld%geom%ny,iz)-q(ix,fld%geom%ny-1,iz)
  enddo
enddo

! Copy 3d field and ensure consistency
if (allocated(x)) then
  fld%x = x
  call qg_fields_complete(fld,'x')
elseif (allocated(q)) then
  fld%q = q
  call qg_fields_complete(fld,'q')
else
  call abor1_ftn('qg_fields_analytic_init: x or q required')
endif

! Check field
call qg_fields_check(fld)

end subroutine qg_fields_analytic_init
! ------------------------------------------------------------------------------
!> Fields statistics
subroutine qg_fields_gpnorm(fld,vpresent,vmin,vmax,vrms)

implicit none

! Passed variables
type(qg_fields),intent(in) :: fld        !< Fields
integer,intent(inout) :: vpresent(6)     !< Variables presence flag
real(kind_real),intent(inout) :: vmin(6) !< Variables minimum
real(kind_real),intent(inout) :: vmax(6) !< Variables maximum
real(kind_real),intent(inout) :: vrms(6) !< Variables RMS

! Check field
call qg_fields_check(fld)

! Initialization
vpresent = 0
vmin = 0.0_kind_real
vmax = 0.0_kind_real
vrms = 0.0_kind_real

! 3d fields
if (allocated(fld%x)) then
  vpresent(1) = 1
  vmin(1) = minval(fld%x)
  vmax(1) = maxval(fld%x)
  vrms(1) = sqrt(sum(fld%x**2)/real(fld%geom%nx*fld%geom%ny*fld%geom%nz,kind_real))
endif
if (allocated(fld%q)) then
  vpresent(2) = 1
  vmin(2) = minval(fld%q)
  vmax(2) = maxval(fld%q)
  vrms(2) = sqrt(sum(fld%q**2)/real(fld%geom%nx*fld%geom%ny*fld%geom%nz,kind_real))
endif
if (allocated(fld%u)) then
  vpresent(3) = 1
  vmin(3) = minval(fld%u)
  vmax(3) = maxval(fld%u)
  vrms(3) = sqrt(sum(fld%u**2)/real(fld%geom%nx*fld%geom%ny*fld%geom%nz,kind_real))
endif
if (allocated(fld%v)) then
  vpresent(4) = 1
  vmin(4) = minval(fld%v)
  vmax(4) = maxval(fld%v)
  vrms(4) = sqrt(sum(fld%v**2)/real(fld%geom%nx*fld%geom%ny*fld%geom%nz,kind_real))
endif

! Boundaries
if (fld%lbc) then
  ! Streamfunction
  vpresent(5) = 1
  vmin(5) = min(minval(fld%x_north),minval(fld%x_south))
  vmax(5) = max(maxval(fld%x_north),maxval(fld%x_south))
  vrms(5) = sqrt(sum(fld%x_north**2+fld%x_south**2)/real(2*fld%geom%nz,kind_real))

  ! Potential vorticity
  vpresent(6) = 1
  vmin(6) = min(minval(fld%q_north),minval(fld%q_south))
  vmax(6) = max(maxval(fld%q_north),maxval(fld%q_south))
  vrms(6) = sqrt(sum(fld%q_north**2+fld%q_south**2)/real(2*fld%geom%nx*fld%geom%nz,kind_real))
endif

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

! Check field
call qg_fields_check(fld)

! Initialization
prms = 0.0_kind_real
norm = 0.0_kind_real

! 3d fields
if (allocated(fld%x)) then
   prms = prms+sum(fld%x**2)
   norm = norm+real(fld%geom%nx*fld%geom%ny*fld%geom%nz,kind_real)
endif
if (allocated(fld%q)) then
   prms = prms+sum(fld%q**2)
   norm = norm+real(fld%geom%nx*fld%geom%ny*fld%geom%nz,kind_real)
endif
if (allocated(fld%u)) then
   prms = prms+sum(fld%u**2)
   norm = norm+real(fld%geom%nx*fld%geom%ny*fld%geom%nz,kind_real)
endif
if (allocated(fld%v)) then
   prms = prms+sum(fld%v**2)
   norm = norm+real(fld%geom%nx*fld%geom%ny*fld%geom%nz,kind_real)
endif

! Boundaries
if (fld%lbc) then
  prms = prms+sum(fld%x_north**2+fld%x_south**2)+sum(fld%q_north**2+fld%q_south**2)
  norm = norm+real(2*(1+fld%geom%nx)*fld%geom%nz,kind_real)
endif

! Normalize and square-root
prms = sqrt(prms/norm)

end subroutine qg_fields_rms
! ------------------------------------------------------------------------------
!> Get fields geometry
subroutine qg_fields_sizes(fld,nx,ny,nz)

implicit none

! Passed variables
type(qg_fields),intent(in) :: fld !< Fields
integer,intent(out) :: nx         !< X size
integer,intent(out) :: ny         !< Y size
integer,intent(out) :: nz         !< Z size

! Copy sizes
nx = fld%geom%nx
ny = fld%geom%ny
nz = fld%geom%nz

end subroutine qg_fields_sizes
! ------------------------------------------------------------------------------
!> Get LBC presence
subroutine qg_fields_lbc(fld,lbc)

implicit none

! Passed variables
type(qg_fields),intent(in) :: fld !< Fields
integer,intent(out) :: lbc        !< LBC presence

! Check field
call qg_fields_check(fld)

if (fld%lbc) then
  lbc = 1
else
  lbc = 0
endif

end subroutine qg_fields_lbc
! ------------------------------------------------------------------------------
!> Convert Fieldset to fields
subroutine qg_fields_to_fieldset(self,afieldset)

implicit none

! Passed variables
type(qg_fields),intent(in) :: self              !< Fields
type(atlas_fieldset),intent(inout) :: afieldset !< FieldSet

! Local variables
integer :: jvar,ix,iy,iz,inode
real(kind_real),pointer :: ptr(:,:)
character(len=1024) :: fieldname
type(atlas_field) :: afield

! Get variable
do jvar=1,self%vars%nvars()
   fieldname = self%vars%variable(jvar)
   if (afieldset%has_field(trim(fieldname))) then
     ! Get afield
     afield = afieldset%field(trim(fieldname))
   else
     ! Create field
     afield = self%geom%afunctionspace%create_field(name=trim(fieldname),kind=atlas_real(kind_real),levels=self%geom%nz)

     ! Add field
     call afieldset%add(afield)
   endif

   ! Copy field
   call afield%data(ptr)
   do iz=1,self%geom%nz
     inode = 0
     do iy=1,self%geom%ny
       do ix=1,self%geom%nx
         inode = inode+1
         select case (trim(fieldname))
         case ('x')
           ptr(iz,inode) = self%x(ix,iy,iz)
         case ('q')
           ptr(iz,inode) = self%q(ix,iy,iz)
         case ('u')
           ptr(iz,inode) = self%u(ix,iy,iz)
         case ('v')
           ptr(iz,inode) = self%v(ix,iy,iz)
         case ('z')
           ptr(iz,inode) = self%geom%z(iz)
         case default
           call abor1_ftn('qg_fields_to_fieldset: wrong variable')
         endselect
       enddo
     enddo
   enddo

   ! Set dirty to false: QG is a serial model so there are no halos
   call afield%set_dirty(.false.)

   ! Release pointer
   call afield%final()
enddo

end subroutine qg_fields_to_fieldset
! ------------------------------------------------------------------------------
!> Convert fields to Fieldset (adjoint)
subroutine qg_fields_to_fieldset_ad(self,afieldset)

implicit none

! Passed variables
type(qg_fields),intent(inout) :: self           !< Fields
type(atlas_fieldset),intent(inout) :: afieldset !< ATLAS fieldset

! Local variables
integer :: jvar,ix,iy,iz,inode
real(kind_real),pointer :: ptr(:,:)
character(len=1024) :: fieldname
type(atlas_field) :: afield

! Get variable
do jvar=1,self%vars%nvars()
   ! Get afield
   fieldname = self%vars%variable(jvar)
   afield = afieldset%field(trim(fieldname))

   ! Copy field
   call afield%data(ptr)
   do iz=1,self%geom%nz
     inode = 0
     do iy=1,self%geom%ny
       do ix=1,self%geom%nx
         inode = inode+1
         select case (trim(fieldname))
         case ('x')
           self%x(ix,iy,iz) = ptr(iz,inode)
         case ('q')
           self%q(ix,iy,iz) = ptr(iz,inode)
         case ('u')
           self%u(ix,iy,iz) = ptr(iz,inode)
         case ('v')
           self%v(ix,iy,iz) = ptr(iz,inode)
         case ('z')
           ! do nothing
         case default
           call abor1_ftn('qg_fields_to_fieldset_ad: wrong variable')
         endselect
       enddo
     enddo
   enddo

   ! Release pointer
   call afield%final()
enddo

end subroutine qg_fields_to_fieldset_ad
! ------------------------------------------------------------------------------
!> Convert Fieldset to fields
subroutine qg_fields_from_fieldset(self,afieldset)

implicit none

! Passed variables
type(qg_fields),intent(inout) :: self           !< Fields
type(atlas_fieldset),intent(inout) :: afieldset !< ATLAS fieldset

! Local variables
integer :: jvar,ix,iy,iz,inode
real(kind_real),pointer :: ptr(:,:)
character(len=1024) :: fieldname
type(atlas_field) :: afield

! Get variable
do jvar=1,self%vars%nvars()
   ! Get afield
   fieldname = self%vars%variable(jvar)
   afield = afieldset%field(trim(fieldname))

   ! Copy field
   call afield%data(ptr)
   do iz=1,self%geom%nz
     inode = 0
     do iy=1,self%geom%ny
       do ix=1,self%geom%nx
         inode = inode+1
         select case (trim(fieldname))
         case ('x')
           self%x(ix,iy,iz) = ptr(iz,inode)
         case ('q')
           self%q(ix,iy,iz) = ptr(iz,inode)
         case ('u')
           self%u(ix,iy,iz) = ptr(iz,inode)
         case ('v')
           self%v(ix,iy,iz) = ptr(iz,inode)
         case ('z')
           ! do nothing
         case default
           call abor1_ftn('qg_fields_from_fieldset: wrong variable')
         endselect
       enddo
     enddo
   enddo

   ! Release pointer
   call afield%final()
enddo

end subroutine qg_fields_from_fieldset
! ------------------------------------------------------------------------------
subroutine qg_fields_getvals(self, vars, lats, lons, vals)

implicit none
type(qg_fields),intent(in)      :: self
type(oops_variables),intent(in) :: vars
real(kind_real), intent(in)     :: lats(:)
real(kind_real), intent(in)     :: lons(:)
real(c_double), intent(inout)   :: vals(:)

integer :: nlocs, levs, jvar, jloc, ii
character(len=1024) :: fname

call qg_fields_check(self)

nlocs = size(lats)
levs = self%geom%nz

ii = 0
do jvar=1,vars%nvars()
  fname = vars%variable(jvar)
  select case (trim(fname))
  case ('x')
    if (.not.allocated(self%x)) call abor1_ftn('qg_fields_getvals: x not allocated')
    do jloc=1,nlocs
      call qg_interp_bilinear(self%geom,lons(jloc),lats(jloc),self%x(:,:,:),vals(ii+1:ii+levs))
      ii = ii + levs
    enddo
  case ('q')
    if (.not.allocated(self%q)) call abor1_ftn('qg_fields_getvals: q not allocated')
    do jloc=1,nlocs
      call qg_interp_bilinear(self%geom,lons(jloc),lats(jloc),self%q(:,:,:),vals(ii+1:ii+levs))
      ii = ii + levs
    enddo
  case ('u')
    if (.not.allocated(self%u)) call abor1_ftn('qg_fields_getvals: u not allocated')
    do jloc=1,nlocs
      call qg_interp_bilinear(self%geom,lons(jloc),lats(jloc),self%u(:,:,:),vals(ii+1:ii+levs))
      ii = ii + levs
    enddo
  case ('v')
    if (.not.allocated(self%v)) call abor1_ftn('qg_fields_getvals: v not allocated')
    do jloc=1,nlocs
      call qg_interp_bilinear(self%geom,lons(jloc),lats(jloc),self%v(:,:,:),vals(ii+1:ii+levs))
      ii = ii + levs
    enddo
  case ('z')
    if (.not.allocated(self%geom%z)) call abor1_ftn('qg_fields_getvals: z not allocated')
    do jloc=1,nlocs
      vals(ii+1:ii+levs) = self%geom%z(:)
      ii = ii + levs
    enddo
  case default
    call abor1_ftn('qg_fields_getvals: wrong input variable')
  endselect
enddo
if (size(vals) /= ii) call abor1_ftn('qg_fields_getvals: error size')

end subroutine qg_fields_getvals
! ------------------------------------------------------------------------------
subroutine qg_fields_getvalsad(self, vars, lats, lons, vals)

implicit none
type(qg_fields),intent(inout)   :: self
type(oops_variables),intent(in) :: vars
real(kind_real), intent(in)     :: lats(:)
real(kind_real), intent(in)     :: lons(:)
real(c_double), intent(in)      :: vals(:)

integer :: nlocs, levs, jvar, jloc, ii
character(len=1024) :: fname

call qg_fields_check(self)

nlocs = size(lats)
levs = self%geom%nz

ii = 0
do jvar=1,vars%nvars()
  fname = vars%variable(jvar)
  select case (trim(fname))
  case ('x')
    if (.not.allocated(self%x)) call abor1_ftn('qg_fields_getvalsad: x not allocated')
    do jloc=1,nlocs
      call qg_interp_bilinear_ad(self%geom,lons(jloc),lats(jloc),vals(ii+1:ii+levs),self%x)
      ii = ii + levs
    enddo
  case ('q')
    if (.not.allocated(self%q)) call abor1_ftn('qg_fields_getvalsad: q not allocated')
    do jloc=1,nlocs
      call qg_interp_bilinear_ad(self%geom,lons(jloc),lats(jloc),vals(ii+1:ii+levs),self%q)
      ii = ii + levs
    enddo
  case ('u')
    if (.not.allocated(self%u)) call abor1_ftn('qg_fields_getvalsad: u not allocated')
    do jloc=1,nlocs
      call qg_interp_bilinear_ad(self%geom,lons(jloc),lats(jloc),vals(ii+1:ii+levs),self%u)
      ii = ii + levs
    enddo
  case ('v')
    if (.not.allocated(self%v)) call abor1_ftn('qg_fields_getvalsad: v not allocated')
    do jloc=1,nlocs
      call qg_interp_bilinear_ad(self%geom,lons(jloc),lats(jloc),vals(ii+1:ii+levs),self%v)
      ii = ii + levs
    enddo
  case ('z')
    ! do nothing
    ii = ii + nlocs * levs
  case default
    call abor1_ftn('qg_fields_getvalsad: wrong input variable')
  endselect
enddo
if (size(vals) /= ii) call abor1_ftn('qg_fields_getvalsad: error size')

end subroutine qg_fields_getvalsad
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
integer :: ii
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

! Initialization
ii = 0

! Get values
if (allocated(fld%x)) then
  vals(ii+1:ii+fld%geom%nz) = fld%x(iter%ilon,iter%ilat,:)
  ii = ii+fld%geom%nz
endif
if (allocated(fld%q)) then
  vals(ii+1:ii+fld%geom%nz) = fld%q(iter%ilon,iter%ilat,:)
  ii = ii+fld%geom%nz
endif
if (allocated(fld%u)) then
  vals(ii+1:ii+fld%geom%nz) = fld%u(iter%ilon,iter%ilat,:)
  ii = ii+fld%geom%nz
endif
if (allocated(fld%v)) then
  vals(ii+1:ii+fld%geom%nz) = fld%v(iter%ilon,iter%ilat,:)
  ii = ii+fld%geom%nz
endif

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
integer :: ii
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

! Initialization
ii = 0

! Get values
if (allocated(fld%x)) then
  fld%x(iter%ilon,iter%ilat,:) = vals(ii+1:ii+fld%geom%nz)
  ii = ii+fld%geom%nz
endif
if (allocated(fld%q)) then
  fld%q(iter%ilon,iter%ilat,:) = vals(ii+1:ii+fld%geom%nz)
  ii = ii+fld%geom%nz
endif
if (allocated(fld%u)) then
  fld%u(iter%ilon,iter%ilat,:) = vals(ii+1:ii+fld%geom%nz)
  ii = ii+fld%geom%nz
endif
if (allocated(fld%v)) then
  fld%v(iter%ilon,iter%ilat,:) = vals(ii+1:ii+fld%geom%nz)
  ii = ii+fld%geom%nz
endif

end subroutine qg_fields_setpoint
! ------------------------------------------------------------------------------
!> Serialize fields
subroutine qg_fields_serialize(fld,vsize,vect_fld)

implicit none

! Passed variables
type(qg_fields),intent(in) :: fld                !< Fields
integer,intent(in) :: vsize                      !< Size
real(kind_real),intent(inout) :: vect_fld(vsize) !< Vector

! Local variables
integer :: ix,iy,iz,ind

! Initialize
ind = 0

! Copy
do iz=1,fld%geom%nz
  do iy=1,fld%geom%ny
    do ix=1,fld%geom%nx
      if (allocated(fld%x)) then
        ind = ind + 1
        vect_fld(ind) = fld%x(ix,iy,iz)
      endif
      if (allocated(fld%q)) then
        ind = ind + 1
        vect_fld(ind) = fld%q(ix,iy,iz)
      endif
      if (allocated(fld%u)) then
        ind = ind + 1
        vect_fld(ind) = fld%u(ix,iy,iz)
      endif
      if (allocated(fld%v)) then
        ind = ind + 1
        vect_fld(ind) = fld%v(ix,iy,iz)
      endif
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
do iz=1,self%geom%nz
  do iy=1,self%geom%ny
    do ix=1,self%geom%nx
      if (allocated(self%x)) then
        self%x(ix,iy,iz) = vect_fld(index)
        index = index+1
      endif
      if (allocated(self%q)) then
        self%q(ix,iy,iz) = vect_fld(index)
        index = index+1
      endif
      if (allocated(self%u)) then
        self%u(ix,iy,iz) = vect_fld(index)
        index = index+1
      endif
      if (allocated(self%v)) then
        self%v(ix,iy,iz) = vect_fld(index)
        index = index+1
      endif
    enddo
  enddo
enddo

! Boundaries
if (self%lbc) then
  do iz=1,self%geom%nz
    self%x_north(iz) = vect_fld(index)
    index = index + 1
    self%x_south(iz) = vect_fld(index)
    index = index + 1
    do ix=1,self%geom%nx
      self%q_north(ix,iz) = vect_fld(index)
      index = index + 1
      self%q_south(ix,iz) = vect_fld(index)
      index = index + 1
    enddo
  enddo
endif

index = index - 1

end subroutine qg_fields_deserialize
! ------------------------------------------------------------------------------
!> Complete missing fields consistently
subroutine qg_fields_complete(self,var)

implicit none

! Passed variables
type(qg_fields),intent(inout) :: self !< Fields
character(len=1),intent(in) :: var    !< Reference variable ('x' or 'q')

! Local variables
real(kind_real) :: x(self%geom%nx,self%geom%ny,self%geom%nz)
real(kind_real) :: q(self%geom%nx,self%geom%ny,self%geom%nz)
real(kind_real) :: u(self%geom%nx,self%geom%ny,self%geom%nz)
real(kind_real) :: v(self%geom%nx,self%geom%ny,self%geom%nz)

select case (var)
case ('x')
  if (allocated(self%q)) then
    if (self%lbc) then
      call convert_x_to_q(self%geom,self%x,self%x_north,self%x_south,self%q)
    else
      call convert_x_to_q_tl(self%geom,self%x,self%q)
    endif
  endif
  if (allocated(self%u)) then
    if (self%lbc) then
      call convert_x_to_u(self%geom,self%x,self%x_north,self%x_south,self%u)
    else
      call convert_x_to_u_tl(self%geom,self%x,self%u)
    endif
  endif
  if (allocated(self%v)) then
    if (self%lbc) then
      call convert_x_to_v(self%geom,self%x,self%v)
    else
      call convert_x_to_v_tl(self%geom,self%x,self%v)
    endif
  endif
case ('q')
  if (allocated(self%x).or.allocated(self%u).or.allocated(self%v)) then
    if (self%lbc) then
      call convert_q_to_x(self%geom,self%q,self%x_north,self%x_south,x)
    else
      call convert_q_to_x_tl(self%geom,self%q,x)
    endif
    if (allocated(self%x)) self%x = x
    if (allocated(self%u)) then
      if (self%lbc) then
        call convert_x_to_u(self%geom,self%x,self%x_north,self%x_south,self%u)
      else
        call convert_x_to_u_tl(self%geom,self%x,self%u)
      endif
    endif
    if (allocated(self%v)) then
      if (self%lbc) then
        call convert_x_to_v(self%geom,self%x,self%v)
      else
        call convert_x_to_v_tl(self%geom,self%x,self%v)
      endif
    endif
  endif
case default
  call abor1_ftn('qg_fields_complete: wrong variable')
end select

end subroutine qg_fields_complete
! ------------------------------------------------------------------------------
!> Check fields
subroutine qg_fields_check(self)

implicit none

! Passed variables
type(qg_fields),intent(in) :: self !< Fields

! Local variables
logical :: bad
character(len=1024) :: record

! Check 3d field
bad = .not.(allocated(self%x).or.allocated(self%q).or.allocated(self%u).or.allocated(self%v))
if (allocated(self%x)) then
  bad = bad.or.(size(self%x,1)/=self%geom%nx)
  bad = bad.or.(size(self%x,2)/=self%geom%ny)
  bad = bad.or.(size(self%x,3)/=self%geom%nz)
endif
if (allocated(self%q)) then
  bad = bad.or.(size(self%q,1)/=self%geom%nx)
  bad = bad.or.(size(self%q,2)/=self%geom%ny)
  bad = bad.or.(size(self%q,3)/=self%geom%nz)
endif
if (allocated(self%u)) then
  bad = bad.or.(size(self%u,1)/=self%geom%nx)
  bad = bad.or.(size(self%u,2)/=self%geom%ny)
  bad = bad.or.(size(self%u,3)/=self%geom%nz)
endif
if (allocated(self%v)) then
  bad = bad.or.(size(self%v,1)/=self%geom%nx)
  bad = bad.or.(size(self%v,2)/=self%geom%ny)
  bad = bad.or.(size(self%v,3)/=self%geom%nz)
endif

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
  if (allocated(self%x)) then
    write(record,*) '  shape(x) = ',shape(self%x)
    call fckit_log%info(record)
  endif
  if (allocated(self%q)) then
    write(record,*) '  shape(q) = ',shape(self%q)
    call fckit_log%info(record)
  endif
  if (allocated(self%u)) then
    write(record,*) '  shape(u) = ',shape(self%u)
    call fckit_log%info(record)
  endif
  if (allocated(self%v)) then
    write(record,*) '  shape(v) = ',shape(self%v)
    call fckit_log%info(record)
  endif
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
end module qg_fields_mod
