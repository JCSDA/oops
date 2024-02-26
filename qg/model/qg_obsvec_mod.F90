! (C) Copyright 2009-2016 ECMWF.
! (C) Copyright 2017-2019 UCAR.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

module qg_obsvec_mod

use iso_c_binding
use kinds
use missing_values_mod
use random_mod

implicit none

private
public :: qg_obsvec
public :: qg_obsvec_registry
public :: qg_obsvec_setup,qg_obsvec_clone,qg_obsvec_delete,qg_obsvec_copy,qg_obsvec_zero, &
        & qg_obsvec_settomissing_ith,qg_obsvec_ones,qg_obsvec_mask,qg_obsvec_mask_with_missing, &
        & qg_obsvec_mul_scal,qg_obsvec_add,qg_obsvec_sub,qg_obsvec_mul,qg_obsvec_div, &
        & qg_obsvec_axpy,qg_obsvec_invert,qg_obsvec_random,qg_obsvec_dotprod,qg_obsvec_stats, &
        & qg_obsvec_size,qg_obsvec_nobs,qg_obsvec_nobs_withmask,qg_obsvec_get_withmask
! ------------------------------------------------------------------------------
interface
  subroutine qg_obsvec_random_i(odb,nn,zz) bind(c,name='qg_obsvec_random_f')
  use iso_c_binding
  implicit none
  type(c_ptr),intent(in) :: odb
  integer(c_int),intent(in) :: nn
  real(c_double),intent(inout) :: zz
  end subroutine qg_obsvec_random_i
end interface
! ------------------------------------------------------------------------------
type qg_obsvec
  integer :: nobs = 0                        !< Number of observations
  integer :: nlev = 0                        !< Number of levels
  real(kind_real),allocatable :: values(:,:) !< Values
  real(kind_real) :: missing                 !< Missing value
end type qg_obsvec

#define LISTED_TYPE qg_obsvec

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry
type(registry_t) :: qg_obsvec_registry
! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
! Public
! ------------------------------------------------------------------------------
!> Linked list implementation
#include "oops/util/linkedList_c.f"
! ------------------------------------------------------------------------------
!> Setup observation vector
subroutine qg_obsvec_setup(self,nlev,nobs)

implicit none

! Passed variables
type(qg_obsvec),intent(inout) :: self !< Observation vector
integer,intent(in) :: nlev            !< Number of levels
integer,intent(in) :: nobs            !< Number of observations

! Set sizes
self%nlev = nlev
self%nobs = nobs

! Release memory
if (allocated(self%values)) deallocate(self%values)

! Allocation
allocate(self%values(self%nlev,self%nobs))

! Initialization
self%values = 0.0_kind_real
self%missing = missing_value(self%missing)

end subroutine qg_obsvec_setup
! ------------------------------------------------------------------------------
!> Clone observation vector
subroutine qg_obsvec_clone(self,other)

implicit none

! Passed variables
type(qg_obsvec),intent(inout) :: self !< Observation vector
type(qg_obsvec),intent(in) :: other   !< Other observation vector

! Set sizes
self%nlev = other%nlev
self%nobs = other%nobs
self%missing = other%missing

! Allocation
allocate(self%values(self%nlev,self%nobs))

end subroutine qg_obsvec_clone
! ------------------------------------------------------------------------------
!> Delete observation vector
subroutine qg_obsvec_delete(self)

implicit none

! Passed variables
type(qg_obsvec),intent(inout) :: self !< Observation vector

! Release memory
deallocate(self%values)

end subroutine qg_obsvec_delete
! ------------------------------------------------------------------------------
!> Copy observation vector
subroutine qg_obsvec_copy(self,other)

implicit none

! Passed variables
type(qg_obsvec),intent(inout) :: self !< Observation vector
type(qg_obsvec),intent(in) :: other   !< Other observation vector

if ((other%nlev/=self%nlev).or.(other%nobs/=self%nobs)) then
  ! Release memory
  deallocate(self%values)

  ! Set sizes
  self%nlev = other%nlev
  self%nobs = other%nobs
  self%missing = other%missing

  ! Allocation
  allocate(self%values(self%nlev,self%nobs))
endif

! Copy data
self%values = other%values

end subroutine qg_obsvec_copy
! ------------------------------------------------------------------------------
!> Set observation vector to zero
subroutine qg_obsvec_zero(self)

implicit none

! Passed variables
type(qg_obsvec),intent(inout) :: self !< Observation vector

! Set observation vector to zero
self%values = 0.0

end subroutine qg_obsvec_zero
! ------------------------------------------------------------------------------
!> Set i-th value of observation vector to missing value
subroutine qg_obsvec_settomissing_ith(self, i)

implicit none

! Passed variables
type(qg_obsvec),intent(inout) :: self !< Observation vector
integer, intent(in) :: i

! Set observation vector to zero
self%values(:,i) = self%missing

end subroutine qg_obsvec_settomissing_ith
! ------------------------------------------------------------------------------
!> Set observation vector to ones
subroutine qg_obsvec_ones(self)

implicit none

! Passed variables
type(qg_obsvec),intent(inout) :: self !< Observation vector

! Set observation vector to ones
self%values = 1.0

end subroutine qg_obsvec_ones
! ------------------------------------------------------------------------------
!> Mask observation vector (set values to missing values where mask == 1)
subroutine qg_obsvec_mask(self,mask)
implicit none

! Passed variables
type(qg_obsvec),intent(inout) :: self !< Observation vector
type(qg_obsvec),intent(in) :: mask    !< mask

if ((self%nobs/=mask%nobs).or.(self%nlev/=mask%nlev)) call abor1_ftn('qg_obsvec_mask: inconsistent sizes')

where(mask%values == 1) self%values = self%missing

end subroutine qg_obsvec_mask
! ------------------------------------------------------------------------------
!> Mask observation vector (set values to missing values where mask == missing value)
subroutine qg_obsvec_mask_with_missing(self,mask)
implicit none

! Passed variables
type(qg_obsvec),intent(inout) :: self !< Observation vector
type(qg_obsvec),intent(in) :: mask    !< mask

if ((self%nobs/=mask%nobs).or.(self%nlev/=mask%nlev)) call abor1_ftn('qg_obsvec_mask: inconsistent sizes')

where(mask%values == mask%missing) self%values = self%missing

end subroutine qg_obsvec_mask_with_missing
! ------------------------------------------------------------------------------
!> Multiply observation vector with a scalar
subroutine qg_obsvec_mul_scal(self,zz)

implicit none

! Passed variables
type(qg_obsvec),intent(inout) :: self !< Observation vector
real(kind_real),intent(in) :: zz      !< Multiplier

! Multiply observation vector with a scalar
where(self%values /= self%missing) self%values = zz*self%values

end subroutine qg_obsvec_mul_scal
! ------------------------------------------------------------------------------
!> Add observation vector
subroutine qg_obsvec_add(self,other)

implicit none

! Passed variables
type(qg_obsvec),intent(inout) :: self !< Observation vector
type(qg_obsvec),intent(in) :: other   !< Other observation vector

if ((self%nobs/=other%nobs).or.(self%nlev/=other%nlev)) call abor1_ftn('qg_obsvec_add: inconsistent sizes')

! Add observation vector
where(self%values /= self%missing .and. other%values /= other%missing)
  self%values = self%values+other%values
elsewhere
  self%values = self%missing
endwhere

end subroutine qg_obsvec_add
! ------------------------------------------------------------------------------
!> Subtract observation vector
subroutine qg_obsvec_sub(self,other)

implicit none

! Passed variables
type(qg_obsvec),intent(inout) :: self !< Observation vector
type(qg_obsvec),intent(in) :: other   !< Other observation vector

if ((self%nobs/=other%nobs).or.(self%nlev/=other%nlev)) call abor1_ftn('qg_obsvec_sub: inconsistent sizes')

! Subtract observation vector
where(self%values /= self%missing .and. other%values /= other%missing)
  self%values = self%values-other%values
elsewhere
  self%values = self%missing
endwhere

end subroutine qg_obsvec_sub
! ------------------------------------------------------------------------------
!> Multiply observation vector
subroutine qg_obsvec_mul(self,other)

implicit none

! Passed variables
type(qg_obsvec),intent(inout) :: self !< Observation vector
type(qg_obsvec),intent(in) :: other   !< Other observation vector

if ((self%nobs/=other%nobs).or.(self%nlev/=other%nlev)) call abor1_ftn('qg_obsvec_mul: inconsistent sizes')

! Multiply observation vector
where(self%values /= self%missing .and. other%values /= other%missing)
  self%values = self%values*other%values
elsewhere
  self%values = self%missing
endwhere

end subroutine qg_obsvec_mul
! ------------------------------------------------------------------------------
!> Divide observation vector
subroutine qg_obsvec_div(self,other)

implicit none

! Passed variables
type(qg_obsvec),intent(inout) :: self !< Observation vector
type(qg_obsvec),intent(in) :: other   !< Other observation vector

if ((self%nobs/=other%nobs).or.(self%nlev/=other%nlev)) call abor1_ftn('qg_obsvec_div: inconsistent sizes')

! Divide observation vector
where(self%values /= self%missing .and. other%values /= other%missing)
  self%values = self%values/other%values
elsewhere
  self%values = self%missing
endwhere

end subroutine qg_obsvec_div
! ------------------------------------------------------------------------------
!> Apply axpy on observation vector
subroutine qg_obsvec_axpy(self,zz,other)

implicit none

! Passed variables
type(qg_obsvec),intent(inout) :: self !< Observation vector
real(kind_real),intent(in) :: zz      !< Multiplier
type(qg_obsvec),intent(in) :: other   !< Other observation vector

if ((self%nobs/=other%nobs).or.(self%nlev/=other%nlev)) call abor1_ftn('qg_obsvec_axpy: inconsistent sizes')

! Apply axpy on observation vector
where(self%values /= self%missing .and. other%values /= other%missing)
  self%values = self%values+zz*other%values
elsewhere
  self%values = self%missing
endwhere

end subroutine qg_obsvec_axpy
! ------------------------------------------------------------------------------
!> Invert observation vector
subroutine qg_obsvec_invert(self)

implicit none

! Passed variables
type(qg_obsvec),intent(inout) :: self !< Observation vector

! Invert observation vector
where(self%values /= self%missing) self%values = 1.0/self%values

end subroutine qg_obsvec_invert
! ------------------------------------------------------------------------------
!> Generate random observation vector
subroutine qg_obsvec_random(c_odb,self)

implicit none

! Passed variables
type(c_ptr),intent(in) :: c_odb       !< Observation data base
type(qg_obsvec),intent(inout) :: self !< Observation vector

! Local variables
integer :: nval

! Compute total size
nval = self%nobs*self%nlev

! Get random values
if (nval.gt.0) call qg_obsvec_random_i(c_odb,nval,self%values(1,1))

end subroutine qg_obsvec_random
! ------------------------------------------------------------------------------
!> Compute dot product between observation vectors
subroutine qg_obsvec_dotprod(obsvec1,obsvec2,zz)

implicit none

! Passed variables
type(qg_obsvec),intent(in) :: obsvec1 !< Observation vector 1
type(qg_obsvec),intent(in) :: obsvec2 !< Observation vector 2
real(kind_real),intent(inout) :: zz   !< Dot product

! Local variables
integer :: jlev,jobs

! Check sizes
if ((obsvec1%nobs/=obsvec2%nobs).or.(obsvec1%nlev/=obsvec2%nlev)) call abor1_ftn('qg_obsvec_dotprod: inconsistent sizes')

! Initalization
zz = 0.0

! Loop over values
do jobs=1,obsvec1%nobs
  do jlev=1,obsvec1%nlev
    if (obsvec1%values(jlev, jobs) /= obsvec1%missing .and. &
        obsvec2%values(jlev, jobs) /= obsvec2%missing)      &
      zz = zz+obsvec1%values(jlev,jobs)*obsvec2%values(jlev,jobs)
  enddo
enddo

end subroutine qg_obsvec_dotprod
! ------------------------------------------------------------------------------
!> Compute observation vector statistics
subroutine qg_obsvec_stats(self,zmin,zmax,zavg)

implicit none

! Passed variables
type(qg_obsvec),intent(in):: self        !< Observation vector
real(kind_real),intent(inout) :: zmin    !< Minimum
real(kind_real),intent(inout) :: zmax    !< Maximum
real(kind_real),intent(inout) :: zavg    !< Average

if (self%nobs*self%nlev>0) then
  ! Compute statistics
  if (.not.allocated(self%values)) call abor1_ftn('qg_obsvec_stats: obs vector not allocated')
  zmin = minval(self%values, mask = (self%values /= self%missing))
  zmax = maxval(self%values, mask = (self%values /= self%missing))
  zavg = sum(self%values, mask = (self%values /= self%missing)) /    &
         count(mask = (self%values /= self%missing))
else
  ! Empty observation vector
  zmin = 0.0
  zmax = 0.0
  zavg = 0.0
endif

end subroutine qg_obsvec_stats

! ------------------------------------------------------------------------------
!> Get observation vector size
subroutine qg_obsvec_size(self,kobs)
implicit none
type(qg_obsvec),intent(in) :: self !< Observation vector
integer,intent(inout) :: kobs      !< Observation vector size
kobs = size(self%values) + 2
end subroutine qg_obsvec_size

! ------------------------------------------------------------------------------
!> Get observation vector size
subroutine qg_obsvec_nobs(self,kobs)

implicit none

! Passed variables
type(qg_obsvec),intent(in) :: self !< Observation vector
integer,intent(inout) :: kobs      !< Observation vector size

! Get observation vector size
kobs = count(mask = (self%values /= self%missing))

end subroutine qg_obsvec_nobs

! ------------------------------------------------------------------------------
!> Get observation vector size (only non-masked observations)
subroutine qg_obsvec_nobs_withmask(self,obsmask,kobs)

implicit none

! Passed variables
type(qg_obsvec),intent(in) :: self    !< Observation vector
type(qg_obsvec),intent(in) :: obsmask !< mask
integer,intent(inout) :: kobs         !< Observation vector size

! Get observation vector size
kobs = count(mask = (self%values /= self%missing) .and.     &
                    (obsmask%values /= obsmask%missing))

end subroutine qg_obsvec_nobs_withmask

! ------------------------------------------------------------------------------
!> Get non-missing values from observation vector into vals array
subroutine qg_obsvec_get_withmask(self,obsmask,vals,nvals)

implicit none

! Passed variables
type(qg_obsvec),intent(in) :: self !< Observation vector
type(qg_obsvec),intent(in) :: obsmask !< mask
integer,intent(in) :: nvals        !< Number of non-missing values
real(kind_real), dimension(nvals), intent(out) :: vals!< returned value

integer :: jobs, jlev, jval

jval = 1
! Loop over values
do jobs=1,self%nobs
  do jlev=1,self%nlev
    if ((self%values(jlev, jobs) /= self%missing) .and.           &
        (obsmask%values(jlev, jobs) /= obsmask%missing)) then
      if (jval > nvals) call abor1_ftn('qg_obsvec_get: inconsistent vector size')
      vals(jval) = self%values(jlev, jobs)
      jval = jval + 1
    endif
  enddo
enddo

end subroutine qg_obsvec_get_withmask

! ------------------------------------------------------------------------------
end module qg_obsvec_mod
