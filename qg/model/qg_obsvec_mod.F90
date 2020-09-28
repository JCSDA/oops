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
use random_mod

implicit none

private
public :: qg_obsvec
public :: qg_obsvec_registry
public :: qg_obsvec_setup,qg_obsvec_clone,qg_obsvec_delete,qg_obsvec_copy,qg_obsvec_zero,qg_obsvec_mul_scal,qg_obsvec_add, &
        & qg_obsvec_sub,qg_obsvec_mul,qg_obsvec_div,qg_obsvec_axpy,qg_obsvec_invert,qg_obsvec_random,qg_obsvec_dotprod, &
        & qg_obsvec_stats,qg_obsvec_nobs,qg_obsvec_copy_local,qg_obsvec_getat
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

  ! Allocation
  allocate(self%values(self%nlev,self%nobs))
endif

! Copy data
self%values = other%values

end subroutine qg_obsvec_copy
! ------------------------------------------------------------------------------
!> Copy a local subset of the observation vector
subroutine qg_obsvec_copy_local(self,other,idx)

implicit none

! Passed variables
type(qg_obsvec),intent(inout) :: self !< Observation vector
type(qg_obsvec),intent(in) :: other   !< Other observation vector
integer,intent(in) :: idx(:)

! local variables
integer :: i

if ((other%nlev/=self%nlev).or.(size(idx)/=self%nobs)) then
  ! Release memory
  deallocate(self%values)

  ! Set sizes
  self%nlev = other%nlev
  self%nobs = size(idx)

  ! Allocation
  allocate(self%values(self%nlev,self%nobs))
endif

! Copy data
do i = 1,self%nobs
  self%values(:,i) = other%values(:,idx(i)+1)
enddo

end subroutine qg_obsvec_copy_local
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
!> Multiply observation vector with a scalar
subroutine qg_obsvec_mul_scal(self,zz)

implicit none

! Passed variables
type(qg_obsvec),intent(inout) :: self !< Observation vector
real(kind_real),intent(in) :: zz      !< Multiplier

! Multiply observation vector with a scalar
self%values = zz*self%values

end subroutine qg_obsvec_mul_scal
! ------------------------------------------------------------------------------
!> Add observation vector
subroutine qg_obsvec_add(self,other)

implicit none

! Passed variables
type(qg_obsvec),intent(inout) :: self !< Observation vector
type(qg_obsvec),intent(in) :: other   !< Other observation vector

! Add observation vector
self%values = self%values+other%values

end subroutine qg_obsvec_add
! ------------------------------------------------------------------------------
!> Subtract observation vector
subroutine qg_obsvec_sub(self,other)

implicit none

! Passed variables
type(qg_obsvec),intent(inout) :: self !< Observation vector
type(qg_obsvec),intent(in) :: other   !< Other observation vector

! Subtract observation vector
self%values = self%values-other%values

end subroutine qg_obsvec_sub
! ------------------------------------------------------------------------------
!> Multiply observation vector
subroutine qg_obsvec_mul(self,other)

implicit none

! Passed variables
type(qg_obsvec),intent(inout) :: self !< Observation vector
type(qg_obsvec),intent(in) :: other   !< Other observation vector

! Multiply observation vector
self%values = self%values*other%values

end subroutine qg_obsvec_mul
! ------------------------------------------------------------------------------
!> Divide observation vector
subroutine qg_obsvec_div(self,other)

implicit none

! Passed variables
type(qg_obsvec),intent(inout) :: self !< Observation vector
type(qg_obsvec),intent(in) :: other   !< Other observation vector

! Divide observation vector
self%values = self%values/other%values

end subroutine qg_obsvec_div
! ------------------------------------------------------------------------------
!> Apply axpy on observation vector
subroutine qg_obsvec_axpy(self,zz,other)

implicit none

! Passed variables
type(qg_obsvec),intent(inout) :: self !< Observation vector
real(kind_real),intent(in) :: zz      !< Multiplier
type(qg_obsvec),intent(in) :: other   !< Other observation vector

! Apply axpy on observation vector
self%values = self%values+zz*other%values

end subroutine qg_obsvec_axpy
! ------------------------------------------------------------------------------
!> Invert observation vector
subroutine qg_obsvec_invert(self)

implicit none

! Passed variables
type(qg_obsvec),intent(inout) :: self !< Observation vector

! Invert observation vector
self%values = 1.0/self%values

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
call qg_obsvec_random_i(c_odb,nval,self%values(1,1))

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
    zz = zz+obsvec1%values(jlev,jobs)*obsvec2%values(jlev,jobs)
  enddo
enddo

end subroutine qg_obsvec_dotprod
! ------------------------------------------------------------------------------
!> Compute observation vector statistics
subroutine qg_obsvec_stats(self,scaling,zmin,zmax,zavg)

implicit none

! Passed variables
type(qg_obsvec),intent(in):: self        !< Observation vector
real(kind_real),intent(inout) :: scaling !< Scaling
real(kind_real),intent(inout) :: zmin    !< Minimum
real(kind_real),intent(inout) :: zmax    !< Maximum
real(kind_real),intent(inout) :: zavg    !< Average

! Local variables
integer :: jlev,jobs
real(kind_real) :: expo

if (self%nobs>0.and.self%nlev>0) then
  ! Compute statistics
  if (.not.allocated(self%values)) call abor1_ftn('qg_obsvec_stats: obs vector not allocated')

  ! Initialization
  zmin = huge(1.0)
  zmax = -huge(1.0)
  zavg = 0.0

  ! Loop over values
  do jobs=1,self%nobs
    do jlev=1,self%nlev
      zmin = min(self%values(jlev,jobs),zmin)
      zmax = max(self%values(jlev,jobs),zmax)
      zavg = zavg+self%values(jlev,jobs)
    enddo
  enddo

  ! Normalization
  zavg = zavg/real(self%nlev*self%nobs,kind_real)
else
  ! Empty observation vector
  zmin = 0.0
  zmax = 0.0
  zavg = 0.0
endif

! Scaling
if (abs(zavg)>0.0) then
  expo = aint(log(abs(zavg))/log(10.0_kind_real))
  scaling = 10.0**expo
else
  scaling = 1.0
endif
zmin = zmin/scaling
zmax = zmax/scaling
zavg = zavg/scaling

end subroutine qg_obsvec_stats
! ------------------------------------------------------------------------------
!> Get observation vector size
subroutine qg_obsvec_nobs(self,kobs)

implicit none

! Passed variables
type(qg_obsvec),intent(in) :: self !< Observation vector
integer,intent(inout) :: kobs      !< Observation vector size

! Get observation vector size
kobs = self%nobs*self%nlev

end subroutine qg_obsvec_nobs

! ------------------------------------------------------------------------------
!> Get value from observation vector at location (iob)
subroutine qg_obsvec_getat(self,iob,val)

implicit none

! Passed variables
type(qg_obsvec),intent(in) :: self !< Observation vector
integer,intent(in) :: iob          !< index into observation vector
real(kind_real), intent(out) :: val!< returned value

integer :: i1, i2

i1 = iob / self%nobs + 1
i2 = iob - self%nobs*(i1-1) + 1
! Retrieve obs. value from vector

if (i1>self%nlev .or. i2>self%nobs) call abor1_ftn ('qg_obsvec_getat: index is out of bounds')

val = self%values(i1,i2)

end subroutine qg_obsvec_getat

! ------------------------------------------------------------------------------
end module qg_obsvec_mod
