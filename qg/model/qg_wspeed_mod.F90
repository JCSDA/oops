! (C) Copyright 2009-2016 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

!> Fortran module to handle wind speed observations for the QG model
module qg_wspeed_mod

use kinds
use iso_c_binding
use missing_values_mod
use qg_gom_mod
use qg_interp_mod
use qg_obsdb_mod
use qg_obsvec_mod
use oops_variables_mod

implicit none

private
public :: qg_wspeed_equiv,qg_wspeed_equiv_tl,qg_wspeed_equiv_ad, &
        & wspeed_traj,qg_wspeed_alloctraj,qg_wspeed_settraj,qg_wspeed_registry

! ------------------------------------------------------------------------------
type :: wspeed_traj
  integer :: nobs
  real(kind_real), allocatable :: u(:)
  real(kind_real), allocatable :: v(:)
  real(kind_real), allocatable :: zsave(:)
end type

#define LISTED_TYPE wspeed_traj

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry
type(registry_t) :: qg_wspeed_registry
! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
! Public
! ------------------------------------------------------------------------------
!> Linked list implementation
#include "oops/util/linkedList_c.f"
! ------------------------------------------------------------------------------
!> Get equivalent for wind speed
subroutine qg_wspeed_equiv(obsdb,gom,hofx,bias)

implicit none

! Passed variables
type(qg_obsdb),intent(in) :: obsdb
type(qg_gom),intent(in) :: gom        !< GOM
type(qg_obsvec),intent(inout) :: hofx !< Observation vector
real(kind_real),intent(in) :: bias    !< Bias

! Local variables
integer :: iobs
real(kind_real) :: valu, valv
type(qg_obsvec) :: zobs

call qg_obsdb_get(obsdb, 'WSpeed', 'Location', zobs)

! Check bias
if (abs(bias)>epsilon(bias)) call abor1_ftn ('qg_wspeed_equiv: bias not implemented')

! Loop over observations
do iobs=1,gom%nobs
  call qg_vert_interp(gom%levs,gom%z(:,iobs),zobs%values(3,iobs),gom%u(:,iobs),valu)
  call qg_vert_interp(gom%levs,gom%z(:,iobs),zobs%values(3,iobs),gom%v(:,iobs),valv)
  hofx%values(1,iobs) = sqrt(valu*valu+valv*valv)
enddo

end subroutine qg_wspeed_equiv
! ------------------------------------------------------------------------------
!> Get equivalent for wind speed - tangent linear
subroutine qg_wspeed_equiv_tl(obsdb,gom,hofx,traj,bias)

implicit none

! Passed variables
type(qg_obsdb),intent(in) :: obsdb
type(qg_gom),intent(in) :: gom        !< GOM
type(qg_obsvec),intent(inout) :: hofx !< Observation vector
type(wspeed_traj),intent(in) :: traj  !< GOM trajectory
real(kind_real),intent(in) :: bias    !< Bias

! Local variables
integer :: iobs
real(kind_real) :: zu,zv,zt,valu,valv
type(qg_obsvec) :: zobs

call qg_obsdb_get(obsdb, 'WSpeed', 'Location', zobs)

! Loop over observations
do iobs=1,gom%nobs
  zu = traj%u(iobs)
  zv = traj%v(iobs)
  zt = sqrt(zu**2+zv**2)
  if (zt>epsilon(zt)) then
    call qg_vert_interp(gom%levs,gom%z(:,iobs),zobs%values(3,iobs),gom%u(:,iobs),valu)
    call qg_vert_interp(gom%levs,gom%z(:,iobs),zobs%values(3,iobs),gom%v(:,iobs),valv)
    hofx%values(1,iobs) = (zu*valu+zv*valv)/zt
  else
    hofx%values(1,iobs) = 0.0
  endif
enddo

end subroutine qg_wspeed_equiv_tl
! ------------------------------------------------------------------------------
!> Get equivalent for wind speed - adjoint
subroutine qg_wspeed_equiv_ad(obsdb,gom,hofx,traj,bias)

implicit none

! Passed variables
type(qg_obsdb),intent(in) :: obsdb
type(qg_gom),intent(inout) :: gom     !< GOM
type(qg_obsvec),intent(in) :: hofx    !< Observation vector
type(wspeed_traj),intent(in) :: traj  !< GOM trajectory
real(kind_real),intent(inout) :: bias !< Bias

! Local variables
integer :: iobs
real(kind_real) :: zu,zv,zt,valu,valv
type(qg_obsvec) :: zobs

call qg_obsdb_get(obsdb, 'WSpeed', 'Location', zobs)
gom%u(:,:) = 0.0_kind_real
gom%v(:,:) = 0.0_kind_real

! Loop over observations
do iobs=1,gom%nobs
  zu = traj%u(iobs)
  zv = traj%v(iobs)
  zt = sqrt(zu**2+zv**2)
  if (zt>epsilon(zt)) then
    valu = zu*hofx%values(1,iobs)/zt
    valv = zv*hofx%values(1,iobs)/zt
    call qg_vert_interp_ad(gom%levs,traj%zsave(:),zobs%values(3,iobs),gom%u(:,iobs),valu)
    call qg_vert_interp_ad(gom%levs,traj%zsave(:),zobs%values(3,iobs),gom%v(:,iobs),valv)
  else
    gom%u(:,iobs) = 0.0
    gom%v(:,iobs) = 0.0
  endif
enddo

end subroutine qg_wspeed_equiv_ad
! ------------------------------------------------------------------------------
subroutine qg_wspeed_alloctraj(traj, nobs)
type(wspeed_traj), intent(inout) :: traj
integer, intent(in) :: nobs
traj%nobs = nobs
allocate(traj%u(nobs))
allocate(traj%v(nobs))
end subroutine qg_wspeed_alloctraj
! ------------------------------------------------------------------------------
!> Set wind speed trajectory
subroutine qg_wspeed_settraj(obsdb,gom,traj)
implicit none
type(qg_obsdb),intent(in) :: obsdb
type(qg_gom),intent(in) :: gom     !< GOM
type(wspeed_traj), intent(inout) :: traj

! Local variables
integer :: iobs
real(kind_real) :: zz
type(qg_obsvec) :: zobs


call qg_obsdb_get(obsdb, 'WSpeed', 'Location', zobs)
if (.not.allocated(traj%zsave)) allocate(traj%zsave(gom%levs))
! Loop over observations
do iobs=1,gom%nobs
  call qg_vert_interp(gom%levs,gom%z(:,iobs),zobs%values(3,iobs),gom%u(:,iobs),traj%u(iobs))
  call qg_vert_interp(gom%levs,gom%z(:,iobs),zobs%values(3,iobs),gom%v(:,iobs),traj%v(iobs))
end do
if (gom%nobs > 0) then
  traj%zsave(:)=gom%z(:,1)
else
  traj%zsave(:)=missing_value(zz)
endif

end subroutine qg_wspeed_settraj
! ------------------------------------------------------------------------------
end module qg_wspeed_mod
