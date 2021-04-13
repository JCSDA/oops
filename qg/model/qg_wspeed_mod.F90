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
use qg_gom_mod
use qg_obsvec_mod
use oops_variables_mod

implicit none

private
public :: qg_wspeed_equiv,qg_wspeed_equiv_tl,qg_wspeed_equiv_ad,qg_wspeed_settraj
! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
! Public
! ------------------------------------------------------------------------------
!> Get equivalent for wind speed
subroutine qg_wspeed_equiv(gom,hofx,bias)

implicit none

! Passed variables
type(qg_gom),intent(in) :: gom        !< GOM
type(qg_obsvec),intent(inout) :: hofx !< Observation vector
real(kind_real),intent(in) :: bias    !< Bias

! Local variables
integer :: iobs

! Check bias
if (abs(bias)>epsilon(bias)) call abor1_ftn ('qg_wspeed_equiv: bias not implemented')

! Loop over observations
do iobs=1,gom%nobs
  hofx%values(1,iobs) = sqrt(gom%u(iobs)*gom%u(iobs)+gom%v(iobs)*gom%v(iobs))
enddo

end subroutine qg_wspeed_equiv
! ------------------------------------------------------------------------------
!> Get equivalent for wind speed - tangent linear
subroutine qg_wspeed_equiv_tl(gom,hofx,traj,bias)

implicit none

! Passed variables
type(qg_gom),intent(in) :: gom        !< GOM
type(qg_obsvec),intent(inout) :: hofx !< Observation vector
type(qg_gom),intent(in) :: traj       !< GOM trajectory
real(kind_real),intent(in) :: bias    !< Bias

! Local variables
integer :: iobs
real(kind_real) :: zu,zv,zt

! Loop over observations
do iobs=1,gom%nobs
  zu = traj%u(iobs)
  zv = traj%v(iobs)
  zt = sqrt(zu**2+zv**2)
  if (zt>epsilon(zt)) then
    hofx%values(1,iobs) = (zu*gom%u(iobs)+zv*gom%v(iobs))/zt
  else
    hofx%values(1,iobs) = 0.0
  endif
enddo

end subroutine qg_wspeed_equiv_tl
! ------------------------------------------------------------------------------
!> Get equivalent for wind speed - adjoint
subroutine qg_wspeed_equiv_ad(gom,hofx,traj,bias)

implicit none

! Passed variables
type(qg_gom),intent(inout) :: gom     !< GOM
type(qg_obsvec),intent(in) :: hofx    !< Observation vector
type(qg_gom),intent(in) :: traj       !< GOM trajectory
real(kind_real),intent(inout) :: bias !< Bias

! Local variables
integer :: iobs
real(kind_real) :: zu,zv,zt

! Loop over observations
do iobs=1,gom%nobs
  zu = traj%u(iobs)
  zv = traj%v(iobs)
  zt = sqrt(zu**2+zv**2)
  if (zt>epsilon(zt)) then
    gom%u(iobs) = zu*hofx%values(1,iobs)/zt
    gom%v(iobs) = zv*hofx%values(1,iobs)/zt
  else
    gom%u(iobs) = 0.0
    gom%v(iobs) = 0.0
  endif
enddo

end subroutine qg_wspeed_equiv_ad
! ------------------------------------------------------------------------------
!> Set wind speed trajectory
subroutine qg_wspeed_settraj(gom,traj)

implicit none

! Passed variables
type(qg_gom),intent(in) :: gom     !< GOM
type(qg_gom),intent(inout) :: traj !< GOM trajectory

! Local variables
integer :: iobs

! Loop over observations
do iobs=1,gom%nobs
  traj%u(iobs) = gom%u(iobs)
  traj%v(iobs) = gom%v(iobs)
end do

end subroutine qg_wspeed_settraj
! ------------------------------------------------------------------------------
end module qg_wspeed_mod
