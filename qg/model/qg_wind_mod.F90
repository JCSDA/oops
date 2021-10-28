! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

module qg_wind_mod

use kinds
use iso_c_binding
use qg_gom_mod
use qg_interp_mod
use qg_obsdb_mod
use qg_obsvec_mod

implicit none

private
public :: qg_wind_equiv,qg_wind_equiv_ad

real(kind_real), allocatable :: zsave(:)

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
! Public
! ------------------------------------------------------------------------------
!> Get equivalent for wind (TL calls this subroutine too)
subroutine qg_wind_equiv(obsdb,gom,hofx,bias)

implicit none

! Passed variables
type(qg_obsdb),intent(in) :: obsdb
type(qg_gom),intent(in) :: gom        !< GOM
type(qg_obsvec),intent(inout) :: hofx !< Observation vector
real(kind_real),intent(in) :: bias(2) !< Bias

! Local variables
integer :: iobs
real(kind_real) :: valu, valv
type(qg_obsvec) :: zobs

call qg_obsdb_get(obsdb, 'Wind', 'Location', zobs)

if (.not.allocated(zsave)) allocate(zsave(gom%levs))

! Loop over observations
do iobs=1,gom%nobs
  zsave(:) = gom%z(:,iobs)
  call qg_vert_interp(gom%levs,gom%z(:,iobs),zobs%values(3,iobs),gom%u(:,iobs),valu)
  call qg_vert_interp(gom%levs,gom%z(:,iobs),zobs%values(3,iobs),gom%v(:,iobs),valv)
  hofx%values(1,iobs) = valu+bias(1)
  hofx%values(2,iobs) = valv+bias(2)
enddo

end subroutine qg_wind_equiv
! ------------------------------------------------------------------------------
!> Get equivalent for wind - adjoint
subroutine qg_wind_equiv_ad(obsdb,gom,hofx,bias)

implicit none

! Passed variables
type(qg_obsdb),intent(in) :: obsdb
type(qg_gom),intent(inout) :: gom        !< GOM
type(qg_obsvec),intent(in) :: hofx       !< Observation vector
real(kind_real),intent(inout) :: bias(2) !< Bias

! Local variables
integer :: iobs
real(kind_real) :: valu, valv
type(qg_obsvec) :: zobs

call qg_obsdb_get(obsdb, 'Wind', 'Location', zobs)
gom%u(:,:) = 0.0_kind_real
gom%v(:,:) = 0.0_kind_real

! Loop over observations
do iobs=1,gom%nobs
  bias(1) = bias(1)+hofx%values(1,iobs)
  bias(2) = bias(2)+hofx%values(2,iobs)
  valu = hofx%values(1,iobs)
  valv = hofx%values(2,iobs)
  call qg_vert_interp_ad(gom%levs,zsave(:),zobs%values(3,iobs),gom%u(:,iobs),valu)
  call qg_vert_interp_ad(gom%levs,zsave(:),zobs%values(3,iobs),gom%v(:,iobs),valv)
enddo

end subroutine qg_wind_equiv_ad
! ------------------------------------------------------------------------------
end module qg_wind_mod
