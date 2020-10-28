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
use qg_obsvec_mod

implicit none

private
public :: qg_wind_equiv,qg_wind_equiv_ad
! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
! Public
! ------------------------------------------------------------------------------
!> Get equivalent for wind (TL calls this subroutine too)
subroutine qg_wind_equiv(gom,hofx,bias)

implicit none

! Passed variables
type(qg_gom),intent(in) :: gom        !< GOM
type(qg_obsvec),intent(inout) :: hofx !< Observation vector
real(kind_real),intent(in) :: bias(2) !< Bias

! Local variables
integer :: iobs

! Loop over observations
do iobs=1,gom%nobs
  hofx%values(1:2,iobs) = gom%values(1:2,iobs)+bias
enddo

end subroutine qg_wind_equiv
! ------------------------------------------------------------------------------
!> Get equivalent for wind - adjoint
subroutine qg_wind_equiv_ad(gom,hofx,bias)

implicit none

! Passed variables
type(qg_gom),intent(inout) :: gom        !< GOM
type(qg_obsvec),intent(in) :: hofx       !< Observation vector
real(kind_real),intent(inout) :: bias(2) !< Bias

! Local variables
integer :: iobs

! Loop over observations
do iobs=1,gom%nobs
  gom%values(1:2,iobs) = hofx%values(1:2,iobs)
  bias(1) = bias(1)+hofx%values(1,iobs)
  bias(2) = bias(2)+hofx%values(2,iobs)
enddo

end subroutine qg_wind_equiv_ad
! ------------------------------------------------------------------------------
end module qg_wind_mod
