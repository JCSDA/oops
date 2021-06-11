! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

module qg_stream_mod

use kinds
use qg_gom_mod
use qg_obsvec_mod

implicit none

private
public :: qg_stream_equiv,qg_stream_equiv_ad
! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
! Public
! ------------------------------------------------------------------------------
!> Get equivalent for streamfunction (TL calls this subroutine too)
subroutine qg_stream_equiv(gom,hofx,bias)

implicit none

! Passed variables
type(qg_gom),intent(in) :: gom        !< GOM
type(qg_obsvec),intent(inout) :: hofx !< Observation vector
real(kind_real),intent(in) :: bias    !< Bias

! Local variables
integer :: iobs

! Loop over observations
do iobs=1,gom%nobs
  hofx%values(1,iobs) = gom%x(iobs)+bias
enddo

end subroutine qg_stream_equiv
! ------------------------------------------------------------------------------
!> Get equivalent for streamfunction - adjoint
subroutine qg_stream_equiv_ad(gom,hofx,bias)

implicit none

! Passed variables
type(qg_gom),intent(inout) :: gom     !< GOM
type(qg_obsvec),intent(in) :: hofx    !< Observation vector
real(kind_real),intent(inout) :: bias !< Bias

! Local variables
integer :: iobs

! Loop over observations
do iobs=1,gom%nobs
  gom%x(iobs) = hofx%values(1,iobs)
  bias = bias+hofx%values(1,iobs)
enddo

end subroutine qg_stream_equiv_ad
! ------------------------------------------------------------------------------
end module qg_stream_mod
