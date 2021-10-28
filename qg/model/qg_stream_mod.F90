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
use qg_interp_mod
use qg_obsdb_mod
use qg_obsvec_mod

implicit none

private
public :: qg_stream_equiv,qg_stream_equiv_ad

real(kind_real), allocatable :: zsave(:)

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
! Public
! ------------------------------------------------------------------------------
!> Get equivalent for streamfunction (TL calls this subroutine too)
subroutine qg_stream_equiv(obsdb,gom,hofx,bias)

implicit none

! Passed variables
type(qg_obsdb),intent(in) :: obsdb
type(qg_gom),intent(in) :: gom        !< GOM
type(qg_obsvec),intent(inout) :: hofx !< Observation vector
real(kind_real),intent(in) :: bias    !< Bias

! Local variables
integer :: iobs
real(kind_real) :: val
type(qg_obsvec) :: zobs

call qg_obsdb_get(obsdb, 'Stream', 'Location', zobs)

if (.not.allocated(zsave)) allocate(zsave(gom%levs))

! Loop over observations
do iobs=1,gom%nobs
  zsave(:) = gom%z(:,iobs)
  call qg_vert_interp(gom%levs,gom%z(:,iobs),zobs%values(3,iobs),gom%x(:,iobs),val)
  hofx%values(1,iobs) = val + bias
enddo

end subroutine qg_stream_equiv
! ------------------------------------------------------------------------------
!> Get equivalent for streamfunction - adjoint
subroutine qg_stream_equiv_ad(obsdb,gom,hofx,bias)

implicit none

! Passed variables
type(qg_obsdb),intent(in) :: obsdb
type(qg_gom),intent(inout) :: gom     !< GOM
type(qg_obsvec),intent(in) :: hofx    !< Observation vector
real(kind_real),intent(inout) :: bias !< Bias

! Local variables
integer :: iobs
real(kind_real) :: val
type(qg_obsvec) :: zobs

call qg_obsdb_get(obsdb, 'Stream', 'Location', zobs)
gom%x(:,:) = 0.0_kind_real

! Loop over observations
do iobs=1,gom%nobs
  bias = bias+hofx%values(1,iobs)
  val = hofx%values(1,iobs)
  call qg_vert_interp_ad(gom%levs,zsave(:)     ,zobs%values(3,iobs),gom%x(:,iobs),val)
enddo

end subroutine qg_stream_equiv_ad
! ------------------------------------------------------------------------------
end module qg_stream_mod
