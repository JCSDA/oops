! (C) Copyright 2020-2021 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module qg_getvalues_mod

use atlas_module, only: atlas_field
use datetime_mod
use fckit_log_module,only: fckit_log
use iso_c_binding
use kinds
!$ use omp_lib
use oops_variables_mod
use qg_fields_mod
use qg_geom_mod
use qg_gom_mod
use qg_interp_mod
use qg_locs_mod

implicit none

private
public :: qg_getvalues_interp, qg_getvalues_interp_tl, qg_getvalues_interp_ad
! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
! Public
! ------------------------------------------------------------------------------
!> Interpolation from fields
subroutine qg_getvalues_interp(locs,fld,t1,t2,gom)

implicit none

! Passed variables
type(qg_locs), intent(in) :: locs      !< Locations
type(qg_fields),intent(in) :: fld      !< Fields
type(datetime),intent(in) :: t1, t2    !< times
type(qg_gom),intent(inout) :: gom      !< Interpolated values

! Local variables
integer :: jloc
real(kind_real), pointer :: lonlat(:,:)
type(atlas_field) :: lonlat_field

! get locations
lonlat_field = locs%lonlat()
call lonlat_field%data(lonlat)

! Check field
call qg_fields_check(fld)

!$omp parallel do schedule(static) private(jloc)
do jloc=1,locs%nlocs()
  ! Check if current obs is in this time frame (t1,t2]
  if (t1 < locs%times(jloc) .and. locs%times(jloc) <= t2) then
    ! Interpolate variables
    if (gom%vars%has('x')) call qg_interp_bilinear(fld%geom,lonlat(1,jloc),lonlat(2,jloc), &
                                                 & fld%x,gom%x(:,jloc))
    if (gom%vars%has('q')) call qg_interp_bilinear(fld%geom,lonlat(1,jloc),lonlat(2,jloc), &
                                                 & fld%q,gom%q(:,jloc))
    if (gom%vars%has('u')) call qg_interp_bilinear(fld%geom,lonlat(1,jloc),lonlat(2,jloc), &
                                                 & fld%u,gom%u(:,jloc))
    if (gom%vars%has('v')) call qg_interp_bilinear(fld%geom,lonlat(1,jloc),lonlat(2,jloc), &
                                                 & fld%v,gom%v(:,jloc))
    if (gom%vars%has('z')) gom%z(:,jloc) = fld%geom%z(:)
  endif
enddo
!$omp end parallel do

call lonlat_field%final()

end subroutine qg_getvalues_interp
! ------------------------------------------------------------------------------
!> Interpolation from fields - tangent linear
subroutine qg_getvalues_interp_tl(locs,fld,t1,t2,gom)

implicit none

! Passed variables
type(qg_locs), intent(in) :: locs      !< Locations
type(qg_fields),intent(in) :: fld      !< Fields
type(datetime),intent(in) :: t1, t2    !< times
type(qg_gom),intent(inout) :: gom      !< Interpolated values

! Local variables
integer :: jloc
real(kind_real), pointer :: lonlat(:,:)
type(atlas_field) :: lonlat_field

! get locations
lonlat_field = locs%lonlat()
call lonlat_field%data(lonlat)

! Check field
call qg_fields_check(fld)

!$omp parallel do schedule(static) private(jloc)
do jloc=1,locs%nlocs()
  ! Check if current obs is in this time frame (t1,t2]
  if (t1 < locs%times(jloc) .and. locs%times(jloc) <= t2) then
    ! Interpolate variables
    if (gom%vars%has('x')) call qg_interp_bilinear(fld%geom,lonlat(1,jloc),lonlat(2,jloc), &
                                                 & fld%x,gom%x(:,jloc))
    if (gom%vars%has('q')) call qg_interp_bilinear(fld%geom,lonlat(1,jloc),lonlat(2,jloc), &
                                                 & fld%q,gom%q(:,jloc))
    if (gom%vars%has('u')) call qg_interp_bilinear(fld%geom,lonlat(1,jloc),lonlat(2,jloc), &
                                                 & fld%u,gom%u(:,jloc))
    if (gom%vars%has('v')) call qg_interp_bilinear(fld%geom,lonlat(1,jloc),lonlat(2,jloc), &
                                                 & fld%v,gom%v(:,jloc))
    if (gom%vars%has('z')) gom%z(:,jloc) = fld%geom%z(:)
  endif
enddo
!$omp end parallel do

call lonlat_field%final()

end subroutine qg_getvalues_interp_tl

! ------------------------------------------------------------------------------
!> Interpolation from fields - adjoint
subroutine qg_getvalues_interp_ad(locs,fld,t1,t2,gom)

implicit none

! Passed variables
type(qg_locs), intent(in) :: locs      !< Locations
type(qg_fields),intent(inout) :: fld   !< Fields
type(datetime),intent(in) :: t1, t2    !< times
type(qg_gom),intent(in) :: gom         !< Interpolated values

! Local variables
integer :: jloc
real(kind_real), pointer :: lonlat(:,:)
type(atlas_field) :: lonlat_field
type(qg_fields) :: fld_tmp

! Get locations
lonlat_field = locs%lonlat()
call lonlat_field%data(lonlat)

! Check field
call qg_fields_check(fld)

! Create field with GOM variables
call qg_fields_create(fld_tmp,fld%geom,gom%vars,.false.)

! Initialization
call qg_fields_zero(fld_tmp)

do jloc=locs%nlocs(),1,-1
  ! Check if current obs is in this time frame (t1,t2]
  if (t1 < locs%times(jloc) .and. locs%times(jloc) <= t2) then
    ! Interpolate variables
    if (gom%vars%has('x')) call qg_interp_bilinear_ad(fld%geom,lonlat(1,jloc),lonlat(2,jloc), &
                                                    & gom%x(:,jloc),fld_tmp%x)
    if (gom%vars%has('q')) call qg_interp_bilinear_ad(fld%geom,lonlat(1,jloc),lonlat(2,jloc), &
                                                    & gom%q(:,jloc),fld_tmp%q)
    if (gom%vars%has('u')) call qg_interp_bilinear_ad(fld%geom,lonlat(1,jloc),lonlat(2,jloc), &
                                                    & gom%u(:,jloc),fld_tmp%u)
    if (gom%vars%has('v')) call qg_interp_bilinear_ad(fld%geom,lonlat(1,jloc),lonlat(2,jloc), &
                                                    & gom%v(:,jloc),fld_tmp%v)
  endif
enddo

call qg_fields_self_add(fld,fld_tmp)

! Release memory
call lonlat_field%final()

end subroutine qg_getvalues_interp_ad
! ------------------------------------------------------------------------------
end module qg_getvalues_mod
