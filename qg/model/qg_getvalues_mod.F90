! (C) Copyright 2020 UCAR
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
use qg_convert_q_to_x_mod
use qg_convert_x_to_q_mod
use qg_convert_x_to_uv_mod
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
real(kind_real),allocatable :: x(:,:,:),q(:,:,:),u(:,:,:),v(:,:,:)
real(kind_real), pointer :: lonlat(:,:), z(:)
type(atlas_field) :: lonlat_field, z_field

! get locations
lonlat_field = locs%lonlat()
call lonlat_field%data(lonlat)

z_field = locs%altitude()
call z_field%data(z)

! Check field
call qg_fields_check(fld)

! Allocation
allocate(x(fld%geom%nx,fld%geom%ny,fld%geom%nz))
allocate(q(fld%geom%nx,fld%geom%ny,fld%geom%nz))
allocate(u(fld%geom%nx,fld%geom%ny,fld%geom%nz))
allocate(v(fld%geom%nx,fld%geom%ny,fld%geom%nz))

! Get variables
if (fld%lq) then
  q = fld%gfld3d
else
  x = fld%gfld3d
endif
if (gom%ix /= 0.or.gom%iu /= 0.or.gom%iv /= 0) then
  if (fld%lq) call convert_q_to_x(fld%geom,q,fld%x_north,fld%x_south,x)
endif
if (gom%iq /= 0) then
  if (.not.fld%lq) call convert_x_to_q(fld%geom,x,fld%x_north,fld%x_south,q)
endif
if (gom%iu /= 0.or.gom%iv /= 0) call convert_x_to_uv(fld%geom,x,fld%x_north,fld%x_south,u,v)

!$omp parallel do schedule(static) private(jloc)
do jloc=1,locs%nlocs()
  ! Check if current obs is in this time frame (t1,t2]
  if (t1 < locs%times(jloc) .and. locs%times(jloc) <= t2) then
    ! Interpolate variables
    if (gom%ix /= 0) call qg_interp_trilinear(fld%geom,lonlat(1,jloc),lonlat(2,jloc), &
    &                                     z(jloc),x,gom%values(gom%ix,jloc))
    if (gom%iq /= 0) call qg_interp_trilinear(fld%geom,lonlat(1,jloc),lonlat(2,jloc), &
    &                                     z(jloc),q,gom%values(gom%iq,jloc))
    if (gom%iu /= 0) call qg_interp_trilinear(fld%geom,lonlat(1,jloc),lonlat(2,jloc), &
    &                                     z(jloc),u,gom%values(gom%iu,jloc))
    if (gom%iv /= 0) call qg_interp_trilinear(fld%geom,lonlat(1,jloc),lonlat(2,jloc), &
    &                                     z(jloc),v,gom%values(gom%iv,jloc))
  endif
enddo
!$omp end parallel do

call lonlat_field%final()
call z_field%final()

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
real(kind_real),allocatable :: x(:,:,:),q(:,:,:),u(:,:,:),v(:,:,:)
real(kind_real), pointer :: lonlat(:,:), z(:)
type(atlas_field) :: lonlat_field, z_field

! get locations
lonlat_field = locs%lonlat()
call lonlat_field%data(lonlat)

z_field = locs%altitude()
call z_field%data(z)

! Check field
call qg_fields_check(fld)

! Allocation
allocate(x(fld%geom%nx,fld%geom%ny,fld%geom%nz))
allocate(q(fld%geom%nx,fld%geom%ny,fld%geom%nz))
allocate(u(fld%geom%nx,fld%geom%ny,fld%geom%nz))
allocate(v(fld%geom%nx,fld%geom%ny,fld%geom%nz))

! Get variables
if (fld%lq) then
  q = fld%gfld3d
else
  x = fld%gfld3d
endif
if (gom%ix /= 0.or.gom%iu /= 0.or.gom%iv /= 0) then
  if (fld%lq) call convert_q_to_x_tl(fld%geom,q,x)
endif
if (gom%iq /= 0) then
  if (.not.fld%lq) call convert_x_to_q_tl(fld%geom,x,q)
endif
if (gom%iu /= 0.or.gom%iv /= 0) call convert_x_to_uv_tl(fld%geom,x,u,v)

!$omp parallel do schedule(static) private(jloc)
do jloc=1,locs%nlocs()
  ! Check if current obs is in this time frame (t1,t2]
  if (t1 < locs%times(jloc) .and. locs%times(jloc) <= t2) then
    ! Interpolate variables
    if (gom%ix /= 0) call qg_interp_trilinear(fld%geom,lonlat(1,jloc),lonlat(2,jloc), &
    &                                     z(jloc),x,gom%values(gom%ix,jloc))
    if (gom%iq /= 0) call qg_interp_trilinear(fld%geom,lonlat(1,jloc),lonlat(2,jloc), &
    &                                     z(jloc),q,gom%values(gom%iq,jloc))
    if (gom%iu /= 0) call qg_interp_trilinear(fld%geom,lonlat(1,jloc),lonlat(2,jloc), &
    &                                     z(jloc),u,gom%values(gom%iu,jloc))
    if (gom%iv /= 0) call qg_interp_trilinear(fld%geom,lonlat(1,jloc),lonlat(2,jloc), &
    &                                     z(jloc),v,gom%values(gom%iv,jloc))
  endif
enddo
!$omp end parallel do

call lonlat_field%final()
call z_field%final()

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
real(kind_real),allocatable :: x(:,:,:),q(:,:,:),u(:,:,:),v(:,:,:)
real(kind_real), pointer :: lonlat(:,:), z(:)
type(atlas_field) :: lonlat_field, z_field

! get locations
lonlat_field = locs%lonlat()
call lonlat_field%data(lonlat)

z_field = locs%altitude()
call z_field%data(z)

! Check field
call qg_fields_check(fld)

! Allocation
allocate(x(fld%geom%nx,fld%geom%ny,fld%geom%nz))
allocate(q(fld%geom%nx,fld%geom%ny,fld%geom%nz))
allocate(u(fld%geom%nx,fld%geom%ny,fld%geom%nz))
allocate(v(fld%geom%nx,fld%geom%ny,fld%geom%nz))

! Initialization
x = 0.0
q = 0.0
u = 0.0
v = 0.0

do jloc=locs%nlocs(),1,-1
  ! Check if current obs is in this time frame (t1,t2]
  if (t1 < locs%times(jloc) .and. locs%times(jloc) <= t2) then
    ! Interpolate variables
    if (gom%ix /= 0) call qg_interp_trilinear_ad(fld%geom,lonlat(1,jloc),lonlat(2,jloc), &
    &                                        z(jloc),gom%values(gom%ix,jloc),x)
    if (gom%iq /= 0) call qg_interp_trilinear_ad(fld%geom,lonlat(1,jloc),lonlat(2,jloc), &
    &                                        z(jloc),gom%values(gom%iq,jloc),q)
    if (gom%iu /= 0) call qg_interp_trilinear_ad(fld%geom,lonlat(1,jloc),lonlat(2,jloc), &
    &                                        z(jloc),gom%values(gom%iu,jloc),u)
    if (gom%iv /= 0) call qg_interp_trilinear_ad(fld%geom,lonlat(1,jloc),lonlat(2,jloc), &
    &                                        z(jloc),gom%values(gom%iv,jloc),v)
  endif
enddo

! Get variables
if (gom%iu /= 0.or.gom%iv /= 0) call convert_x_to_uv_ad(fld%geom,u,v,x)
if (gom%iq /= 0) then
  if (.not.fld%lq) call convert_x_to_q_ad(fld%geom,q,x)
endif
if (gom%ix /= 0.or.gom%iu /= 0.or.gom%iv /= 0) then
  if (fld%lq) call convert_q_to_x_ad(fld%geom,x,q)
endif
if (fld%lq) then
  fld%gfld3d = fld%gfld3d+q
else
  fld%gfld3d = fld%gfld3d+x
endif

call lonlat_field%final()
call z_field%final()

end subroutine qg_getvalues_interp_ad
! ------------------------------------------------------------------------------
end module qg_getvalues_mod
