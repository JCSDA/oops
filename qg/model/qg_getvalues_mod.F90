! (C) Copyright 2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module qg_getvalues_mod

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
public :: qg_getvalues
public :: qg_getvalues_registry
public :: qg_getvalues_create,qg_getvalues_delete, &
        & qg_getvalues_interp,qg_getvalues_interp_tl,qg_getvalues_interp_ad

!> type for GetValues and GetValuesTLAD
type :: qg_getvalues
  type(qg_locs)         :: locs   !< all locations in the window
end type qg_getvalues

#define LISTED_TYPE qg_getvalues

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry
type(registry_t) :: qg_getvalues_registry
! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
! Public
! ------------------------------------------------------------------------------
!> Linked list implementation
#include "oops/util/linkedList_c.f"
! ------------------------------------------------------------------------------
!> Create GetValues from geometry and locations (only saves locs currently)
subroutine qg_getvalues_create(self,geom,locs)

implicit none

! Passed variables
type(qg_getvalues),intent(inout) :: self   !< GetValues
type(qg_geom),target,intent(in) :: geom    !< Geometry
type(qg_locs),intent(in) :: locs           !< Locations

! Copy locations
call qg_locs_copy(self%locs, locs)

end subroutine qg_getvalues_create
! ------------------------------------------------------------------------------
!> Delete fields
subroutine qg_getvalues_delete(self)

implicit none

! Passed variables
type(qg_getvalues),intent(inout) :: self !< GetValues

! Release memory
call qg_locs_delete(self%locs)

end subroutine qg_getvalues_delete
! ------------------------------------------------------------------------------
!> Interpolation from fields
subroutine qg_getvalues_interp(self,fld,t1,t2,gom)

implicit none

! Passed variables
type(qg_getvalues),intent(in) :: self  !< GetValues
type(qg_fields),intent(in) :: fld      !< Fields
type(datetime),intent(in) :: t1, t2    !< times
type(qg_gom),intent(inout) :: gom      !< Interpolated values

! Local variables
integer :: jloc
real(kind_real),allocatable :: x(:,:,:),q(:,:,:),u(:,:,:),v(:,:,:)

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
do jloc=1,self%locs%nlocs
  ! Check if current obs is in this time frame (t1,t2]
  if (t1 < self%locs%times(jloc) .and. self%locs%times(jloc) <= t2) then
    ! Interpolate variables
    if (gom%ix /= 0) call qg_interp_trilinear(fld%geom,self%locs%lon(jloc),self%locs%lat(jloc), &
    &                                     self%locs%z(jloc),x,gom%values(gom%ix,jloc))
    if (gom%iq /= 0) call qg_interp_trilinear(fld%geom,self%locs%lon(jloc),self%locs%lat(jloc), &
    &                                     self%locs%z(jloc),q,gom%values(gom%iq,jloc))
    if (gom%iu /= 0) call qg_interp_trilinear(fld%geom,self%locs%lon(jloc),self%locs%lat(jloc), &
    &                                     self%locs%z(jloc),u,gom%values(gom%iu,jloc))
    if (gom%iv /= 0) call qg_interp_trilinear(fld%geom,self%locs%lon(jloc),self%locs%lat(jloc), &
    &                                     self%locs%z(jloc),v,gom%values(gom%iv,jloc))
  endif
enddo
!$omp end parallel do

end subroutine qg_getvalues_interp
! ------------------------------------------------------------------------------
!> Interpolation from fields - tangent linear
subroutine qg_getvalues_interp_tl(self,fld,t1,t2,gom)

implicit none

! Passed variables
type(qg_getvalues),intent(in) :: self  !< GetValues
type(qg_fields),intent(in) :: fld      !< Fields
type(datetime),intent(in) :: t1, t2    !< times
type(qg_gom),intent(inout) :: gom      !< Interpolated values

! Local variables
integer :: jloc
real(kind_real),allocatable :: x(:,:,:),q(:,:,:),u(:,:,:),v(:,:,:)

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
do jloc=1,self%locs%nlocs
  ! Check if current obs is in this time frame (t1,t2]
  if (t1 < self%locs%times(jloc) .and. self%locs%times(jloc) <= t2) then
    ! Interpolate variables
    if (gom%ix /= 0) call qg_interp_trilinear(fld%geom,self%locs%lon(jloc),self%locs%lat(jloc), &
    &                                     self%locs%z(jloc),x,gom%values(gom%ix,jloc))
    if (gom%iq /= 0) call qg_interp_trilinear(fld%geom,self%locs%lon(jloc),self%locs%lat(jloc), &
    &                                     self%locs%z(jloc),q,gom%values(gom%iq,jloc))
    if (gom%iu /= 0) call qg_interp_trilinear(fld%geom,self%locs%lon(jloc),self%locs%lat(jloc), &
    &                                     self%locs%z(jloc),u,gom%values(gom%iu,jloc))
    if (gom%iv /= 0) call qg_interp_trilinear(fld%geom,self%locs%lon(jloc),self%locs%lat(jloc), &
    &                                     self%locs%z(jloc),v,gom%values(gom%iv,jloc))
  endif
enddo
!$omp end parallel do

end subroutine qg_getvalues_interp_tl
! ------------------------------------------------------------------------------
!> Interpolation from fields - adjoint
subroutine qg_getvalues_interp_ad(self,fld,t1,t2,gom)

implicit none

! Passed variables
type(qg_getvalues),intent(in) :: self  !< GetValues
type(qg_fields),intent(inout) :: fld   !< Fields
type(datetime),intent(in) :: t1, t2    !< times
type(qg_gom),intent(in) :: gom         !< Interpolated values

! Local variables
integer :: jloc
real(kind_real),allocatable :: x(:,:,:),q(:,:,:),u(:,:,:),v(:,:,:)

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

do jloc=self%locs%nlocs,1,-1
  ! Check if current obs is in this time frame (t1,t2]
  if (t1 < self%locs%times(jloc) .and. self%locs%times(jloc) <= t2) then
    ! Interpolate variables
    if (gom%ix /= 0) call qg_interp_trilinear_ad(fld%geom,self%locs%lon(jloc),self%locs%lat(jloc), &
    &                                        self%locs%z(jloc),gom%values(gom%ix,jloc),x)
    if (gom%iq /= 0) call qg_interp_trilinear_ad(fld%geom,self%locs%lon(jloc),self%locs%lat(jloc), &
    &                                        self%locs%z(jloc),gom%values(gom%iq,jloc),q)
    if (gom%iu /= 0) call qg_interp_trilinear_ad(fld%geom,self%locs%lon(jloc),self%locs%lat(jloc), &
    &                                        self%locs%z(jloc),gom%values(gom%iu,jloc),u)
    if (gom%iv /= 0) call qg_interp_trilinear_ad(fld%geom,self%locs%lon(jloc),self%locs%lat(jloc), &
    &                                        self%locs%z(jloc),gom%values(gom%iv,jloc),v)
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

end subroutine qg_getvalues_interp_ad
! ------------------------------------------------------------------------------
end module qg_getvalues_mod
