! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

module unstructured_grid_mod

use atlas_module
use datetime_mod
use iso_c_binding
use kinds
use oops_variables_mod

implicit none

private
public :: unstructured_grid
public :: unstructured_grid_registry
public :: create_ug,delete_ug,ug_set_atlas_lonlat,ug_fill_atlas_fieldset,ug_set_atlas,ug_to_atlas,ug_from_atlas
public :: allocate_unstructured_grid_coord,allocate_unstructured_grid_field
! ------------------------------------------------------------------------------
type grid_type
  integer :: igrid                                       !> Index of the grid
  integer :: nmga                                        !> Number of gridpoints (on a given MPI task)
  integer :: nl0                                         !> Number of levels
  integer :: nv                                          !> Number of variables
  integer :: nts                                         !> Number of timeslots
  real(kind=kind_real),allocatable :: lon(:)             !> Longitude (in degrees: -180 to 180)
  real(kind=kind_real),allocatable :: lat(:)             !> Latitude (in degrees: -90 to 90)
  real(kind=kind_real),allocatable :: area(:)            !> Area (in m^2)
  real(kind=kind_real),allocatable :: vunit(:,:)         !> Vertical unit
  logical,allocatable :: lmask(:,:)                      !> Mask
  real(kind=kind_real),allocatable :: fld(:,:,:,:)       !> Data
end type grid_type

type unstructured_grid
  integer :: colocated                                   !> Colocation flag
  integer :: nts                                         !> Number of timeslots
  integer :: ngrid                                       !> Number of different grids
  type(grid_type),allocatable :: grid(:)                 !> Grid instance
  type(atlas_functionspace_pointcloud) :: afunctionspace !< ATLAS function space
  type(atlas_fieldset) :: afieldset                      !< ATLAS fieldset
end type unstructured_grid

#define LISTED_TYPE unstructured_grid

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry
type(registry_t) :: unstructured_grid_registry
!-------------------------------------------------------------------------------
contains
!-------------------------------------------------------------------------------
!> Linked list implementation
#include "oops/util/linkedList_c.f"
!-------------------------------------------------------------------------------
!> Create unstructured grid
subroutine create_ug(self,colocated,nts)

implicit none

! Passed variables
type(unstructured_grid),intent(inout) :: self !< Unstructured grid
integer,intent(in) :: colocated               !< Colocated flag
integer,intent(in) :: nts                     !< Number of timeslots

! Set colocated flag
self%colocated = colocated

! Set number of timeslots
self%nts = nts

end subroutine create_ug
!-------------------------------------------------------------------------------
!> Delete unstructured grid
subroutine delete_ug(self)

implicit none

! Passed variables
type(unstructured_grid),intent(inout) :: self !< Unstructured grid

! Local variables
integer :: igrid

! Release memory
if (allocated(self%grid)) then
  do igrid=1,self%ngrid
    if (allocated(self%grid(igrid)%lon)) deallocate(self%grid(igrid)%lon)
    if (allocated(self%grid(igrid)%lat)) deallocate(self%grid(igrid)%lat)
    if (allocated(self%grid(igrid)%area)) deallocate(self%grid(igrid)%area)
    if (allocated(self%grid(igrid)%vunit)) deallocate(self%grid(igrid)%vunit)
    if (allocated(self%grid(igrid)%lmask)) deallocate(self%grid(igrid)%lmask)
    if (allocated(self%grid(igrid)%fld)) deallocate(self%grid(igrid)%fld)
  enddo
  deallocate(self%grid)
endif
call self%afunctionspace%final()
call self%afieldset%final()

end subroutine delete_ug
! ------------------------------------------------------------------------------
!> Set ATLAS grid lon/lat in fieldset
subroutine ug_set_atlas_lonlat(self,afieldset)

! Passed variables
type(unstructured_grid),intent(inout) :: self   !< Unstructured grid
type(atlas_fieldset),intent(inout) :: afieldset !< ATLAS fieldset

! Local variables
integer :: imga
real(kind_real),parameter :: pi = 4.d0*atan(1.0)
real(kind_real),parameter :: rad2deg = 180.0/pi
real(kind_real),pointer :: real_ptr(:,:)
type(atlas_field) :: afield

! Create lon/lat field
afield = atlas_field(name="lonlat",kind=atlas_real(kind_real),shape=(/2,self%grid(1)%nmga/))
call afield%data(real_ptr)
do imga=1,self%grid(1)%nmga
   real_ptr(1,imga) = self%grid(1)%lon(imga)
   real_ptr(2,imga) = self%grid(1)%lat(imga)
end do
call afieldset%add(afield)

end subroutine ug_set_atlas_lonlat
! ------------------------------------------------------------------------------
!> Fill ATLAS fieldset from unstructured grid
subroutine ug_fill_atlas_fieldset(self,afieldset)

! Passed variables
type(unstructured_grid),intent(inout) :: self   !< Unstructured grid
type(atlas_fieldset),intent(inout) :: afieldset !< ATLAS fieldset

! Local variables
integer :: imga,il0
integer,pointer :: int_ptr_2(:,:)
real(kind_real),pointer :: real_ptr_1(:),real_ptr_2(:,:)
type(atlas_field) :: afield

! Add area
afield = self%afunctionspace%create_field(name='area',kind=atlas_real(kind_real),levels=0)
call afield%data(real_ptr_1)
real_ptr_1(1:self%grid(1)%nmga) = self%grid(1)%area
call afieldset%add(afield)
call afield%final()

! Add vertical unit
afield = self%afunctionspace%create_field(name='vunit',kind=atlas_real(kind_real),levels=self%grid(1)%nl0)
call afield%data(real_ptr_2)
real_ptr_2(1:self%grid(1)%nl0,1:self%grid(1)%nmga) = transpose(self%grid(1)%vunit)
call afieldset%add(afield)
call afield%final()

! Add geometry mask
afield = self%afunctionspace%create_field(name='gmask',kind=atlas_integer(kind_int),levels=self%grid(1)%nl0)
call afield%data(int_ptr_2)
do il0=1,self%grid(1)%nl0
  do imga=1,self%grid(1)%nmga
    if (self%grid(1)%lmask(imga,il0)) then
      int_ptr_2(il0,imga) = 1
    else
      int_ptr_2(il0,imga) = 0
    endif
  enddo
enddo
call afieldset%add(afield)
call afield%final()

end subroutine ug_fill_atlas_fieldset
! ------------------------------------------------------------------------------
!> Set ATLAS fieldset
subroutine ug_set_atlas(self,afieldset)

implicit none

! Passed variables
type(unstructured_grid),intent(in) :: self      !< Unstructured grid
type(atlas_fieldset),intent(inout) :: afieldset !< ATLAS fieldset

! Local variables
integer :: iv,its,igrid
character(len=1024) :: fieldname
type(atlas_field) :: afield

! Copy fields
do igrid=1,self%ngrid
  do its=1,self%nts
    do iv=1,self%grid(igrid)%nv
      ! Get or create field
      write(fieldname,'(a,i2.2,a,i2.2)') 'var_',iv,'_',its
      if (afieldset%has_field(trim(fieldname))) then
        ! Get field
        afield = afieldset%field(trim(fieldname))
      else
        ! Create field
        afield = self%afunctionspace%create_field(name=trim(fieldname),kind=atlas_real(kind_real),levels=self%grid(igrid)%nl0)

        ! Add field
        call afieldset%add(afield)
      endif

      ! Release pointer
      call afield%final()
    enddo
  enddo
enddo

end subroutine ug_set_atlas
! ------------------------------------------------------------------------------
!> Unstructured grid to ATLAS fieldset
subroutine ug_to_atlas(self,afieldset)

implicit none

! Passed variables
type(unstructured_grid),intent(in) :: self      !< Unstructured grid
type(atlas_fieldset),intent(inout) :: afieldset !< ATLAS fieldset

! Local variables
integer :: iv,its,igrid
real(kind_real),pointer :: real_ptr_2(:,:)
character(len=1024) :: fieldname
type(atlas_field) :: afield

! Copy fields
do igrid=1,self%ngrid
  do its=1,self%nts
    do iv=1,self%grid(igrid)%nv
      ! Get or create field
      write(fieldname,'(a,i2.2,a,i2.2)') 'var_',iv,'_',its
      if (afieldset%has_field(trim(fieldname))) then
        ! Get field
        afield = afieldset%field(trim(fieldname))
      else
        ! Create field
        afield = self%afunctionspace%create_field(name=trim(fieldname),kind=atlas_real(kind_real),levels=self%grid(igrid)%nl0)

        ! Add field
        call afieldset%add(afield)
      endif

      ! Copy data
      call afield%data(real_ptr_2)
      real_ptr_2(1:self%grid(igrid)%nl0,1:self%grid(igrid)%nmga) = transpose(self%grid(igrid)%fld(:,:,iv,its))

      ! Release pointer
      call afield%final()
    enddo
  enddo
enddo

end subroutine ug_to_atlas
! ------------------------------------------------------------------------------
!> Convert ATLAS fieldset to unstructured grid
subroutine ug_from_atlas(self,afieldset)

implicit none

! Passed variables
type(unstructured_grid),intent(inout) :: self   !< Unstructured grid
type(atlas_fieldset),intent(inout) :: afieldset !< ATLAS fieldset

! Local variables
integer :: iv,its,igrid
real(kind_real),pointer :: real_ptr_2(:,:)
character(len=1024) :: fieldname
type(atlas_field) :: afield

! Allocation
call allocate_unstructured_grid_field(self)

! Copy fields
do igrid=1,self%ngrid
  do its=1,self%nts
    do iv=1,self%grid(igrid)%nv
      ! Get field
      write(fieldname,'(a,i2.2,a,i2.2)') 'var_',iv,'_',its
      afield = afieldset%field(trim(fieldname))

      ! Copy data
      call afield%data(real_ptr_2)
      self%grid(igrid)%fld(:,:,iv,its) = transpose(real_ptr_2(1:self%grid(igrid)%nl0,1:self%grid(igrid)%nmga))

      ! Release pointer
      call afield%final()
    enddo
  enddo
enddo

end subroutine ug_from_atlas
! ------------------------------------------------------------------------------
!> Allocate unstructured grid coordinates
subroutine allocate_unstructured_grid_coord(self)

implicit none

! Passed variables
type(unstructured_grid),intent(inout) :: self !< Unstructured grid

! Local variables
integer :: igrid

! Allocation
do igrid=1,self%ngrid
   if (.not.allocated(self%grid(igrid)%lon)) allocate(self%grid(igrid)%lon(self%grid(igrid)%nmga))
   if (.not.allocated(self%grid(igrid)%lat)) allocate(self%grid(igrid)%lat(self%grid(igrid)%nmga))
   if (.not.allocated(self%grid(igrid)%area)) allocate(self%grid(igrid)%area(self%grid(igrid)%nmga))
   if (.not.allocated(self%grid(igrid)%vunit)) allocate(self%grid(igrid)%vunit(self%grid(igrid)%nmga,self%grid(igrid)%nl0))
   if (.not.allocated(self%grid(igrid)%lmask)) allocate(self%grid(igrid)%lmask(self%grid(igrid)%nmga,self%grid(igrid)%nl0))
enddo

end subroutine allocate_unstructured_grid_coord
! ------------------------------------------------------------------------------
!> Allocate unstructured grid fields
subroutine allocate_unstructured_grid_field(self)

implicit none

! Passed variables
type(unstructured_grid),intent(inout) :: self !< Unstructured grid

! Local variables
integer :: igrid

! Allocation
do igrid=1,self%ngrid
   if (.not.allocated(self%grid(igrid)%fld)) allocate(self%grid(igrid)%fld(self%grid(igrid)%nmga,self%grid(igrid)%nl0,&
 & self%grid(igrid)%nv,self%grid(igrid)%nts))
enddo

end subroutine allocate_unstructured_grid_field
!-------------------------------------------------------------------------------
end module unstructured_grid_mod
