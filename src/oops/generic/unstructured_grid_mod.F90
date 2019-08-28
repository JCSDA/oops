! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

module unstructured_grid_mod

use iso_c_binding
use kinds

implicit none

private
public :: unstructured_grid
public :: unstructured_grid_registry
public :: create_ug,delete_ug,allocate_unstructured_grid_coord,allocate_unstructured_grid_field
! ------------------------------------------------------------------------------
type grid_type
  integer :: igrid                                 !> Index of the grid
  integer :: nmga                                  !> Number of gridpoints (on a given MPI task)
  integer :: nl0                                   !> Number of levels
  integer :: nv                                    !> Number of variables
  integer :: nts                                   !> Number of timeslots
  real(kind=kind_real),allocatable :: lon(:)       !> Longitude (in degrees: -180 to 180)
  real(kind=kind_real),allocatable :: lat(:)       !> Latitude (in degrees: -90 to 90)
  real(kind=kind_real),allocatable :: area(:)      !> Area (in m^2)
  real(kind=kind_real),allocatable :: vunit(:,:)   !> Vertical unit
  logical,allocatable :: lmask(:,:)                !> Mask
  real(kind=kind_real),allocatable :: fld(:,:,:,:) !> Data
end type grid_type

type unstructured_grid
  integer :: colocated                             !> Colocation flag
  integer :: nts                                   !> Number of timeslots
  integer :: ngrid                                 !> Number of different grids
  type(grid_type),allocatable :: grid(:)           !> Grid instance
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
integer,intent(in) :: colocated               !< Colocation flag
integer,intent(in) :: nts                     !< Number of timeslots

! Set parameters
if ((colocated==0).or.(colocated==1)) then
   self%colocated = colocated
else
   call abor1_ftn('create_ug: wrong colocation flag')
endif
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

end subroutine delete_ug
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
type(unstructured_grid),intent(inout) :: self

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
