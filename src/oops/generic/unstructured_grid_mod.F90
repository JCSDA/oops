! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!>  Fortran module for handling generic unstructured grid

module unstructured_grid_mod

use iso_c_binding
use config_mod
use kinds

implicit none
private
public unstructured_grid, allocate_unstructured_grid_coord, allocate_unstructured_grid_field, &
     & delete_unstructured_grid, unstructured_grid_registry

! ------------------------------------------------------------------------------

!>  Derived type containing the data

type grid_type
  integer :: igrid                                 !> Index of the grid
  integer :: nmga                                  !> Number of gridpoints (on a given MPI task)
  integer :: nl0                                   !> Number of levels
  integer :: nv                                    !> Number of variables
  integer :: nts                                   !> Number of timeslots
  real(kind=kind_real),allocatable :: lon(:)       !> Longitude (in degrees: -180 to 180)
  real(kind=kind_real),allocatable :: lat(:)       !> Latitude (in degrees: -90 to 90)
  real(kind=kind_real),allocatable :: area(:)      !> Area (in m^2)
  real(kind=kind_real), allocatable :: vunit(:,:)  !> Vertical unit
  logical,allocatable :: lmask(:,:)                !> Mask
  real(kind=kind_real),allocatable :: fld(:,:,:,:) !> Data
end type grid_type

type unstructured_grid
  integer :: colocated                            !> Colocation flag
  integer :: ngrid                                !> Number of different grids
  type(grid_type),allocatable :: grid(:)          !> Grid instance
end type unstructured_grid

! ------------------------------------------------------------------------------

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

! ------------------------------------------------------------------------------
!  C++ interfaces
! ------------------------------------------------------------------------------

subroutine create_ug_c(key) bind(c, name='create_ug_f90')
implicit none
integer(c_int), intent(inout) :: key

call unstructured_grid_registry%init()
call unstructured_grid_registry%add(key)

end subroutine

! ------------------------------------------------------------------------------

subroutine delete_ug_c(key) bind(c, name='delete_ug_f90')
implicit none
integer(c_int), intent(inout) :: key

type(unstructured_grid), pointer :: self

call unstructured_grid_registry%get(key,self)
call delete_unstructured_grid(self)
call unstructured_grid_registry%remove(key)

end subroutine

! ------------------------------------------------------------------------------

subroutine allocate_unstructured_grid_coord(self)
implicit none
type(unstructured_grid), intent(inout) :: self
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

subroutine allocate_unstructured_grid_field(self)
implicit none
type(unstructured_grid), intent(inout) :: self
integer :: igrid

! Allocation
do igrid=1,self%ngrid
   if (.not.allocated(self%grid(igrid)%fld)) allocate(self%grid(igrid)%fld(self%grid(igrid)%nmga,self%grid(igrid)%nl0, &
 & self%grid(igrid)%nv,self%grid(igrid)%nts))
enddo

end subroutine allocate_unstructured_grid_field

!-------------------------------------------------------------------------------

subroutine delete_unstructured_grid(self)
implicit none
type(unstructured_grid), intent(inout) :: self
integer :: igrid

! Release memory
do igrid=1,self%ngrid
   if (allocated(self%grid(igrid)%lon)) deallocate(self%grid(igrid)%lon)
   if (allocated(self%grid(igrid)%lat)) deallocate(self%grid(igrid)%lat)
   if (allocated(self%grid(igrid)%area)) deallocate(self%grid(igrid)%area)
   if (allocated(self%grid(igrid)%vunit)) deallocate(self%grid(igrid)%vunit)
   if (allocated(self%grid(igrid)%lmask)) deallocate(self%grid(igrid)%lmask)
   if (allocated(self%grid(igrid)%fld)) deallocate(self%grid(igrid)%fld)
enddo

end subroutine delete_unstructured_grid

!-------------------------------------------------------------------------------

end module unstructured_grid_mod
