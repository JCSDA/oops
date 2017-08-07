! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!>  Fortran module for handling generic unstructured grid

module unstructured_grid_mod

use iso_c_binding
use kinds
use column_data_mod

implicit none
private
public unstructured_grid, create_unstructured_grid, delete_unstructured_grid, &
     & unstructured_grid_registry, add_column, column_element

! ------------------------------------------------------------------------------

!>  Derived types containing the data

type :: column_element
  type(column_data) :: column
  type(column_element), pointer :: next => null()
end type column_element

type unstructured_grid
  integer :: nlevs                              !> Number of levels
  real(kind=kind_real), allocatable :: levs(:)  !> Definition of vertical coordinate
  type(column_element), pointer :: head         !> Linked list containing the columns
  type(column_element), pointer :: last         !> Last element of linked list
end type unstructured_grid

! ------------------------------------------------------------------------------

#define LISTED_TYPE unstructured_grid

!> Linked list interface - defines registry_t type
#include "linkedList_i.f"

!> Global registry
type(registry_t) :: unstructured_grid_registry

!-------------------------------------------------------------------------------
contains
!-------------------------------------------------------------------------------

!> Linked list implementation
#include "linkedList_c.f"

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

subroutine create_unstructured_grid(self, klevs, plevs)
implicit none
type(unstructured_grid), intent(inout) :: self
integer, intent(in) :: klevs
real(kind=kind_real), intent(in) :: plevs(klevs)

self%nlevs = klevs
allocate(self%levs(self%nlevs))
self%levs(:) = plevs(:)
self%head => null()
self%last => null()

end subroutine create_unstructured_grid

!-------------------------------------------------------------------------------

subroutine delete_unstructured_grid(self)
implicit none
type(unstructured_grid), intent(inout) :: self
type(column_element), pointer :: current, prev

do while (associated(self%head))
  call delete_column_data(self%head%column)
  prev => self%head
  self%head => self%head%next
  deallocate(prev)
enddo
self%head => null()
self%last => null()

deallocate(self%levs)

end subroutine delete_unstructured_grid

!-------------------------------------------------------------------------------

subroutine add_column(self, plat, plon, klevs, kvars, ksurf)
implicit none
type(unstructured_grid), intent(inout) :: self
real(kind=kind_real), intent(in) :: plat
real(kind=kind_real), intent(in) :: plon
integer, intent(in) :: klevs
integer, intent(in) :: kvars
integer, intent(in) :: ksurf

if (associated(self%last)) then
  allocate(self%last%next)
  self%last => self%last%next
else
  allocate(self%head)
  self%last => self%head
endif
call create_column_data(self%last%column, plat, plon, klevs, kvars, ksurf)

end subroutine add_column

!-------------------------------------------------------------------------------

end module unstructured_grid_mod
