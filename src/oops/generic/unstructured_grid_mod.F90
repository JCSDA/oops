! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!>  Fortran module for handling generic unstructured grid

module unstructured_grid_mod

use iso_c_binding
use config_mod
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
  integer :: ncols
  integer :: nlevs 
  real(kind=kind_real), allocatable :: levs(:)  !> Definition of vertical coordinate
  type(column_element), pointer :: head         !> Linked list containing the columns
  type(column_element), pointer :: last         !> Last element of linked list
end type unstructured_grid

! ------------------------------------------------------------------------------

#define LISTED_TYPE unstructured_grid

!> Linked list interface - defines registry_t type
#include "util/linkedList_i.f"

!> Global registry
type(registry_t) :: unstructured_grid_registry

!-------------------------------------------------------------------------------
contains
!-------------------------------------------------------------------------------

!> Linked list implementation
#include "util/linkedList_c.f"

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

subroutine get_nlevs_c(key, nlevs) bind(c, name='get_nlevs_f90')
implicit none
integer(c_int), intent(inout) :: key
integer,intent(out) :: nlevs
type(unstructured_grid), pointer :: self
call unstructured_grid_registry%get(key,self)
nlevs = self%nlevs
end subroutine get_nlevs_c

! ------------------------------------------------------------------------------

subroutine get_ncols_c(key, ncols) bind(c, name='get_ncols_f90')
implicit none
integer(c_int), intent(inout) :: key
integer,intent(out) :: ncols
type(unstructured_grid), pointer :: self
call unstructured_grid_registry%get(key,self)
ncols = self%ncols
end subroutine get_ncols_c

!-------------------------------------------------------------------------------

subroutine get_lats_c(key, ncols, lats) bind(c, name='get_lats_f90')
implicit none
integer(c_int), intent(inout) :: key
integer(c_int), intent(in) :: ncols
real(kind=kind_real),intent(out) :: lats(ncols)
integer :: icols
type(column_element), pointer :: current
type(unstructured_grid), pointer :: self

! Get self
call unstructured_grid_registry%get(key,self)

! Get lats
icols = 0
current => self%head
do while (associated(current))
   icols = icols+1
   lats(icols) = current%column%lat
   current => current%next
end do

end subroutine get_lats_c

!-------------------------------------------------------------------------------

subroutine get_lons_c(key, ncols, lons) bind(c, name='get_lons_f90')
implicit none
integer(c_int), intent(inout) :: key
integer(c_int), intent(in) :: ncols
real(kind=kind_real),intent(out) :: lons(ncols)
integer :: icols
type(column_element), pointer :: current
type(unstructured_grid), pointer :: self

! Get self
call unstructured_grid_registry%get(key,self)

! Get lons
icols = 0
current => self%head
do while (associated(current))
   icols = icols+1
   lons(icols) = current%column%lon
   current => current%next
end do

end subroutine get_lons_c

!-------------------------------------------------------------------------------

subroutine get_areas_c(key, ncols, areas) bind(c, name='get_areas_f90')
implicit none
integer(c_int), intent(inout) :: key
integer(c_int), intent(in) :: ncols
real(kind=kind_real),intent(out) :: areas(ncols)
integer :: icols
type(column_element), pointer :: current
type(unstructured_grid), pointer :: self

! Get self
call unstructured_grid_registry%get(key,self)

! Get area
icols = 0
current => self%head
do while (associated(current))
   icols = icols+1
   areas(icols) = current%column%area
   current => current%next
end do

end subroutine get_areas_c

!-------------------------------------------------------------------------------

subroutine get_levs_c(key, nlevs, levs) bind(c, name='get_levs_f90')
implicit none
integer(c_int), intent(inout) :: key
integer(c_int), intent(in) :: nlevs
real(kind=kind_real),intent(out) :: levs(nlevs)
type(unstructured_grid), pointer :: self

! Get self
call unstructured_grid_registry%get(key,self)

! Get levs
levs = self%levs

end subroutine get_levs_c

!-------------------------------------------------------------------------------

subroutine get_mask3d_c(key, ncols, ilev, mask3d) bind(c, name='get_mask3d_f90')
implicit none
integer(c_int), intent(inout) :: key
integer,intent(in) :: ncols
integer,intent(in) :: ilev
integer,intent(out) :: mask3d(ncols)
integer :: icols
type(column_element), pointer :: current
type(unstructured_grid), pointer :: self

! Get self
call unstructured_grid_registry%get(key,self)

! Get mask
icols = 0
current => self%head
do while (associated(current))
   icols = icols+1
   mask3d(icols) = current%column%mask3d(ilev+1) ! +1 for arrays offset from C++ to Fortran
   current => current%next
end do

end subroutine get_mask3d_c

!-------------------------------------------------------------------------------

subroutine get_mask2d_c(key, ncols, mask2d) bind(c, name='get_mask2d_f90')
implicit none
integer(c_int), intent(inout) :: key
integer,intent(in) :: ncols
integer,intent(out) :: mask2d(ncols)
integer :: icols
type(column_element), pointer :: current
type(unstructured_grid), pointer :: self

! Get self
call unstructured_grid_registry%get(key,self)

! Get mask
icols = 0
current => self%head
do while (associated(current))
   icols = icols+1
   mask2d(icols) = current%column%mask2d
   current => current%next
end do

end subroutine get_mask2d_c

!-------------------------------------------------------------------------------

subroutine get_glbind_c(key, ncols, glbind) bind(c, name='get_glbind_f90')
implicit none
integer(c_int), intent(inout) :: key
integer,intent(in) :: ncols
integer,intent(out) :: glbind(ncols)
integer :: icols
type(column_element), pointer :: current
type(unstructured_grid), pointer :: self

! Get self
call unstructured_grid_registry%get(key,self)

! Get mask
icols = 0
current => self%head
do while (associated(current))
   icols = icols+1
   glbind(icols) = current%column%glbind
   current => current%next
end do

end subroutine get_glbind_c

! ------------------------------------------------------------------------------

subroutine create_unstructured_grid(self, klevs, plevs)
implicit none
type(unstructured_grid), intent(inout) :: self
integer, intent(in) :: klevs
real(kind=kind_real), intent(in) :: plevs(klevs)

self%ncols = 0
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

subroutine add_column(self, plat, plon, parea, klevs, kvar3d, kvar2d, kmask3d, kmask2d, kglbind)
implicit none
type(unstructured_grid), intent(inout) :: self
real(kind=kind_real), intent(in) :: plat
real(kind=kind_real), intent(in) :: plon
real(kind=kind_real), intent(in) :: parea
integer, intent(in) :: klevs
integer, intent(in) :: kvar3d
integer, intent(in) :: kvar2d
integer, intent(in) :: kmask3d(klevs)
integer, intent(in) :: kmask2d
integer, optional, intent(in) :: kglbind
integer :: glbind

! Update pointer
if (associated(self%last)) then
  allocate(self%last%next)
  self%last => self%last%next
else
  allocate(self%head)
  self%last => self%head
endif

! Global index
glbind = -1
if (present(kglbind)) glbind = kglbind

! Create column
call create_column_data(self%last%column, plat, plon, parea, klevs, kvar3d, kvar2d, kmask3d, kmask2d, glbind)
self%ncols = self%ncols+1

end subroutine add_column

!-------------------------------------------------------------------------------

end module unstructured_grid_mod
