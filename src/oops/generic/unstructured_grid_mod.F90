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
  real(kind=kind_real), allocatable :: vunit(:) !> Definition of the vertical unit
  type(column_element), pointer :: head         !> Linked list containing the columns
  type(column_element), pointer :: last         !> Last element of linked list
end type unstructured_grid

! ------------------------------------------------------------------------------

#define LISTED_TYPE unstructured_grid

!> Linked list interface - defines registry_t type
#include "util/linkedList_i.f"
!#include "linkedList.intf.h"

!> Global registry
type(registry_t) :: unstructured_grid_registry

!-------------------------------------------------------------------------------
contains
!-------------------------------------------------------------------------------

!> Linked list implementation
#include "util/linkedList_c.f"
!#include "linkedList.h"

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

subroutine get_vunit_c(key, nlevs, vunit) bind(c, name='get_vunit_f90')
implicit none
integer(c_int), intent(inout) :: key
integer(c_int), intent(in) :: nlevs
real(kind=kind_real),intent(out) :: vunit(nlevs)
type(unstructured_grid), pointer :: self

! Get self
call unstructured_grid_registry%get(key,self)

! Get vunit
vunit = self%vunit

end subroutine get_vunit_c

!-------------------------------------------------------------------------------

subroutine get_mask_c(key, ncols, ilev, mask) bind(c, name='get_mask_f90')
implicit none
integer(c_int), intent(inout) :: key
integer,intent(in) :: ncols
integer,intent(in) :: ilev
integer,intent(out) :: mask(ncols)
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
   mask(icols) = current%column%mask(ilev+1) ! +1 for arrays offset from C++ to Fortran
   current => current%next
end do

end subroutine get_mask_c

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

! Get global index
icols = 0
current => self%head
do while (associated(current))
   icols = icols+1
   glbind(icols) = current%column%glbind
   current => current%next
end do

end subroutine get_glbind_c

!-------------------------------------------------------------------------------

subroutine get_nvar_c(key, nvar) bind(c, name='get_nvar_f90')
implicit none
integer(c_int), intent(inout) :: key
integer,intent(out) :: nvar
class(unstructured_grid), pointer :: self

! Get self
call unstructured_grid_registry%get(key,self)

! Get nvar
nvar = self%head%column%nvar

end subroutine get_nvar_c

!-------------------------------------------------------------------------------

subroutine get_data_c(key, ntot, fld) bind(c, name='get_data_f90')
implicit none
integer(c_int), intent(inout) :: key
integer(c_int), intent(in) :: ntot
real(kind=kind_real),intent(out) :: fld(ntot)
integer :: offset,ivar
type(column_element), pointer :: current
type(unstructured_grid), pointer :: self

! Get self
call unstructured_grid_registry%get(key,self)

! Loop over columns
offset = 0
current => self%head
do while (associated(current))
  fld(offset+1:offset+self%head%column%nvar*self%head%column%nlevs) = pack(current%column%fld,mask=.true.)
  offset = offset+self%head%column%nvar*self%head%column%nlevs
  current => current%next
end do

end subroutine get_data_c

! ------------------------------------------------------------------------------

subroutine create_unstructured_grid(self, nlevs, vunit)
implicit none
type(unstructured_grid), intent(inout) :: self
integer, intent(in) :: nlevs
real(kind=kind_real), intent(in) :: vunit(nlevs)

self%ncols = 0
self%nlevs = nlevs
allocate(self%vunit(self%nlevs))
self%vunit(:) = vunit(:)
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

deallocate(self%vunit)

end subroutine delete_unstructured_grid

!-------------------------------------------------------------------------------

subroutine add_column(self, plat, plon, parea, klevs, kvar, kmask, kglbind)
implicit none
type(unstructured_grid), intent(inout) :: self
real(kind=kind_real), intent(in) :: plat
real(kind=kind_real), intent(in) :: plon
real(kind=kind_real), intent(in) :: parea
integer, intent(in) :: klevs
integer, intent(in) :: kvar
integer, intent(in) :: kmask(klevs)
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
call create_column_data(self%last%column, plat, plon, parea, klevs, kvar, kmask, glbind)
self%ncols = self%ncols+1

end subroutine add_column

!-------------------------------------------------------------------------------

end module unstructured_grid_mod
