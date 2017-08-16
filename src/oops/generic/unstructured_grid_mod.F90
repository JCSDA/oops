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
use tools_dirac, only: setup_dirac_points

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

subroutine setup_dirac_c(key, c_conf) bind(c, name='setup_dirac_f90')
implicit none
integer(c_int), intent(inout) :: key
type(c_ptr), intent(in) :: c_conf
type(unstructured_grid), pointer :: self
call unstructured_grid_registry%get(key,self)
call setup_dirac(self, c_conf)
end subroutine setup_dirac_c

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

subroutine setup_dirac(self, c_conf)
implicit none
type(unstructured_grid), intent(inout) :: self
type(c_ptr), intent(in) :: c_conf
integer :: ndir,idir,levdir,ivardir
integer,allocatable :: ic0dir(:) 
real(kind_real),allocatable :: lon(:),lat(:),londir(:),latdir(:)
integer :: nc0a,ic0a
logical,allocatable :: lmask(:)
type(column_element), pointer :: current
character(len=3) :: idirchar

! Get nc0a
nc0a = 0
current => self%head
do while (associated(current))
   nc0a = nc0a+1
   current => current%next
end do

! Allocation
allocate(lon(nc0a))
allocate(lat(nc0a))
allocate(lmask(nc0a))

! Get lon/lat/mask
ic0a = 0
current => self%head
do while (associated(current))
   ic0a = ic0a+1
   lon(ic0a) = current%column%lon
   lat(ic0a) = current%column%lat
   lmask(ic0a) = .true. ! TODO: add mask in column data
   current => current%next
end do

! Get Dirac properties from JSON
ndir = config_get_int(c_conf,"ndir")
allocate(londir(ndir))
allocate(latdir(ndir))
do idir=1,ndir
   write(idirchar,'(i3)') idir
   londir(idir) = config_get_real(c_conf,"londir("//trim(adjustl(idirchar))//")")
   latdir(idir) = config_get_real(c_conf,"latdir("//trim(adjustl(idirchar))//")")
end do
levdir = config_get_int(c_conf,"levdir")
ivardir = config_get_int(c_conf,"ivardir")

! Setup dirac field
allocate(ic0dir(ndir))
call setup_dirac_points(nc0a,lon,lat,lmask,ndir,londir,latdir,ic0dir)

! Set all columns to zero
current => self%head
do while (associated(current))
   current%column%cols = 0.0
   current => current%next
end do

! Set Dirac
ic0a = 0
current => self%head
do while (associated(current))
   ic0a = ic0a+1
   if (any((latdir==current%column%lat).and.(londir==current%column%lon))) & 
 & current%column%cols((ivardir-1)*self%head%column%nlevs+levdir) = 1.0 
   current => current%next
end do

end subroutine setup_dirac

!-------------------------------------------------------------------------------

end module unstructured_grid_mod
