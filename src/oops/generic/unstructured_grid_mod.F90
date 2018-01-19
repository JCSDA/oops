! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!>  Fortran module for handling generic unstructured grid

module unstructured_grid_mod

use iso_c_binding
use config_mod
use kinds
use tools_display, only: msgerror

implicit none
private
public unstructured_grid, create_unstructured_grid, delete_unstructured_grid, &
     & unstructured_grid_registry

! ------------------------------------------------------------------------------

!>  Derived type containing the data

type unstructured_grid
  integer :: nc0a                                  !> Number of gridpoints (on a given MPI task)
  integer :: nl0                                   !> Number of levels
  integer :: nv                                    !> Number of variables
  integer :: nts                                   !> Number of timeslots
  real(kind=kind_real),allocatable :: lon(:)       !> Longitude
  real(kind=kind_real),allocatable :: lat(:)       !> Latitude
  real(kind=kind_real),allocatable :: area(:)      !> Area
  real(kind=kind_real), allocatable :: vunit(:)    !> Vertical unit
  integer,allocatable :: imask(:,:)                !> Mask
  real(kind=kind_real),allocatable :: fld(:,:,:,:) !> Data
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

subroutine get_size_c(key, ind, isize) bind(c, name='get_size_f90')
implicit none
integer(c_int), intent(inout) :: key
<<<<<<< HEAD
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
=======
integer(c_int), intent(in) :: ind
integer,intent(out) :: isize

class(unstructured_grid), pointer :: self

>>>>>>> feature/nicas_latest
call unstructured_grid_registry%get(key,self)
select case (ind)
case (1)
   ! Number of gridpoints
   isize = self%nc0a
case (2)
   ! Number of levels
   isize = self%nl0
case (3)
   ! Number of variables
   isize = self%nv
case (4)
   ! Number of timeslots
   isize= self%nts
case default
   call msgerror('wrong index in get_size_c')
end select

end subroutine get_size_c

!-------------------------------------------------------------------------------

subroutine get_lon_c(key, nc0a, lon) bind(c, name='get_lon_f90')
implicit none
integer(c_int), intent(inout) :: key
<<<<<<< HEAD
integer(c_int), intent(in) :: ncols
real(kind=kind_real),intent(out) :: lats(ncols)
integer :: icols
type(column_element), pointer :: current
type(unstructured_grid), pointer :: self
=======
integer(c_int), intent(in) :: nc0a
real(kind=kind_real),intent(out) :: lon(nc0a)

class(unstructured_grid), pointer :: self
>>>>>>> feature/nicas_latest

call unstructured_grid_registry%get(key,self)
lon = self%lon

end subroutine get_lon_c

!-------------------------------------------------------------------------------

subroutine get_lat_c(key, nc0a, lat) bind(c, name='get_lat_f90')
implicit none
integer(c_int), intent(inout) :: key
<<<<<<< HEAD
integer(c_int), intent(in) :: ncols
real(kind=kind_real),intent(out) :: lons(ncols)
integer :: icols
type(column_element), pointer :: current
type(unstructured_grid), pointer :: self
=======
integer(c_int), intent(in) :: nc0a
real(kind=kind_real),intent(out) :: lat(nc0a)

class(unstructured_grid), pointer :: self
>>>>>>> feature/nicas_latest

call unstructured_grid_registry%get(key,self)
lat = self%lat

end subroutine get_lat_c

!-------------------------------------------------------------------------------

subroutine get_area_c(key, nc0a, area) bind(c, name='get_area_f90')
implicit none
integer(c_int), intent(inout) :: key
<<<<<<< HEAD
integer(c_int), intent(in) :: ncols
real(kind=kind_real),intent(out) :: areas(ncols)
integer :: icols
type(column_element), pointer :: current
type(unstructured_grid), pointer :: self
=======
integer(c_int), intent(in) :: nc0a
real(kind=kind_real),intent(out) :: area(nc0a)

class(unstructured_grid), pointer :: self
>>>>>>> feature/nicas_latest

call unstructured_grid_registry%get(key,self)
area = self%area

end subroutine get_area_c

!-------------------------------------------------------------------------------

subroutine get_vunit_c(key, nl0, vunit) bind(c, name='get_vunit_f90')
implicit none
integer(c_int), intent(inout) :: key
<<<<<<< HEAD
integer(c_int), intent(in) :: nlevs
real(kind=kind_real),intent(out) :: vunit(nlevs)
type(unstructured_grid), pointer :: self
=======
integer(c_int), intent(in) :: nl0
real(kind=kind_real),intent(out) :: vunit(nl0)

class(unstructured_grid), pointer :: self
>>>>>>> feature/nicas_latest

call unstructured_grid_registry%get(key,self)
vunit = self%vunit

end subroutine get_vunit_c

!-------------------------------------------------------------------------------

<<<<<<< HEAD
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
type(unstructured_grid), pointer :: self
=======
subroutine get_imask_c(key, nc0a, nl0, imask) bind(c, name='get_imask_f90')
implicit none
integer(c_int), intent(inout) :: key
integer,intent(in) :: nc0a
integer,intent(in) :: nl0
integer,intent(out) :: imask(nc0a*nl0)

class(unstructured_grid), pointer :: self
>>>>>>> feature/nicas_latest

call unstructured_grid_registry%get(key,self)
imask = pack(self%imask,mask=.true.)

end subroutine get_imask_c

!-------------------------------------------------------------------------------

subroutine get_data_c(key, ntot, fld) bind(c, name='get_data_f90')
implicit none
integer(c_int), intent(inout) :: key
integer(c_int), intent(in) :: ntot
real(kind=kind_real),intent(out) :: fld(ntot)
<<<<<<< HEAD
integer :: offset,ivar
type(column_element), pointer :: current
type(unstructured_grid), pointer :: self
=======

class(unstructured_grid), pointer :: self
>>>>>>> feature/nicas_latest

call unstructured_grid_registry%get(key,self)
fld = pack(self%fld,mask=.true.)

end subroutine get_data_c

! ------------------------------------------------------------------------------

subroutine create_unstructured_grid(self, nc0a, nl0, nv, nts, lon, lat, area, vunit, imask)
implicit none
<<<<<<< HEAD
type(unstructured_grid), intent(inout) :: self
integer, intent(in) :: nlevs
real(kind=kind_real), intent(in) :: vunit(nlevs)

self%ncols = 0
self%nlevs = nlevs
allocate(self%vunit(self%nlevs))
self%vunit(:) = vunit(:)
self%head => null()
self%last => null()
=======
class(unstructured_grid), intent(inout) :: self
integer, intent(in) :: nc0a
integer, intent(in) :: nl0
integer, intent(in) :: nv
integer, intent(in) :: nts
real(kind=kind_real), intent(in) :: lon(nc0a)
real(kind=kind_real), intent(in) :: lat(nc0a)
real(kind=kind_real), intent(in) :: area(nc0a)
real(kind=kind_real), intent(in) :: vunit(nl0)
integer, intent(in) :: imask(nc0a,nl0)

! Copy sizes
self%nc0a = nc0a
self%nl0 = nl0
self%nv = nv
self%nts = nts
if (nts>1) call msgerror('not ready yet for nts>1')

! Allocation
allocate(self%lon(nc0a))
allocate(self%lat(nc0a))
allocate(self%area(nc0a))
allocate(self%vunit(nl0))
allocate(self%imask(nc0a,nl0))
allocate(self%fld(nc0a,nl0,nv,nts))

! Copy coordinates
self%lon = lon
self%lat = lat
self%area = area
self%vunit = vunit
self%imask = imask
>>>>>>> feature/nicas_latest

end subroutine create_unstructured_grid

!-------------------------------------------------------------------------------

subroutine delete_unstructured_grid(self)
implicit none
<<<<<<< HEAD
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
=======
class(unstructured_grid), intent(inout) :: self
>>>>>>> feature/nicas_latest

! Release memory 
deallocate(self%lon)
deallocate(self%lat)
deallocate(self%area)
deallocate(self%vunit)
deallocate(self%imask)
deallocate(self%fld)

end subroutine delete_unstructured_grid

!-------------------------------------------------------------------------------

<<<<<<< HEAD
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

=======
>>>>>>> feature/nicas_latest
end module unstructured_grid_mod
