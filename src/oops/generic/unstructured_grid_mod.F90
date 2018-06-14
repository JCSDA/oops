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
  real(kind=kind_real),allocatable :: lon(:)       !> Longitude (in degrees: -180 to 180)
  real(kind=kind_real),allocatable :: lat(:)       !> Latitude (in degrees: -90 to 90)
  real(kind=kind_real),allocatable :: area(:)      !> Area (in m^2)
  real(kind=kind_real), allocatable :: vunit(:,:)  !> Vertical unit
  integer,allocatable :: imask(:,:)                !> Mask
  real(kind=kind_real),allocatable :: fld(:,:,:,:) !> Data
end type unstructured_grid

! ------------------------------------------------------------------------------

#define LISTED_TYPE unstructured_grid

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"
!#include "linkedList.intf.h"

!> Global registry
type(registry_t) :: unstructured_grid_registry

!-------------------------------------------------------------------------------
contains
!-------------------------------------------------------------------------------

!> Linked list implementation
#include "oops/util/linkedList_c.f"
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
integer(c_int), intent(in) :: ind
integer,intent(out) :: isize

type(unstructured_grid), pointer :: self

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

subroutine get_lon_c(key, n, lon) bind(c, name='get_lon_f90')
implicit none
integer(c_int), intent(inout) :: key
integer(c_int), intent(in) :: n
real(kind=kind_real),intent(out) :: lon(n)

type(unstructured_grid), pointer :: self

call unstructured_grid_registry%get(key,self)
lon = self%lon

end subroutine get_lon_c

!-------------------------------------------------------------------------------

subroutine get_lat_c(key, n, lat) bind(c, name='get_lat_f90')
implicit none
integer(c_int), intent(inout) :: key
integer(c_int), intent(in) :: n
real(kind=kind_real),intent(out) :: lat(n)

type(unstructured_grid), pointer :: self

call unstructured_grid_registry%get(key,self)
lat = self%lat

end subroutine get_lat_c

!-------------------------------------------------------------------------------

subroutine get_area_c(key, n, area) bind(c, name='get_area_f90')
implicit none
integer(c_int), intent(inout) :: key
integer(c_int), intent(in) :: n
real(kind=kind_real),intent(out) :: area(n)

type(unstructured_grid), pointer :: self

call unstructured_grid_registry%get(key,self)
area = self%area

end subroutine get_area_c

!-------------------------------------------------------------------------------

subroutine get_vunit_c(key, n, vunit) bind(c, name='get_vunit_f90')
implicit none
integer(c_int), intent(inout) :: key
integer(c_int), intent(in) :: n
real(kind=kind_real),intent(out) :: vunit(n)

type(unstructured_grid), pointer :: self

call unstructured_grid_registry%get(key,self)
vunit = pack(self%vunit,mask=.true.)

end subroutine get_vunit_c

!-------------------------------------------------------------------------------

subroutine get_imask_c(key, n, imask) bind(c, name='get_imask_f90')
implicit none
integer(c_int), intent(inout) :: key
integer,intent(in) :: n
integer,intent(out) :: imask(n)

type(unstructured_grid), pointer :: self

call unstructured_grid_registry%get(key,self)
imask = pack(self%imask,mask=.true.)

end subroutine get_imask_c

!-------------------------------------------------------------------------------

subroutine get_data_c(key, n, fld) bind(c, name='get_data_f90')
implicit none
integer(c_int), intent(inout) :: key
integer(c_int), intent(in) :: n
real(kind=kind_real),intent(out) :: fld(n)

type(unstructured_grid), pointer :: self

call unstructured_grid_registry%get(key,self)
fld = pack(self%fld,mask=.true.)

end subroutine get_data_c

! ------------------------------------------------------------------------------

subroutine create_unstructured_grid(self, nc0a, nl0, nv, nts, lon, lat, area, vunit, imask)
implicit none
type(unstructured_grid), intent(inout) :: self
integer, intent(in) :: nc0a
integer, intent(in) :: nl0
integer, intent(in) :: nv
integer, intent(in) :: nts
real(kind=kind_real), intent(in) :: lon(nc0a)
real(kind=kind_real), intent(in) :: lat(nc0a)
real(kind=kind_real), intent(in) :: area(nc0a)
real(kind=kind_real), intent(in) :: vunit(nc0a,nl0)
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
allocate(self%vunit(nc0a,nl0))
allocate(self%imask(nc0a,nl0))
allocate(self%fld(nc0a,nl0,nv,nts))

! Copy coordinates
self%lon = lon
self%lat = lat
self%area = area
self%vunit = vunit
self%imask = imask

end subroutine create_unstructured_grid

!-------------------------------------------------------------------------------

subroutine delete_unstructured_grid(self)
implicit none
type(unstructured_grid), intent(inout) :: self

! Release memory 
deallocate(self%lon)
deallocate(self%lat)
deallocate(self%area)
deallocate(self%vunit)
deallocate(self%imask)
deallocate(self%fld)

end subroutine delete_unstructured_grid

!-------------------------------------------------------------------------------

end module unstructured_grid_mod
