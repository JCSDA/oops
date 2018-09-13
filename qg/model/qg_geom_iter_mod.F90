! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

!> Fortran module handling geometry iterator for the QG model

module qg_geom_iter_mod

use iso_c_binding
use config_mod
use kinds
use qg_geom_mod

implicit none
private
public :: qg_geom_iter
public :: qg_geom_iter_registry

! ------------------------------------------------------------------------------

!> Fortran derived type to hold geom_iteretry data for the QG model
type :: qg_geom_iter
  type(qg_geom), pointer :: geom => null()
  integer :: ilat = 1
  integer :: ilon = 1
end type qg_geom_iter

#define LISTED_TYPE qg_geom_iter

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry
type(registry_t) :: qg_geom_iter_registry

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

!> Linked list implementation
#include "oops/util/linkedList_c.f"

! ------------------------------------------------------------------------------

subroutine c_qg_geo_iter_setup(c_key_self, c_key_geom, c_index) bind(c,name='qg_geo_iter_setup_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self
integer(c_int), intent(in)    :: c_key_geom
integer(c_int), intent(in)    :: c_index

type(qg_geom_iter), pointer :: self
type(qg_geom), pointer :: geom

call qg_geom_iter_registry%init()
call qg_geom_iter_registry%add(c_key_self)
call qg_geom_iter_registry%get(c_key_self,self)

call qg_geom_registry%get(c_key_geom, geom)

self%geom => geom
self%ilat = (c_index-1)/geom%nx + 1
self%ilon = c_index - (self%ilat-1)*geom%nx

end subroutine c_qg_geo_iter_setup

! ------------------------------------------------------------------------------

subroutine c_qg_geo_iter_clone(c_key_self, c_key_other) bind(c,name='qg_geo_iter_clone_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self
integer(c_int), intent(in)    :: c_key_other

type(qg_geom_iter), pointer :: self, other

call qg_geom_iter_registry%init()
call qg_geom_iter_registry%add(c_key_self)
call qg_geom_iter_registry%get(c_key_self,self)

call qg_geom_iter_registry%get(c_key_other, other)

self%geom => other%geom
self%ilon = other%ilon
self%ilat = other%ilat

end subroutine c_qg_geo_iter_clone

! ------------------------------------------------------------------------------

subroutine c_qg_geo_iter_delete(c_key_self) bind(c,name='qg_geo_iter_delete_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self     

call qg_geom_iter_registry%remove(c_key_self)

end subroutine c_qg_geo_iter_delete

! ------------------------------------------------------------------------------

subroutine c_qg_geo_iter_equals(c_key_self, c_key_other, c_equals) bind(c,name='qg_geo_iter_equals_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self
integer(c_int), intent(in)    :: c_key_other
integer(c_int), intent(inout) :: c_equals

type(qg_geom_iter), pointer :: self, other

call qg_geom_iter_registry%get(c_key_self, self)
call qg_geom_iter_registry%get(c_key_other, other)

c_equals = 0
if (associated(self%geom, other%geom) .and. &
             (self%ilon == other%ilon) .and. (self%ilat == other%ilat)) then 
  c_equals = 1
endif

end subroutine c_qg_geo_iter_equals

! ------------------------------------------------------------------------------

subroutine c_qg_geo_iter_current(c_key_self, c_lat, c_lon) bind(c,name='qg_geo_iter_current_f90')
implicit none
integer(c_int), intent(in   ) :: c_key_self
real(c_double), intent(inout) :: c_lat
real(c_double), intent(inout) :: c_lon
type(qg_geom_iter), pointer :: self

call qg_geom_iter_registry%get(c_key_self , self )

if (self%ilon*self%ilat > self%geom%nx * self%geom%ny) then
  print *, 'qg_geo_iter_current: iterator out of bounds'
  stop
endif

c_lon = self%geom%lon(self%ilon)
c_lat = self%geom%lat(self%ilat)

end subroutine c_qg_geo_iter_current

! ------------------------------------------------------------------------------

subroutine c_qg_geo_iter_next(c_key_self) bind(c,name='qg_geo_iter_next_f90')
implicit none
integer(c_int), intent(in   ) :: c_key_self
type(qg_geom_iter), pointer :: self

call qg_geom_iter_registry%get(c_key_self , self )

if (self%ilon == self%geom%nx) then
  self%ilon = 1; self%ilat = self%ilat + 1
else
  self%ilon = self%ilon + 1
endif

end subroutine c_qg_geo_iter_next


end module qg_geom_iter_mod
