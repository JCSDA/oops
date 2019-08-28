! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

module qg_geom_iter_mod

use iso_c_binding
use kinds
use qg_geom_mod

implicit none

private
public :: qg_geom_iter
public :: qg_geom_iter_registry
public :: qg_geom_iter_setup,qg_geom_iter_clone,qg_geom_iter_equals,qg_geom_iter_current,qg_geom_iter_next
! ------------------------------------------------------------------------------
type :: qg_geom_iter
  type(qg_geom),pointer :: geom => null() !< Geometry
  integer :: ilon = 1                     !< Longitude index
  integer :: ilat = 1                     !< Latitude index
end type qg_geom_iter

#define LISTED_TYPE qg_geom_iter

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry
type(registry_t) :: qg_geom_iter_registry
! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
! Public
! ------------------------------------------------------------------------------
!> Linked list implementation
#include "oops/util/linkedList_c.f"
! ------------------------------------------------------------------------------
!> Setup for the QG model's geometry iterator
subroutine qg_geom_iter_setup(self,geom,ind)

! Passed variables
type(qg_geom_iter),intent(inout) :: self !< Geometry iterator
type(qg_geom),pointer,intent(in) :: geom !< Geometry
integer,intent(in) :: ind                !< Index

! Associate geometry
self%geom => geom

! Define ilon/ilat
self%ilat = (ind-1)/geom%nx+1
self%ilon = ind-(self%ilat-1)*geom%nx

end subroutine qg_geom_iter_setup
! ------------------------------------------------------------------------------
!> Clone for the QG model's geometry iterator
subroutine qg_geom_iter_clone(self,other)

! Passed variables
type(qg_geom_iter),intent(inout) :: self !< Geometry iterator
type(qg_geom_iter),intent(in) :: other   !< Other geometry iterator

! Associate geometry
self%geom => other%geom

! Copy ilon/ilat
self%ilon = other%ilon
self%ilat = other%ilat

end subroutine qg_geom_iter_clone
! ------------------------------------------------------------------------------
!> Check for the QG model's geometry iterator equality
subroutine qg_geom_iter_equals(self,other,equals)

! Passed variables
type(qg_geom_iter),intent(in) :: self  !< Geometry iterator
type(qg_geom_iter),intent(in) :: other !< Other geometry iterator
integer,intent(out) :: equals          !< Equality flag

! Initialization
equals = 0

! Check equality
if (associated(self%geom,other%geom).and.(self%ilon==other%ilon).and.(self%ilat==other%ilat)) equals = 1

end subroutine qg_geom_iter_equals
! ------------------------------------------------------------------------------
!> Get geometry iterator current lat/lon
subroutine qg_geom_iter_current(self,lat,lon)

! Passed variables
type(qg_geom_iter),intent(in) :: self !< Geometry iterator
real(kind_real),intent(out) :: lat    !< Latitude
real(kind_real),intent(out) :: lon    !< Longitude

! Check ilon/ilat
if (self%ilon*self%ilat>self%geom%nx*self%geom%ny) call abor1_ftn('qg_geom_iter_current: iterator out of bounds')

! Get lat/lon
lat = self%geom%lat(self%ilon,self%ilat)
lon = self%geom%lon(self%ilon,self%ilat)

end subroutine qg_geom_iter_current
! ------------------------------------------------------------------------------
!> Update geometry iterator to next point
subroutine qg_geom_iter_next(self)

! Passed variables
type(qg_geom_iter),intent(inout) :: self !< Geometry iterator

! Update ilon/ilat
if (self%ilon==self%geom%nx) then
  self%ilon = 1
  self%ilat = self%ilat+1
else
  self%ilon = self%ilon+1
endif

end subroutine qg_geom_iter_next
! ------------------------------------------------------------------------------
end module qg_geom_iter_mod
