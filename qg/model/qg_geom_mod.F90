! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

!> Fortran module handling geometry for the QG model

module qg_geom_mod

use iso_c_binding
use config_mod
use kinds

implicit none
private
public :: qg_geom
public :: qg_geom_registry

! ------------------------------------------------------------------------------

!> Fortran derived type to hold geometry data for the QG model
type :: qg_geom
  integer :: nx
  integer :: ny
  real(kind=kind_real),allocatable :: lats(:)
  real(kind=kind_real),allocatable :: lons(:)
  real(kind=kind_real),allocatable :: areas(:,:)
end type qg_geom

#define LISTED_TYPE qg_geom

!> Linked list interface - defines registry_t type
#include "util/linkedList_i.f"

!> Global registry
type(registry_t) :: qg_geom_registry

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Linked list implementation
#include "util/linkedList_c.f"

! ------------------------------------------------------------------------------

subroutine c_qg_geo_setup(c_key_self, c_conf) bind(c,name='qg_geo_setup_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self
type(c_ptr), intent(in)    :: c_conf

integer :: ix,iy
real(kind=kind_real) :: dx,dytot,dy
type(qg_geom), pointer :: self

call qg_geom_registry%init()
call qg_geom_registry%add(c_key_self)
call qg_geom_registry%get(c_key_self,self)

self%nx = config_get_int(c_conf, "nx")
self%ny = config_get_int(c_conf, "ny")

allocate(self%lons(self%nx))
allocate(self%lats(self%ny))
allocate(self%areas(self%nx,self%ny))

dx = 360.0 / real(self%nx,kind=kind_real);
dytot = 360.0 * real(self%ny,kind=kind_real) / real(self%nx,kind=kind_real);
dy = dytot / real(self%ny,kind=kind_real);
do ix=1,self%nx
   self%lons(ix) = -180.0+(real(ix,kind=kind_real)-0.5)*dx
end do
do iy=1,self%ny
   self%lats(iy) = -0.5*dytot+(real(iy,kind=kind_real)-0.5)*dy;
end do
do iy=1,self%ny
   self%areas(:,iy) = 6.371e6**2*cos(self%lats(iy)*acos(-1.0)/180.0)*dx*dy
end do

end subroutine c_qg_geo_setup

! ------------------------------------------------------------------------------

subroutine c_qg_geo_clone(c_key_self, c_key_other) bind(c,name='qg_geo_clone_f90')
implicit none
integer(c_int), intent(in   ) :: c_key_self
integer(c_int), intent(inout) :: c_key_other

type(qg_geom), pointer :: self, other

call qg_geom_registry%add(c_key_other)
call qg_geom_registry%get(c_key_other, other)
call qg_geom_registry%get(c_key_self , self )
other%nx = self%nx
other%ny = self%ny
allocate(other%lons(other%nx))
allocate(other%lats(other%ny))
allocate(other%areas(other%nx,other%ny))
other%lons = self%lons
other%lats = self%lats
other%areas = self%areas

end subroutine c_qg_geo_clone

! ------------------------------------------------------------------------------

subroutine c_qg_geo_delete(c_key_self) bind(c,name='qg_geo_delete_f90')

implicit none
integer(c_int), intent(inout) :: c_key_self     

call qg_geom_registry%remove(c_key_self)

end subroutine c_qg_geo_delete

! ------------------------------------------------------------------------------

subroutine c_qg_geo_info(c_key_self, c_nx, c_ny) bind(c,name='qg_geo_info_f90')
implicit none
integer(c_int), intent(in   ) :: c_key_self
integer(c_int), intent(inout) :: c_nx
integer(c_int), intent(inout) :: c_ny
type(qg_geom), pointer :: self

call qg_geom_registry%get(c_key_self , self )
c_nx = self%nx
c_ny = self%ny

end subroutine c_qg_geo_info

! ------------------------------------------------------------------------------

end module qg_geom_mod
