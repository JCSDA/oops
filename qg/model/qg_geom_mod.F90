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
use mpi

implicit none
private
public :: qg_geom
public :: qg_geom_registry

! ------------------------------------------------------------------------------

!> Fortran derived type to hold geometry data for the QG model
type :: qg_geom
  integer :: nx
  integer :: ny
  real(kind=kind_real),allocatable :: lat(:)
  real(kind=kind_real),allocatable :: lon(:)
  real(kind=kind_real),allocatable :: area(:,:)
  integer,allocatable :: myproc(:,:)
end type qg_geom

#define LISTED_TYPE qg_geom

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry
type(registry_t) :: qg_geom_registry

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

!> Linked list implementation
#include "oops/util/linkedList_c.f"

! ------------------------------------------------------------------------------

subroutine c_qg_geo_setup(c_key_self, c_conf) bind(c,name='qg_geo_setup_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self
type(c_ptr), intent(in)    :: c_conf

integer :: ix,iy,nproc,info,myproc
real(kind=kind_real) :: dx,dy
real(kind=kind_real),parameter :: pi = acos(-1.0)
real(kind=kind_real),parameter :: req = 6371229.0
type(qg_geom), pointer :: self

call qg_geom_registry%init()
call qg_geom_registry%add(c_key_self)
call qg_geom_registry%get(c_key_self,self)

self%nx = config_get_int(c_conf, "nx")
self%ny = config_get_int(c_conf, "ny")

! Allocate
allocate(self%lon(self%nx))
allocate(self%lat(self%ny))
allocate(self%area(self%nx,self%ny))
allocate(self%myproc(self%nx,self%ny))

! Define longitude/latitude
dx = 2.0*pi/real(self%nx,kind=kind_real)
dy = pi/real(self%ny,kind=kind_real)
do ix=1,self%nx
   self%lon(ix) = -pi+(real(ix,kind=kind_real)-0.5_kind_real)*dx
end do
do iy=1,self%ny
   self%lat(iy) = -0.5*pi+(real(iy,kind=kind_real)-0.5_kind_real)*dy
end do

! Define area
do iy=1,self%ny
   self%area(:,iy) = 2.0*pi*req**2*(sin((self%lat(iy)+0.5_kind_real*dy))-sin((self%lat(iy)-0.5_kind_real*dy))) &
            & /real(self%nx,kind=kind_real)
end do

! Define processor
call mpi_comm_size(mpi_comm_world,nproc,info)
myproc = 1
do ix=1,self%nx
   self%myproc(ix,:) = myproc
   if (ix>myproc*self%nx/nproc) myproc = myproc+1
end do

! Convert longitude/latitude to degrees
self%lon = self%lon*180_kind_real/pi
self%lat = self%lat*180_kind_real/pi

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
allocate(other%lon(other%nx))
allocate(other%lat(other%ny))
allocate(other%area(other%nx,other%ny))
allocate(other%myproc(other%nx,other%ny))
other%lon = self%lon
other%lat = self%lat
other%area = self%area
other%myproc = self%myproc

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
