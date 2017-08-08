! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!>  Fortran module for handling generic nicas localization

module nicas_mod

use iso_c_binding
use kinds
use config_mod
use unstructured_grid_mod

implicit none
private
public nicas, create_nicas, delete_nicas, nicas_multiply

! ------------------------------------------------------------------------------

!>  Derived type containing the data

type nicas
  real(kind=kind_real) :: length
end type nicas

! ------------------------------------------------------------------------------

#define LISTED_TYPE nicas

!> Linked list interface - defines registry_t type
#include "util/linkedList_i.f"

!> Global registry
type(registry_t) :: nicas_registry

!-------------------------------------------------------------------------------
contains
!-------------------------------------------------------------------------------

!> Linked list implementation
#include "util/linkedList_c.f"

! ------------------------------------------------------------------------------
!  C++ interfaces
! ------------------------------------------------------------------------------

subroutine create_nicas_c(key, c_conf, cnh, clats, clons, cnv, clevs, cmask) bind(c, name='create_nicas_f90')
implicit none
integer(c_int), intent(inout) :: key
type(c_ptr), intent(in) :: c_conf
integer(c_int), intent(in) :: cnh, cnv
real(c_double), intent(in) :: clats(cnh), clons(cnh), clevs(cnv)
integer(c_int), intent(in) :: cmask(cnh*cnv)
type(nicas), pointer :: self
real(kind=kind_real) :: lats(cnh), lons(cnh), levs(cnv)
integer :: mask(cnh*cnv)
call nicas_registry%init()
call nicas_registry%add(key)
call nicas_registry%get(key,self)
lats(:)=clats(:)
lons(:)=clons(:)
levs(:)=clevs(:)
mask(:)=cmask(:)
call create_nicas(self, c_conf, lats, lons, levs, mask)
end subroutine create_nicas_c

! ------------------------------------------------------------------------------

subroutine delete_nicas_c(key) bind(c, name='delete_nicas_f90')
implicit none
integer(c_int), intent(inout) :: key
type(nicas), pointer :: self
call nicas_registry%get(key,self)
call delete_nicas(self)
call nicas_registry%remove(key)
end subroutine delete_nicas_c

! ------------------------------------------------------------------------------

subroutine nicas_multiply_c(key, idx) bind(c, name='nicas_multiply_f90')
implicit none
integer(c_int), intent(in) :: key
integer(c_int), intent(in) :: idx
type(nicas), pointer :: self
type(unstructured_grid), pointer :: udx
call nicas_registry%get(key,self)
call unstructured_grid_registry%get(idx, udx)
call nicas_multiply(self, udx)
end subroutine nicas_multiply_c

! ------------------------------------------------------------------------------
!  End C++ interfaces
! ------------------------------------------------------------------------------

subroutine create_nicas(self, c_conf, lats, lons, levs, mask)
implicit none
type(nicas), intent(inout) :: self
type(c_ptr), intent(in) :: c_conf
real(kind=kind_real), intent(in) :: lats(:), lons(:), levs(:)
integer, intent(in) :: mask(:)

!self%length = config_get_real(c_conf, "length_scale")
!write(*,*)'lats : ',lats(:)
!write(*,*)'lons : ',lons(:)
!write(*,*)'levs : ',levs(:)
!write(*,*)'mask : ',mask(:)

end subroutine create_nicas

!-------------------------------------------------------------------------------

subroutine delete_nicas(self)
implicit none
type(nicas), intent(inout) :: self

end subroutine delete_nicas

!-------------------------------------------------------------------------------

subroutine nicas_multiply(self,dx)
implicit none
type(nicas), intent(in) :: self
type(unstructured_grid), intent(inout) :: dx

end subroutine nicas_multiply

!-------------------------------------------------------------------------------

end module nicas_mod
