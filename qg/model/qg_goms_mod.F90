! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

!> Fortran module handling interpolated (to obs locations) model variables

module qg_goms_mod

use iso_c_binding
use qg_vars_mod
use kinds

implicit none
private
public :: qg_goms, gom_setup
public :: qg_goms_registry

! ------------------------------------------------------------------------------

!> Fortran derived type to hold interpolated fields required by the obs operators
type :: qg_goms
  integer :: nobs
  integer :: nvar
  integer :: used
  integer, allocatable :: indx(:)
  real(kind=kind_real), allocatable :: values(:,:)
  character(len=1), allocatable :: variables(:)
  logical :: lalloc
end type qg_goms

#define LISTED_TYPE qg_goms

!> Linked list interface - defines registry_t type
#include "linkedList_i.f"

!> Global registry
type(registry_t) :: qg_goms_registry

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Linked list implementation
#include "linkedList_c.f"

! ------------------------------------------------------------------------------

subroutine c_qg_gom_create(c_key_self) bind(c,name='qg_gom_create_f90')

implicit none
integer(c_int), intent(inout) :: c_key_self

type(qg_goms), pointer :: self
call qg_goms_registry%init()
call qg_goms_registry%add(c_key_self)
call qg_goms_registry%get(c_key_self, self)

self%lalloc = .false.

end subroutine c_qg_gom_create

! ------------------------------------------------------------------------------

subroutine gom_setup(self, vars, kobs)
implicit none
type(qg_goms), intent(inout) :: self
type(qg_vars), intent(in) :: vars
integer, intent(in) :: kobs(:)

self%nobs=size(kobs)
self%nvar=vars%nv
self%used=0

allocate(self%indx(self%nobs))
self%indx(:)=kobs(:)

allocate(self%variables(self%nvar))
self%variables(:)=vars%fldnames(:)

allocate(self%values(self%nvar,self%nobs))

self%lalloc = .true.

end subroutine gom_setup

! ------------------------------------------------------------------------------

subroutine c_qg_gom_delete(c_key_self) bind(c,name='qg_gom_delete_f90')

implicit none
integer(c_int), intent(inout) :: c_key_self

type(qg_goms), pointer :: self

call qg_goms_registry%get(c_key_self, self)
if (self%lalloc) then
  deallocate(self%values)
  deallocate(self%indx)
  deallocate(self%variables)
endif
call qg_goms_registry%remove(c_key_self)

end subroutine c_qg_gom_delete

! ------------------------------------------------------------------------------

subroutine c_qg_gom_zero(c_key_self) bind(c,name='qg_gom_zero_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
type(qg_goms), pointer :: self
call qg_goms_registry%get(c_key_self, self)
self%values(:,:)=0.0_kind_real
end subroutine c_qg_gom_zero

! ------------------------------------------------------------------------------

subroutine c_qg_gom_dotprod(c_key_self, c_key_other, prod) bind(c,name='qg_gom_dotprod_f90')
implicit none
integer(c_int), intent(in) :: c_key_self, c_key_other
real(c_double), intent(inout) :: prod
type(qg_goms), pointer :: self, other
integer :: jo, jv

call qg_goms_registry%get(c_key_self, self)
call qg_goms_registry%get(c_key_other, other)
prod=0.0_kind_real
do jo=1,self%nobs
  do jv=1,self%nvar
    prod=prod+self%values(jv,jo)*other%values(jv,jo)
  enddo
enddo

end subroutine c_qg_gom_dotprod

! ------------------------------------------------------------------------------

subroutine c_qg_gom_minmaxavg(c_key_self, kobs, pmin, pmax, prms) bind(c,name='qg_gom_minmaxavg_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(inout) :: kobs
real(c_double), intent(inout) :: pmin, pmax, prms
type(qg_goms), pointer :: self

call qg_goms_registry%get(c_key_self, self)

kobs = self%nobs
pmin=minval(self%values(:,:))
pmax=maxval(self%values(:,:))
prms=sqrt(sum(self%values(:,:)**2)/real(self%nobs*self%nvar,kind_real))

end subroutine c_qg_gom_minmaxavg

! ------------------------------------------------------------------------------

end module qg_goms_mod
