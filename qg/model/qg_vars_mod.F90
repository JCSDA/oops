! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

!> Fortran module to handle variables for the QG model

module qg_vars_mod

use iso_c_binding
use config_mod

implicit none
private
public :: qg_vars, qg_vars_create

! ------------------------------------------------------------------------------
!> Fortran derived type to represent QG model variables

type :: qg_vars
  integer :: nv
  character(len=1), allocatable :: fldnames(:) !< Variable identifiers
  logical :: lbc
end type qg_vars

#define LISTED_TYPE qg_vars

!> Linked list interface - defines registry_t type
#include "util/linkedList_i.f"

!> Global registry
type(registry_t) :: qg_vars_registry

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Linked list implementation
#include "util/linkedList_c.f"

! ------------------------------------------------------------------------------

subroutine qg_vars_create(self, kvars)
implicit none
type(qg_vars), intent(inout) :: self
integer(c_int), dimension(*), intent(in) :: kvars
integer :: ii, jj

if (kvars(1)<1 .or. kvars(1)>5) call abor1_ftn ("qg_vars_create: error variables")
if (kvars(kvars(1)+2)/=999) call abor1_ftn ("qg_vars_create: error check")

self%lbc = .false.
self%nv = 0

do jj=1,kvars(1)
  ii=jj+1
  if (kvars(ii)<1 .or. kvars(ii)>5) call abor1_ftn ("qg_vars_create: unknown index")
  if (kvars(ii)==5) then
    self%lbc = .true.
  else
    self%nv=self%nv+1
  endif
enddo

allocate(self%fldnames(self%nv))

ii = 0
do jj=1,kvars(1)
  if (kvars(jj+1)/=5) ii=ii+1
  if (kvars(jj+1)==1) self%fldnames(ii)="x"
  if (kvars(jj+1)==2) self%fldnames(ii)="q"
  if (kvars(jj+1)==3) self%fldnames(ii)="u"
  if (kvars(jj+1)==4) self%fldnames(ii)="v"
enddo

end subroutine qg_vars_create

! ------------------------------------------------------------------------------

end module qg_vars_mod
