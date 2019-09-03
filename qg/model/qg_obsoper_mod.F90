! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

module qg_obsoper_mod

use fckit_configuration_module, only: fckit_configuration
use iso_c_binding
use kinds

implicit none

private
public :: qg_obsoper
public :: qg_obsoper_registry
public :: qg_oper_setup
! ------------------------------------------------------------------------------
type :: qg_obsoper
  integer :: ncol              !< Number of columns
end type qg_obsoper

#define LISTED_TYPE qg_obsoper

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry
type(registry_t) :: qg_obsoper_registry
! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Linked list implementation
#include "oops/util/linkedList_c.f"
! ------------------------------------------------------------------------------
subroutine qg_oper_setup(self,f_conf,svars,ncol)

implicit none

! Passed variables
type(qg_obsoper),intent(inout) :: self         !< Observations
type(fckit_configuration),intent(in) :: f_conf !< FCKIT configuration
character(len=*),intent(in) :: svars(:)        !< Variables
integer :: ncol                                !< Number of columns

! Set number of columns
self%ncol = ncol

end subroutine qg_oper_setup
! ------------------------------------------------------------------------------
end module qg_obsoper_mod
