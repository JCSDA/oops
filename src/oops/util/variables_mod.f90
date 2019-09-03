! (C) Copyright 2019  UCAR. 
!  
! This software is licensed under the terms of the Apache Licence Version 2.0 
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!>  Fortran module handling the variables list

module variables_mod

use fckit_configuration_module, only: fckit_configuration
use iso_c_binding

implicit none
private
public :: oops_vars, oops_vars_create, oops_vars_delete

!> Derived type holding the variable names and count
type :: oops_vars
  integer :: nv                                  !< Number of variables
  character(len=100), allocatable :: fldnames(:) !< Variable identifiers
end type oops_vars

! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

subroutine oops_vars_create(f_conf,self)
implicit none

type(fckit_configuration),intent(in) :: f_conf !< FCKIT configuration
type(oops_vars), intent(inout) :: self         !< OOPS variables

character(len=1023) :: varlist
character(len=:),allocatable :: str

! Get long comma seperated list of variables
call f_conf%get_or_die("variables",str)
varlist = str

! Count number of commas
self%nv = 1 + count(transfer(varlist, 'a', len(varlist)) == ",")

! Allocate array to hold variable names
allocate(self%fldnames(self%nv))

! Place list into fldnames
read(varlist,*) self%fldnames

end subroutine oops_vars_create

! ------------------------------------------------------------------------------

subroutine oops_vars_delete(self)
implicit none

type(oops_vars), intent(inout) :: self

if (allocated(self%fldnames)) deallocate(self%fldnames)
self%nv = 0

end subroutine oops_vars_delete

! ------------------------------------------------------------------------------

end module variables_mod
