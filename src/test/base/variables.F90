!
! (C) Copyright 2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!

!> Test interface for oops variables

module test_oops_variables

use fckit_configuration_module, only: fckit_configuration
use, intrinsic :: iso_c_binding
use kinds
use oops_variables_mod
use fckit_log_module, only : fckit_log

implicit none
private

integer, parameter :: var_length = 10
integer, parameter :: max_string = 800

!-------------------------------------------------------------------------------
contains
!-------------------------------------------------------------------------------
!> Test the Fortran interface to Variables
!
subroutine c_test_vars_interface(c_conf, c_vars) bind(c,name='test_vars_interface_f')
implicit none

type(c_ptr), intent(in), value :: c_conf
type(c_ptr), intent(in), value :: c_vars

character(len=*), parameter :: myname_="test_vars_interface"
type(fckit_configuration) :: f_conf
type(oops_variables) vars
character(len=:), allocatable :: test_vars(:), varlist(:)
character(var_length) :: varname
character(max_string) :: err_msg
integer :: jvar

f_conf = fckit_configuration(c_conf)
vars = oops_variables(c_vars)

! Get variable list from config file and push to C++
call f_conf%get_or_die("test variables",test_vars)
call vars%push_back(test_vars)

! add another variable to check single name interface
varname = "newvar"
call vars%push_back(varname)

! check the clear method
call vars%clear()
if (vars%nvars() /= 0) then
  write(err_msg,*) myname_ // " clear did not remove all variables"
  call abor1_ftn(err_msg)
endif

! add the variables again
call vars%push_back(test_vars)
call vars%push_back(varname)

! check varlist method
varlist = vars%varlist()

do jvar = 1, size(test_vars)
    if (trim(test_vars(jvar)) /= trim(varlist(jvar))) then
        write(err_msg,*) myname_ // " varlist incorrect: ", jvar, &
        & trim(test_vars(jvar)) // " /= " // trim(varlist(jvar))
        call abor1_ftn(err_msg)
    endif
enddo

jvar = size(varlist)
if (trim(varname) /= trim(varlist(jvar))) then
   write(err_msg,*) myname_ // " varlist incorrect: ", jvar, &
       & trim(varname) // " /= " // trim(varlist(jvar))
    call abor1_ftn(err_msg)
endif

! check the find method
jvar = vars%find(varname)
if (jvar < 0 .or. jvar > vars%nvars()) then
  write(err_msg,*) myname_ // " find returned an out-of-bound index ", jvar
  call abor1_ftn(err_msg)
endif
if (trim(vars%variable(jvar)) /= varname) then
  write(err_msg,*) myname_ // " find returned an incorrect index ", jvar
  call abor1_ftn(err_msg)
endif

end subroutine c_test_vars_interface

end module test_oops_variables
