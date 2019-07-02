!
! (C) Copyright 2019 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
!

!> Test interface for Fortran/C++ interface tools for strings

module test_f_c_string

use, intrinsic :: iso_c_binding
use config_mod
use kinds
use string_f_c_mod
use fckit_log_module, only : fckit_log

implicit none
private

integer, parameter :: max_string = 60

!-------------------------------------------------------------------------------
contains
!-------------------------------------------------------------------------------
!> Pass a string vector from Fortran to C++
!
  
subroutine c_test_push_string_vector(c_conf,vec,var) bind(c,name='test_push_string_vector_f')
implicit none
type(c_ptr), intent(in) :: c_conf
type(c_ptr), intent(in), value :: vec
type(c_ptr), intent(in), value :: var

character(len=max_string), allocatable :: fort_vec(:)

fort_vec = config_get_string_vector(c_conf, max_string, "string_vec")

call f_c_push_string_vector(vec, fort_vec)

call f_c_push_string_varlist(var, fort_vec)

end subroutine c_test_push_string_vector

!-------------------------------------------------------------------------------

end module test_f_c_string
