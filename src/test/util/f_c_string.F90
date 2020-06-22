!
! (C) Copyright 2019 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
!

!> Test interface for Fortran/C++ interface tools for strings

module test_f_c_string

use fckit_configuration_module, only: fckit_configuration
use, intrinsic :: iso_c_binding
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
  
subroutine c_test_push_string_vector(c_conf,vec) bind(c,name='test_push_string_vector_f')
implicit none
type(c_ptr), intent(in) :: c_conf
type(c_ptr), intent(in), value :: vec

character(len=:), allocatable :: fort_vec(:)
type(fckit_configuration) :: f_conf

f_conf = fckit_configuration(c_conf)
call f_conf%get_or_die("string_vec",fort_vec)

call f_c_push_string_vector(vec, fort_vec)

end subroutine c_test_push_string_vector

!-------------------------------------------------------------------------------

end module test_f_c_string
