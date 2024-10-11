! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

!> Abort function

!> Prints an error message and stops the execution.

subroutine abor1_ftn(cderror)
use, intrinsic :: iso_c_binding, only: c_char
use string_f_c_mod
#ifdef NAG
use f90_unix_proc, only: exit
#endif
implicit none

character(len=*), intent(in) :: cderror !< The error message

character(kind=c_char,len=1), allocatable :: c_string(:)

interface
  subroutine abor1_cpp(str) BIND(C,NAME='abor1_cpp_')
    use, intrinsic :: iso_c_binding
    implicit none
    character(kind=c_char, len=1), intent(in) :: str(*)
  end subroutine abor1_cpp
end interface

call f_c_string(cderror, c_string)

call abor1_cpp(c_string)
call exit(1) ! just in case!!!

end subroutine abor1_ftn
