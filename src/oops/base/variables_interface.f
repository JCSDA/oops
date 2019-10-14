!
! (C) Copyright 2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!>  Define interface for C++ Variables code called from Fortran

!-------------------------------------------------------------------------------
interface
!-------------------------------------------------------------------------------

subroutine c_variables_push_back(vars, str) bind(C, name='variables_push_back_f')
  use, intrinsic :: iso_c_binding
  implicit none

  type(c_ptr), value :: vars
  character(kind=c_char,len=1), intent(in) :: str(*)
end subroutine c_variables_push_back

!-------------------------------------------------------------------------------

integer(kind=c_size_t) function c_variables_size(vars) bind(C,name='variables_size_f')
  use, intrinsic :: iso_c_binding
  implicit none

  type(c_ptr), value :: vars
end function c_variables_size

!-------------------------------------------------------------------------------

subroutine c_variables_getvariable(vars, jj, lcname, cname) bind (C,name='variables_getvariable_f')
  use, intrinsic :: iso_c_binding
  implicit none
  type(c_ptr), value :: vars
  integer(c_size_t),intent(in) :: jj
  integer(c_size_t),intent(in) :: lcname
  character(kind=c_char,len=1), intent(inout) :: cname(*)
end subroutine c_variables_getvariable

!-------------------------------------------------------------------------------
end interface
!-------------------------------------------------------------------------------

