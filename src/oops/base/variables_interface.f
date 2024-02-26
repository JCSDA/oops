!
! (C) Copyright 2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!>  Define interface for C++ Variables code called from Fortran

!-------------------------------------------------------------------------------
interface
!-------------------------------------------------------------------------------

type(c_ptr) function c_variables_empty_ctor() bind(C, name='variables_empty_ctor_f')
  use, intrinsic :: iso_c_binding
  implicit none
end function c_variables_empty_ctor

!-------------------------------------------------------------------------------

subroutine c_variables_destruct(vars) bind(C, name='variables_destruct_f')
  use, intrinsic :: iso_c_binding
  implicit none

  type(c_ptr), value :: vars
end subroutine c_variables_destruct

!-------------------------------------------------------------------------------

subroutine c_variables_push_back(vars, str) bind(C, name='variables_push_back_f')
  use, intrinsic :: iso_c_binding
  implicit none

  type(c_ptr), value :: vars
  character(kind=c_char,len=1), intent(in) :: str(*)
end subroutine c_variables_push_back

!-------------------------------------------------------------------------------

subroutine c_variables_clear(vars) bind(C, name='variables_clear_f')
  use, intrinsic :: iso_c_binding
  implicit none

  type(c_ptr), value :: vars
end subroutine c_variables_clear

!-------------------------------------------------------------------------------

integer(kind=c_size_t) function c_variables_size(vars) bind(C,name='variables_size_f')
  use, intrinsic :: iso_c_binding
  implicit none

  type(c_ptr), value :: vars
end function c_variables_size

!-------------------------------------------------------------------------------

subroutine c_variables_getvariablelength(vars, jj, lcname) bind (C,name='variables_getvariablelength_f')
  use, intrinsic :: iso_c_binding
  implicit none
  type(c_ptr), value :: vars
  integer(c_size_t),intent(in) :: jj
  integer(c_size_t),intent(inout) :: lcname
end subroutine c_variables_getvariablelength

!-------------------------------------------------------------------------------

subroutine c_variables_getvariable(vars, jj, lcname, lfname, cname) bind (C,name='variables_getvariable_f')
  use, intrinsic :: iso_c_binding
  implicit none
  type(c_ptr), value :: vars
  integer(c_size_t),intent(in) :: jj
  integer(c_size_t),intent(inout) :: lcname
  integer(c_size_t),intent(in) :: lfname    
  character(kind=c_char,len=1), intent(inout) :: cname(*)
end subroutine c_variables_getvariable

!-------------------------------------------------------------------------------

logical*1 function c_variables_has(vars, str) bind(C, name='variables_has_f')
  use, intrinsic :: iso_c_binding
  implicit none

  type(c_ptr), value :: vars
  character(kind=c_char,len=1), intent(in) :: str(*)
end function c_variables_has

!-------------------------------------------------------------------------------

integer(c_int) function c_variables_find(vars, str) bind(C, name='variables_find_f')
  use, intrinsic :: iso_c_binding
  implicit none

  type(c_ptr), value :: vars
  character(kind=c_char,len=1), intent(in) :: str(*)
end function c_variables_find
!-------------------------------------------------------------------------------
end interface
!-------------------------------------------------------------------------------

