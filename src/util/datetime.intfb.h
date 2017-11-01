! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

!>  Define interface for C++ datetime code called from Fortran

!-------------------------------------------------------------------------------
interface
!-------------------------------------------------------------------------------
type(c_ptr) function c_datetime_construct(str) bind(C,name='datetime_construct_f')
  use, intrinsic :: iso_c_binding
  implicit none
  character(kind=c_char, len=1), intent(in) :: str(21)
end function c_datetime_construct
!-------------------------------------------------------------------------------
subroutine c_datetime_set(str, ptr) bind(C,name='datetime_set_f')
  use, intrinsic :: iso_c_binding
  implicit none
  character(kind=c_char, len=1), intent(in) :: str(21)
  type(c_ptr), value :: ptr
end subroutine c_datetime_set
!-------------------------------------------------------------------------------
subroutine c_datetime_destruct(ptr) bind(C,name='datetime_destruct_f')
  use, intrinsic :: iso_c_binding
  implicit none
  type(c_ptr), value :: ptr
end subroutine c_datetime_destruct
!-------------------------------------------------------------------------------
subroutine c_datetime_string(ptr, c_string) bind(C,name='datetime_string_f')
  use, intrinsic :: iso_c_binding
  implicit none
  type(c_ptr), value :: ptr
  character(kind=c_char,len=1), intent(inout) :: c_string(21)
end subroutine c_datetime_string
!-------------------------------------------------------------------------------
subroutine c_datetime_getints(ptr, kdate, ksecs) bind(C,name='datetime_getints_f')
  use, intrinsic :: iso_c_binding
  implicit none
  type(c_ptr), value :: ptr
  integer(kind=c_int64_t), intent(inout) :: kdate
  integer(kind=c_int), intent(inout) :: ksecs
end subroutine c_datetime_getints
!-------------------------------------------------------------------------------
subroutine c_datetime_setints(ptr, kdate, ksecs) bind(C,name='datetime_setints_f')
  use, intrinsic :: iso_c_binding
  implicit none
  type(c_ptr), value :: ptr
  integer(kind=c_int64_t), intent(in) :: kdate
  integer(kind=c_int), intent(in) :: ksecs
end subroutine c_datetime_setints
!-------------------------------------------------------------------------------
integer(kind=c_int64_t) function c_datetime_diff(d1, d2) bind(C,name='datetime_diff_f')
  use, intrinsic :: iso_c_binding
  implicit none
  type(c_ptr), value :: d1, d2
end function c_datetime_diff
!-------------------------------------------------------------------------------
subroutine c_datetime_update(ptt,kdt) bind(C,name='datetime_update_f')
  use, intrinsic :: iso_c_binding
  implicit none
  type(c_ptr), value :: ptt
  integer(c_int64_t), intent(in) :: kdt
end subroutine c_datetime_update
!-------------------------------------------------------------------------------
end interface
!-------------------------------------------------------------------------------
