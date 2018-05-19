! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

!>  Define interface for C++ duration code called from Fortran

!-------------------------------------------------------------------------------
interface
!-------------------------------------------------------------------------------
type(c_ptr) function c_duration_construct(str) bind(C,name='duration_construct_f')
  use, intrinsic :: iso_c_binding
  implicit none
  character(kind=c_char, len=1), intent(in) :: str(21)
end function c_duration_construct
!-------------------------------------------------------------------------------
subroutine c_duration_destruct(ptr) bind(C,name='duration_destruct_f')
  use, intrinsic :: iso_c_binding
  implicit none
  type(c_ptr), value :: ptr
end subroutine c_duration_destruct
!-------------------------------------------------------------------------------
integer(kind=c_int64_t) function c_duration_str_int(str) bind(C,name='duration_str_int_f')
  use, intrinsic :: iso_c_binding
  implicit none
  character(kind=c_char, len=1), intent(in) :: str(*)
end function c_duration_str_int
!-------------------------------------------------------------------------------
integer(kind=c_int64_t) function c_duration_int(ptr) bind(C,name='duration_int_f')
  use, intrinsic :: iso_c_binding
  implicit none
  type(c_ptr), value :: ptr
end function c_duration_int
!-------------------------------------------------------------------------------
subroutine c_duration_int_str(ksecs, str) bind(C,name='duration_int_str_f')
  use, intrinsic :: iso_c_binding
  implicit none
  integer(kind=c_int64_t), intent(in) :: ksecs
  character(kind=c_char, len=1), intent(out) :: str(21)
end subroutine c_duration_int_str
!-------------------------------------------------------------------------------
end interface
!-------------------------------------------------------------------------------
