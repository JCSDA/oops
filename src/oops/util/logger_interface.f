! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

!>  Define interface for C++ config code called from Fortran

!-------------------------------------------------------------------------------
interface
!-------------------------------------------------------------------------------

subroutine c_log_info(msg,newl,flush) bind(c, name='log_info_f')
  use, intrinsic :: iso_c_binding
  implicit none
  character(kind=c_char), dimension(*) :: msg
  integer(c_int32_t), value :: newl
  integer(c_int32_t), value :: flush
end subroutine c_log_info

!-------------------------------------------------------------------------------

subroutine c_log_error(msg,newl,flush) bind(c, name='log_error_f')
  use, intrinsic :: iso_c_binding
  implicit none
  character(kind=c_char), dimension(*) :: msg
  integer(c_int32_t), value :: newl
  integer(c_int32_t), value :: flush
end subroutine c_log_error

!-------------------------------------------------------------------------------

subroutine c_log_warning(msg,newl,flush) bind(c, name='log_warning_f')
  use, intrinsic :: iso_c_binding
  implicit none
  character(kind=c_char), dimension(*) :: msg
  integer(c_int32_t), value :: newl
  integer(c_int32_t), value :: flush
end subroutine c_log_warning

!-------------------------------------------------------------------------------

subroutine c_log_debug(msg,newl,flush) bind(c, name='log_debug_f')
  use, intrinsic :: iso_c_binding
  implicit none
  character(kind=c_char), dimension(*) :: msg
  integer(c_int32_t), value :: newl
  integer(c_int32_t), value :: flush
end subroutine c_log_debug

!-------------------------------------------------------------------------------

subroutine c_log_trace(msg,newl,flush) bind(c, name='log_trace_f')
  use, intrinsic :: iso_c_binding
  implicit none
  character(kind=c_char), dimension(*) :: msg
  integer(c_int32_t), value :: newl
  integer(c_int32_t), value :: flush
end subroutine c_log_trace

!-------------------------------------------------------------------------------

subroutine c_log_stats(msg,newl,flush) bind(c, name='log_stats_f')
  use, intrinsic :: iso_c_binding
  implicit none
  character(kind=c_char), dimension(*) :: msg
  integer(c_int32_t), value :: newl
  integer(c_int32_t), value :: flush
end subroutine c_log_stats

!-------------------------------------------------------------------------------

subroutine c_log_test(msg,newl,flush) bind(c, name='log_test_f')
  use, intrinsic :: iso_c_binding
  implicit none
  character(kind=c_char), dimension(*) :: msg
  integer(c_int32_t), value :: newl
  integer(c_int32_t), value :: flush
end subroutine c_log_test

!-------------------------------------------------------------------------------
end interface
!-------------------------------------------------------------------------------
