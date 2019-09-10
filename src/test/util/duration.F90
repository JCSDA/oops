! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

!> Test interface for C++ config code called from Fortran

#include <fckit/fctest.h>

module DurationTestFixture
  use, intrinsic :: iso_c_binding
  
  implicit none
  
  type(c_ptr) :: ptr = c_null_ptr
  
  include 'oops/util/duration.intfb.h'
end module

TESTSUITE_WITH_FIXTURE(duration_intfb, DurationTestFixture)

TESTSUITE_INIT
END_TESTSUITE_INIT

TESTSUITE_FINALIZE
END_TESTSUITE_FINALIZE

!> Test the c_duration_int_str
! Pass in the total number of seconds that equates to 4 days, 3 hours, 2 minutes
! and 1 second and check that the string returned is correct.
TEST(test_c_duration_int_str)
  implicit none
  character(kind=c_char, len=21) :: returned_string, dur_string
  integer(kind=8) time_in_secs
  dur_string="P4DT3H2M1S"//c_null_char
  returned_string = c_null_char
  time_in_secs=356521
  call c_duration_int_str(time_in_secs, returned_string)
  CHECK_EQUAL(dur_string, returned_string)
END_TEST

!> Test the c_duration_int method
TEST(test_c_duration_int)
  implicit none
  character(kind=c_char, len=21) :: dur_string
  integer(kind=8) returned_time
  dur_string="P4DT3H2M1S"//c_null_char
  ptr=c_duration_construct(dur_string)
  returned_time=c_duration_int(ptr)
  CHECK(356521 == returned_time)
  call c_duration_destruct(ptr)
END_TEST

!> Test the c_duration_str_int method.
!! Create a valid ISO8601 duration in a null-terminated string and pass this 
!! into the method, check that it returns the correct number of seconds.
TEST(test_c_duration_str_int_1)
  implicit none
  character(kind=c_char, len=21) :: dur_string
  integer(kind=8) time_in_secs, returned_time
  dur_string="P4DT3H2M1S"//c_null_char
  time_in_secs=356521
  returned_time=c_duration_str_int(dur_string)
  CHECK(time_in_secs == returned_time)
END_TEST

!> Test the c_duration_str_int method.
!! Create a valid ISO8601 duration in a null-terminated string and pass this 
!! into the method, check that it returns the correct number of seconds.
TEST(test_c_duration_str_int_2)
  implicit none
  character(kind=c_char, len=21) :: dur_string
  integer(kind=8) time_in_secs, returned_time
  dur_string="P4DT3H121S"//c_null_char
  time_in_secs=356521
  returned_time=c_duration_str_int(dur_string)
  CHECK(time_in_secs == returned_time)
END_TEST

!> Test the c_duration_str_int method.
!! Create a valid ISO8601 duration in a null-terminated string and pass this 
!! into the method, check that it returns the correct number of seconds.
TEST(test_c_duration_str_int_3)
  implicit none
  character(kind=c_char, len=21) :: dur_string
  integer(kind=8) time_in_secs, returned_time
  dur_string="PT99H121S"//c_null_char
  time_in_secs=356521
  returned_time=c_duration_str_int(dur_string)
  CHECK(time_in_secs == returned_time)
END_TEST

!> Test the c_duration_str_int method.
!! Create a valid ISO8601 duration in a null-terminated string and pass this 
!! into the method, check that it returns the correct number of seconds.
!! This is one second less than all the other tests - we needed a whole number
!! of minutes for this one.
TEST(test_c_duration_str_int_4)
  implicit none
  character(kind=c_char, len=21) :: dur_string
  integer(kind=8) time_in_secs, returned_time
  dur_string="PT5942M"//c_null_char
  time_in_secs=356520
  returned_time=c_duration_str_int(dur_string)
  CHECK(time_in_secs == returned_time)
END_TEST

!> Test the c_duration_str_int method.
!! Create a valid ISO8601 duration in a null-terminated string and pass this 
!! into the method, check that it returns the correct number of seconds.
TEST(test_c_duration_str_int_5)
  implicit none
  character(kind=c_char, len=21) :: dur_string
  integer(kind=8) time_in_secs, returned_time
  dur_string="PT356521S"//c_null_char
  time_in_secs=356521
  returned_time=c_duration_str_int(dur_string)
  CHECK(time_in_secs == returned_time)
END_TEST

END_TESTSUITE
