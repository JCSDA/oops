! (C) Copyright 2009-2016 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

!> Test interface for C++ config code called from Fortran

#include <fckit/fctest.h>

module DateTimeTestFixture
  use, intrinsic :: iso_c_binding
  use kinds

  implicit none

  type(c_ptr) :: ptr1 = c_null_ptr
  type(c_ptr) :: ptr2 = c_null_ptr

  include 'oops/util/datetime.intfb.h'
end module

TESTSUITE_WITH_FIXTURE(datetime_intfb, DateTimeTestFixture)

TESTSUITE_INIT
END_TESTSUITE_INIT

TESTSUITE_FINALIZE
END_TESTSUITE_FINALIZE

!> Test the c_datetime_construct and the c_datetime_destruct methods.
!! Create a valid ISO8601 date in a null-terminated string and pass this into
!! the method, check that it returns a C pointer associated with the
!! constructed datetime.
!! No need to check invalid months, days or dates, this is already covered in
!! dateFunctions.cc
!!
!! Then check that the C pointer is deleted using the destructor.
TEST(test_c_datetime_construct_destruct)
  implicit none
  character(21) :: dtstring, returnstring
  dtstring="2014-08-11T12:34:56Z"//c_null_char
  ptr1 = c_datetime_construct(dtstring)
  CHECK(c_associated(ptr1)) ! check datetime is constructed correctly
  call c_datetime_destruct(ptr1)
!  CHECK(.not.c_associated(ptr1)) ! check datetime is destructed correctly
END_TEST

!> Test the c_datetime_string method
!! Constructs a new datetime pointer and then confirms that the method creates
!! the same string that the datetime was originally constructed from.
TEST(test_c_datetime_string)
  implicit none
  character(21) :: dtstring, returnstring
  dtstring="2014-08-11T12:34:56Z"//c_null_char
  ptr1 = c_datetime_construct(dtstring)
  call c_datetime_string(ptr1, returnstring)
  CHECK_EQUAL(dtstring, returnstring) !
  call c_datetime_destruct(ptr1)
END_TEST

!> Test the c_datetime_set method.
!! Constructs two new datetime pointers (ptr1 and ptr) from strings
!! (dtstring1 and dtstring2). Passes one of the strings (dtstring1) and one
!! of the pointers (ptr2) to c_datetime_set method. The pointer that is
!! returned (ptr2) is then converted to a string (returnstring) and compared
!! to the string that the pointer value should be set to.
TEST(test_c_datetime_set)
  implicit none
  character(21) :: dtstring1, dtstring2, returnstring
  dtstring1="2013-07-10T01:23:45Z"//c_null_char
  dtstring2="2014-08-11T12:34:56Z"//c_null_char
  ptr1 = c_datetime_construct(dtstring1)
  ptr2 = c_datetime_construct(dtstring2)
  call c_datetime_set(dtstring1, ptr2)
  call c_datetime_string(ptr2, returnstring)
  CHECK_EQUAL(dtstring1, returnstring)
  call c_datetime_destruct(ptr1)
  call c_datetime_destruct(ptr2)
END_TEST

!> Test the c_datetime_diff method.
!! Constructs two datetime pointers and compares the difference in seconds
!! between the two values.
TEST(test_c_datetime_diff)
  implicit none
  character(21) :: dtstring1, dtstring2
  integer timediff
  dtstring1="2014-08-11T12:34:46Z"//c_null_char
  dtstring2="2014-08-11T12:34:56Z"//c_null_char
  ptr1 = c_datetime_construct(dtstring1)
  ptr2 = c_datetime_construct(dtstring2)
  timediff = c_datetime_diff(ptr1, ptr2)
  CHECK_EQUAL(abs(timediff), 10)
  call c_datetime_destruct(ptr1)
  call c_datetime_destruct(ptr2)
END_TEST

!> Test the c_datetime_diff method
!! Constructs two datetime pointers and compares the difference in seconds
!! between the two values. Ensures the method still works if they are the same
!! value.
TEST(test_c_datetime_diff_same_date)
  implicit none
  character(21) :: dtstring1, dtstring2
  integer timediff
  dtstring1="2014-08-11T12:34:46Z"//c_null_char
  dtstring2="2014-08-11T12:34:46Z"//c_null_char
  ptr1 = c_datetime_construct(dtstring1)
  ptr2 = c_datetime_construct(dtstring2)
  timediff = c_datetime_diff(ptr1, ptr2)
  CHECK_EQUAL(abs(timediff), 0)
  call c_datetime_destruct(ptr1)
  call c_datetime_destruct(ptr2)
END_TEST

!> Test the c_datetime_diff method
!! Constructs two datetime pointers and compares the difference in seconds
!! between the two values. Ensures the method still works regardless of which
!! datetime value is entered first, either the later or earlier datetime.
TEST(test_c_datetime_diff_later_date_first)
  implicit none
  character(21) :: dtstring1, dtstring2
  integer timediff
  dtstring1="2014-08-11T12:34:56Z"//c_null_char
  dtstring2="2014-08-11T12:34:46Z"//c_null_char
  ptr1 = c_datetime_construct(dtstring1)
  ptr2 = c_datetime_construct(dtstring2)
  timediff = c_datetime_diff(ptr1, ptr2)
  CHECK_EQUAL(abs(timediff), 10)
  call c_datetime_destruct(ptr1)
  call c_datetime_destruct(ptr2)
END_TEST

!> Test the c_datetime_diff method
!! Constructs two datetime pointers and compares the difference in seconds
!! between the two values. Ensures the method still works when using a datetime
!! of one year difference.
TEST(test_c_datetime_diff_future)
  implicit none
  character(21) :: dtstring1, dtstring2
  integer timediff
  dtstring1="2014-08-11T12:34:56Z"//c_null_char
  dtstring2="2015-08-11T12:34:56Z"//c_null_char
  ptr1 = c_datetime_construct(dtstring1)
  ptr2 = c_datetime_construct(dtstring2)
  timediff = c_datetime_diff(ptr1, ptr2)
  CHECK_EQUAL(abs(timediff), 31536000)
  call c_datetime_destruct(ptr1)
  call c_datetime_destruct(ptr2)
END_TEST

!> Test the c_datetime_update.
!! Construct a datetime pointer and pass this and a time difference in seconds
!! to the c_datetime_update method. Take the returned pointer and convert to a
!! string. The returned string should be equal to the original string + time
!! difference.
TEST(test_c_datetime_update)
  implicit none
  character(21) :: dtstring, timediffstring, returnstring
  integer(kind=8) timediff
  dtstring="2014-08-11T12:34:56Z"//c_null_char
  timediffstring="2014-08-11T12:35:56Z"//c_null_char
  timediff=60
  ptr1 = c_datetime_construct(dtstring)
  call c_datetime_update(ptr1, timediff)
  call c_datetime_string(ptr1, returnstring)
  CHECK_EQUAL(returnstring, timediffstring)
  call c_datetime_destruct(ptr1)
END_TEST

!> Test the c_datetime_update.
!! Construct a datetime pointer and pass this and a time difference in seconds
!! to the c_datetime_update method. Take the returned pointer and convert to a
!! string. The returned string should be equal to the original string + time
!! difference. Ensure this method works even if a difference value of zero is
!! used.
TEST(test_c_datetime_update_zero)
  implicit none
  character(21) :: dtstring, timediffstring, returnstring
  integer(kind=8) timediff
  dtstring="2014-08-11T12:34:56Z"//c_null_char
  timediffstring="2014-08-11T12:34:56Z"//c_null_char
  timediff=0
  ptr1 = c_datetime_construct(dtstring)
  call c_datetime_update(ptr1, timediff)
  call c_datetime_string(ptr1, returnstring)
  CHECK_EQUAL(returnstring, timediffstring)
  call c_datetime_destruct(ptr1)
END_TEST

!> Test the c_datetime_update.
!! Construct a datetime pointer and pass this and a time difference in seconds
!! to the c_datetime_update method. Take the returned pointer and convert to a
!! string. The returned string should be equal to the original string + time
!! difference. Ensure this method works even if a negative difference value is
!! used.
TEST(test_c_datetime_update_negative)
  implicit none
  character(21) :: dtstring, timediffstring, returnstring
  integer(kind=8) timediff
  dtstring="2014-08-11T12:34:56Z"//c_null_char
  timediffstring="2014-08-11T12:33:56Z"//c_null_char
  timediff=-60
  ptr1 = c_datetime_construct(dtstring)
  call c_datetime_update(ptr1, timediff)
  call c_datetime_string(ptr1, returnstring)
  CHECK_EQUAL(returnstring, timediffstring)
  call c_datetime_destruct(ptr1)
END_TEST

!> Test datetime_to_ifs.
!! Conversion from DateTime to IFS format
TEST(test_datetime_to_ifs)
  use, intrinsic :: iso_c_binding
  use datetime_mod
  implicit none
  character(len=20) :: fstring
  type(datetime) :: fdt
  integer(kind=c_int) :: idate, isecs

  fstring="2016-08-02T12:34:56Z"
  call datetime_create(fstring, fdt)
  call datetime_to_ifs(fdt, idate, isecs)

  CHECK_EQUAL(idate, 20160802)
  CHECK_EQUAL(isecs, 45296)
END_TEST

!> Test datetime_from_ifs.
!! Conversion from IFS format to DateTime
TEST(test_datetime_from_ifs)
  use, intrinsic :: iso_c_binding
  use datetime_mod
  implicit none
  character(len=20) :: fstring
  type(datetime) :: fdt
  integer(kind=c_int) :: idate, isecs

  fstring="2016-08-02T12:34:56Z"
  call datetime_create(fstring, fdt)

  idate=20160702
  isecs=77696

  call datetime_from_ifs(fdt, idate, isecs)

  call datetime_to_string(fdt, fstring)
  CHECK_EQUAL(fstring, "2016-07-02T21:34:56Z")
END_TEST

!> Test datetime_to_yyyymmddhhmmss.
!! Extraction of individual date and time components from DateTime
TEST(test_datetime_to_yyyymmddhhmmss)
  use, intrinsic :: iso_c_binding
  use datetime_mod
  implicit none
  character(len=20) :: fstring
  type(datetime) :: fdt
  integer(kind=c_int) :: year, month, day, hour, minute, second

  fstring="2016-08-02T12:34:56Z"
  call datetime_create(fstring, fdt)
  call datetime_to_yyyymmddhhmmss(fdt, year, month, day, hour, minute, second)

  CHECK_EQUAL(year, 2016)
  CHECK_EQUAL(month, 8)
  CHECK_EQUAL(day, 2)
  CHECK_EQUAL(hour, 12)
  CHECK_EQUAL(minute, 34)
  CHECK_EQUAL(second, 56)
END_TEST

!> Test datetime_seconds_since_jan1.
TEST(test_datetime_seconds_since_jan1)
  use datetime_mod
  implicit none
  character(len=20) :: fstring
  type(datetime) :: fdt
  integer :: seconds_since_jan1

  fstring="2020-05-18T03:27:45Z"
  call datetime_create(fstring, fdt)
  seconds_since_jan1 = datetime_seconds_since_jan1(fdt)

  ! Below, the number 138 comes from the number of elapsed days before May 18th,
  ! which was the 139th day of 2020
  CHECK_EQUAL(seconds_since_jan1, 45 + 27 * 60 + 3 * 3600 + 138 * 86400)
END_TEST

END_TESTSUITE
