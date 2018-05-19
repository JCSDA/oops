! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

!> Test interface for C++ config code called from Fortran

#include <fckit/fctest.h>

!-------------------------------------------------------------------------------
module ConfigTestFixture
  use, intrinsic :: iso_c_binding
  use kinds
  use config_mod
  implicit none
  include 'util/config.intfb.h'
end module
!-------------------------------------------------------------------------------

TESTSUITE_WITH_FIXTURE(config_intfb, ConfigTestFixture)

TESTSUITE_INIT
END_TESTSUITE_INIT

TESTSUITE_FINALIZE
END_TESTSUITE_FINALIZE

!> Test the c_config_construct and c_config_destruct methods
!! Calling the config_construct generates a structure that provides enough data
!! to test the functionality implemented within the interface.
TEST(test_c_config_construct_destruct)
  implicit none
  type(c_ptr) :: cptr
  cptr=c_config_get_testdata()
  CHECK(c_associated(cptr))
  CHECK(config_element_exists(cptr, "grandparent"))
  CHECK(config_element_exists(cptr, "grandparent.son"))
  CHECK(config_element_exists(cptr, "grandparent.son.maxdouble"))
  CHECK(config_element_exists(cptr, "grandparent.daughter"))
  CHECK(config_element_exists(cptr, "grandparent.daughter.grandson"))
  CHECK(.not.config_element_exists(cptr, "grandparent.daughter.hamsters"))
END_TEST

!> Test the c_config_element_exists method
!! Check that an existing path is found
TEST(test_c_config_element_exists_pass)
  implicit none
  type(c_ptr) :: ptr
  ptr=c_config_get_testdata()
  CHECK(config_element_exists(ptr, "grandparent.daughter.grandson"))
END_TEST

!> Test the c_config_element_exists method
!! Check that an non-existent path is not found
TEST(test_c_config_element_exists_fail)
  implicit none
  type(c_ptr) :: ptr
  ptr=c_config_get_testdata()
  CHECK(config_element_exists(ptr, "grandparent.daughter.granddaughter") .eqv. .false.)
END_TEST

!> Test the c_config_get_data_as_int method
!! Ensure that data stored as an integer can be retrieved
TEST(test_c_config_get_data_as_int)
  implicit none
  type(c_ptr) :: ptr
  integer(kind=4) :: my_int, return_int
  my_int=21
  ptr=c_config_get_testdata()
  return_int=config_get_int(ptr, "grandparent.daughter.grandson")
  CHECK_EQUAL(return_int, my_int)
END_TEST

!> Test the c_config_get_data_as_int method
!! Ensure that data stored as the largest possible integer can be retrieved
TEST(test_c_config_get_data_as_max_int)
  implicit none
  type(c_ptr) :: ptr
  integer(kind=4) :: max_int, return_int
  max_int=2147483647
  ptr=c_config_get_testdata()
  return_int=config_get_int(ptr, "grandparent.son.maxint")
  CHECK_EQUAL(return_int, max_int)
END_TEST

!> Test the c_config_get_data_as_int method
!! Ensure that data stored as the smallest possible integer can be
!! retrieved
TEST(test_c_config_get_data_as_min_int)
  implicit none
  type(c_ptr) :: ptr
  integer(kind=8) :: min_int, return_int
  min_int=-2147483647
  min_int=min_int-1
  ptr=c_config_get_testdata()
  return_int=config_get_int(ptr, "grandparent.son.minint")
  CHECK(return_int == min_int)
END_TEST

!> Test the c_config_get_data_as_double method
!! Ensure that data stored as a double can be retrieved
TEST(test_c_config_get_data_as_double)
  implicit none
  type(c_ptr) :: ptr
  real(kind_real) :: return_double
  ptr=c_config_get_testdata()
  return_double=config_get_real(ptr, "grandparent.son.granddaughter")
  CHECK_EQUAL(return_double, 3.14159_kind_real)
END_TEST

!> Test the c_config_get_data_as_double method
!! Ensure that data stored as a largest possible double can be retrieved
TEST(test_c_config_get_data_as_max_double)
  implicit none
  type(c_ptr) :: ptr
  double precision return_double
  ptr=c_config_get_testdata()
  return_double=config_get_real(ptr, "grandparent.son.maxdouble")
  CHECK(return_double == 1.79769d308)
END_TEST

!> Test the c_config_get_data_as_double method
!! Ensure that data stored as a smallest possible double can be retrieved
TEST(test_c_config_get_data_as_min_double)
  implicit none
  type(c_ptr) :: ptr
  double precision return_double
  ptr=c_config_get_testdata()
  return_double=config_get_real(ptr, "grandparent.son.mindouble")
  CHECK(return_double == -1.79769d308)
END_TEST

!> Test the c_config_get_data method
!! Ensure that data stored as a string can be retrieved
TEST(test_c_config_get_data)
  implicit none
  type(c_ptr) :: ptr
  character(12) :: rtnstr
  ptr=c_config_get_testdata()
  rtnstr = config_get_string(ptr, len(rtnstr), "grandparent.daughter.hamster")
  CHECK_EQUAL(trim(rtnstr), "Errol")
END_TEST

END_TESTSUITE
