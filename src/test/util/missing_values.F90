! (C) Copyright 2020 MEt Office UK

! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> \file Test the Fortran interface for the missingValue() C++ functions.

#include <fckit/fctest.h>

TESTSUITE(missing_values_f)

TESTSUITE_INIT
END_TESTSUITE_INIT

TESTSUITE_FINALIZE
END_TESTSUITE_FINALIZE

!> Test the single- and double-precision overloads of missing_value().
TEST(test_real_missing_values)
  use iso_c_binding
  use missing_values_mod
  implicit none
  real(c_float) :: missing_float
  real(c_double) :: missing_double

  missing_float = missing_value(missing_float)
  missing_double = missing_value(missing_double)

  CHECK(missing_float /= 0.0)
  CHECK(missing_double /= 0.0)
  CHECK(missing_float < -1e30)
  CHECK(missing_double < -1d30)
  CHECK(missing_float /= missing_double)
END_TEST

!> Test the int32_t and int64_t overloads of missing_value().
TEST(test_int_missing_values)
  use iso_c_binding
  use missing_values_mod
  implicit none
  integer(c_int32_t) :: missing_int32_t
  integer(c_int64_t) :: missing_int64_t

  missing_int32_t = missing_value(missing_int32_t)
  missing_int64_t = missing_value(missing_int64_t)

  CHECK(missing_int32_t /= 0)
  CHECK(missing_int64_t /= 0)
  CHECK(missing_int32_t < -1000000)
  CHECK(missing_int64_t < -1000000)
  CHECK(missing_int64_t < missing_int32_t)
END_TEST

END_TESTSUITE
