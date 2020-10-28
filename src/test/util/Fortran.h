/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 *
 */

#ifndef TEST_UTIL_FORTRAN_H_
#define TEST_UTIL_FORTRAN_H_

#include <cstdint>
#include <string>
#include <vector>
#include "eckit/config/Configuration.h"

namespace test {

// -----------------------------------------------------------------------------
/*! This is intended as a general interface for testing objects in the **util**
 *  namespace from Fortran
 */
// -----------------------------------------------------------------------------

extern "C" {

  std::int32_t test_uniform_real_f(const eckit::Configuration * const *);
  std::int32_t test_uniform_double_f(const eckit::Configuration * const *);
  std::int32_t test_uniform_int_f(const eckit::Configuration * const *);
  std::int32_t test_uniform_long_f(const eckit::Configuration * const *);
  std::int32_t test_normal_real_f(const eckit::Configuration * const *);
  std::int32_t test_normal_double_f(const eckit::Configuration * const *);

  void test_push_string_vector_f(const eckit::Configuration * const *,
                                 std::vector<std::string> &);

}

#endif  // TEST_UTIL_FORTRAN_H_

}  // namespace test
