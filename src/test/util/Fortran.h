/*
 * (C) Copyright 2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 *
 */

#ifndef TEST_UTIL_F_H_
#define TEST_UTIL_F_H_

#include <cstdint>
#include "eckit/config/Configuration.h"

namespace test {

// -----------------------------------------------------------------------------
/*! This is intended as a general interface for testing objects in the **util**
 *  namespace from Fortran
 */
// -----------------------------------------------------------------------------

extern "C" {

  std::int32_t test_uniform_real_f(const eckit::Configuration * const *);

}

#endif  // TEST_UTIL_F_H_

} // namespace test 
