/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 *
 */

#ifndef TEST_BASE_FORTRAN_H_
#define TEST_BASE_FORTRAN_H_

namespace eckit {
  class Configuration;
}

namespace oops {
  class Variables;
}

namespace test {

// -----------------------------------------------------------------------------
/*! This is intended as a general interface for testing objects in the **util**
 *  namespace from Fortran
 */
// -----------------------------------------------------------------------------

extern "C" {

  void test_vars_interface_f(const eckit::Configuration &, oops::Variables &);

}

#endif  // TEST_BASE_FORTRAN_H_

}  // namespace test
