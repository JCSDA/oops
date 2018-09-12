/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef TEST_TESTFIXTURE_H_
#define TEST_TESTFIXTURE_H_

#pragma once

#include "eckit/runtime/Main.h"
#include "oops/util/LibOOPS.h"

namespace test {

// -----------------------------------------------------------------------------

struct TestFixture {
  TestFixture() {
    // This method of initializing boost is currently only used for the Lorentz model
    // (as of August, 2018)
    eckit::Main::initialise(boost::unit_test::framework::master_test_suite().argc,
                            boost::unit_test::framework::master_test_suite().argv);
  }

  ~TestFixture() {
    oops::LibOOPS::instance().finalise();
  }
};

// -----------------------------------------------------------------------------

}  // namespace test

#endif  // TEST_TESTFIXTURE_H_
