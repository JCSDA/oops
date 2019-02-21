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

#include "Setup.h"
#include "oops/util/LibOOPS.h"

namespace test {

// -----------------------------------------------------------------------------
template <bool useBoost>
struct TestFixtureBase : public testing::Setup<useBoost> {
  TestFixtureBase() : testing::Setup<useBoost>() {
    // Common setup for every unit-test goes here
  }

  ~TestFixtureBase() {
    oops::LibOOPS::instance().finalise();
  }
};

#ifndef OOPS_TEST_NO_BOOST
using TestFixture = TestFixtureBase<true>;
#endif
// -----------------------------------------------------------------------------

}  // namespace test

#endif  // TEST_TESTFIXTURE_H_
