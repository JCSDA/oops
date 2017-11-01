/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE QG_TEST_OP_OBS

#include <iostream>
#include <iomanip>
#include <cmath>

#include <boost/test/unit_test.hpp>
#include "oops/runs/Run.h"
#include "model/QgTraits.h"
#include "model/instantiateCovarFactory.h"
#include "model/instantiateQgObsFactory.h"

#include "test/base/TestSuiteOpObsFixture.h"

namespace qg {

/*!
 *  TestEnv initializes the test suite for the QG model as an OOPS application
 */
class TestEnv: public oops::Run {
  public:
    TestEnv() :
        oops::Run(
            (const int) boost::unit_test::framework::master_test_suite().argc,
            (const char**) boost::unit_test::framework::master_test_suite().argv) {
      setup(config());
    }

    virtual ~TestEnv() {
    }

  private:
    // Initialize datas for the whole test suite
    void setup(const eckit::Configuration & fullConfig) {
      instantiateQgObsFactory();
      instantiateCovarFactory();
      test::TestSuiteOpObsFixture<qg::QgTraits>::getInstance().setup(
          fullConfig);
    }
};
}

// -----------------------------------------------------------------------------
typedef qg::TestEnv TestEnv_;
BOOST_GLOBAL_FIXTURE(TestEnv_);

/*!
 *  Run the generic tests for the QG model
 */

// -----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_SUITE(tl_ad, test::TestSuiteOpObsFixture<qg::QgTraits>)
#include "test/base/TestSuiteOpObsTLAD.h"
BOOST_AUTO_TEST_SUITE_END()
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_SUITE(tl, test::TestSuiteOpObsFixture<qg::QgTraits>)
#include "test/base/TestSuiteOpObsTL.h"
BOOST_AUTO_TEST_SUITE_END()
// -----------------------------------------------------------------------------
