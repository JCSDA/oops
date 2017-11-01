/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <boost/scoped_ptr.hpp>
#include <boost/test/unit_test.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "./TestConfig.h"
#include "lorenz95/GomL95.h"
#include "lorenz95/NoVariables.h"
#include "lorenz95/ObservationL95.h"
#include "lorenz95/ObsTable.h"
#include "lorenz95/Resolution.h"
#include "util/DateTime.h"
#include "test/TestFixture.h"

namespace test {

// -----------------------------------------------------------------------------
class GomTestFixture : TestFixture {
 public:
  GomTestFixture() {
    const eckit::LocalConfiguration conf(TestConfig::config(), "Observations");
    const util::DateTime bgn(conf.getString("window_begin"));
    const util::DateTime end(conf.getString("window_end"));
    const eckit::LocalConfiguration otconf(conf, "Observation");
    ot_.reset(new lorenz95::ObsTable(otconf, bgn, end));
    t1_.reset(new util::DateTime("2010-01-01T03:00:00Z"));
    t2_.reset(new util::DateTime("2010-01-02T06:00:00Z"));
    novar_.reset(new lorenz95::NoVariables());
  }
  ~GomTestFixture() {}
  lorenz95::Resolution * res_;  // Not actually used
  boost::scoped_ptr<lorenz95::ObsTable> ot_;
  boost::scoped_ptr<util::DateTime> t1_;
  boost::scoped_ptr<util::DateTime> t2_;
  boost::scoped_ptr<lorenz95::NoVariables> novar_;
};
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_SUITE(test_GomL95, GomTestFixture)
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_GomL95_constructor) {
    boost::scoped_ptr<lorenz95::GomL95> gom(new lorenz95::GomL95(*ot_, *novar_, *t1_, *t2_, *res_));
    BOOST_CHECK(gom.get() != NULL);
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_GomL95_nobs) {
    boost::scoped_ptr<lorenz95::GomL95> gom(new lorenz95::GomL95(*ot_, *novar_, *t1_, *t2_, *res_));
    BOOST_CHECK_EQUAL(gom->nobs(), 160);
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_gomL95_classname) {
    boost::scoped_ptr<lorenz95::GomL95> gom(new lorenz95::GomL95(*ot_, *novar_, *t1_, *t2_, *res_));
    BOOST_CHECK_EQUAL(gom->classname(), "lorenz95::GomL95");
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_gomL95_zero) {
    boost::scoped_ptr<lorenz95::GomL95> gom(new lorenz95::GomL95(*ot_, *novar_, *t1_, *t2_, *res_));
    gom->zero();
    for (int i = 0; i < gom->nobs(); ++i) {
      BOOST_CHECK_EQUAL((*gom)[i], 0.0);
    }
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_gomL95_dot_product_with) {
    boost::scoped_ptr<lorenz95::GomL95> gom1(new lorenz95::GomL95(*ot_, *novar_, *t1_, *t2_, *res_));
    gom1->zero();
    boost::scoped_ptr<lorenz95::GomL95> gom2(new lorenz95::GomL95(*ot_, *novar_, *t1_, *t2_, *res_));
    gom2->zero();

    double zz = gom1->dot_product_with(*gom2);
    BOOST_CHECK_EQUAL(zz, 0.0);
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_gomL95_operator) {
    boost::scoped_ptr<lorenz95::GomL95> gom1(new lorenz95::GomL95(*ot_, *novar_, *t1_, *t2_, *res_));
    gom1->zero();
    boost::scoped_ptr<lorenz95::GomL95> gom2(new lorenz95::GomL95(*ot_, *novar_, *t1_, *t2_, *res_));
    gom2->zero();

    (*gom1)[1] = 1.0;
    (*gom2)[3] = 1.0;
    double zz = gom1->dot_product_with(*gom2);
    BOOST_CHECK_EQUAL(zz, 0.0);

    (*gom1)[3] = 1.0;
    zz = gom1->dot_product_with(*gom2);
    BOOST_CHECK_EQUAL(zz, 1.0);
  }
// -----------------------------------------------------------------------------
BOOST_AUTO_TEST_SUITE_END()
// -----------------------------------------------------------------------------

}  // namespace test
