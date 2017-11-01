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
#include "lorenz95/LocsL95.h"
#include "lorenz95/ObsTable.h"
#include "util/DateTime.h"
#include "test/TestFixture.h"

namespace test {

// -----------------------------------------------------------------------------
class LocsTestFixture : TestFixture {
 public:
  LocsTestFixture() {
    const eckit::LocalConfiguration conf(TestConfig::config(), "Observations");
    const util::DateTime bgn(conf.getString("window_begin"));
    const util::DateTime end(conf.getString("window_end"));
    const eckit::LocalConfiguration otconf(conf, "Observation");
    ot_.reset(new lorenz95::ObsTable(otconf, bgn, end));
    t1_.reset(new util::DateTime("2010-01-01T12:00:00Z"));
    t2_.reset(new util::DateTime("2010-01-02T00:00:00Z"));
  }
  ~LocsTestFixture() {}
  boost::scoped_ptr<lorenz95::ObsTable> ot_;
  boost::scoped_ptr<util::DateTime> t1_;
  boost::scoped_ptr<util::DateTime> t2_;
};
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_SUITE(test_LocsL95, LocsTestFixture)
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_LocsL95_constructor) {
    boost::scoped_ptr<lorenz95::LocsL95> locs(new lorenz95::LocsL95(*ot_, *t1_, *t2_));
    BOOST_CHECK(locs.get() != NULL);
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_LocsL95_nobs) {
    boost::scoped_ptr<lorenz95::LocsL95> locs(new lorenz95::LocsL95(*ot_, *t1_, *t2_));
    BOOST_CHECK_EQUAL(locs->nobs(), 80);
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_LocsL95_operator) {
    boost::scoped_ptr<lorenz95::LocsL95> locs(new lorenz95::LocsL95(*ot_, *t1_, *t2_));
    double pos = 0.0;
    for (int i = 0; i < locs->nobs(); ++i) {
      BOOST_CHECK_CLOSE((*locs)[i], pos, 0.000001);
      pos += 0.05;
      if (pos >= 1.0) pos=0.0;
    }
  }
// -----------------------------------------------------------------------------
BOOST_AUTO_TEST_SUITE_END()
// -----------------------------------------------------------------------------

}  // namespace test
