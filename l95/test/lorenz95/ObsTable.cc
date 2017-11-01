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

#include "./TestConfig.h"
#include "eckit/config/LocalConfiguration.h"
#include "lorenz95/ObsTable.h"
#include "test/TestFixture.h"

namespace test {

// -----------------------------------------------------------------------------
class ObsTableTestFixture : public TestFixture {
 public:
  ObsTableTestFixture() {
    obsconf_.reset(new eckit::LocalConfiguration(TestConfig::config(), "Observations"));
    bgn_.reset(new util::DateTime(obsconf_->getString("window_begin")));
    end_.reset(new util::DateTime(obsconf_->getString("window_end")));
    testconf_.reset(new eckit::LocalConfiguration(*obsconf_, "Observation"));
  }
  ~ObsTableTestFixture() {}
  boost::scoped_ptr<const eckit::LocalConfiguration> obsconf_;
  boost::scoped_ptr<const eckit::LocalConfiguration> testconf_;
  boost::scoped_ptr<const util::DateTime> bgn_;
  boost::scoped_ptr<const util::DateTime> end_;
};
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_SUITE(test_ObsTable, ObsTableTestFixture)
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_ObsTable_constructor) {
    boost::scoped_ptr<lorenz95::ObsTable> ot(new lorenz95::ObsTable(*testconf_, *bgn_, *end_));
    BOOST_CHECK(ot.get() != NULL);
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_ObsTable_nobs) {
    boost::scoped_ptr<lorenz95::ObsTable> ot(new lorenz95::ObsTable(*testconf_, *bgn_, *end_));
    const unsigned int nobs = 160;
    BOOST_CHECK_EQUAL(ot->nobs(), nobs);
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_observationL95_put_get) {
    boost::scoped_ptr<lorenz95::ObsTable> ot(new lorenz95::ObsTable(*testconf_, *bgn_, *end_));

    unsigned int nn = ot->nobs();
    std::vector<double> v1(nn);
    for (unsigned int jj = 0; jj < nn; ++jj) {
      v1[jj] = 3.14*jj;
    }
    ot->putdb("ObsTest", v1);

    std::vector<double> v2;
    ot->getdb("ObsTest", v2);

    BOOST_CHECK_EQUAL(v2.size(), nn);
    for (unsigned int jj = 0; jj < nn; ++jj) {
      BOOST_CHECK_EQUAL(v2[jj], v1[jj]);
    }
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_observationL95_timeSelect) {
    boost::scoped_ptr<lorenz95::ObsTable> ot(new lorenz95::ObsTable(*testconf_, *bgn_, *end_));
    util::DateTime t1("2010-01-01T09:00:00Z");
    util::DateTime t2("2010-01-01T21:00:00Z");
    std::vector<int> mask = ot->timeSelect(t1, t2);  // t1 not includede, t2 is
    const unsigned int size = 80;
    BOOST_CHECK_EQUAL(mask.size(), size);
    BOOST_CHECK_EQUAL(mask[0], 40);
    BOOST_CHECK_EQUAL(mask[79], 119);
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_observationL95_locations) {
    boost::scoped_ptr<lorenz95::ObsTable> ot(new lorenz95::ObsTable(*testconf_, *bgn_, *end_));
    util::DateTime t1("2010-01-01T09:00:00Z");
    util::DateTime t2("2010-01-01T21:00:00Z");
    std::vector<double> locs = ot->locations(t1, t2);
    const unsigned int size = 80;
    BOOST_CHECK_EQUAL(locs.size(), size);
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_observationL95_distribute) {
    eckit::LocalConfiguration otconf;
    util::DateTime t1("2010-01-01T00:00:00Z");
    util::DateTime t2("2010-01-02T23:59:59Z");
    boost::scoped_ptr<lorenz95::ObsTable> ot(new lorenz95::ObsTable(otconf, t1, t2));

    // More complete test in makeobs* tests.
    eckit::LocalConfiguration genconf(*obsconf_, "Generate");

    ot->generateDistribution(genconf);
  }
// -----------------------------------------------------------------------------
BOOST_AUTO_TEST_SUITE_END()
// -----------------------------------------------------------------------------
}  // namespace test
