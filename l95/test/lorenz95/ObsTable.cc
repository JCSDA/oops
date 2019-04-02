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

#include "./TestConfig.h"
#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "lorenz95/LocsL95.h"
#include "lorenz95/ObsTable.h"
#include "oops/runs/Test.h"
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
CASE("test_ObsTable") {
  ObsTableTestFixture fix;
// -----------------------------------------------------------------------------
  SECTION("test_ObsTable_constructor") {
    boost::scoped_ptr<lorenz95::ObsTable>
      ot(new lorenz95::ObsTable(*fix.testconf_, *fix.bgn_, *fix.end_));
    EXPECT(ot.get() != NULL);
  }
// -----------------------------------------------------------------------------
  SECTION("test_ObsTable_nobs") {
    boost::scoped_ptr<lorenz95::ObsTable>
      ot(new lorenz95::ObsTable(*fix.testconf_, *fix.bgn_, *fix.end_));
    const unsigned int nobs = 160;
    EXPECT(ot->nobs() == nobs);
  }
// -----------------------------------------------------------------------------
  SECTION("test_observationL95_put_get") {
    boost::scoped_ptr<lorenz95::ObsTable>
      ot(new lorenz95::ObsTable(*fix.testconf_, *fix.bgn_, *fix.end_));

    unsigned int nn = ot->nobs();
    std::vector<double> v1(nn);
    for (unsigned int jj = 0; jj < nn; ++jj) {
      v1[jj] = 3.14*jj;
    }
    ot->putdb("ObsTest", v1);

    std::vector<double> v2;
    ot->getdb("ObsTest", v2);

    EXPECT(v2.size() == nn);
    for (unsigned int jj = 0; jj < nn; ++jj) {
      EXPECT(v2[jj] == v1[jj]);
    }
  }
// -----------------------------------------------------------------------------
  SECTION("test_observationL95_timeSelect") {
    boost::scoped_ptr<lorenz95::ObsTable>
      ot(new lorenz95::ObsTable(*fix.testconf_, *fix.bgn_, *fix.end_));
    util::DateTime t1("2010-01-01T09:00:00Z");
    util::DateTime t2("2010-01-01T21:00:00Z");
    std::vector<int> mask = ot->timeSelect(t1, t2);  // t1 not includede, t2 is
    const unsigned int size = 80;
    EXPECT(mask.size() == size);
    EXPECT(mask[0] == 40);
    EXPECT(mask[79] == 119);
  }
// -----------------------------------------------------------------------------
  SECTION("test_observationL95_locations") {
    boost::scoped_ptr<lorenz95::ObsTable>
      ot(new lorenz95::ObsTable(*fix.testconf_, *fix.bgn_, *fix.end_));
    util::DateTime t1("2010-01-01T09:00:00Z");
    util::DateTime t2("2010-01-01T21:00:00Z");
    lorenz95::LocsL95 * locs = ot->locations(t1, t2);
    const size_t size = 80;
    EXPECT(locs->size() == size);
  }
// -----------------------------------------------------------------------------
  SECTION("test_observationL95_distribute") {
    eckit::LocalConfiguration otconf;
    util::DateTime t1("2010-01-01T00:00:00Z");
    util::DateTime t2("2010-01-02T23:59:59Z");
    boost::scoped_ptr<lorenz95::ObsTable>
      ot(new lorenz95::ObsTable(otconf, t1, t2));

    // More complete test in makeobs* tests.
    eckit::LocalConfiguration genconf(*fix.obsconf_, "Generate");

    ot->generateDistribution(genconf);
  }
// -----------------------------------------------------------------------------
}  //  CASE
// -----------------------------------------------------------------------------
}  // namespace test
int main(int argc, char **argv)
{
    return eckit::testing::run_tests ( argc, argv );
}
