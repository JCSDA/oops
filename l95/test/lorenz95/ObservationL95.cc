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
#include "lorenz95/ObservationL95.h"
#include "lorenz95/ObsTable.h"
#include "test/TestFixture.h"

namespace test {

// -----------------------------------------------------------------------------
class ObsTestFixture : TestFixture {
 public:
  ObsTestFixture() {
    const eckit::LocalConfiguration conf(TestConfig::config(), "Observations");
    const util::DateTime bgn(conf.getString("window_begin"));
    const util::DateTime end(conf.getString("window_end"));
    const eckit::LocalConfiguration otconf(conf, "Observation");
    ot_.reset(new lorenz95::ObsTable(otconf, bgn, end));
  }
  ~ObsTestFixture() {}
  boost::scoped_ptr<lorenz95::ObsTable> ot_;
};
// -----------------------------------------------------------------------------
CASE("test_ObsL95") {
  ObsTestFixture f;
// -----------------------------------------------------------------------------
  SECTION("test_ObsL95_constructor") {
    boost::scoped_ptr<lorenz95::ObservationL95>
      obs(new lorenz95::ObservationL95(*f.ot_, TestConfig::config()));
    EXPECT(obs.get() != NULL);
  }
// -----------------------------------------------------------------------------
  SECTION("test_observationL95_classname") {
    boost::scoped_ptr<lorenz95::ObservationL95>
      obs(new lorenz95::ObservationL95(*f.ot_, TestConfig::config()));
    EXPECT(obs->classname() == "lorenz95::ObservationL95");
  }
// -----------------------------------------------------------------------------
}  //  CASE
// -----------------------------------------------------------------------------
}  // namespace test
int main(int argc, char **argv)
{
    return eckit::testing::run_tests ( argc, argv );
}
