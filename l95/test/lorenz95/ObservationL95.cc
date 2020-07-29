/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <memory>

#include "./TestConfig.h"
#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "lorenz95/ObservationL95.h"
#include "lorenz95/ObsTableView.h"
#include "oops/parallel/mpi/mpi.h"
#include "oops/runs/Test.h"
#include "test/TestFixture.h"

namespace test {

// -----------------------------------------------------------------------------
class ObsTestFixture : TestFixture {
 public:
  ObsTestFixture() {
    const eckit::LocalConfiguration conf(TestConfig::config(), "observations");
    const util::DateTime bgn(conf.getString("window begin"));
    const util::DateTime end(conf.getString("window end"));
    const eckit::LocalConfiguration otconf(conf, "observation");
    ot_.reset(new lorenz95::ObsTableView(otconf, oops::mpi::comm(), bgn, end));
  }
  ~ObsTestFixture() {}
  std::unique_ptr<lorenz95::ObsTableView> ot_;
};
// -----------------------------------------------------------------------------
CASE("test_ObsL95") {
  ObsTestFixture fix;
// -----------------------------------------------------------------------------
  SECTION("test_ObsL95_constructor") {
    std::unique_ptr<lorenz95::ObservationL95>
      obs(new lorenz95::ObservationL95(*fix.ot_, TestConfig::config()));
    EXPECT(obs.get() != NULL);
  }
// -----------------------------------------------------------------------------
  SECTION("test_observationL95_classname") {
    std::unique_ptr<lorenz95::ObservationL95>
      obs(new lorenz95::ObservationL95(*fix.ot_, TestConfig::config()));
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
