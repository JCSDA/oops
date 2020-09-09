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
#include "lorenz95/LocalizationMatrixL95.h"
#include "lorenz95/Resolution.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Test.h"
#include "test/TestFixture.h"

namespace test {
class LocalizationMatrixFixture : TestFixture {
 public:
  LocalizationMatrixFixture() {
    eckit::LocalConfiguration res(TestConfig::config(), "geometry");
    resol_.reset(new lorenz95::Resolution(res, oops::mpi::world()));

    eckit::LocalConfiguration cfg(TestConfig::config(), "background error");
    cfg_.reset(new eckit::LocalConfiguration(cfg));
  }
  ~LocalizationMatrixFixture() {}
  std::unique_ptr<lorenz95::Resolution> resol_;
  std::unique_ptr<const eckit::LocalConfiguration> cfg_;
};
// -----------------------------------------------------------------------------
CASE("test_localizationMatrixL95") {
  LocalizationMatrixFixture fix;
// -----------------------------------------------------------------------------
  SECTION("test_localizationMatrixL95_constructor") {
    std::unique_ptr<lorenz95::LocalizationMatrixL95> locmat(
        new lorenz95::LocalizationMatrixL95(*fix.resol_, *fix.cfg_));

    EXPECT(locmat.get() != NULL);
  }
// -----------------------------------------------------------------------------
}  //  CASE
// -----------------------------------------------------------------------------
}  // namespace test
int main(int argc, char **argv)
{
    return eckit::testing::run_tests ( argc, argv );
}
