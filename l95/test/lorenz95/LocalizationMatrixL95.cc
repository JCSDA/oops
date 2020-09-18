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
#include "oops/util/parameters/IgnoreOtherParameters.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "test/TestFixture.h"

namespace test {
// -----------------------------------------------------------------------------
class LocalizationMatrixTestParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(LocalizationMatrixTestParameters, Parameters)

 public:
  oops::RequiredParameter<lorenz95::ResolutionParameters> resol{"geometry", this};
  oops::RequiredParameter<eckit::LocalConfiguration> backgroundError{"background error", this};
  /// \brief Don't treat the presence of other parameter groups as an error (this makes it
  /// possible to reuse a single YAML file in tests of implementations of multiple oops interfaces).
  oops::IgnoreOtherParameters ignoreOthers{this};
};
// -----------------------------------------------------------------------------
class LocalizationMatrixFixture : TestFixture {
 public:
  LocalizationMatrixFixture() {
    LocalizationMatrixTestParameters parameters;
    parameters.validateAndDeserialize(TestConfig::config());
    resol_.reset(new lorenz95::Resolution(parameters.resol, oops::mpi::world()));

    cfg_.reset(new eckit::LocalConfiguration(parameters.backgroundError));
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
