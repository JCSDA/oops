/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <fstream>
#include <iostream>
#include <memory>


#include "./TestConfig.h"
#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "lorenz95/ModelBias.h"
#include "lorenz95/ModelBiasCorrection.h"
#include "lorenz95/ModelL95.h"
#include "lorenz95/ModelTrajectory.h"
#include "lorenz95/Resolution.h"
#include "lorenz95/StateL95.h"
#include "lorenz95/TLML95.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Test.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/IgnoreOtherParameters.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "test/TestFixture.h"

namespace test {

// -----------------------------------------------------------------------------
class TlmTestParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(TlmTestParameters, Parameters)

 public:
  oops::RequiredParameter<lorenz95::ResolutionParameters> resol{"geometry", this};
  oops::RequiredParameter<lorenz95::ModelL95Parameters> model{"model", this};
  oops::RequiredParameter<eckit::LocalConfiguration> linearModel{"linear model", this};
  /// \brief Don't treat the presence of other parameter groups as an error (this makes it
  /// possible to reuse a single YAML file in tests of implementations of multiple oops interfaces).
  oops::IgnoreOtherParameters ignoreOthers{this};
};

// -----------------------------------------------------------------------------
class TlmTestFixture : TestFixture {
 public:
  TlmTestFixture() {
    TlmTestParameters parameters;
    parameters.validateAndDeserialize(TestConfig::config());

    resol_.reset(new lorenz95::Resolution(parameters.resol, oops::mpi::world()));
    model_.reset(new lorenz95::ModelL95(*resol_, parameters.model));

    tlconf_.reset(new eckit::LocalConfiguration(parameters.linearModel));
  }
  ~TlmTestFixture() {}

  std::unique_ptr<lorenz95::ModelL95>   model_;
  std::unique_ptr<lorenz95::Resolution> resol_;
  std::unique_ptr<const eckit::LocalConfiguration> tlconf_;
};
// -----------------------------------------------------------------------------
CASE("test_tlmL95") {
  TlmTestFixture fix;
// -----------------------------------------------------------------------------
  SECTION("test_tlmL95_constructor") {
    std::unique_ptr<lorenz95::TLML95> tlm(new lorenz95::TLML95(*fix.resol_, *fix.tlconf_));
    EXPECT(tlm.get() != NULL);
  }
// -----------------------------------------------------------------------------
  SECTION("test_tlmL95_get_classname") {
    lorenz95::TLML95 tlm(*fix.resol_, *fix.tlconf_);
    EXPECT(tlm.classname() == "lorenz95::TLML95");
  }
// -----------------------------------------------------------------------------
  SECTION("test_tlmL95_get_timestep") {
    lorenz95::TLML95 tlm(*fix.resol_, *fix.tlconf_);
    util::Duration dt(fix.tlconf_->getString("tstep"));
    EXPECT(tlm.timeResolution().toSeconds() == dt.toSeconds());
  }
// -----------------------------------------------------------------------------
  SECTION("test_tlmL95_get_resolution") {
    eckit::LocalConfiguration rescf(TestConfig::config(), "geometry");
    lorenz95::ResolutionParameters resolParameters;
    resolParameters.validateAndDeserialize(rescf);
    lorenz95::Resolution resol(resolParameters, oops::mpi::world());

    lorenz95::TLML95 tlm(*fix.resol_, *fix.tlconf_);
    EXPECT(tlm.resolution().npoints() == resol.npoints());
  }
// -----------------------------------------------------------------------------
  SECTION("test_tlmL95_stream_output") {
    lorenz95::TLML95 tlm(*fix.resol_, *fix.tlconf_);

    // use the operator<< method to write the value to a file
    std::filebuf fb;
    std::string filename("TLML95Test.txt");
    fb.open(filename.c_str(), std::ios::out);
    std::ostream os(&fb);
    os << tlm;
    fb.close();

    // then read the value that was written to the file
    std::string inputString;
    std::string inputStartsWith;
    std::string testString("TLML95: resol = ");
    int endPos = testString.size();
    std::ifstream inputFile(filename.c_str());
    if (inputFile.is_open()) {
      getline(inputFile, inputString);

      inputStartsWith = inputString.substr(0, endPos);

      EXPECT(inputStartsWith == testString);
    } else {
      // if we can't open the file then we can't
      // verify that the value was correctly written
      oops::Log::error() <<"operator<< functionality cannot be determined" << std::endl;
    }
    inputFile.close();
  }
// -----------------------------------------------------------------------------
}  //  CASE
// -----------------------------------------------------------------------------
}  // namespace test
int main(int argc, char **argv)
{
    return eckit::testing::run_tests ( argc, argv );
}
