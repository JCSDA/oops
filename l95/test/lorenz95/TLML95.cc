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

#include <boost/scoped_ptr.hpp>

#include "./TestConfig.h"
#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "lorenz95/ModelBias.h"
#include "lorenz95/ModelBiasCorrection.h"
#include "lorenz95/ModelBiasCovariance.h"
#include "lorenz95/ModelL95.h"
#include "lorenz95/ModelTrajectory.h"
#include "lorenz95/Resolution.h"
#include "lorenz95/StateL95.h"
#include "lorenz95/TLML95.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"
#include "test/TestFixture.h"

namespace test {

// -----------------------------------------------------------------------------
class TlmTestFixture : TestFixture {
 public:
  TlmTestFixture() {
    eckit::LocalConfiguration res(TestConfig::config(), "resolution");
    resol_.reset(new lorenz95::Resolution(res));
    eckit::LocalConfiguration mod(TestConfig::config(), "model");
    model_.reset(new lorenz95::ModelL95(*resol_, mod));
    tlconf_.reset(new eckit::LocalConfiguration(TestConfig::config(), "linearmodel"));
  }
  ~TlmTestFixture() {}
  boost::scoped_ptr<lorenz95::ModelL95>   model_;
  boost::scoped_ptr<lorenz95::Resolution> resol_;
  boost::scoped_ptr<const eckit::LocalConfiguration> tlconf_;
};
// -----------------------------------------------------------------------------
CASE("test_tlmL95") {
  TlmTestFixture fix;
// -----------------------------------------------------------------------------
  SECTION("test_tlmL95_constructor") {
    boost::scoped_ptr<lorenz95::TLML95> tlm(new lorenz95::TLML95(*fix.resol_, *fix.tlconf_));
    EXPECT(tlm.get() != NULL);
  }
// -----------------------------------------------------------------------------
/*
  SECTION("test_tlmL95_set_get_Trajectory") {
    lorenz95::TLML95 tlm(*fix.resol_, *fix.tlconf_, *fix.model_);

    // construct a StateL95 object
    std::string date_string("2014-09-12T09:35:00Z");
    util::DateTime dt(date_string);
    oops::Variables vars();
    boost::scoped_ptr<lorenz95::StateL95> stateL95(new lorenz95::StateL95(*fix.resol_, vars, dt));

    // construct a ModelBias object
    boost::shared_ptr<const eckit::LocalConfiguration> biasCfg(new eckit::LocalConfiguration(TestConfig::config(), "ModelBias"));
    boost::scoped_ptr<lorenz95::ModelBias> modelBias(new lorenz95::ModelBias(*biasCfg, *fix.resol_));

    tlm.setTrajectory(*stateL95, *modelBias);

    boost::scoped_ptr<lorenz95::ModelTrajectory> modelTraj(new lorenz95::ModelTrajectory(true));
    oops::Log::error() << "PMC: getTraj result <" << tlm.getTrajectory(dt)->get(1) << ">" << std::endl;
    modelTraj->set(tlm.getTrajectory(dt)->get(1));
    //oops::Log::error() << "PMC: modelTraj 1st " << modelTraj->get(1)
    oops::Log::error() << "PMC: modelTraj resol  <" << modelTraj->get(1).resol() << ">" << std::endl;
    for(int i = 0; i < modelTraj->get(1).resol(); ++i) {
      oops::Log::error() << "PMC: modelTraj FieldL95 elements " << modelTraj->get(1)[i] << std::endl;
    }
  }
*/
// -----------------------------------------------------------------------------
/*
  SECTION("test_tlmL95_stepTL") {
    lorenz95::TLML95 tlm(*fix.resol_, *fix.tlconf_, *fix.model_);

    // construct a FieldL95 object
    boost::scoped_ptr<lorenz95::FieldL95> fieldL95(new lorenz95::FieldL95(*fix.resol_));

    // construct a datetime object
    std::string dateString("2014-10-08T11:25:45Z");
    util::DateTime dt(dateString);

    //construct a ModelBiasCorrection object
    boost::scoped_ptr<const eckit::LocalConfiguration> covarCfg(new eckit::LocalConfiguration(TestConfig::config(), "Covariance"));
    boost::scoped_ptr<lorenz95::ModelBiasCovariance> modelBiasCovariance(
        new lorenz95::ModelBiasCovariance(*covarCfg, *fix.resol_));
    boost::scoped_ptr<lorenz95::ModelBiasCorrection> modelBiasCorrection(
        new lorenz95::ModelBiasCorrection(*modelBiasCovariance));

    tlm.stepTL(*fieldL95, dt, *modelBiasCorrection);

    oops::Log::error() << "PMC: original dt was " << dateString << std::endl;
    oops::Log::error() << "PMC: modified dt is  " << dt.toString() << std::endl;
  }
// -----------------------------------------------------------------------------
  SECTION("test_tlmL95_stepAD") {
  }
*/
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
    eckit::LocalConfiguration rescf(TestConfig::config(), "resolution");
    lorenz95::Resolution resol(rescf);
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
