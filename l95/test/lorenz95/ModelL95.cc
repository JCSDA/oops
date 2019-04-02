/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <iostream>

#include <boost/scoped_ptr.hpp>

#include "./TestConfig.h"
#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "lorenz95/FieldL95.h"
#include "lorenz95/ModelBias.h"
#include "lorenz95/ModelL95.h"
#include "lorenz95/ModelTrajectory.h"
#include "lorenz95/Resolution.h"
#include "oops/runs/Test.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"
#include "test/TestFixture.h"

namespace test {

// -----------------------------------------------------------------------------
class ModelTestFixture : TestFixture {
 public:
  ModelTestFixture() {
    eckit::LocalConfiguration res(TestConfig::config(), "resolution");
    resol_.reset(new lorenz95::Resolution(res));
    nlconf_.reset(new eckit::LocalConfiguration(TestConfig::config(), "model"));
  }
  ~ModelTestFixture() {}
  boost::scoped_ptr<lorenz95::Resolution> resol_;
  boost::scoped_ptr<const eckit::LocalConfiguration> nlconf_;
};
// -----------------------------------------------------------------------------
CASE("test_modelL95") {
ModelTestFixture fix;
// -----------------------------------------------------------------------------
  SECTION("test_modelL95_constructor") {
    boost::scoped_ptr<lorenz95::ModelL95> model(new lorenz95::ModelL95(*fix.resol_, *fix.nlconf_));
    EXPECT(model != NULL);
  }
// -----------------------------------------------------------------------------
  SECTION("test_modelL95_get_classname") {
    lorenz95::ModelL95 model(*fix.resol_, *fix.nlconf_);
    EXPECT(model.classname() == "lorenz95::ModelL95");
  }
// -----------------------------------------------------------------------------
  SECTION("test_modelL95_get_timestep") {
    lorenz95::ModelL95 model(*fix.resol_, *fix.nlconf_);
    util::Duration dt(fix.nlconf_->getString("tstep"));
    EXPECT(model.timeResolution().toSeconds() == dt.toSeconds());
  }
// -----------------------------------------------------------------------------
  SECTION("test_modelL95_stepRk") {
    lorenz95::ModelL95 model(*fix.resol_, *fix.nlconf_);

    // construct a FieldL95 object
    lorenz95::FieldL95 fieldL95(*fix.resol_);

    // construct a ModelBias object
    eckit::LocalConfiguration biasCfg(TestConfig::config(), "ModelBias");
    lorenz95::ModelBias modelBias(*fix.resol_, biasCfg);

    // construct a ModelTrajectory object
    lorenz95::ModelTrajectory modelTraj(true);
    modelTraj.set(fieldL95);

//    for(int i = 0; i < fieldL95.resol(); ++i) {
//      oops::Log::error() << "PMC: pre  stepRK " << fieldL95[i] << std::endl;
//    }

    model.stepRK(fieldL95, modelBias, modelTraj);

//    for(int i = 0; i < fieldL95.resol(); ++i) {
//      oops::Log::error() << "PMC: post stepRK " << fieldL95[i] << std::endl;
//    }
  }
// -----------------------------------------------------------------------------
}  //  CASE
// -----------------------------------------------------------------------------
}  // namespace test
int main(int argc, char **argv)
{
    return eckit::testing::run_tests ( argc, argv );
}
