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
#include <memory>


#include "./TestConfig.h"
#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "lorenz95/FieldL95.h"
#include "lorenz95/ModelBias.h"
#include "lorenz95/ModelL95.h"
#include "lorenz95/ModelTrajectory.h"
#include "lorenz95/Resolution.h"
#include "oops/mpi/mpi.h"
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
    eckit::LocalConfiguration res(TestConfig::config(), "geometry");
    resol_.reset(new lorenz95::Resolution(res, oops::mpi::world()));

    eckit::LocalConfiguration model(TestConfig::config(), "model");
    nlparams_.validateAndDeserialize(model);
  }
  ~ModelTestFixture() {}
  std::unique_ptr<lorenz95::Resolution> resol_;
  lorenz95::ModelL95Parameters nlparams_;
};
// -----------------------------------------------------------------------------
CASE("test_modelL95") {
ModelTestFixture fix;
// -----------------------------------------------------------------------------
  SECTION("test_modelL95_constructor") {
    std::unique_ptr<lorenz95::ModelL95> model(new lorenz95::ModelL95(*fix.resol_, fix.nlparams_));
    EXPECT(model != NULL);
  }
// -----------------------------------------------------------------------------
  SECTION("test_modelL95_get_classname") {
    lorenz95::ModelL95 model(*fix.resol_, fix.nlparams_);
    EXPECT(model.classname() == "lorenz95::ModelL95");
  }
// -----------------------------------------------------------------------------
  SECTION("test_modelL95_get_timestep") {
    lorenz95::ModelL95 model(*fix.resol_, fix.nlparams_);
    util::Duration dt(fix.nlparams_.tstep);
    EXPECT(model.timeResolution().toSeconds() == dt.toSeconds());
  }
// -----------------------------------------------------------------------------
  SECTION("test_modelL95_stepRk") {
    lorenz95::ModelL95 model(*fix.resol_, fix.nlparams_);

    // construct a FieldL95 object
    lorenz95::FieldL95 fieldL95(*fix.resol_);

    // construct a ModelBias object
    eckit::LocalConfiguration biasCfg(TestConfig::config(), "model aux control");
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
