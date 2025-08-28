/*
 * (C) Copyright 2022 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_RUNS_ADJOINTFORECAST_H_
#define OOPS_RUNS_ADJOINTFORECAST_H_

#include <memory>
#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "oops/base/Model.h"
#include "oops/base/PostProcessor.h"
#include "oops/base/PostProcessorTLAD.h"
#include "oops/base/StateInfo.h"
#include "oops/base/StateWriter.h"
#include "oops/base/TrajectorySaver.h"
#include "oops/base/Variables.h"
#include "oops/generic/instantiateLinearModelFactory.h"
#include "oops/generic/instantiateModelFactory.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
#include "oops/interface/ModelAuxControl.h"
#include "oops/interface/ModelAuxIncrement.h"
#include "oops/interface/NormGradient.h"
#include "oops/interface/State.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

namespace oops {

// -----------------------------------------------------------------------------

template <typename MODEL> class AdjointForecast : public Application {
  typedef Geometry<MODEL>                  Geometry_;
  typedef Increment<MODEL>                 Increment_;
  typedef Model<MODEL>                     Model_;
  typedef ModelAuxControl<MODEL>           ModelAux_;
  typedef ModelAuxIncrement<MODEL>         ModelAuxIncr_;
  typedef NormGradient<MODEL>              NormGradient_;
  typedef State<MODEL>                     State_;
  typedef oops::LinearModel<MODEL>         LinearModel_;

 public:
// -----------------------------------------------------------------------------
  explicit AdjointForecast(const eckit::mpi::Comm & comm = oops::mpi::world()) : Application(comm) {
    oops::instantiateModelFactory<MODEL>();
    oops::instantiateLinearModelFactory<MODEL>();
  }
// -----------------------------------------------------------------------------
  virtual ~AdjointForecast() {}
// -----------------------------------------------------------------------------
  int execute(const eckit::Configuration & fullConfig) const override {
    eckit::LocalConfiguration tlConf(fullConfig, "linear forecast");
    eckit::LocalConfiguration adConf(fullConfig, "adjoint forecast");
    eckit::LocalConfiguration aspectConf(fullConfig, "forecast aspect");

    // Create the linear model
    const Geometry_ geomAD(eckit::LocalConfiguration(tlConf, "geometry"), this->getComm());
    oops::PostProcessor<State_> post;

    const eckit::LocalConfiguration fcConf(fullConfig, "forecast");
    if (fcConf.has("prints")) {
      const eckit::LocalConfiguration prtConfig(fcConf, "prints");
      post.enrollProcessor(new StateInfo<State_>("fc", prtConfig));
    }
    const eckit::LocalConfiguration outConfig(fcConf, "output");
    post.enrollProcessor(new StateWriter<State_>(outConfig));

    oops::PostProcessorTLAD<MODEL> pptraj;
    std::shared_ptr<LinearModel_> linearmodel_;
    linearmodel_.reset(new LinearModel_(geomAD, eckit::LocalConfiguration(tlConf, "linear model")));

    // Generate the model trajectory
    // -----------------------------

    // Setup resolution
    const Geometry_ geomNL(eckit::LocalConfiguration(fcConf, "geometry"), this->getComm());

    // Setup Model
    const Model_ model(geomNL, eckit::LocalConfiguration(fcConf, "model"));

    // Setup initial state
    State_ xxf(geomNL, eckit::LocalConfiguration(fcConf, "initial condition"));
    Log::test() << "Initial state: " << xxf << std::endl;

    // Setup augmented state
    const eckit::LocalConfiguration auxConf = fcConf.getSubConfiguration("model aux control");
    const ModelAux_ moderr(geomNL, auxConf);

    // Forecast length
    const util::Duration fclength(fcConf.getString("forecast length"));

    // Run forecast to get the trajectory
    post.enrollProcessor(new oops::TrajectorySaver<MODEL>
      (eckit::LocalConfiguration(tlConf, "linear model"), geomAD, moderr, linearmodel_, pptraj));
    model.forecast(xxf, moderr, fclength, post);
    Log::test() << "Forecast state: " << xxf << std::endl;

    // Run linear forecast
    // -------------------

    // Setup verification resolution
    const Geometry_ verifGeom(eckit::LocalConfiguration(aspectConf, "verification resolution"),
                              this->getComm());

    // Setup verification state
    State_ xxv(verifGeom, eckit::LocalConfiguration(aspectConf, "verification state"));
    Log::test() << "Verifying state: " << xxv << std::endl;

    // Set datetime for the increment
    util::DateTime incdate(xxf.validTime());

    // Increment variables
    const Variables vars(xxv.variables());

    // Create increment and diff two states for the increment
    Increment_ dx(verifGeom, vars, incdate);
    dx.diff(xxf, xxv);
    Log::test() << "Created perturbation from states diff: " << dx << std::endl;

    // Initialization type for increment
    NormGradient_ normgrad(xxf.geometry(), xxf, eckit::LocalConfiguration(aspectConf, "norm"));
    normgrad.apply(dx);

    // Write ADM initial conditions
    if (adConf.has("adjoint initial condition output")) {
      dx.write(eckit::LocalConfiguration(adConf, "adjoint initial condition output"));
    }

    // Setup augmented state for TLAD
    ModelAuxIncr_ admaux(verifGeom, auxConf);

    // Run ADM forecast
    linearmodel_->forecastAD(dx, admaux, fclength);

    // Write ADM final conditions
    dx.write(eckit::LocalConfiguration(adConf, "adjoint forecast output"));

    Log::test() << "Final increment state: " << dx << std::endl;

    return 0;
  }
// -----------------------------------------------------------------------------
 private:
  std::string appname() const override {
    return "oops::AdjointForecast<" + MODEL::name() + ">";
  }
// -----------------------------------------------------------------------------
};

}  // namespace oops
#endif  // OOPS_RUNS_ADJOINTFORECAST_H_
