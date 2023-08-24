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

#include "oops/base/ForecastParameters.h"
#include "oops/base/Model.h"
#include "oops/base/PostProcessor.h"
#include "oops/base/PostProcessorTLAD.h"
#include "oops/base/StateInfo.h"
#include "oops/base/StateWriter.h"
#include "oops/base/TrajectorySaver.h"
#include "oops/base/Variables.h"
#include "oops/generic/instantiateLinearModelFactory.h"
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
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

namespace oops {

// -----------------------------------------------------------------------------

/// \brief Top-level options taken by the AdjointForecast application.
template <typename MODEL>
class AdjointForecastParameters : public ApplicationParameters {
  OOPS_CONCRETE_PARAMETERS(AdjointForecastParameters, ApplicationParameters)

 public:
  typedef typename Geometry<MODEL>::Parameters_ GeometryParameters_;
  typedef ForecastParameters<MODEL> ForecastParameters_;
  typedef ForecastAspectParameters<MODEL> ForecastAspectParameters_;
  typedef ForecastTLParameters<MODEL> ForecastTLParameters_;
  typedef ForecastADParameters<MODEL> ForecastADParameters_;

  /// Forecast parameters.
  RequiredParameter<ForecastParameters_> fcstConf{"forecast", this};

  /// Forecast aspect parameters.
  RequiredParameter<ForecastAspectParameters_> fcstAspectConf{"forecast aspect", this};

  /// Linear forecast parameters.
  RequiredParameter<ForecastTLParameters_> linearFcstConf{"linear forecast", this};

  /// Adjoint forecast parameters.
  RequiredParameter<ForecastADParameters_> adjointForecast{"adjoint forecast", this};
};

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
  typedef AdjointForecastParameters<MODEL> AdjointForecastParameters_;

 public:
// -----------------------------------------------------------------------------
  explicit AdjointForecast(const eckit::mpi::Comm & comm = oops::mpi::world()) : Application(comm) {
  }
// -----------------------------------------------------------------------------
  virtual ~AdjointForecast() {}
// -----------------------------------------------------------------------------
  int execute(const eckit::Configuration & fullConfig, bool validate) const override {
    // Deserialize parameters
    AdjointForecastParameters_ params;
    if (validate) params.validate(fullConfig);
    params.deserialize(fullConfig);

    // Create the linear model
    const Geometry_ adjointForecastModelGeometry
      (params.linearFcstConf.value().geometry, this->getComm());
    oops::instantiateLinearModelFactory<MODEL>();
    oops::PostProcessor<State_> post;
    post.enrollProcessor(new StateInfo<State_>("fc", params.fcstConf.value().prints));
    post.enrollProcessor(new StateWriter<State_>(params.fcstConf.value().output));
    oops::PostProcessorTLAD<MODEL> pptraj;
    std::shared_ptr<LinearModel_> linearmodel_;
    linearmodel_.reset(
      new LinearModel_(adjointForecastModelGeometry, params.linearFcstConf.value().model));

    // Generate the model trajectory
    // -----------------------------

    // Setup resolution
    const Geometry_ fcstModelGeom(params.fcstConf.value().geometry, this->getComm());

    // Setup Model
    eckit::LocalConfiguration fconf(fullConfig, "forecast");
    const Model_ model(fcstModelGeom,  eckit::LocalConfiguration(fconf, "model"));

    // Setup initial state
    State_ xxf(fcstModelGeom, params.fcstConf.value().initialCondition);
    Log::test() << "Initial state: " << xxf << std::endl;

    // Setup augmented state
    const ModelAux_ moderr(fcstModelGeom, params.fcstConf.value().modelAuxControl);

    // Forecast length
    const util::Duration fclength = params.fcstConf.value().forecastLength;

    // Run forecast to get the trajectory
    post.enrollProcessor(new oops::TrajectorySaver<MODEL>
      (params.linearFcstConf.value().model, adjointForecastModelGeometry,
      moderr, linearmodel_, pptraj));
    model.forecast(xxf, moderr, fclength, post);
    Log::test() << "Forecast state: " << xxf << std::endl;

    // Run linear forecast
    // -------------------

    // Setup verification resolution
    const Geometry_ verificationGeom(params.fcstAspectConf.value().verificationResolution,
      this->getComm());

    // Setup verification state
    State_ xxv(verificationGeom, params.fcstAspectConf.value().verificationConfig);
    Log::test() << "Verifying state: " << xxv << std::endl;

    // Set datetime for the increment
    util::DateTime incdate(xxf.validTime());

    // Increment variables
    const Variables vars(xxv.variables());

    // Create increment and diff two states for the increment
    Increment_ dx(verificationGeom, vars, incdate);
    dx.diff(xxf, xxv);
    Log::test() << "Created perturbation from states diff: " << dx << std::endl;

    // Initialization type for increment
    NormGradient_ normgrad(xxf.geometry(), xxf,
      params.fcstAspectConf.value().normConfig);
    normgrad.apply(dx);

    // Write ADM initial conditions
    if (params.adjointForecast.value().adjointInitialConditionOutput.value() != boost::none) {
      dx.write(params.adjointForecast.value().adjointInitialConditionOutput.value().value());
    }

    // Setup augmented state for TLAD
    ModelAuxIncr_ admaux(verificationGeom, params.adjointForecast.value().modelAuxIncrement);

    // Run ADM forecast
    linearmodel_->forecastAD(dx, admaux, fclength);

    // Write ADM final conditions
    dx.write(params.adjointForecast.value().adjointForecastOutput);

    Log::test() << "Final increment state: " << dx << std::endl;

    return 0;
  }
// -----------------------------------------------------------------------------
  void outputSchema(const std::string & outputPath) const override {
    AdjointForecastParameters_ params;
    params.outputSchema(outputPath);
  }
// -----------------------------------------------------------------------------
  void validateConfig(const eckit::Configuration & fullConfig) const override {
    AdjointForecastParameters_ params;
    params.validate(fullConfig);
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
