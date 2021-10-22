/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Copyright 2017-2021 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_RUNS_FORECAST_H_
#define OOPS_RUNS_FORECAST_H_

#include <string>

#include "eckit/config/LocalConfiguration.h"
#include "oops/base/Geometry.h"
#include "oops/base/Model.h"
#include "oops/base/PostProcessor.h"
#include "oops/base/State.h"
#include "oops/base/StateInfo.h"
#include "oops/base/StateWriter.h"
#include "oops/interface/ModelAuxControl.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

namespace oops {

/// Options taken by the Forecast application.
template <typename MODEL> class ForecastParameters : public ApplicationParameters {
  OOPS_CONCRETE_PARAMETERS(ForecastParameters, ApplicationParameters);

 public:
  typedef typename Geometry<MODEL>::Parameters_                 GeometryParameters_;
  typedef ModelParametersWrapper<MODEL>                         ModelParameters_;
  typedef State<MODEL>                                          State_;
  typedef typename StateParametersND<MODEL>::StateParameters3D_ StateParameters_;
  typedef StateWriterParameters<State_>                         StateWriterParameters_;

  /// Geometry parameters.
  RequiredParameter<GeometryParameters_> geometry{"geometry", this};

  /// Model parameters.
  RequiredParameter<ModelParameters_> model{"model", this};

  /// Initial state parameters.
  RequiredParameter<StateParameters_> initialCondition{"initial condition", this};

  /// Augmented model state.
  Parameter<eckit::LocalConfiguration> modelAuxControl{
    "model aux control", eckit::LocalConfiguration(), this};

  /// Forecast length.
  RequiredParameter<util::Duration> forecastLength{"forecast length", this};

  /// Where to write the output.
  RequiredParameter<StateWriterParameters_> output{"output", this};

  /// Options passed to the object writing out forecast fields.
  Parameter<PostTimerParameters> prints{"prints", {}, this};
};

// -----------------------------------------------------------------------------

/// Application that runs a forecast from a model and initial condition
template <typename MODEL> class Forecast : public Application {
  typedef Geometry<MODEL>           Geometry_;
  typedef Model<MODEL>              Model_;
  typedef ModelAuxControl<MODEL>    ModelAux_;
  typedef State<MODEL>              State_;

  typedef ForecastParameters<MODEL> ForecastParameters_;

 public:
// -----------------------------------------------------------------------------
  explicit Forecast(const eckit::mpi::Comm & comm = oops::mpi::world()) : Application(comm) {}
// -----------------------------------------------------------------------------
  virtual ~Forecast() {}
// -----------------------------------------------------------------------------
  int execute(const eckit::Configuration & fullConfig) const {
//  Deserialize parameters
    ForecastParameters_ params;
    params.validateAndDeserialize(fullConfig);

//  Setup resolution
    const Geometry_ resol(params.geometry, this->getComm(), oops::mpi::myself());

//  Setup Model
    const Model_ model(resol, params.model.value().modelParameters);

//  Setup initial state
    State_ xx(resol, params.initialCondition);
    Log::test() << "Initial state: " << xx << std::endl;

//  Setup augmented state
    const ModelAux_ moderr(resol, params.modelAuxControl);

//  Setup times
    const util::Duration fclength = params.forecastLength;
    const util::DateTime bgndate(xx.validTime());
    const util::DateTime enddate(bgndate + fclength);
    Log::info() << "Running forecast from " << bgndate << " to " << enddate << std::endl;

//  Setup forecast outputs
    PostProcessor<State_> post;
    post.enrollProcessor(new StateInfo<State_>("fc", params.prints));
    post.enrollProcessor(new StateWriter<State_>(params.output));

//  Run forecast
    model.forecast(xx, moderr, fclength, post);

    Log::test() << "Final state: " << xx << std::endl;

    return 0;
  }
// -----------------------------------------------------------------------------
 private:
  std::string appname() const {
    return "oops::Forecast<" + MODEL::name() + ">";
  }
// -----------------------------------------------------------------------------
};

}  // namespace oops
#endif  // OOPS_RUNS_FORECAST_H_
