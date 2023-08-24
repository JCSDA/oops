/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Copyright 2017-2022 UCAR.
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

#include "oops/base/ForecastParameters.h"
#include "oops/base/Geometry.h"
#include "oops/base/LatLonGridPostProcessor.h"
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
template <typename MODEL> class ForecastAppParameters : public ApplicationParameters {
  OOPS_CONCRETE_PARAMETERS(ForecastAppParameters, ApplicationParameters);

 public:
  typedef ForecastParameters<MODEL> ForecastParameters_;

  /// Forecast parameters.
  ForecastParameters_ fcstConf{this};
};

// -----------------------------------------------------------------------------

/// Application that runs a forecast from a model and initial condition
template <typename MODEL> class Forecast : public Application {
  typedef Geometry<MODEL>              Geometry_;
  typedef Model<MODEL>                 Model_;
  typedef ModelAuxControl<MODEL>       ModelAux_;
  typedef State<MODEL>                 State_;
  typedef ForecastAppParameters<MODEL> ForecastAppParameters_;

 public:
// -----------------------------------------------------------------------------
  explicit Forecast(const eckit::mpi::Comm & comm = oops::mpi::world()) : Application(comm) {}
// -----------------------------------------------------------------------------
  virtual ~Forecast() {}
// -----------------------------------------------------------------------------
  int execute(const eckit::Configuration & fullConfig, bool validate) const override {
//  Deserialize parameters
    ForecastAppParameters_ params;
    if (validate) params.validate(fullConfig);
    params.deserialize(fullConfig);

//  Setup resolution
    const Geometry_ resol(params.fcstConf.geometry, this->getComm());

//  Setup Model
    const Model_ model(resol, eckit::LocalConfiguration(fullConfig, "model"));

//  Setup initial state
    State_ xx(resol, params.fcstConf.initialCondition);
    Log::test() << "Initial state: " << xx << std::endl;

//  Setup augmented state
    const ModelAux_ moderr(resol, params.fcstConf.modelAuxControl);

//  Setup times
    const util::Duration fclength = params.fcstConf.forecastLength;
    const util::DateTime bgndate(xx.validTime());
    const util::DateTime enddate(bgndate + fclength);
    Log::info() << "Running forecast from " << bgndate << " to " << enddate << std::endl;

//  Setup forecast outputs
    PostProcessor<State_> post;
    post.enrollProcessor(new StateInfo<State_>("fc", params.fcstConf.prints));
//    params.output.date = bgndate;     DATE SHOULD BE SET HERE, NOT IN YAML
    post.enrollProcessor(new StateWriter<State_>(params.fcstConf.output));
    if (params.fcstConf.latlonGridOutput.value() != boost::none) {
      post.enrollProcessor(new LatLonGridPostProcessor<MODEL, State_>(
            params.fcstConf.latlonGridOutput.value().value(), resol));
    }

//  Run forecast
    model.forecast(xx, moderr, fclength, post);

    Log::test() << "Final state: " << xx << std::endl;

    return 0;
  }
// -----------------------------------------------------------------------------
  void outputSchema(const std::string & outputPath) const override {
    ForecastAppParameters_ params;
    params.outputSchema(outputPath);
  }
// -----------------------------------------------------------------------------
  void validateConfig(const eckit::Configuration & fullConfig) const override {
    ForecastAppParameters_ params;
    params.validate(fullConfig);
  }
// -----------------------------------------------------------------------------
 private:
  std::string appname() const override {
    return "oops::Forecast<" + MODEL::name() + ">";
  }
// -----------------------------------------------------------------------------
};

}  // namespace oops
#endif  // OOPS_RUNS_FORECAST_H_
