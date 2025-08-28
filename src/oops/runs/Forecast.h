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

#include "eckit/config/LocalConfiguration.h"
#include "oops/base/Geometry.h"
#include "oops/base/Model.h"
#include "oops/base/PostProcessor.h"
#include "oops/base/State.h"
#include "oops/base/StateInfo.h"
#include "oops/base/StateWriter.h"
#include "oops/base/StructuredGridPostProcessor.h"
#include "oops/generic/instantiateModelFactory.h"
#include "oops/interface/ModelAuxControl.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"
#include "oops/util/WorkflowUpdater.h"

namespace oops {

// -----------------------------------------------------------------------------

/// Application that runs a forecast from a model and initial condition
template <typename MODEL> class Forecast : public Application {
  typedef Geometry<MODEL>              Geometry_;
  typedef Model<MODEL>                 Model_;
  typedef ModelAuxControl<MODEL>       ModelAux_;
  typedef State<MODEL>                 State_;

 public:
// -----------------------------------------------------------------------------
  explicit Forecast(const eckit::mpi::Comm & comm = oops::mpi::world()) : Application(comm) {
    instantiateModelFactory<MODEL>();
  }
// -----------------------------------------------------------------------------
  virtual ~Forecast() {}
// -----------------------------------------------------------------------------
  int execute(const eckit::Configuration & fullConfig) const override {
//  Setup resolution
    const Geometry_ resol(eckit::LocalConfiguration(fullConfig, "geometry"), this->getComm());

//  Setup Model
    const Model_ model(resol, eckit::LocalConfiguration(fullConfig, "model"));

//  Setup initial state
    const eckit::LocalConfiguration initialConfig(fullConfig, "initial condition");
    State_ xx(resol, initialConfig);
    Log::test() << "Initial state: " << xx << std::endl;

//  Setup augmented state
    const ModelAux_ moderr(resol, fullConfig.getSubConfiguration("model aux control"));

//  Setup times
    const util::Duration fclength(fullConfig.getString("forecast length"));
    const util::DateTime bgndate(xx.validTime());
    const util::DateTime enddate(bgndate + fclength);
    Log::info() << "Running forecast from " << bgndate << " to " << enddate << std::endl;

//  Setup forecast outputs
    PostProcessor<State_> post;

    eckit::LocalConfiguration prtConfig = fullConfig.getSubConfiguration("prints");
    post.enrollProcessor(new StateInfo<State_>("fc", prtConfig));

    eckit::LocalConfiguration outConfig(fullConfig, "output");
    outConfig.set("date", bgndate.toString());
    post.enrollProcessor(new StateWriter<State_>(outConfig));

    if (fullConfig.has("forecast to structured grid")) {
      eckit::LocalConfiguration structConfig(fullConfig, "forecast to structured grid");
      structConfig.set("date", bgndate.toString());
      post.enrollProcessor(new StructuredGridPostProcessor<MODEL, State_>(structConfig, resol));
    }

    eckit::LocalConfiguration wflow;
    wflow.set("frequency", "PT3H");
    post.enrollProcessor(new WorkflowUpdater<State_>(wflow));

//  Run forecast
    model.forecast(xx, moderr, fclength, post);

    Log::test() << "Final state: " << xx << std::endl;

    return 0;
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
