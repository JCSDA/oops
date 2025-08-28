/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Copyright 2020-2021 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_RUNS_HOFX4D_H_
#define OOPS_RUNS_HOFX4D_H_

#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"
#include "oops/base/Departures.h"
#include "oops/base/Geometry.h"
#include "oops/base/Model.h"
#include "oops/base/ObsAuxControls.h"
#include "oops/base/ObsErrors.h"
#include "oops/base/Observations.h"
#include "oops/base/Observers.h"
#include "oops/base/ObsSpaces.h"
#include "oops/base/PostProcessor.h"
#include "oops/base/State.h"
#include "oops/base/StateInfo.h"
#include "oops/base/StateWriter.h"
#include "oops/base/StructuredGridPostProcessor.h"
#include "oops/generic/instantiateModelFactory.h"
#include "oops/generic/instantiateObsErrorFactory.h"
#include "oops/interface/ModelAuxControl.h"
#include "oops/interface/ObsDataVector.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"
#include "oops/util/TimeWindow.h"

namespace oops {

// -----------------------------------------------------------------------------

/// Application runs model forecast from "initial condition" for the "forecast length"
/// and computes H(x) on the run.
template <typename MODEL, typename OBS> class HofX4D : public Application {
  typedef Departures<OBS>            Departures_;
  typedef Geometry<MODEL>            Geometry_;
  typedef Model<MODEL>               Model_;
  typedef ModelAuxControl<MODEL>     ModelAux_;
  typedef ObsAuxControls<OBS>        ObsAux_;
  typedef Observations<OBS>          Observations_;
  typedef ObsDataVector<OBS, int>    ObsDataInt_;
  typedef ObsErrors<OBS>             ObsErrors_;
  typedef Observers<MODEL, OBS>      Observers_;
  typedef ObsSpaces<OBS>             ObsSpaces_;
  typedef State<MODEL>               State_;

 public:
// -----------------------------------------------------------------------------
  explicit HofX4D(const eckit::mpi::Comm & comm = oops::mpi::world()) : Application(comm) {
    instantiateModelFactory<MODEL>();
    instantiateObsErrorFactory<OBS>();
  }
// -----------------------------------------------------------------------------
  virtual ~HofX4D() = default;
// -----------------------------------------------------------------------------
  int execute(const eckit::Configuration & fullConfig) const override {
//  Setup observation window
    const util::TimeWindow timeWindow(fullConfig.getSubConfiguration("time window"));
    Log::info() << "HofX4D observation window: " << timeWindow << std::endl;

//  Setup geometry
    const eckit::LocalConfiguration resolConfig(fullConfig, "geometry");
    const Geometry_ geometry(resolConfig, this->getComm(), mpi::myself());

//  Setup initial state
    const eckit::LocalConfiguration initialConfig(fullConfig, "initial condition");
    State_ xx(geometry, initialConfig);
    Log::test() << "Initial state: " << xx << std::endl;

//  Check that window specified for forecast is at least the same as obs window
    const util::Duration fclength(fullConfig.getString("forecast length"));
    const util::DateTime bgndate(xx.validTime());

    if (timeWindow.start() < bgndate ||
        timeWindow.end() > bgndate + fclength) {
        Log::error() << "Observation window can not be outside of forecast window." << std::endl;
        Log::error() << "Obs window: " << timeWindow.start() << " to "
                     << timeWindow.end() << std::endl;
        Log::error() << "Forecast runs from: " << bgndate << " for "
                     << fclength << std::endl;
        throw eckit::BadValue("Observation window can not be outside of forecast window.");
    }

//  Setup observations
    const eckit::LocalConfiguration oConfig(fullConfig, "observations");
    const eckit::LocalConfiguration obsConfig(oConfig, "observers");
    ObsSpaces_ obspaces(obsConfig, this->getComm(), timeWindow);
    ObsAux_ obsaux(obspaces, obsConfig);
    ObsErrors_ Rmat(obsConfig, obspaces);

//  Setup and initialize observer
    PostProcessor<State_> post;
    Observers_ hop(obspaces, oConfig);
    hop.initialize(geometry, obsaux, Rmat, post);

//  Setup Model
    const Model_ model(geometry, eckit::LocalConfiguration(fullConfig, "model"));
    ModelAux_ moderr(geometry, fullConfig.getSubConfiguration("model aux control"));

    const eckit::LocalConfiguration prtConfig = fullConfig.getSubConfiguration("prints");
    post.enrollProcessor(new StateInfo<State_>("fc", prtConfig));

//  Setup forecast outputs if requested
    if (fullConfig.has("output")) {
      eckit::LocalConfiguration outConfig(fullConfig, "output");
      outConfig.set("date", bgndate.toString());
      post.enrollProcessor(new StateWriter<State_>(outConfig));
    }
    if (fullConfig.has("forecast to structured grid")) {
      eckit::LocalConfiguration structConfig(fullConfig, "forecast to structured grid");
      structConfig.set("date", bgndate.toString());
      post.enrollProcessor(new StructuredGridPostProcessor<MODEL, State_>(structConfig, geometry));
    }

//  Run the model and compute H(x)
    model.forecast(xx, moderr, fclength, post);
    Log::test() << "Final state: " << xx << std::endl;

//  Get observations from observer
    Observations_ hofx(obspaces);
    std::vector<ObsDataInt_> qcflags;
    for (size_t jj = 0; jj < obspaces.size(); ++jj) {
      ObsDataInt_ qc(obspaces[jj], obspaces[jj].obsvariables());
      qcflags.push_back(qc);
    }
    hop.finalize(hofx, qcflags);
    Log::info() << "H(x): " << hofx.info("H(x): ") << std::endl;
    Log::test() << "H(x): " << hofx << std::endl << "End H(x)" << std::endl;

//  Perturb H(x) if needed
    if (oConfig.getBool("obs perturbations", false)) {
      hofx.perturb(Rmat);
      Log::info() << "Perturbed H(x): " << hofx.info("Perturbed H(x): ") << std::endl;
      Log::test() << "Perturbed H(x): " << hofx << std::endl << "End Perturbed H(x)" << std::endl;
    }

//  O-B diagnostics if obs available
    if (obspaces.has("ObsValue")) {
      Observations_ yobs(obspaces, "ObsValue");
      Departures_ ydep(hofx - yobs);
      Log::info() << "O-B :" << ydep.info("O-B") << std::endl;
    }

//  Save H(x) as observations (if "make obs" == true)
    if (fullConfig.getBool("make obs", false)) hofx.save("ObsValue");
    obspaces.save();

    return 0;
  }
// -----------------------------------------------------------------------------
 private:
  std::string appname() const override {
    return "oops::HofX4D<" + MODEL::name() + ", " + OBS::name() + ">";
  }
// -----------------------------------------------------------------------------
};

}  // namespace oops

#endif  // OOPS_RUNS_HOFX4D_H_
