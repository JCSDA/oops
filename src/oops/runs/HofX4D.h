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
#include "oops/generic/instantiateObsErrorFactory.h"
#include "oops/interface/ModelAuxControl.h"
#include "oops/interface/ObsDataVector.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "oops/util/TimeWindow.h"

namespace oops {

// -----------------------------------------------------------------------------

/// \brief Top-level options taken by the HofX4D application.
template <typename MODEL, typename OBS>
class HofX4DParameters : public ApplicationParameters {
  OOPS_CONCRETE_PARAMETERS(HofX4DParameters, ApplicationParameters)

  typedef Geometry<MODEL> Geometry_;
  typedef State<MODEL> State_;
  typedef ModelAuxControl<MODEL>     ModelAux_;

 public:
  typedef typename Geometry_::Parameters_ GeometryParameters_;

  /// Options describing the assimilation time window.
  RequiredParameter<eckit::LocalConfiguration> timeWindow{"time window", this};

  /// Options describing the observations and their treatment
  RequiredParameter<eckit::LocalConfiguration> observations{"observations", this};

  /// Geometry parameters.
  RequiredParameter<GeometryParameters_> geometry{"geometry", this};

  /// Options passed to the object writing out forecast fields.
  Parameter<PostTimerParameters> prints{"prints", {}, this};

  /// Whether to save the H(x) vector as ObsValues.
  Parameter<bool> makeObs{"make obs", false, this};

  /// Forecast length.
  RequiredParameter<util::Duration> forecastLength{"forecast length", this};

  /// Model parameters.
  RequiredParameter<eckit::LocalConfiguration> model{"model", this};

  /// Initial state parameters.
  RequiredParameter<eckit::LocalConfiguration> initialCondition{"initial condition", this};

  /// Augmented model state.
  Parameter<eckit::LocalConfiguration> modelAuxControl{"model aux control",
                                                       eckit::LocalConfiguration(), this};
};

// -----------------------------------------------------------------------------

/// Application runs model forecast from "initial condition" for the "forecast length"
/// and computes H(x) on the run.
template <typename MODEL, typename OBS> class HofX4D : public Application {
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

  typedef HofX4DParameters<MODEL, OBS> HofX4DParameters_;

 public:
// -----------------------------------------------------------------------------
  explicit HofX4D(const eckit::mpi::Comm & comm = oops::mpi::world()) : Application(comm) {
    instantiateObsErrorFactory<OBS>();
  }
// -----------------------------------------------------------------------------
  virtual ~HofX4D() = default;
// -----------------------------------------------------------------------------
  int execute(const eckit::Configuration & fullConfig, bool validate) const override {
//  Deserialize parameters
    HofX4DParameters_ params;
    if (validate) params.validate(fullConfig);
    params.deserialize(fullConfig);

//  Setup observation window
    const util::TimeWindow timeWindow(fullConfig.getSubConfiguration("time window"));
    Log::info() << "HofX4D observation window: " << timeWindow << std::endl;

//  Setup geometry
    const Geometry_ geometry(params.geometry, this->getComm(), mpi::myself());

//  Setup initial state
    const eckit::LocalConfiguration initialConfig(fullConfig, "initial condition");
    State_ xx(geometry, initialConfig);
    Log::test() << "Initial state: " << xx << std::endl;

//  Check that window specified for forecast is at least the same as obs window
    const util::Duration fclength(fullConfig.getString("forecast length"));

    if (timeWindow.start() < xx.validTime() ||
        timeWindow.end() > xx.validTime() + fclength) {
        Log::error() << "Observation window can not be outside of forecast window." << std::endl;
        Log::error() << "Obs window: " << timeWindow.start() << " to "
                     << timeWindow.end() << std::endl;
        Log::error() << "Forecast runs from: " << xx.validTime() << " for "
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
    Observers_ hofx(obspaces, oConfig);
    hofx.initialize(geometry, obsaux, Rmat, post);

//  Setup Model
    const Model_ model(geometry, eckit::LocalConfiguration(fullConfig, "model"));
    ModelAux_ moderr(geometry, fullConfig.getSubConfiguration("model aux control"));

    const eckit::LocalConfiguration prtConfig = fullConfig.getSubConfiguration("prints");
    post.enrollProcessor(new StateInfo<State_>("fc", prtConfig));

//  Run the model and compute H(x)
    model.forecast(xx, moderr, fclength, post);
    Log::test() << "Final state: " << xx << std::endl;

//  Get observations from observer
    Observations_ yobs(obspaces);
    std::vector<ObsDataInt_> qcflags;
    for (size_t jj = 0; jj < obspaces.size(); ++jj) {
      ObsDataInt_ qc(obspaces[jj], obspaces[jj].obsvariables());
      qcflags.push_back(qc);
    }
    hofx.finalize(yobs, qcflags);
    Log::info() << "H(x): " << std::endl << yobs << "End H(x)" << std::endl;
    Log::test() << "H(x): " << std::endl << yobs << "End H(x)" << std::endl;

//  Perturb H(x) if needed
    if (oConfig.getBool("obs perturbations", false)) {
      yobs.perturb(Rmat);
      Log::test() << "Perturbed H(x): " << std::endl << yobs << "End Perturbed H(x)" << std::endl;
    }

//  Save H(x) as observations (if "make obs" == true)
    if (fullConfig.getBool("make obs", false)) yobs.save("ObsValue");
    obspaces.save();

    return 0;
  }
// -----------------------------------------------------------------------------
  void outputSchema(const std::string & outputPath) const override {
    HofX4DParameters_ params;
    params.outputSchema(outputPath);
  }
// -----------------------------------------------------------------------------
  void validateConfig(const eckit::Configuration & fullConfig) const override {
    HofX4DParameters_ params;
    params.validate(fullConfig);
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
