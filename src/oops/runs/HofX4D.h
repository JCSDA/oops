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

#include "eckit/exception/Exceptions.h"
#include "oops/base/Geometry.h"
#include "oops/base/instantiateObsFilterFactory.h"
#include "oops/base/Model.h"
#include "oops/base/ObsAuxControls.h"
#include "oops/base/ObsErrors.h"
#include "oops/base/Observations.h"
#include "oops/base/Observers.h"
#include "oops/base/ObsSpaces.h"
#include "oops/base/ObsTypeParameters.h"
#include "oops/base/PostProcessor.h"
#include "oops/base/State.h"
#include "oops/base/StateInfo.h"
#include "oops/generic/instantiateObsErrorFactory.h"
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
  typedef ModelParametersWrapper<MODEL> ModelParameters_;
  typedef typename State_::Parameters_ StateParameters_;
  typedef typename ModelAux_::Parameters_ ModelAuxParameters_;

  /// Only observations taken at times lying in the (`window begin`, `window begin` + `window
  /// length`] interval will be included in observation spaces.
  RequiredParameter<util::DateTime> windowBegin{"window begin", this};
  RequiredParameter<util::Duration> windowLength{"window length", this};

  /// A list whose elements determine treatment of observations from individual observation spaces.
  Parameter<std::vector<ObsTypeParameters<OBS>>> observations{"observations", {}, this};

  /// Geometry parameters.
  RequiredParameter<GeometryParameters_> geometry{"geometry", this};

  /// Options passed to the object writing out forecast fields.
  Parameter<PostTimerParameters> prints{"prints", {}, this};

  /// Whether to perturb the H(x) vector before saving.
  Parameter<bool> obsPerturbations{"obs perturbations", false, this};

  /// Whether to save the H(x) vector as ObsValues.
  Parameter<bool> makeObs{"make obs", false, this};

  /// Forecast length.
  RequiredParameter<util::Duration> forecastLength{"forecast length", this};

  /// Model parameters.
  RequiredParameter<ModelParameters_> model{"model", this};

  /// Initial state parameters.
  RequiredParameter<StateParameters_> initialCondition{"initial condition", this};

  /// Augmented model state.
  Parameter<ModelAuxParameters_> modelAuxControl{"model aux control", {}, this};
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
  typedef ObsErrors<OBS>             ObsErrors_;
  typedef Observers<MODEL, OBS>      Observers_;
  typedef ObsSpaces<OBS>             ObsSpaces_;
  typedef State<MODEL>               State_;

  typedef HofX4DParameters<MODEL, OBS> HofX4DParameters_;

 public:
// -----------------------------------------------------------------------------
  explicit HofX4D(const eckit::mpi::Comm & comm = oops::mpi::world()) : Application(comm) {
    instantiateObsErrorFactory<OBS>();
    instantiateObsFilterFactory<OBS>();
  }
// -----------------------------------------------------------------------------
  virtual ~HofX4D() = default;
// -----------------------------------------------------------------------------
  int execute(const eckit::Configuration & fullConfig) const override {
//  Deserialize parameters
    HofX4DParameters_ params;
    params.validateAndDeserialize(fullConfig);

//  Setup observation window
    const util::Duration winlen = params.windowLength;
    const util::DateTime winbgn = params.windowBegin;
    const util::DateTime winend(winbgn + winlen);
    Log::info() << "Observation window from " << winbgn << " to " << winend << std::endl;

//  Setup geometry
    const Geometry_ geometry(params.geometry, this->getComm(), mpi::myself());

//  Setup initial state
    State_ xx(geometry, params.initialCondition);
    Log::test() << "Initial state: " << xx << std::endl;

//  Check that window specified for forecast is at least the same as obs window
    const util::Duration fclength = params.forecastLength;
    if (winbgn < xx.validTime() || winend > xx.validTime() + fclength) {
      Log::error() << "Observation window can not be outside of forecast window." << std::endl;
      Log::error() << "Obs window: " << winbgn << " to " << winend << std::endl;
      Log::error() << "Forecast runs from: " << xx.validTime() << " for " << fclength << std::endl;
      throw eckit::BadValue("Observation window can not be outside of forecast window.");
    }

//  Setup observations
    ObsSpaces_ obspaces(obsSpaceParameters(params.observations.value()),
                        this->getComm(), winbgn, winend);
    ObsAux_ obsaux(obspaces, obsAuxParameters(params.observations.value()));
    ObsErrors_ Rmat(obsErrorParameters(params.observations.value()), obspaces);

//  Setup and initialize observer
    PostProcessor<State_> post;
    Observers_ hofx(obspaces, observerParameters(params.observations.value()));
    hofx.initialize(geometry, obsaux, Rmat, post);

//  Setup Model
    const Model_ model(geometry, params.model.value().modelParameters);
    ModelAux_ moderr(geometry, params.modelAuxControl);

    post.enrollProcessor(new StateInfo<State_>("fc", params.prints));

//  Run the model and compute H(x)
    model.forecast(xx, moderr, fclength, post);
    Log::test() << "Final state: " << xx << std::endl;

//  Get observations from observer
    Observations_ yobs(obspaces);
    hofx.finalize(yobs);
    Log::test() << "H(x): " << std::endl << yobs << "End H(x)" << std::endl;

//  Perturb H(x) if needed
    if (params.obsPerturbations) {
      yobs.perturb(Rmat);
      Log::test() << "Perturbed H(x): " << std::endl << yobs << "End Perturbed H(x)" << std::endl;
    }

//  Save H(x) as observations (if "make obs" == true)
    if (params.makeObs) yobs.save("ObsValue");
    obspaces.save();

    return 0;
  }
// -----------------------------------------------------------------------------
  void outputSchema(const std::string & outputPath) const override {
    HofX4DParameters_ params;
    params.outputSchema(outputPath);
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
