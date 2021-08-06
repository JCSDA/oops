/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Copyright 2020-2020 UCAR.
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
#include "oops/base/instantiateObsFilterFactory.h"
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
#include "oops/generic/instantiateVariableChangeFactory.h"
#include "oops/interface/ModelAuxControl.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/algorithms.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

namespace oops {

// -----------------------------------------------------------------------------

/// \brief Options controlling the processing of observations from a single obs space.
template <typename OBS>
class ObsTypeParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ObsTypeParameters, Parameters)

 public:
  typedef typename ObsAuxControl<OBS>::Parameters_ ObsAuxControlParameters_;

  /// Options used to configure the observation space.
  oops::RequiredParameter<eckit::LocalConfiguration> obsSpace{"obs space", this};

  /// Options used to configure the observation operator, observation filters and GetValues.
  ObserverParameters<OBS> observer{this};

  /// Options used to configure the observation error model.
  oops::Parameter<eckit::LocalConfiguration> obsError{
    "obs error", eckit::LocalConfiguration(), this};

  /// Options used to configure bias correction.
  oops::Parameter<ObsAuxControlParameters_> obsBias{"obs bias", {}, this};
};

// -----------------------------------------------------------------------------

/// \brief Top-level options taken by the HofX4D application.
template <typename MODEL, typename OBS>
class HofX4DParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(HofX4DParameters, Parameters)

 public:
  typedef typename Geometry<MODEL>::Parameters_ GeometryParameters_;
  typedef ModelParametersWrapper<MODEL> ModelParameters_;

  /// Only observations taken at times lying in the (`window begin`, `window begin` + `window
  /// length`] interval will be included in observation spaces.
  oops::RequiredParameter<util::DateTime> windowBegin{"window begin", this};
  oops::RequiredParameter<util::Duration> windowLength{"window length", this};

  /// Forecast length.
  oops::RequiredParameter<util::Duration> forecastLength{"forecast length", this};

  /// A list whose elements determine treatment of observations from individual observation spaces.
  oops::Parameter<std::vector<ObsTypeParameters<OBS>>> observations{"observations", {}, this};

  /// Geometry parameters.
  oops::RequiredParameter<GeometryParameters_> geometry{"geometry", this};

  /// Model parameters.
  oops::RequiredParameter<ModelParameters_> model{"model", this};

  /// Initial state parameters.
  oops::RequiredParameter<eckit::LocalConfiguration> initialCondition{"initial condition", this};

  /// Options passed to the object writing out forecast fields.
  oops::Parameter<eckit::LocalConfiguration> prints{"prints", eckit::LocalConfiguration(), this};

  /// Whether to perturb the H(x) vector before saving.
  oops::Parameter<bool> obsPerturbations{"obs perturbations", false, this};

  /// Whether to save the H(x) vector as ObsValues.
  oops::Parameter<bool> makeObs{"make obs", false, this};

  /// Parameters used by regression tests comparing results produced by the application against
  /// known good outputs.
  oops::Parameter<eckit::LocalConfiguration> test{"test", eckit::LocalConfiguration(), this};
};

// -----------------------------------------------------------------------------

/// Application runs model forecast from "initial condition" for the "forecast length"
/// and computes H(x) on the run. If "obspert" is specified in the config, the resulting
/// H(x) is perturbed. It is saved as "hofx" by default, or as specified "hofx group name"
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
  typedef typename ObsAuxControl<OBS>::Parameters_ ObsAuxControlParameters_;
  typedef ObserverParameters<OBS> ObserverParameters_;
  typedef ObsTypeParameters<OBS> ObsTypeParameters_;

 public:
// -----------------------------------------------------------------------------
  explicit HofX4D(const eckit::mpi::Comm & comm = oops::mpi::world()) : Application(comm) {
    instantiateObsErrorFactory<OBS>();
    instantiateObsFilterFactory<OBS>();
    instantiateVariableChangeFactory<MODEL>();
  }
// -----------------------------------------------------------------------------
  virtual ~HofX4D() = default;
// -----------------------------------------------------------------------------
  int execute(const eckit::Configuration & fullConfig) const {
// Deserialize parameters
    HofX4DParameters_ params;
    params.validateAndDeserialize(fullConfig);

//  Setup observation window
    const util::Duration &winlen = params.windowLength;
    const util::DateTime &winbgn = params.windowBegin;
    const util::DateTime winend(winbgn + winlen);
    Log::info() << "Observation window from " << winbgn << " to " << winend << std::endl;

//  Setup geometry
    const Geometry_ geometry(params.geometry, this->getComm(), oops::mpi::myself());

//  Setup Model
    const Model_ model(geometry, params.model.value().modelParameters);

//  Setup initial state
    State_ xx(geometry, params.initialCondition);
    ModelAux_ moderr(geometry, params.initialCondition);
    const util::Duration &flength = params.forecastLength;
    Log::test() << "Initial state: " << xx << std::endl;

//  Check that window specified for forecast is at least the same as obs window
    if (winbgn < xx.validTime() || winend > xx.validTime() + flength) {
      Log::error() << "Observation window can not be outside of forecast window." << std::endl;
      Log::error() << "Obs window: " << winbgn << " to " << winend << std::endl;
      Log::error() << "Forecast runs from: " << xx.validTime() << " for " << flength << std::endl;
      throw eckit::BadValue("Observation window can not be outside of forecast window.");
    }

//  Setup forecast outputs
    PostProcessor<State_> post;

    post.enrollProcessor(new StateInfo<State_>("fc", params.prints));

//  Setup observations

    const eckit::LocalConfiguration obsConfig(fullConfig, "observations");
    ObsSpaces_ obspaces(obsConfig, this->getComm(), winbgn, winend);

    const std::vector<ObsAuxControlParameters_> obsauxParams = util::transformVector(
          params.observations.value(),
          [](const ObsTypeParameters_ & obsTypeParams) { return obsTypeParams.obsBias.value(); });
    ObsAux_ obsaux(obspaces, obsauxParams);

    ObsErrors_ Rmat(obsConfig, obspaces);

//  Setup and initialize observer
    std::vector<ObserverParameters_> observerParams = util::transformVector(
          params.observations.value(),
          [](const ObsTypeParameters_ & obsTypeParams) { return obsTypeParams.observer; });
    Observers_ hofx(obspaces, observerParams);
    hofx.initialize(geometry, obsaux, Rmat, post);

//  run the model and compute H(x)
    model.forecast(xx, moderr, flength, post);
    Log::test() << "Final state: " << xx << std::endl;

//  Get observations from observer
    Observations_ yobs(obspaces);
    hofx.finalize(yobs);
    Log::test() << "H(x): " << std::endl << yobs << "End H(x)" << std::endl;

//  Perturb H(x) if needed (can be used for generating obs in OSSE: perturbed H(x) could be saved
//  as ObsValue if "hofx group name" == ObsValue.
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
 private:
  std::string appname() const {
    return "oops::HofX4D<" + MODEL::name() + ", " + OBS::name() + ">";
  }
// -----------------------------------------------------------------------------
};

}  // namespace oops

#endif  // OOPS_RUNS_HOFX4D_H_
