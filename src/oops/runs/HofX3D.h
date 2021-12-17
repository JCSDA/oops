/*
 * (C) Copyright 2018-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_RUNS_HOFX3D_H_
#define OOPS_RUNS_HOFX3D_H_

#include <string>
#include <vector>

#include "eckit/exception/Exceptions.h"
#include "oops/base/Geometry.h"
#include "oops/base/instantiateObsFilterFactory.h"
#include "oops/base/ObsAuxControls.h"
#include "oops/base/ObsErrors.h"
#include "oops/base/Observations.h"
#include "oops/base/Observers.h"
#include "oops/base/ObsSpaces.h"
#include "oops/base/ObsTypeParameters.h"
#include "oops/base/PostProcessor.h"
#include "oops/base/State.h"
#include "oops/generic/instantiateObsErrorFactory.h"
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

/// \brief Top-level options taken by the HofX3D application.
template <typename MODEL, typename OBS>
class HofX3DParameters : public ApplicationParameters {
  OOPS_CONCRETE_PARAMETERS(HofX3DParameters, ApplicationParameters)

  typedef Geometry<MODEL> Geometry_;
  typedef State<MODEL> State_;

 public:
  typedef typename Geometry_::Parameters_ GeometryParameters_;
  typedef typename State_::Parameters_ StateParameters_;

  /// Only observations taken at times lying in the (`window begin`, `window begin` + `window
  /// length`] interval will be included in observation spaces.
  RequiredParameter<util::DateTime> windowBegin{"window begin", this};
  RequiredParameter<util::Duration> windowLength{"window length", this};

  /// A list whose elements determine treatment of observations from individual observation spaces.
  Parameter<std::vector<ObsTypeParameters<OBS>>> observations{"observations", {}, this};

  /// Geometry parameters.
  RequiredParameter<GeometryParameters_> geometry{"geometry", this};

  /// Whether to perturb the H(x) vector before saving.
  Parameter<bool> obsPerturbations{"obs perturbations", false, this};

  /// Whether to save the H(x) vector as ObsValues.
  Parameter<bool> makeObs{"make obs", false, this};

  /// Initial state parameters.
  RequiredParameter<StateParameters_> initialCondition{"state", this};
};

// -----------------------------------------------------------------------------

/// Application computes H(x) in 3D mode from the state that is read in.
template <typename MODEL, typename OBS> class HofX3D : public Application {
  typedef Geometry<MODEL>            Geometry_;
  typedef ObsAuxControls<OBS>        ObsAux_;
  typedef Observations<OBS>          Observations_;
  typedef ObsErrors<OBS>             ObsErrors_;
  typedef Observers<MODEL, OBS>      Observers_;
  typedef ObsSpaces<OBS>             ObsSpaces_;
  typedef State<MODEL>               State_;

  typedef HofX3DParameters<MODEL, OBS> HofX3DParameters_;

 public:
// -----------------------------------------------------------------------------
  explicit HofX3D(const eckit::mpi::Comm & comm = oops::mpi::world()) : Application(comm) {
    instantiateObsErrorFactory<OBS>();
    instantiateObsFilterFactory<OBS>();
  }
// -----------------------------------------------------------------------------
  virtual ~HofX3D() = default;
// -----------------------------------------------------------------------------
  int execute(const eckit::Configuration & fullConfig) const override {
//  Deserialize parameters
    HofX3DParameters_ params;
    params.validateAndDeserialize(fullConfig);

//  Setup observation window
    const util::Duration winlen = params.windowLength;
    const util::DateTime winbgn = params.windowBegin;
    const util::DateTime winend(winbgn + winlen);
    const util::DateTime winHalf = winbgn + winlen/2;
    Log::info() << "Observation window from " << winbgn << " to " << winend << std::endl;

//  Setup geometry
    const Geometry_ geometry(params.geometry, this->getComm(), oops::mpi::myself());

//  Setup state
    State_ xx(geometry, params.initialCondition);
    Log::test() << "State: " << xx << std::endl;

//  Check that state is inside the obs window
    if (xx.validTime() < winbgn || xx.validTime() > winend) {
      Log::error() << "State time: " << xx.validTime() << std::endl;
      Log::error() << "Obs window: " << winbgn << " to " << winend << std::endl;
      throw eckit::BadValue("State cannot be outside observation window.");
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

//  Compute H(x)
    post.initialize(xx, winHalf, winlen);
    post.process(xx);
    post.finalize(xx);

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
    HofX3DParameters_ params;
    params.outputSchema(outputPath);
  }
// -----------------------------------------------------------------------------
 private:
  std::string appname() const override {
    return "oops::HofX3D<" + MODEL::name() + ", " + OBS::name() + ">";
  }
// -----------------------------------------------------------------------------
};

}  // namespace oops

#endif  // OOPS_RUNS_HOFX3D_H_
