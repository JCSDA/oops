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

#include <memory>
#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"
#include "oops/base/instantiateObsFilterFactory.h"
#include "oops/base/ObsAuxControls.h"
#include "oops/base/ObsErrors.h"
#include "oops/base/Observations.h"
#include "oops/base/Observers.h"
#include "oops/base/ObsSpaces.h"
#include "oops/base/PostProcessor.h"
#include "oops/base/StateInfo.h"
#include "oops/generic/instantiateObsErrorFactory.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Model.h"
#include "oops/interface/ModelAuxControl.h"
#include "oops/interface/ObsDataVector.h"
#include "oops/interface/ObsVector.h"
#include "oops/interface/State.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

namespace oops {

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
  typedef ObsVector<OBS>             ObsVector_;
  typedef State<MODEL>               State_;
  typedef ObsDataVector<OBS, float>  ObsData_;

 public:
// -----------------------------------------------------------------------------
  explicit HofX4D(const eckit::mpi::Comm & comm = oops::mpi::world()) : Application(comm) {
    instantiateObsErrorFactory<OBS>();
    instantiateObsFilterFactory<OBS>();
  }
// -----------------------------------------------------------------------------
  virtual ~HofX4D() = default;
// -----------------------------------------------------------------------------
  int execute(const eckit::Configuration & fullConfig) const {
//  Setup observation window
    const util::Duration winlen(fullConfig.getString("window length"));
    const util::DateTime winbgn(fullConfig.getString("window begin"));
    const util::DateTime winend(winbgn + winlen);
    Log::info() << "Observation window from " << winbgn << " to " << winend << std::endl;

//  Setup geometry
    const eckit::LocalConfiguration geometryConfig(fullConfig, "geometry");
    const Geometry_ geometry(geometryConfig, this->getComm());

//  Setup Model
    const eckit::LocalConfiguration modelConfig(fullConfig, "model");
    const Model_ model(geometry, modelConfig);

//  Setup initial state
    const eckit::LocalConfiguration initialConfig(fullConfig, "initial condition");
    State_ xx(geometry, initialConfig);
    ModelAux_ moderr(geometry, initialConfig);
    const util::Duration flength(fullConfig.getString("forecast length"));
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

    eckit::LocalConfiguration prtConf;
    fullConfig.get("prints", prtConf);
    post.enrollProcessor(new StateInfo<State_>("fc", prtConf));

//  Setup observations
    const eckit::LocalConfiguration obsConfig(fullConfig, "observations");
    ObsSpaces_ obspaces(obsConfig, this->getComm(), winbgn, winend);
    ObsAux_ obsaux(obspaces, obsConfig);

//  Setup and run observer
    Observers_ hofx(obspaces, obsConfig);

    std::vector<std::shared_ptr<ObsVector_>> obserrs;
    for (size_t jj = 0; jj < obspaces.size(); ++jj) {
      obserrs.emplace_back(std::make_shared<ObsVector_>(obspaces[jj], "ObsError"));
    }

//  run the model and compute H(x)
    std::shared_ptr<PostBase<State<MODEL> > > getvals(hofx.initialize(obsaux, obserrs));
    post.enrollProcessor(getvals);
    model.forecast(xx, moderr, flength, post);

    Observations_ yobs = hofx.finalize();
    for (size_t jj = 0; jj < obserrs.size(); ++jj) {
      obserrs[jj]->save("EffectiveError");  // Obs error covariance is looking for that for now
    }

    Log::test() << "Final state: " << xx << std::endl;
    Log::test() << "H(x): " << std::endl << yobs << "End H(x)" << std::endl;

//  Perturb H(x) if needed (can be used for generating obs in OSSE: perturbed H(x) could be saved
//  as ObsValue if "hofx group name" == ObsValue.
    bool obspert = fullConfig.getBool("obs perturbations", false);
    if (obspert) {
      ObsErrors_ matR(obsConfig, obspaces);
      yobs.perturb(matR);
      Log::test() << "Perturbed H(x): " << std::endl << yobs << "End Perturbed H(x)" << std::endl;
    }

//  Save H(x) either as observations (if "make obs" == true) or as "hofx"
    const bool makeobs = fullConfig.getBool("make obs", false);
    if (makeobs) {
      yobs.save("ObsValue");
    } else {
      yobs.save("hofx");
    }

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
