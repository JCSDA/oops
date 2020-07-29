/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_RUNS_MAKEOBS_H_
#define OOPS_RUNS_MAKEOBS_H_

#include <string>

#include "eckit/config/LocalConfiguration.h"
#include "oops/assimilation/CalcHofX.h"
#include "oops/base/Departures.h"
#include "oops/base/instantiateObsFilterFactory.h"
#include "oops/base/ObsErrors.h"
#include "oops/base/Observations.h"
#include "oops/base/ObsSpaces.h"
#include "oops/base/PostProcessor.h"
#include "oops/base/StateInfo.h"
#include "oops/generic/instantiateObsErrorFactory.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Model.h"
#include "oops/interface/State.h"
#include "oops/parallel/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

namespace oops {

template <typename MODEL, typename OBS> class MakeObs : public Application {
  typedef Departures<OBS>            Departures_;
  typedef Geometry<MODEL>            Geometry_;
  typedef Model<MODEL>               Model_;
  typedef Observations<OBS>          Observations_;
  typedef ObsErrors<OBS>             ObsErrors_;
  typedef ObsSpaces<OBS>             ObsSpaces_;
  typedef State<MODEL>               State_;

 public:
// -----------------------------------------------------------------------------
  explicit MakeObs(const eckit::mpi::Comm & comm = oops::mpi::comm()) : Application(comm) {
    instantiateObsErrorFactory<OBS>();
    instantiateObsFilterFactory<OBS>();
  }
// -----------------------------------------------------------------------------
  virtual ~MakeObs() {}
// -----------------------------------------------------------------------------
  int execute(const eckit::Configuration & fullConfig) const {
//  Setup observation window
    const util::DateTime winbgn(fullConfig.getString("window begin"));
    const util::Duration winlen(fullConfig.getString("window length"));
    const util::DateTime winend(winbgn + winlen);
    Log::info() << "Observation window from " << winbgn << " to " << winend << std::endl;

//  Setup geometry
    const eckit::LocalConfiguration geometryConfig(fullConfig, "geometry");
    const Geometry_ geometry(geometryConfig, this->getComm());

//  Setup Model
    const eckit::LocalConfiguration modelConfig(fullConfig, "model");
    const Model_ model(geometry, modelConfig);

//  Setup initial "true" state
    const eckit::LocalConfiguration initialConfig(fullConfig, "initial condition");
    State_ xx(geometry, initialConfig);
    Log::test() << "Initial state: " << xx << std::endl;

//  Setup forecast outputs
    PostProcessor<State_> post;

    eckit::LocalConfiguration prtConf;
    fullConfig.get("prints", prtConf);
    post.enrollProcessor(new StateInfo<State_>("fc", prtConf));

//  Setup observations
    const eckit::LocalConfiguration obsConfig(fullConfig);  // to be replaced
    ObsSpaces_ obspace(obsConfig, this->getComm(), winbgn, winend);

    CalcHofX<MODEL, OBS> hofx(obspace, geometry, fullConfig);
    const util::Duration flength(fullConfig.getString("forecast length"));
    Observations_ yobs = hofx.compute(model, xx, post, flength);

    Log::test() << "Final state: " << xx << std::endl;

    Log::info() << "Generated observation: " << yobs << std::endl;

//  Perturb observations
    if (fullConfig.has("obspert")) {
      Departures_ ypert(obspace);
      ObsErrors_ matR(obsConfig, obspace);
      matR.randomize(ypert);
      double opert = fullConfig.getDouble("obspert");
      ypert *= opert;
      yobs += ypert;
      Log::info() << "Perturbed observation: " << yobs << std::endl;
    }

//  Save observations
    for (std::size_t jj = 0; jj < yobs.size(); ++jj) {
      Log::test() << "Generated observation: " << yobs[jj] << std::endl;
    }
    yobs.save("ObsValue");

    return 0;
  }
// -----------------------------------------------------------------------------
 private:
  std::string appname() const {
    return "oops::MakeObs<" + MODEL::name() + ", " + OBS::name() + ">";
  }
// -----------------------------------------------------------------------------
};

}  // namespace oops
#endif  // OOPS_RUNS_MAKEOBS_H_
