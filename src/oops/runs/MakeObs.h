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

#include <memory>
#include <string>
#include <vector>

#include <boost/shared_ptr.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "oops/base/Departures.h"
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
#include "oops/interface/State.h"
#include "oops/parallel/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

namespace oops {

template <typename MODEL> class MakeObs : public Application {
  typedef Departures<MODEL>          Departures_;
  typedef Geometry<MODEL>            Geometry_;
  typedef Model<MODEL>               Model_;
  typedef ModelAuxControl<MODEL>     ModelAux_;
  typedef ObsAuxControls<MODEL>      ObsAuxCtrls_;
  typedef Observations<MODEL>        Observations_;
  typedef ObsSpaces<MODEL>           ObsSpaces_;
  typedef State<MODEL>               State_;

 public:
// -----------------------------------------------------------------------------
  explicit MakeObs(const eckit::mpi::Comm & comm = oops::mpi::comm()) : Application(comm) {
    instantiateObsErrorFactory<MODEL>();
    instantiateObsFilterFactory<MODEL>();
  }
// -----------------------------------------------------------------------------
  virtual ~MakeObs() {}
// -----------------------------------------------------------------------------
  int execute(const eckit::Configuration & fullConfig) const {
//  Setup observation window
    const eckit::LocalConfiguration windowConf(fullConfig, "Assimilation Window");
    const util::DateTime bgn(windowConf.getString("Begin"));
    const util::DateTime end(windowConf.getString("End"));
    const util::Duration fclen(end - bgn);
    Log::info() << "Observation window is:" << windowConf << std::endl;

//  Setup resolution
    const eckit::LocalConfiguration resolConfig(fullConfig, "Geometry");
    const Geometry_ resol(resolConfig, this->getComm());

//  Setup Model
    const eckit::LocalConfiguration modelConfig(fullConfig, "Model");
    const Model_ model(resol, modelConfig);

//  Setup initial "true" state
    const eckit::LocalConfiguration initialConfig(fullConfig, "Initial Condition");
    Log::info() << "Initial configuration is:" << initialConfig << std::endl;
    State_ xx(resol, model.variables(), initialConfig);
    Log::test() << "Initial state: " << xx << std::endl;

//  Setup augmented state
    ModelAux_ moderr(resol, initialConfig);

//  Setup forecast outputs
    PostProcessor<State_> post;

    eckit::LocalConfiguration prtConf;
    fullConfig.get("prints", prtConf);
    post.enrollProcessor(new StateInfo<State_>("fc", prtConf));

//  Setup observations
    const eckit::LocalConfiguration obsconf(fullConfig, "Observations");
    Log::info() << "Observation configuration is:" << obsconf << std::endl;
    ObsSpaces_ obspace(obsconf, this->getComm(), bgn, end);

//  Setup observations bias
    ObsAuxCtrls_ ybias(obspace, obsconf);

//  Setup Observer
    boost::shared_ptr<Observers<MODEL, State_> >
      pobs(new Observers<MODEL, State_>(obsconf, obspace, ybias));
    post.enrollProcessor(pobs);

//  Run forecast and generate observations
    model.forecast(xx, moderr, fclen, post);
    Log::test() << "Final state: " << xx << std::endl;
    Log::info() << "MakeObs: Finished observation generation." << std::endl;

    std::unique_ptr<Observations_> yobs(pobs->release());
    Log::info() << "Generated observation: " << *yobs << std::endl;

//  Perturb observations
    if (obsconf.has("obspert")) {
      Departures_ ypert(obspace);
      ObsErrors<MODEL> matR(obsconf, obspace);
      matR.randomize(ypert);
      double opert = obsconf.getDouble("obspert");
      ypert *= opert;
      *yobs += ypert;
      Log::info() << "Perturbed observation: " << *yobs << std::endl;
    }

//  Save observations
    for (std::size_t jj = 0; jj < yobs->size(); ++jj) {
      Log::test() << "Generated observation: " << (*yobs)[jj] << std::endl;
    }
    yobs->save("ObsValue");

    return 0;
  }
// -----------------------------------------------------------------------------
 private:
  std::string appname() const {
    return "oops::MakeObs<" + MODEL::name() + ">";
  }
// -----------------------------------------------------------------------------
};

}  // namespace oops
#endif  // OOPS_RUNS_MAKEOBS_H_
