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
#include <vector>

#include <boost/scoped_ptr.hpp>
#include <boost/shared_ptr.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "oops/base/Departures.h"
#include "oops/base/instantiateObsFilterFactory.h"
#include "oops/base/ObsErrors.h"
#include "oops/base/Observations.h"
#include "oops/base/Observer.h"
#include "oops/base/ObsFilters.h"
#include "oops/base/ObsOperators.h"
#include "oops/base/ObsSpaces.h"
#include "oops/base/PostProcessor.h"
#include "oops/base/StateInfo.h"
#include "oops/generic/instantiateObsErrorFactory.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Model.h"
#include "oops/interface/ModelAuxControl.h"
#include "oops/interface/ObsAuxControl.h"
#include "oops/interface/State.h"
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
  typedef ObsAuxControl<MODEL>       ObsAuxCtrl_;
  typedef Observations<MODEL>        Observations_;
  typedef ObsFilters<MODEL>          ObsFilters_;
  typedef ObsSpaces<MODEL>           ObsSpaces_;
  typedef ObsOperators<MODEL>        ObsOperators_;
  typedef State<MODEL>               State_;

 public:
// -----------------------------------------------------------------------------
  MakeObs() {
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
    const Geometry_ resol(resolConfig);

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
    eckit::LocalConfiguration biasConf;
    fullConfig.get("ObsBias", biasConf);
    ObsAuxCtrl_ ybias(biasConf);

//  Setup observations
    const eckit::LocalConfiguration obsconf(fullConfig, "Observations");
    Log::info() << "Observation configuration is:" << obsconf << std::endl;
    ObsSpaces_ obspace(obsconf, bgn, end);
    ObsOperators_ hop(obspace, obsconf);

//  Setup QC filters
    std::vector<eckit::LocalConfiguration> typeconfs;
    obsconf.get("ObsTypes", typeconfs);
    std::vector<ObsFilters_> filters;
    for (size_t jj = 0; jj < obspace.size(); ++jj) {
      filters.push_back(ObsFilters_(obspace[jj], typeconfs[jj]));
    }

//  Setup Observer
    boost::shared_ptr<Observer<MODEL, State_> >
      pobs(new Observer<MODEL, State_>(obspace, hop, ybias, filters));
    post.enrollProcessor(pobs);

//  Run forecast and generate observations
    model.forecast(xx, moderr, fclen, post);
    Log::test() << "Final state: " << xx << std::endl;
    Log::info() << "MakeObs: Finished observation generation." << std::endl;

    boost::scoped_ptr<Observations_> yobs(pobs->release());
    Log::info() << "Generated observation: " << *yobs << std::endl;

//  Perturb observations
    if (obsconf.has("obspert")) {
      Departures_ ypert(obspace, hop);
      ObsErrors<MODEL> matR(obsconf, obspace, hop);
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
