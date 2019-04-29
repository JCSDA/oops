/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_RUNS_HOFX_H_
#define OOPS_RUNS_HOFX_H_

#include <string>
#include <vector>

#include <boost/scoped_ptr.hpp>
#include <boost/shared_ptr.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "oops/base/instantiateObsFilterFactory.h"
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

template <typename MODEL> class HofX : public Application {
  typedef Geometry<MODEL>            Geometry_;
  typedef Model<MODEL>               Model_;
  typedef ModelAuxControl<MODEL>     ModelAux_;
  typedef ObsAuxControl<MODEL>       ObsAuxCtrl_;
  typedef Observations<MODEL>        Observations_;
  typedef ObsFilters<MODEL>          ObsFilters_;
  typedef ObsSpaces<MODEL>           ObsSpaces_;
  typedef ObsOperators<MODEL>        ObsOperators_;
  typedef State<MODEL>               State_;
  typedef boost::shared_ptr<ObsFilters_> PtrFilters_;

 public:
// -----------------------------------------------------------------------------
  HofX() {
    instantiateObsFilterFactory<MODEL>();
  }
// -----------------------------------------------------------------------------
  virtual ~HofX() {}
// -----------------------------------------------------------------------------
  int execute(const eckit::Configuration & fullConfig) const {
//  Setup observation window
    const eckit::LocalConfiguration windowConf(fullConfig, "Assimilation Window");
    const util::Duration winlen(windowConf.getString("Length"));
    const util::DateTime winbgn(windowConf.getString("Begin"));
    const util::DateTime winend(winbgn + winlen);
    Log::info() << "Observation window is:" << windowConf << std::endl;

//  Setup resolution
    const eckit::LocalConfiguration resolConfig(fullConfig, "Geometry");
    const Geometry_ resol(resolConfig);

//  Setup Model
    const eckit::LocalConfiguration modelConfig(fullConfig, "Model");
    const Model_ model(resol, modelConfig);

//  Setup initial state
    const eckit::LocalConfiguration initialConfig(fullConfig, "Initial Condition");
    Log::info() << "Initial configuration is:" << initialConfig << std::endl;
    State_ xx(resol, model.variables(), initialConfig);
    Log::test() << "Initial state: " << xx << std::endl;

//  Setup augmented state
    ModelAux_ moderr(resol, initialConfig);

//  Setup forecast outputs
    PostProcessor<State_> post;

    eckit::LocalConfiguration prtConf;
    fullConfig.get("Prints", prtConf);
    post.enrollProcessor(new StateInfo<State_>("fc", prtConf));

//  Setup observations
    eckit::LocalConfiguration biasConf;
    fullConfig.get("ObsBias", biasConf);
    ObsAuxCtrl_ ybias(biasConf);

//  Setup observations
    const eckit::LocalConfiguration obsconf(fullConfig, "Observations");
    Log::info() << "Observation configuration is:" << obsconf << std::endl;
    ObsSpaces_ obspace(obsconf, winbgn, winend);
    ObsOperators_ hop(obspace, obsconf);

//  Setup QC filters
    std::vector<eckit::LocalConfiguration> typeconfs;
    obsconf.get("ObsTypes", typeconfs);
    std::vector<PtrFilters_> filters;
    for (size_t jj = 0; jj < obspace.size(); ++jj) {
      PtrFilters_ tmp(new ObsFilters_(obspace[jj], typeconfs[jj], hop[jj].observed()));
      filters.push_back(tmp);
    }

//  Setup Observer
    boost::shared_ptr<Observer<MODEL, State_> >
      pobs(new Observer<MODEL, State_>(obspace, hop, ybias, filters));
    post.enrollProcessor(pobs);

//  Compute H(x)
    model.forecast(xx, moderr, winlen, post);
    Log::info() << "HofX: Finished observation computation." << std::endl;
    Log::test() << "Final state: " << xx << std::endl;

//  Save H(x)
    boost::scoped_ptr<Observations_> yobs(pobs->release());
    Log::test() << "H(x): " << std::endl << *yobs << "End H(x)" << std::endl;
    yobs->save("hofx");

    return 0;
  }
// -----------------------------------------------------------------------------
 private:
  std::string appname() const {
    return "oops::HofX<" + MODEL::name() + ">";
  }
// -----------------------------------------------------------------------------
};

}  // namespace oops

#endif  // OOPS_RUNS_HOFX_H_
