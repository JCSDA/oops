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
  typedef ObsOperators<MODEL>        ObsOperators_;
  typedef ObsSpaces<MODEL>           ObsSpaces_;
  typedef State<MODEL>               State_;

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
    eckit::LocalConfiguration obsconf(fullConfig, "Observations");
    Log::debug() << "Observations configuration is:" << obsconf << std::endl;
    ObsSpaces_ obsdb(obsconf, winbgn, winend);
    ObsOperators_ hop(obsdb, obsconf);

//  Setup QC filters
    std::vector<ObsFilters_> filters(obsdb.size());

//  Setup Observer
    boost::shared_ptr<Observer<MODEL, State_> >
      pobs(new Observer<MODEL, State_>(obsdb, hop, ybias, filters));
    post.enrollProcessor(pobs);

//  Compute H(x)
    model.forecast(xx, moderr, winlen, post);
    Log::info() << "HofX: Finished observation computation." << std::endl;
    Log::test() << "Final state: " << xx << std::endl;

//  Save H(x)
    boost::scoped_ptr<Observations_> yobs(pobs->release());
    Log::test() << "H(x): " << *yobs << std::endl;
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
