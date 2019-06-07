/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_RUNS_FINDLOCALOBS_H_
#define OOPS_RUNS_FINDLOCALOBS_H_

#include <string>
#include <vector>

#include <boost/scoped_ptr.hpp>
#include <boost/shared_ptr.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "oops/base/instantiateObsFilterFactory.h"
#include "oops/base/ObsAuxControls.h"
#include "oops/base/Observations.h"
#include "oops/base/Observers.h"
#include "oops/base/ObsOperators.h"
#include "oops/base/ObsSpaces.h"
#include "oops/base/PostProcessor.h"
#include "oops/base/StateInfo.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/GeometryIterator.h"
#include "oops/interface/Model.h"
#include "oops/interface/ModelAuxControl.h"
#include "oops/interface/State.h"
#include "oops/runs/Application.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

namespace oops {

template <typename MODEL> class FindLocalObs : public Application {
  typedef Geometry<MODEL>            Geometry_;
  typedef GeometryIterator<MODEL>    GeometryIterator_;
  typedef Model<MODEL>               Model_;
  typedef ModelAuxControl<MODEL>     ModelAux_;
  typedef ObsAuxControls<MODEL>      ObsAuxCtrls_;
  typedef Observations<MODEL>        Observations_;
  typedef ObsOperators<MODEL>        ObsOperators_;
  typedef ObsSpaces<MODEL>           ObsSpaces_;
  typedef State<MODEL>               State_;

 public:
// -----------------------------------------------------------------------------
  FindLocalObs() {
    instantiateObsFilterFactory<MODEL>();
  }
// -----------------------------------------------------------------------------
  virtual ~FindLocalObs() {}
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
    eckit::LocalConfiguration obsconf(fullConfig, "Observations");
    Log::debug() << "Observations configuration is:" << obsconf << std::endl;
    ObsSpaces_ obsdb(obsconf, winbgn, winend);
    ObsOperators_ hop(obsdb, obsconf);

//  Setup observations bias
    ObsAuxCtrls_ ybias(obsconf);

//  Setup Observers
    boost::shared_ptr<Observers<MODEL, State_> >
      pobs(new Observers<MODEL, State_>(obsdb, hop, ybias));
    post.enrollProcessor(pobs);

//  Compute H(x)
    model.forecast(xx, moderr, winlen, post);
    Log::info() << "FindLocalObs: Finished observation computation." << std::endl;
    Log::test() << "Final state: " << xx << std::endl;

//  Save H(x)
    boost::scoped_ptr<Observations_> yobs(pobs->release());
    Log::test() << "H(x): " << *yobs << std::endl;
    yobs->save("hofx");

//  Iterate over all gridpoints and find local observations (do nothing else for now):
    eckit::LocalConfiguration localconfig(fullConfig, "Localization");
    double dist = localconfig.getDouble("distance");
    int max_nobs = localconfig.getInt("max_nobs");
    for (GeometryIterator_ i = resol.begin(); i != resol.end(); ++i) {
       // find all local observations around current gridpoint (*i) with dist and maxnum from config
       ObsSpaces_ localobs(obsdb, *i, dist, max_nobs);
       Log::test() << *i << localobs << std::endl;
    }


    return 0;
  }
// -----------------------------------------------------------------------------
 private:
  std::string appname() const {
    return "oops::FindLocalObs<" + MODEL::name() + ">";
  }
// -----------------------------------------------------------------------------
};

}  // namespace oops

#endif  // OOPS_RUNS_FINDLOCALOBS_H_
