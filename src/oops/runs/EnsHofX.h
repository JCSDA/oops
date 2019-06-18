/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_RUNS_ENSHOFX_H_
#define OOPS_RUNS_ENSHOFX_H_

#include <string>
#include <vector>

#include <boost/scoped_ptr.hpp>
#include <boost/shared_ptr.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "oops/base/instantiateObsFilterFactory.h"
#include "oops/base/ObsAuxControls.h"
#include "oops/base/ObsEnsemble.h"
#include "oops/base/Observations.h"
#include "oops/base/Observers.h"
#include "oops/base/ObsSpaces.h"
#include "oops/base/PostProcessor.h"
#include "oops/base/StateInfo.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Model.h"
#include "oops/interface/ModelAuxControl.h"
#include "oops/interface/State.h"
#include "oops/runs/Application.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

namespace oops {

template <typename MODEL> class EnsHofX : public Application {
  typedef Geometry<MODEL>            Geometry_;
  typedef Model<MODEL>               Model_;
  typedef ModelAuxControl<MODEL>     ModelAux_;
  typedef ObsAuxControls<MODEL>      ObsAuxCtrls_;
  typedef Observations<MODEL>        Observations_;
  typedef ObsEnsemble<MODEL>         ObsEnsemble_;
  typedef ObsSpaces<MODEL>           ObsSpace_;
  typedef State<MODEL>               State_;

 public:
// -----------------------------------------------------------------------------
  EnsHofX() {
    instantiateObsFilterFactory<MODEL>();
  }
// -----------------------------------------------------------------------------
  virtual ~EnsHofX() {}
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

//  Setup observations
    eckit::LocalConfiguration obsconf(fullConfig, "Observations");
    Log::debug() << "Observations configuration is:" << obsconf << std::endl;
    ObsSpace_ obsdb(obsconf, winbgn, winend);

//  Setup observation bias
    ObsAuxCtrls_ ybias(obsconf);

//  Setup initial states
    const eckit::LocalConfiguration initialConfig(fullConfig, "Initial Condition");
    std::vector<eckit::LocalConfiguration> members;
    initialConfig.get("state", members);
    Log::debug() << "EnsHofX: using " << members.size() << " states." << std::endl;

//  Setup ObsEnsemble
    ObsEnsemble_ obsens(obsdb, members.size());
//  Loop on all ensemble members
    for (unsigned jj = 0; jj < members.size(); ++jj) {
//    Setup initial state for jj-th member
      Log::info() << "Initial configuration for member " << jj << " is :" << members[jj]
                  << std::endl;
      State_ xx(resol, model.variables(), members[jj]);
      Log::test() << "Initial state for member " << jj << ":" << xx << std::endl;

//    Setup augmented state
      ModelAux_ moderr(resol, members[jj]);

//    Setup postprocessor: forecast outputs
      PostProcessor<State_> post;
      eckit::LocalConfiguration prtConf;
      fullConfig.get("Prints", prtConf);
      post.enrollProcessor(new StateInfo<State_>("fc", prtConf));

//    Setup postprocessor: Observers
      boost::shared_ptr<Observers<MODEL, State_> >
      pobs(new Observers<MODEL, State_>(obsconf, obsdb, ybias));
      post.enrollProcessor(pobs);

//    Compute H(x)
      model.forecast(xx, moderr, winlen, post);
      Log::info() << "Finished observation computation for member " << jj << ":" << std::endl;
      Log::test() << "Final state for member " << jj << ":" << xx << std::endl;

//    Save H(x)
      boost::scoped_ptr<Observations_> yobs(pobs->release());
      Log::test() << "H(x) for member " << jj << ":" << std::endl << *yobs
                  << "End H(x) for member " << jj << std::endl;;
      obsens[jj] = *yobs;
      obsens[jj].save("hofx_"+std::to_string(jj+1));
    }

    return 0;
  }
// -----------------------------------------------------------------------------
 private:
  std::string appname() const {
    return "oops::EnsHofX<" + MODEL::name() + ">";
  }
// -----------------------------------------------------------------------------
};

}  // namespace oops

#endif  // OOPS_RUNS_ENSHOFX_H_
