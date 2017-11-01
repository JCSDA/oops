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
#include "util/Logger.h"
#include "oops/base/Observations.h"
#include "oops/base/Observer.h"
#include "oops/base/PostProcessor.h"
#include "oops/base/StateInfo.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Model.h"
#include "oops/interface/ModelAuxControl.h"
#include "oops/interface/ObsAuxControl.h"
#include "oops/interface/ObservationSpace.h"
#include "oops/interface/ObsOperator.h"
#include "oops/interface/State.h"
#include "oops/runs/Application.h"
#include "util/DateTime.h"
#include "util/Duration.h"

namespace oops {

template <typename MODEL> class HofX : public Application {
  typedef Geometry<MODEL>            Geometry_;
  typedef Model<MODEL>               Model_;
  typedef ModelAuxControl<MODEL>     ModelAux_;
  typedef ObsAuxControl<MODEL>       ObsAuxCtrl_;
  typedef Observations<MODEL>        Observations_;
  typedef ObsOperator<MODEL>         ObsOperator_;
  typedef ObservationSpace<MODEL>    ObsSpace_;
  typedef State<MODEL>               State_;

 public:
// -----------------------------------------------------------------------------
  HofX() {}
// -----------------------------------------------------------------------------
  virtual ~HofX() {}
// -----------------------------------------------------------------------------
  int execute(const eckit::Configuration & fullConfig) const {
//  Setup observation window
    const eckit::LocalConfiguration windowConf(fullConfig, "assimilation_window");
    const util::Duration winlen(windowConf.getString("window_length"));
    const util::DateTime winbgn(windowConf.getString("window_begin"));
    const util::DateTime winend(winbgn + winlen);
    Log::info() << "Observation window is:" << windowConf << std::endl;

//  Setup resolution
    const eckit::LocalConfiguration resolConfig(fullConfig, "resolution");
    const Geometry_ resol(resolConfig);

//  Setup Model
    const eckit::LocalConfiguration modelConfig(fullConfig, "model");
    const Model_ model(resol, modelConfig);

//  Setup initial state
    const eckit::LocalConfiguration initialConfig(fullConfig, "initial");
    Log::info() << "Initial configuration is:" << initialConfig << std::endl;
    State_ xx(resol, initialConfig);
    Log::test() << "Initial state: " << xx.norm() << std::endl;

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
    std::vector<boost::shared_ptr<Observer<MODEL, State_> > > pobs;

    std::vector<eckit::LocalConfiguration> obsconf;
    fullConfig.get("Observation", obsconf);
    size_t nobs = obsconf.size();

    for (size_t jobs = 0; jobs < nobs; ++jobs) {
      Log::debug() << "Observation configuration is:" << obsconf[jobs] << std::endl;
      ObsSpace_ obsdb(obsconf[jobs], winbgn, winend);
      ObsOperator_ hop(obsdb, obsconf[jobs]);
      Observations_ yy(obsdb);

      boost::shared_ptr<Observer<MODEL, State_> >
        pp(new Observer<MODEL, State_>(obsdb, hop, yy, ybias));
      post.enrollProcessor(pp);
      pobs.push_back(pp);
    }

//  Compute H(x)
    model.forecast(xx, moderr, winlen, post);
    Log::info() << "HofX: Finished observation computation." << std::endl;
    Log::test() << "Final state: " << xx.norm() << std::endl;

//  Save H(x)
    for (size_t jobs = 0; jobs < nobs; ++jobs) {
      boost::scoped_ptr<Observations_> yobs(pobs[jobs]->release());
      Log::test() << "H(x): " << *yobs << std::endl;
      yobs->save("hofx");
    }

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
