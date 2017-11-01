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

#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/shared_ptr.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "util/Logger.h"
#include "oops/base/Departures.h"
#include "oops/base/Observations.h"
#include "oops/base/Observer.h"
#include "oops/base/PostProcessor.h"
#include "oops/base/StateInfo.h"
#include "oops/generic/instantiateObsErrorFactory.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Model.h"
#include "oops/interface/ModelAuxControl.h"
#include "oops/interface/ObsAuxControl.h"
#include "oops/interface/ObservationSpace.h"
#include "oops/interface/ObsErrorCovariance.h"
#include "oops/interface/ObsOperator.h"
#include "oops/interface/State.h"
#include "oops/runs/Application.h"
#include "util/DateTime.h"
#include "util/Duration.h"

namespace oops {

template <typename MODEL> class MakeObs : public Application {
  typedef Departures<MODEL>          Departures_;
  typedef Geometry<MODEL>            Geometry_;
  typedef Model<MODEL>               Model_;
  typedef ModelAuxControl<MODEL>     ModelAux_;
  typedef ObsAuxControl<MODEL>       ObsAuxCtrl_;
  typedef Observations<MODEL>        Observations_;
  typedef ObservationSpace<MODEL>    ObsSpace_;
  typedef ObsOperator<MODEL>         ObsOperator_;
  typedef State<MODEL>               State_;

 public:
// -----------------------------------------------------------------------------
  MakeObs() {
    instantiateObsErrorFactory<MODEL>();
  }
// -----------------------------------------------------------------------------
  virtual ~MakeObs() {}
// -----------------------------------------------------------------------------
  int execute(const eckit::Configuration & fullConfig) const {
//  Setup observation window
    const eckit::LocalConfiguration windowConf(fullConfig, "assimilation_window");
    const util::DateTime bgn(windowConf.getString("begin"));
    const util::DateTime end(windowConf.getString("end"));
    const util::Duration fclen(end - bgn);
    Log::info() << "Observation window is:" << windowConf << std::endl;

//  Setup resolution
    const eckit::LocalConfiguration resolConfig(fullConfig, "resolution");
    const Geometry_ resol(resolConfig);

//  Setup Model
    const eckit::LocalConfiguration modelConfig(fullConfig, "model");
    const Model_ model(resol, modelConfig);

//  Setup initial "true" state
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
    boost::ptr_vector<ObsSpace_> obspaces;

    std::vector<eckit::LocalConfiguration> obsconf;
    fullConfig.get("Observation", obsconf);
    size_t nobs = obsconf.size();

    for (size_t jobs = 0; jobs < nobs; ++jobs) {
      Log::info() << "Observation configuration is:" << obsconf[jobs] << std::endl;
      obspaces.push_back(new ObsSpace_(obsconf[jobs], bgn, end));
      const eckit::LocalConfiguration genConf(obsconf[jobs], "Generate");
      obspaces.at(jobs).generateDistribution(genConf);

      ObsOperator_ hop(obspaces.at(jobs), obsconf[jobs]);
      hop.generateObsError(genConf);

      Observations_ yy(obspaces.at(jobs));

      boost::shared_ptr<Observer<MODEL, State_> >
        pp(new Observer<MODEL, State_>(obspaces.at(jobs), hop, yy, ybias));
      post.enrollProcessor(pp);
      pobs.push_back(pp);
    }

//  Run forecast and generate observations
    model.forecast(xx, moderr, fclen, post);
    Log::info() << "MakeObs: Finished observation generation." << std::endl;
    Log::test() << "Final state: " << xx.norm() << std::endl;

    for (size_t jobs = 0; jobs < nobs; ++jobs) {
      boost::scoped_ptr<Observations_> yobs(pobs[jobs]->release());

//    Perturb observations
      if (obsconf[jobs].has("obspert")) {
        ObsErrorCovariance<MODEL> matR(obspaces.at(jobs), obsconf[jobs]);
        Departures_ ypert(obspaces.at(jobs));
        matR.randomize(ypert);
        double opert = obsconf[jobs].getDouble("obspert");
        ypert *= opert;
        *yobs += ypert;
      }

//    Save observations
      Log::test() << "Generated observation: " << *yobs << std::endl;
      yobs->save("ObsVal");
    }

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
