/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_RUNS_ENSFORECASTS_H_
#define OOPS_RUNS_ENSFORECASTS_H_

#include <string>
#include <vector>

#include "eckit/config/Configuration.h"
#include "oops/base/PostProcessor.h"
#include "oops/base/StateInfo.h"
#include "oops/base/StateWriter.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Model.h"
#include "oops/interface/ModelAuxControl.h"
#include "oops/interface/State.h"
#include "oops/runs/Application.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

namespace oops {

template <typename MODEL> class EnsForecast : public Application {
  typedef Geometry<MODEL>         Geometry_;
  typedef Model<MODEL>            Model_;
  typedef ModelAuxControl<MODEL>  ModelAux_;
  typedef State<MODEL>            State_;

 public:
// -----------------------------------------------------------------------------
  EnsForecast() {}
// -----------------------------------------------------------------------------
  virtual ~EnsForecast() {}
// -----------------------------------------------------------------------------
  int execute(const eckit::Configuration & fullConfig) const {
//  Setup resolution
    const eckit::LocalConfiguration resolConfig(fullConfig, "resolution");
    const Geometry_ resol(resolConfig);

//  Setup Model
    const eckit::LocalConfiguration modelConfig(fullConfig, "model");
    const Model_ model(resol, modelConfig);

    std::vector<eckit::LocalConfiguration> memberConf;
    fullConfig.get("members", memberConf);
    int members = memberConf.size();

    for (int jm = 0; jm < members; ++jm) {
//    Setup initial state
      State_ xx(resol, model.variables(), memberConf[jm]);
      Log::test() << "Initial state " << jm << " : " << xx << std::endl;

//    Setup augmented state
      const ModelAux_ moderr(resol, memberConf[jm]);

//    Setup times
      const util::Duration fclength(fullConfig.getString("forecast_length"));
      const util::DateTime bgndate(xx.validTime());
      const util::DateTime enddate(bgndate + fclength);
      Log::info() << "Running forecast " << jm << " from " << bgndate
                  << " to " << enddate << std::endl;

//    Setup forecast outputs
      PostProcessor<State_> post;

      eckit::LocalConfiguration prtConfig;
      if (fullConfig.has("prints")) {
        prtConfig = eckit::LocalConfiguration(fullConfig, "prints");
      }
      post.enrollProcessor(new StateInfo<State_>("fc", prtConfig));

      eckit::LocalConfiguration outConfig(fullConfig, "output");
      outConfig.set("member", jm+1);

      post.enrollProcessor(new StateWriter<State_>(bgndate, outConfig));

//    Run forecast
      model.forecast(xx, moderr, fclength, post);

      Log::test() << "Final state " << jm << " : " << xx << std::endl;
    }
    return 0;
  }
// -----------------------------------------------------------------------------
 private:
  std::string appname() const {
    return "oops::EnsForecast<" + MODEL::name() + ">";
  }
// -----------------------------------------------------------------------------
};

}  // namespace oops
#endif  // OOPS_RUNS_ENSFORECASTS_H_
