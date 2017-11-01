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

#include "util/Logger.h"
#include "oops/base/PostProcessor.h"
#include "oops/base/StateInfo.h"
#include "oops/base/StateWriter.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Model.h"
#include "oops/interface/ModelAuxControl.h"
#include "oops/interface/State.h"
#include "oops/runs/Application.h"
#include "eckit/config/Configuration.h"
#include "util/DateTime.h"
#include "util/Duration.h"

namespace oops {

template <typename MODEL> class EnsForecast : public Application {
  typedef Geometry<MODEL>            Geometry_;
  typedef Model<MODEL>               Model_;
  typedef ModelAuxControl<MODEL>      ModelAux_;
  typedef State<MODEL>               State_;

 public:
// -----------------------------------------------------------------------------
  EnsForecast() {}
// -----------------------------------------------------------------------------
  virtual ~EnsForecast() {}
// -----------------------------------------------------------------------------
  int execute(const eckit::Configuration & fullConfig) const {
//  Setup resolution
    const eckit::Configuration resolConfig(fullConfig, "resolution");
    const Geometry_ resol(resolConfig);

//  Setup Model
    const eckit::Configuration modelConfig(fullConfig, "model");
    const Model_ model(resol, modelConfig);

    unsigned nm = config.getElementSize("member");
    for (unsigned jj = 0; jj < nm; ++jj) {
//    Setup initial state
      const eckit::Configuration initialConfig(fullConfig, "member", jj);
      State_ xx(resol, initialConfig);
      Log::test() << "Initial state: " << xx.norm() << std::endl;

//    Setup augmented state
      const ModelAux_ moderr(initialConfig);

//    Setup times
      const util::Duration fclength(fullConfig.getString("forecast_length"));
      const util::DateTime bgndate(xx.validTime());
      const util::DateTime enddate(bgndate + fclength);
      Log::info() << "Running forecast " << jj << " from " << bgndate << " to " << enddate << std::endl;

//    Setup forecast outputs
      PostProcessor<State_> post;

      const eckit::Configuration prtConfig(fullConfig, "prints", true);
      post.enrollProcessor(new StateInfo<State_>("fc", prtConfig));

      const eckit::Configuration outConfig(fullConfig, "output");
      post.enrollProcessor(new StateWriter<State_>(bgndate, outConfig));

//    Run forecast
      model.forecast(xx, moderr, fclength, post);

      Log::test() << "Final state " << jj << " : " << xx.norm() << std::endl;
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
