/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Crown Copyright 2023, the Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_ASSIMILATION_INCREMENTALASSIMILATION_H_
#define OOPS_ASSIMILATION_INCREMENTALASSIMILATION_H_

#include <memory>
#include <vector>

#include "eckit/config/Configuration.h"
#include "oops/assimilation/ControlIncrement.h"
#include "oops/assimilation/ControlVariable.h"
#include "oops/assimilation/CostFunction.h"
#include "oops/assimilation/Minimizer.h"
#include "oops/base/PostProcessor.h"
#include "oops/base/State.h"
#include "oops/base/StateInfo.h"
#include "oops/util/Logger.h"
#include "oops/util/printRunStats.h"

namespace oops {

template<typename MODEL, typename OBS>
int IncrementalAssimilation(ControlVariable<MODEL, OBS> & xx, CostFunction<MODEL, OBS> & J,
                            const eckit::Configuration & config) {
  typedef ControlIncrement<MODEL, OBS>    CtrlInc_;
  typedef Minimizer<MODEL, OBS>           Minimizer_;
  typedef State<MODEL>                    State_;

  util::printRunStats("IncrementalAssimilation start");

// Setup outer loop
  std::vector<eckit::LocalConfiguration> iterconfs;
  config.get("iterations", iterconfs);
  const unsigned int nouter = iterconfs.size();
  Log::info() << "Running incremental assimilation with " << nouter
              << " outer iterations." << std::endl;

// Setup minimizer
  eckit::LocalConfiguration minConf(config, "minimizer");
  minConf.set("nouter", static_cast<const int>(nouter));
  std::unique_ptr<Minimizer_> minim(MinFactory<MODEL, OBS>::create(minConf, J));

  for (unsigned jouter = 0; jouter < nouter; ++jouter) {
    iterconfs[jouter].set("iteration", static_cast<int>(jouter));
    iterconfs[jouter].set("total iterations", static_cast<int>(nouter));

//  Get configuration for current outer iteration
    Log::info() << "IncrementalAssimilation: Configuration for outer iteration "
                << jouter << ":" << std::endl << iterconfs[jouter] << std::endl;
    util::printRunStats("IncrementalAssimilation iteration " + std::to_string(jouter));

//  Append any new obs if they are in the config
    if (iterconfs[jouter].has("obs append directory")) {
      J.appendObs(iterconfs[jouter]);
    }

//  Setup for the trajectory run
    PostProcessor<State_> post;
    if (iterconfs[jouter].has("prints")) {
      const eckit::LocalConfiguration prtConfig(iterconfs[jouter], "prints");
      post.enrollProcessor(new StateInfo<State_>("traj", prtConfig));
    }

//  Evaluate cost function and setup quadratic problem
    iterconfs[jouter].set("linearize", true);
    J.evaluate(xx, iterconfs[jouter], post);
    util::printRunStats("IncrementalAssimilation linearize " + std::to_string(jouter));

//  Minimization
    std::unique_ptr<CtrlInc_> dx(minim->minimize(iterconfs[jouter]));

//  Compute analysis in physical space
    J.addIncrement(xx, *dx);

//  Clean-up trajectory, etc...
    J.resetLinearization();
  }
  util::printRunStats("IncrementalAssimilation end");
  return nouter;
}

}  // namespace oops
#endif  // OOPS_ASSIMILATION_INCREMENTALASSIMILATION_H_
