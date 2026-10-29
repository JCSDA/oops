/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Copyright 2019-2026 UCAR
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

#include "eckit/config/LocalConfiguration.h"
#include "oops/assimilation/ControlIncrement.h"
#include "oops/assimilation/ControlVariable.h"
#include "oops/assimilation/CostFunction.h"
#include "oops/assimilation/Minimizer.h"
#include "oops/base/IncrementWriter.h"
#include "oops/base/PostProcessor.h"
#include "oops/base/State.h"
#include "oops/base/StateInfo.h"
#include "oops/base/StateWriter.h"
#include "oops/base/StructuredIncrementWriter.h"
#include "oops/base/StructuredStateWriter.h"
#include "oops/util/ConfigFunctions.h"
#include "oops/util/Logger.h"
#include "oops/util/printRunStats.h"

namespace oops {

template<typename MODEL, typename OBS>
int IncrementalAssimilation(ControlVariable<MODEL, OBS> & xx, CostFunction<MODEL, OBS> & J,
                            const eckit::Configuration & config,
                            const eckit::Configuration & outConf,
                            const bool doFinalEvaluate = true) {
  typedef ControlIncrement<MODEL, OBS>     CtrlInc_;
  typedef Minimizer<MODEL, OBS>            Minimizer_;
  typedef State<MODEL>                     State_;
  typedef IncrementWriter<MODEL>           IncrementWriter_;
  typedef StructuredIncrementWriter<MODEL> StructuredIncrementWriter_;

  util::printRunStats("IncrementalAssimilation start");

// Setup outer loop
  std::vector<eckit::LocalConfiguration> iterconfs;
  config.get("iterations", iterconfs);
  const unsigned int nouter = iterconfs.size();
  Log::info() << "Running incremental assimilation with " << nouter
              << " outer iterations." << std::endl;

// Setup minimizer
  eckit::LocalConfiguration minConf(config, "minimizer");
  minConf.set("nouter", static_cast<int>(nouter));
  std::unique_ptr<Minimizer_> minim(MinFactory<MODEL, OBS>::create(minConf, J));

// Setup increment outputs
  PostProcessor<State_> ppInc;
  ASSERT(!outConf.has("increment to structured grid"));  // TODO(yt): remove
  if (doFinalEvaluate && outConf.has("increment")) {
    eckit::LocalConfiguration incConfig(outConf, "increment");
    std::vector<eckit::LocalConfiguration> incouts;
    if (util::isVector(incConfig)) {
      incouts = incConfig.getSubConfigurations();
    } else {
      incouts.push_back(incConfig);
    }
    for (eckit::LocalConfiguration & incout : incouts) {
      if (incout.getString("grid type", "native") == "native") {
        if (!incout.has("variables")) {
          incout.set("variables", J.jb().jbState().variables().variables());
        }
        ppInc.enrollProcessor(new IncrementWriter_(incout));
      } else {
        if (!incout.has("geometry")) {  // will be removed when we move to mist
          incout.set("geometry", eckit::LocalConfiguration(iterconfs.back(), "geometry"));
        }
        ppInc.enrollProcessor(new StructuredIncrementWriter_(incout));
      }
    }
  }
  const eckit::LocalConfiguration prtConfig = config.getSubConfiguration("prints");

// Outer Loop
  for (unsigned jouter = 0; jouter < nouter; ++jouter) {
    iterconfs[jouter].set("iteration", static_cast<int>(jouter));
    iterconfs[jouter].set("total iterations", static_cast<int>(nouter));

//  Get configuration for current outer iteration
    Log::info() << "IncrementalAssimilation: Configuration for outer iteration "
                << jouter << ":" << std::endl << iterconfs[jouter] << std::endl;
    util::printRunStats("IncrementalAssimilation iteration " + std::to_string(jouter));

//  Setup for the trajectory run
    PostProcessor<State_> post;
    if (!prtConfig.empty()) post.enrollProcessor(new StateInfo<State_>("traj", prtConfig));
    if (jouter == 0) post.append(ppInc);

//  Append any new obs if they are in the config
    if (iterconfs[jouter].has("continuous da")) {
      eckit::LocalConfiguration CdaConf(iterconfs[jouter], "continuous da");
      util::DateTime xxTimePrev(xx.state(0).validTime());
      J.applyContDaUpdate(xx, CdaConf);
      if (xx.state(0).validTime() != xxTimePrev) {
        minim.reset(MinFactory<MODEL, OBS>::create(minConf, J));
      }
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
  util::printRunStats("IncrementalAssimilation end outer loop");

// Outputs and Diagnostics
  if (doFinalEvaluate) {
    PostProcessor<State_> post;
    if (!prtConfig.empty()) post.enrollProcessor(new StateInfo<State_>("final", prtConfig));
    post.append(ppInc);

    ASSERT(!outConf.has("analysis to structured grid"));  // TODO(yt): remove
//  Configure writing of analysis
    if (outConf.has("analysis")) {
      const eckit::LocalConfiguration anConf(outConf, "analysis");
      std::vector<eckit::LocalConfiguration> anouts;
      if (util::isVector(anConf)) {
        anouts = anConf.getSubConfigurations();
      } else {
        anouts.push_back(anConf);
      }
      for (const eckit::Configuration & anout : anouts) {
        if (anout.getString("grid type", "native") == "native") {
          post.enrollProcessor(new StateWriter<State_>(anout));
        } else {
          post.enrollProcessor(new StructuredStateWriter<MODEL>(anout));
        }
      }
    }

    if (outConf.getBool("do not evaluate o minus a", false)) {
//    Hacky way to save analysis and increments (at initial time if 4D)
//    *** O-A will not be produced, not recommended ***
      Log::warning() << "NOT evaluating O-A" << std::endl;
      post.initialize(xx.state(0), xx.state(0).validTime(), util::Duration(0));
      post.process(xx.state(0));
    } else {
      eckit::LocalConfiguration finalConfig;
      finalConfig.set("iteration", static_cast<int>(nouter));
      finalConfig.set("total iterations", static_cast<int>(nouter));

//    Produce outputs and diagnostics (including o-a)
      J.evaluate(xx, finalConfig, post);
    }
  }

  util::printRunStats("IncrementalAssimilation end");
  return nouter;
}

}  // namespace oops
#endif  // OOPS_ASSIMILATION_INCREMENTALASSIMILATION_H_
