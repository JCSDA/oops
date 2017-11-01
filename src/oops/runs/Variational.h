/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_RUNS_VARIATIONAL_H_
#define OOPS_RUNS_VARIATIONAL_H_

#include <string>

#include <boost/scoped_ptr.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "util/Logger.h"
#include "oops/assimilation/ControlVariable.h"
#include "oops/assimilation/CostFunction.h"
#include "oops/assimilation/IncrementalAssimilation.h"
#include "oops/assimilation/instantiateCostFactory.h"
#include "oops/assimilation/instantiateMinFactory.h"
#include "oops/base/Observations.h"
#include "oops/base/PostProcessor.h"
#include "oops/base/StateInfo.h"
#include "oops/base/StateWriter.h"
#include "oops/base/instantiateCovarFactory.h"
#include "oops/generic/instantiateObsErrorFactory.h"
#include "oops/generic/instantiateTlmFactory.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Model.h"
#include "oops/interface/State.h"
#include "oops/runs/Application.h"
#include "util/DateTime.h"
#include "util/Duration.h"

namespace oops {

template <typename MODEL> class Variational : public Application {
  typedef Geometry<MODEL>            Geometry_;
  typedef Model<MODEL>               Model_;
  typedef State<MODEL>               State_;

 public:
// -----------------------------------------------------------------------------
  Variational() {
    instantiateCostFactory<MODEL>();
    instantiateCovarFactory<MODEL>();
    instantiateMinFactory<MODEL>();
    instantiateObsErrorFactory<MODEL>();
    instantiateTlmFactory<MODEL>();
  }
// -----------------------------------------------------------------------------
  virtual ~Variational() {}
// -----------------------------------------------------------------------------
  int execute(const eckit::Configuration & fullConfig) const {
//  Setup resolution
    const eckit::LocalConfiguration resolConfig(fullConfig, "resolution");
    const Geometry_ resol(resolConfig);

//  Setup Model
    const eckit::LocalConfiguration modelConfig(fullConfig, "model");
    const Model_ model(resol, modelConfig);
    Log::trace() << "Variational: model has been set up" << std::endl;

/// The background is constructed inside the cost function because its valid
/// time within the assimilation window can be different (3D-Var vs. 4D-Var),
/// it can be 3D or 4D (strong vs weak constraint), etc...

//  Setup cost function
    const eckit::LocalConfiguration cfConf(fullConfig, "cost_function");
    boost::scoped_ptr< CostFunction<MODEL> > J(CostFactory<MODEL>::create(cfConf, resol, model));
    Log::trace() << "Variational: cost function has been set up" << std::endl;

//  Initialize first guess from background
    ControlVariable<MODEL> xx(J->jb().getBackground());
    Log::trace() << "Variational: first guess has been set up" << std::endl;

//  Perform Incremental Variational Assimilation
    IncrementalAssimilation<MODEL>(xx, *J, fullConfig);
    Log::trace() << "Variational: incremantal assimilation done" << std::endl;

//  Save analysis and final diagnostics
    PostProcessor<State_> post;
    const util::DateTime winbgn(cfConf.getString("window_begin"));
    const eckit::LocalConfiguration outConfig(fullConfig, "output");
    post.enrollProcessor(new StateWriter<State_>(winbgn, outConfig));

    const eckit::LocalConfiguration finalConfig(fullConfig, "final");
    if (finalConfig.has("prints")) {
      const eckit::LocalConfiguration prtConfig(finalConfig, "prints");
      post.enrollProcessor(new StateInfo<State_>("final", prtConfig));
    }

    J->evaluate(xx, finalConfig, post);

    return 0;
  }
// -----------------------------------------------------------------------------
 private:
  std::string appname() const {
    return "oops::Variational<" + MODEL::name() + ">";
  }
// -----------------------------------------------------------------------------
};

}  // namespace oops
#endif  // OOPS_RUNS_VARIATIONAL_H_
