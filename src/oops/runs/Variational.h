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

#include <memory>
#include <string>


#include "eckit/config/LocalConfiguration.h"
#include "oops/assimilation/ControlVariable.h"
#include "oops/assimilation/CostFunction.h"
#include "oops/assimilation/IncrementalAssimilation.h"
#include "oops/assimilation/instantiateCostFactory.h"
#include "oops/assimilation/instantiateMinFactory.h"
#include "oops/base/instantiateCovarFactory.h"
#include "oops/base/instantiateObsFilterFactory.h"
#include "oops/base/PostProcessor.h"
#include "oops/base/State.h"
#include "oops/base/StateInfo.h"
#include "oops/base/StateWriter.h"
#include "oops/generic/instantiateLinearModelFactory.h"
#include "oops/generic/instantiateObsErrorFactory.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"
#include "oops/util/printRunStats.h"

namespace oops {

template <typename MODEL, typename OBS> class Variational : public Application {
  typedef State<MODEL>               State_;

 public:
// -----------------------------------------------------------------------------
  explicit Variational(const eckit::mpi::Comm & comm = oops::mpi::world()) : Application(comm) {
    instantiateCostFactory<MODEL, OBS>();
    instantiateCovarFactory<MODEL>();
    instantiateMinFactory<MODEL, OBS>();
    instantiateObsErrorFactory<OBS>();
    instantiateObsFilterFactory<OBS>();
    instantiateLinearModelFactory<MODEL>();
  }
// -----------------------------------------------------------------------------
  virtual ~Variational() {}
// -----------------------------------------------------------------------------
  int execute(const eckit::Configuration & fullConfig) const override {
/// The background is constructed inside the cost function because its valid
/// time within the assimilation window can be different (3D-Var vs. 4D-Var),
/// it can be 3D or 4D (strong vs weak constraint), etc...
    util::printRunStats("Variational start");

//  Setup cost function
    eckit::LocalConfiguration cfConf(fullConfig, "cost function");
    CostFunctionParametersWrapper<MODEL, OBS> cfParams;
    cfParams.validateAndDeserialize(cfConf);
    std::unique_ptr<CostFunction<MODEL, OBS>>
      J(CostFactory<MODEL, OBS>::create(cfParams.costTypeParameters, this->getComm()));
    Log::trace() << "Variational: cost function has been set up" << std::endl;

//  Initialize first guess from background
    ControlVariable<MODEL, OBS> xx(J->jb().getBackground());
    Log::trace() << "Variational: first guess has been set up" << std::endl;

//  Perform Incremental Variational Assimilation
    eckit::LocalConfiguration varConf(fullConfig, "variational");
    int iouter = IncrementalAssimilation<MODEL, OBS>(xx, *J, varConf);
    Log::info() << "Variational: incremental assimilation done "
                << iouter << " iterations." << std::endl;

//  Save analysis and final diagnostics
    PostProcessor<State_> post;
    const util::DateTime winbgn(cfConf.getString("window begin"));
    const eckit::LocalConfiguration outConfig(fullConfig, "output");
    post.enrollProcessor(new StateWriter<State_>(outConfig));

    eckit::LocalConfiguration finalConfig(fullConfig, "final");
    finalConfig.set("iteration", iouter);
    if (finalConfig.has("prints")) {
      const eckit::LocalConfiguration prtConfig(finalConfig, "prints");
      post.enrollProcessor(new StateInfo<State_>("final", prtConfig));
    }

    J->evaluate(xx, finalConfig, post);

//  Save ObsAux
    xx.obsVar().write(cfConf);

    util::printRunStats("Variational end");
    return 0;
  }
// -----------------------------------------------------------------------------
 private:
  std::string appname() const override {
    return "oops::Variational<" + MODEL::name() + ", " + OBS::name() + ">";
  }
// -----------------------------------------------------------------------------
};

}  // namespace oops
#endif  // OOPS_RUNS_VARIATIONAL_H_
