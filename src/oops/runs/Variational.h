/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Crown Copyright 2023, the Met Office
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
#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/instantiateCovarFactory.h"
#include "oops/base/instantiateObsFilterFactory.h"
#include "oops/base/LatLonGridPostProcessor.h"
#include "oops/base/LatLonGridWriter.h"
#include "oops/base/PostProcessor.h"
#include "oops/base/State.h"
#include "oops/base/StateInfo.h"
#include "oops/base/StateWriter.h"
#include "oops/generic/instantiateLinearModelFactory.h"
#include "oops/generic/instantiateNormFactory.h"
#include "oops/generic/instantiateObsErrorFactory.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/printRunStats.h"

namespace oops {

template <typename MODEL, typename OBS> class Variational : public Application {
  typedef Increment<MODEL>           Increment_;
  typedef State<MODEL>               State_;

 public:
// -----------------------------------------------------------------------------
  explicit Variational(const eckit::mpi::Comm & comm = oops::mpi::world()) : Application(comm) {
    instantiateCostFactory<MODEL, OBS>();
    instantiateCovarFactory<MODEL>();
    instantiateMinFactory<MODEL, OBS>();
    instantiateNormFactory<MODEL>();
    instantiateObsErrorFactory<OBS>();
    instantiateObsFilterFactory<OBS>();
    instantiateLinearModelFactory<MODEL>();
  }
// -----------------------------------------------------------------------------
  virtual ~Variational() {}
// -----------------------------------------------------------------------------
  int execute(const eckit::Configuration & fullConfig, bool validate) const override {
    Log::trace() << "Variational: execute start" << std::endl;
    util::printRunStats("Variational start");
/// The background is constructed inside the cost function because its valid
/// time within the assimilation window can be different (3D-Var vs. 4D-Var),
/// it can be 3D or 4D (strong vs weak constraint), etc...

//  Setup cost function
    eckit::LocalConfiguration cfConf(fullConfig, "cost function");
    std::unique_ptr<CostFunction<MODEL, OBS>>
      J(CostFactory<MODEL, OBS>::create(cfConf, this->getComm()));

//  Initialize first guess from background
    ControlVariable<MODEL, OBS> xx(J->jb().getBackground());

//  Perform Incremental Variational Assimilation
    eckit::LocalConfiguration varConf(fullConfig, "variational");
    int iouter = IncrementalAssimilation<MODEL, OBS>(xx, *J, varConf);
    Log::info() << "Variational: incremental assimilation done "
                << iouter << " iterations." << std::endl;

//  Save analysis and final diagnostics
    PostProcessor<State_> post;
    if (fullConfig.has("output")) {
      const eckit::LocalConfiguration outConfig(fullConfig, "output");
      post.enrollProcessor(new StateWriter<State_>(outConfig));
    }

    eckit::LocalConfiguration finalConfig(fullConfig, "final");
    finalConfig.set("iteration", iouter);

//  Save increment if desired
    if (finalConfig.has("increment")) {
      const eckit::LocalConfiguration incConfig(finalConfig, "increment");
      ControlVariable<MODEL, OBS> x_b(J->jb().getBackground());
      const eckit::LocalConfiguration incGeomConfig(incConfig, "geometry");
      Geometry<MODEL> incGeom(incGeomConfig,
                              xx.states().geometry().getComm(),
                              xx.states().commTime());
      ControlIncrement<MODEL, OBS> dx_tmp(J->jb());
      ControlIncrement<MODEL, OBS> dx(incGeom, dx_tmp);
      dx.diff(xx, x_b);
      const eckit::LocalConfiguration incOutConfig(incConfig, "output");
      dx.write(incOutConfig);
    }

    if (finalConfig.has("increment to latlon")) {
      const eckit::LocalConfiguration incLatlonConf(finalConfig, "increment to latlon");

      ControlVariable<MODEL, OBS> x_b(J->jb().getBackground());
      ControlIncrement<MODEL, OBS> dx(J->jb());
      dx.diff(xx, x_b);

      const LatLonGridWriter<MODEL> latlon(incLatlonConf, dx.states().geometry());
      for (size_t jtime = 0; jtime < dx.states().size(); ++jtime) {
        latlon.interpolateAndWrite(dx.states()[jtime], xx.states()[jtime]);
      }
    }

    if (finalConfig.has("prints")) {
      const eckit::LocalConfiguration prtConfig(finalConfig, "prints");
      post.enrollProcessor(new StateInfo<State_>("final", prtConfig));
    }

    if (finalConfig.has("analysis to latlon")) {
      const eckit::LocalConfiguration anLatlonConf(finalConfig, "analysis to latlon");
      post.enrollProcessor(new LatLonGridPostProcessor<MODEL, State_>(
            anLatlonConf, xx.state().geometry() ));
    }

    J->evaluate(xx, finalConfig, post);

//  Save ObsAux
    xx.obsVar().write(cfConf);

    util::printRunStats("Variational end");
    Log::trace() << "Variational: execute done" << std::endl;
    return 0;
  }
// -----------------------------------------------------------------------------
  void validateConfig(const eckit::Configuration & fullConfig) const override {
    // Note: Variational app doesn't have application level Parameters yet;
    // not validating anything.
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
