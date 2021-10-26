/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_RUNS_HYBRIDGAIN_H_
#define OOPS_RUNS_HYBRIDGAIN_H_

#include <memory>
#include <string>
#include <vector>


#include "eckit/config/LocalConfiguration.h"
#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/State.h"
#include "oops/base/Variables.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/DateTime.h"

namespace oops {

template <typename MODEL> class HybridGain : public Application {
  typedef Geometry<MODEL>   Geometry_;
  typedef Increment<MODEL>  Increment_;
  typedef State<MODEL>      State_;

 public:
  // -----------------------------------------------------------------------------
  explicit HybridGain(const eckit::mpi::Comm & comm = oops::mpi::world()) : Application(comm) {}
  // -----------------------------------------------------------------------------
  virtual ~HybridGain() {}
  // -----------------------------------------------------------------------------
  int execute(const eckit::Configuration & fullConfig) const {
    // Setup Geometry
    const eckit::LocalConfiguration resolConfig(fullConfig, "geometry");
    const Geometry_ resol(resolConfig, this->getComm());

    // Read averaging weights
    double alphaControl = fullConfig.getDouble("hybrid weights.control");
    double alphaEnsemble = fullConfig.getDouble("hybrid weights.ensemble");

    // Read hybrid type
    std::string hybridType = fullConfig.getString("hybrid type", "average analysis");

    // Get control state
    const eckit::LocalConfiguration bkgConfig(fullConfig, "control");
    State_ xaControl(resol, bkgConfig);
    Log::test() << "Control: " << std::endl << xaControl << std::endl;
    const Variables vars = xaControl.variables();

    // Get posterior ens mean
    const eckit::LocalConfiguration emeanConfig(fullConfig, "ensemble mean posterior");
    State_ xaEmeanPost(resol, emeanConfig);
    Log::test() << "Ensemble mean posterior: " << std::endl << xaEmeanPost << std::endl;

    // Compute new center
    State_ xNewCenter(resol, vars, xaControl.validTime());
    if (hybridType == "average analysis") {
        // using average of analysis (following Bonavita)
        // xa_hybrid = a1*xa1 + a2*xa2
        // a1+a2 have to equal to 1

        ASSERT(alphaControl + alphaEnsemble == 1.0);
        xNewCenter.zero();
        xNewCenter.accumul(alphaControl, xaControl);
        xNewCenter.accumul(alphaEnsemble, xaEmeanPost);
    } else if (hybridType == "average increment") {
        // using average of analysis incerments (following Whitaker)
        // xa_hybrid = xf_prior + a1*xinc1 + a2*xinc2
        // Note: a1+a2 no longer need to add to one

        // Get prior ens mean
        const eckit::LocalConfiguration emeanConfigPrior(fullConfig, "ensemble mean prior");
        State_ xfEmeanPrior(resol, emeanConfigPrior);
        Log::test() << "Ensemble mean prior: " << std::endl << xfEmeanPrior << std::endl;
        // compute ensemble mean increment
        Increment_ pertEns(resol, vars, xaControl.validTime());
        pertEns.diff(xaEmeanPost, xfEmeanPrior);
        pertEns *= alphaEnsemble;
        // compute control increment
        Increment_ pertControl(resol, vars, xaControl.validTime());
        pertControl.diff(xaControl, xfEmeanPrior);
        pertControl *= alphaControl;
        // compute hybrid posterior
        xNewCenter = xfEmeanPrior;
        xNewCenter += pertEns;
        xNewCenter += pertControl;
    } else {
      ABORT("Unknown hybrid gain type: " + hybridType);
    }
    // Output new center
    eckit::LocalConfiguration centerOut(fullConfig, "recentered output");
    centerOut.set("member", static_cast<int>(0) );
    xNewCenter.write(centerOut);
    Log::test() << "new center : " << xNewCenter << std::endl;

    // Get ensemble configuration
    std::vector<eckit::LocalConfiguration> ensConfig;
    fullConfig.get("ensemble", ensConfig);
    unsigned nens = ensConfig.size();

    // Recenter ensemble around new center and save
    for (unsigned jj = 0; jj < nens; ++jj) {
      State_ x(resol, ensConfig[jj]);
      Increment_ pert(resol, vars, x.validTime());
      pert.diff(x, xaEmeanPost);
      x = xNewCenter;
      x += pert;

      // Save recentered member
      eckit::LocalConfiguration recenterout(fullConfig, "recentered output");
      recenterout.set("member", static_cast<int>(jj+1) );
      x.write(recenterout);
      Log::test() << "Recentered member " << jj << " : " << x << std::endl;
    }

    return 0;
  }
  // -----------------------------------------------------------------------------
 private:
  std::string appname() const {
    return "oops::HybridGain<" + MODEL::name() + ">";
  }
  // -----------------------------------------------------------------------------
};

}  // namespace oops

#endif  // OOPS_RUNS_HYBRIDGAIN_H_
