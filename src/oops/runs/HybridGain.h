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
#include "oops/base/Variables.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
#include "oops/interface/State.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
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

    // Get control state
    const eckit::LocalConfiguration bkgConfig(fullConfig, "control");
    State_ x_control(resol, bkgConfig);
    Log::test() << "Control prior: " << std::endl << x_control << std::endl;
    const Variables vars = x_control.variables();

    // Get ens mean
    const eckit::LocalConfiguration emeanConfig(fullConfig, "ensemble mean");
    State_ x_emean(resol, emeanConfig);
    Log::test() << "Ensemble mean: " << std::endl << x_emean << std::endl;

    // Get ensemble configuration
    std::vector<eckit::LocalConfiguration> ensConfig;
    fullConfig.get("ensemble", ensConfig);
    unsigned nens = ensConfig.size();

    // Compute new center and save
    State_ x_new_center(resol, vars, x_control.validTime());
    x_new_center.zero();
    x_new_center.accumul(alphaControl, x_control);
    x_new_center.accumul(alphaEnsemble, x_emean);
    eckit::LocalConfiguration centerOut(fullConfig, "recentered output");
    centerOut.set("member", static_cast<int>(0) );
    x_new_center.write(centerOut);
    Log::test() << "new center : " << x_new_center << std::endl;

    // Recenter ensemble around new centr and save
    for (unsigned jj = 0; jj < nens; ++jj) {
      State_ x(resol, ensConfig[jj]);
      Increment_ pert(resol, vars, x.validTime());
      pert.diff(x, x_emean);
      x = x_new_center;
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
