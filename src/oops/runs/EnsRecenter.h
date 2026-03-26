/*
 * (C) Copyright 2019-2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_RUNS_ENSRECENTER_H_
#define OOPS_RUNS_ENSRECENTER_H_

#include <memory>
#include <string>
#include <vector>

#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/State.h"
#include "oops/base/StateSet.h"
#include "oops/base/Variables.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/ConfigHelpers.h"
#include "oops/util/DateTime.h"

namespace oops {

template <typename MODEL> class EnsRecenter : public Application {
  typedef Geometry<MODEL>   Geometry_;
  typedef Increment<MODEL>  Increment_;
  typedef State<MODEL>      State_;
  typedef StateSet<MODEL>      StateSet_;

 public:
  // -----------------------------------------------------------------------------
  explicit EnsRecenter(const eckit::mpi::Comm & comm = oops::mpi::world()) : Application(comm) {}
  // -----------------------------------------------------------------------------
  virtual ~EnsRecenter() {}
  // -----------------------------------------------------------------------------
  int execute(const eckit::Configuration & fullConfig) const override {
    // Setup Geometry
    const Geometry_ resol(eckit::LocalConfiguration(fullConfig, "geometry"), this->getComm());

    // Get central state
    State_ x_center(resol, eckit::LocalConfiguration(fullConfig, "center"));

    // Optionally zero the center
    if (fullConfig.getBool("zero center", false)) {
      x_center.zero();
    }

    // Compute ensemble mean using StateSet
    eckit::LocalConfiguration ensConf(fullConfig, "ensemble");
    const StateSet_ ensemble(resol, ensConf);
    const StateSet_ ensmean = ensemble.ens_mean();
    Log::test() << "Ensemble mean: " << std::endl << ensmean << std::endl;

    // Optionally write the mean out
    if (fullConfig.has("ensemble meanoutput")) {
      ensmean.write(eckit::LocalConfiguration(fullConfig, "ensemble meanoutput"));
    }

    // Recenter ensemble around central and save
    Variables vars(fullConfig, "recenter variables");
    for (unsigned jj = 0; jj < ensemble.ens_size(); ++jj) {
      State_ x(resol, ensemble(0, jj));
      Increment_ pert(resol, vars, x.validTime());
      pert.diff(x, ensmean[0]);
      x = x_center;
      x += pert;

      // Save recentered member
      eckit::LocalConfiguration recenteredOutput(fullConfig, "recentered output");
      util::setMember(recenteredOutput, jj+1);
      x.write(recenteredOutput);
      Log::test() << "Recentered member " << jj << " : " << x << std::endl;
    }

    return 0;
  }
// -----------------------------------------------------------------------------
 private:
  std::string appname() const override {
    return "oops::EnsRecenter<" + MODEL::name() + ">";
  }
  // -----------------------------------------------------------------------------
};

}  // namespace oops

#endif  // OOPS_RUNS_ENSRECENTER_H_
