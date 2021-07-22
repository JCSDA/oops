/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_RUNS_ENSRECENTER_H_
#define OOPS_RUNS_ENSRECENTER_H_

#include <memory>
#include <string>
#include <vector>


#include "eckit/config/LocalConfiguration.h"
#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/Variables.h"
#include "oops/interface/State.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/DateTime.h"

namespace oops {

template <typename MODEL> class EnsRecenter : public Application {
  typedef Geometry<MODEL>   Geometry_;
  typedef Increment<MODEL>  Increment_;
  typedef State<MODEL>      State_;

 public:
  // -----------------------------------------------------------------------------
  explicit EnsRecenter(const eckit::mpi::Comm & comm = oops::mpi::world()) : Application(comm) {}
  // -----------------------------------------------------------------------------
  virtual ~EnsRecenter() {}
  // -----------------------------------------------------------------------------
  int execute(const eckit::Configuration & fullConfig) const {
    // Setup Geometry
    const eckit::LocalConfiguration resolConfig(fullConfig, "geometry");
    const Geometry_ resol(resolConfig, this->getComm());

    // Get central state
    const eckit::LocalConfiguration bkgConfig(fullConfig, "center");
    State_ x_center(resol, bkgConfig);

    // Optionally zero the center
    bool zeroCenter = fullConfig.getBool("zero center", false);
    if (zeroCenter) {
      x_center.zero();
    }

    // Get ensemble configuration
    std::vector<eckit::LocalConfiguration> ensConfig;
    fullConfig.get("ensemble", ensConfig);

    // Get ensemble size
    unsigned nm = ensConfig.size();

    // Compute ensemble mean
    State_ ensmean(x_center);
    ensmean.zero();
    const double rk = 1.0/(static_cast<double>(nm));
    for (unsigned jj = 0; jj < nm; ++jj) {
      State_ x(resol, ensConfig[jj]);
      ensmean.accumul(rk, x);
      Log::test() << "Original member " << jj << " : " << x << std::endl;
    }
    Log::test() << "Ensemble mean: " << std::endl << ensmean << std::endl;

    // Optionally write the mean out
    if (fullConfig.has("ensemble mean output")) {
      eckit::LocalConfiguration meanout(fullConfig, "ensemble mean output");
      ensmean.write(meanout);
    }

    // Setup variables
    const Variables vars(fullConfig, "recenter variables");

    // Recenter ensemble around central and save
    for (unsigned jj = 0; jj < nm; ++jj) {
      State_ x(resol, ensConfig[jj]);
      Increment_ pert(resol, vars, x.validTime());
      pert.diff(x, ensmean);
      x = x_center;
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
    return "oops::EnsRecenter<" + MODEL::name() + ">";
  }
  // -----------------------------------------------------------------------------
};

}  // namespace oops

#endif  // OOPS_RUNS_ENSRECENTER_H_
