/*
 * (C) Copyright 2022 UCAR
 * (C) Crown copyright 2024, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_RUNS_ENSMEANANDVARIANCE_H_
#define OOPS_RUNS_ENSMEANANDVARIANCE_H_

#include <string>

#include "eckit/config/LocalConfiguration.h"
#include "oops/base/Geometry.h"
#include "oops/base/IncrementSet.h"
#include "oops/base/StateSet.h"
#include "oops/base/StructuredGridWriter.h"
#include "oops/base/Variables.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"

namespace oops {

// -----------------------------------------------------------------------------

template <typename MODEL> class EnsMeanAndVariance : public Application {
  typedef Geometry<MODEL>                          Geometry_;
  typedef IncrementSet<MODEL>                      IncrementSet_;
  typedef StateSet<MODEL>                          StateSet_;

 public:
  // -----------------------------------------------------------------------------
  explicit EnsMeanAndVariance(const eckit::mpi::Comm & comm = oops::mpi::world())
    : Application(comm) {}
  // -----------------------------------------------------------------------------
  virtual ~EnsMeanAndVariance() {}
  // -----------------------------------------------------------------------------
  int execute(const eckit::Configuration & fullConfig) const override {
//  Setup Geometry
    const Geometry_ resol(eckit::LocalConfiguration(fullConfig, "geometry"), this->getComm());

//  Setup ensemble of states
    eckit::LocalConfiguration ensConf(fullConfig, "ensemble");
    const StateSet_ stateEnsemble(resol, ensConf);
    const StateSet_ ensmean = stateEnsemble.ens_mean();
    // Convert StateSet to IncrementSet for variance calculation and free stateEnsemble memory
    const IncrementSet_ ensemble(resol, stateEnsemble.variables(), stateEnsemble);
    const IncrementSet_ sigb2 = ensemble.ens_var();

//  Write mean to file
    if (fullConfig.has("mean output")) {
      ensmean.write(eckit::LocalConfiguration(fullConfig, "mean output"));
    }

    if (fullConfig.has("ensmean to structured grid")) {
      const eckit::LocalConfiguration latlonConf(fullConfig, "ensmean to structured grid");
      const StructuredGridWriter<MODEL> latlon(latlonConf, resol);
      for (size_t jt = 0; jt < ensmean.time_size(); ++jt) {
        latlon.interpolateAndWrite(ensmean[jt]);
      }
    }
    Log::test() << "Mean: " << std::endl << ensmean << std::endl;

//  Write variance to file
    if (fullConfig.has("variance output")) {
      sigb2.write(eckit::LocalConfiguration(fullConfig, "variance output"));
    }

    if (fullConfig.has("ensvariance to structured grid")) {
      const eckit::LocalConfiguration latlonConf(fullConfig, "ensvariance to structured grid");
      const StructuredGridWriter<MODEL> latlon(latlonConf, resol);
      for (size_t jt = 0; jt < ensmean.time_size(); ++jt) {
        latlon.interpolateAndWrite(sigb2[jt], ensmean[jt]);
      }
    }
    Log::test() << "Variance: " << std::endl << sigb2 << std::endl;

//  Compute and write standard deviation to file, if it was requested
//  The std dev is computed via a generic atlas algorithm; computing it only when requested allows
//  this executable to work with model interfaces with non-conforming atlas interfaces.
    if (fullConfig.has("standard deviation output")
        || fullConfig.has("standard deviation to structured grid")) {
      const IncrementSet_ sigb = ensemble.ens_stddev();

      if (fullConfig.has("standard deviation output")) {
        sigb.write(eckit::LocalConfiguration(fullConfig, "standard deviation output"));
      }
      if (fullConfig.has("standard deviation to structured grid")) {
        const eckit::LocalConfiguration latlonConf(fullConfig,
                                                   "standard deviation to structured grid");
        const StructuredGridWriter<MODEL> latlon(latlonConf, resol);
        for (size_t jt = 0; jt < sigb.time_size(); ++jt) {
          latlon.interpolateAndWrite(sigb[jt], ensmean[jt]);
        }
      }
      Log::test() << "Standard Deviation: " << std::endl << sigb << std::endl;
    }

    return 0;
  }
// -----------------------------------------------------------------------------
 private:
  std::string appname() const override {
    return "oops::EnsMeanAndVariance<" + MODEL::name() + ">";
  }
// -----------------------------------------------------------------------------
};

}  // namespace oops

#endif  // OOPS_RUNS_ENSMEANANDVARIANCE_H_
