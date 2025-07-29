/*
 * (C) Copyright 2022 UCAR
 * (C) Crown copyright 2024, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_RUNS_ENSMEANANDVARIANCE_H_
#define OOPS_RUNS_ENSMEANANDVARIANCE_H_

#include <memory>
#include <string>
#include <vector>


#include "eckit/config/LocalConfiguration.h"
#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/PostProcessor.h"
#include "oops/base/State.h"
#include "oops/base/StateEnsemble.h"
#include "oops/base/StructuredGridWriter.h"
#include "oops/base/Variables.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"

namespace oops {

// -----------------------------------------------------------------------------

template <typename MODEL> class EnsMeanAndVariance : public Application {
  typedef Geometry<MODEL>                          Geometry_;
  typedef Increment<MODEL>                         Increment_;
  typedef State<MODEL>                             State_;
  typedef StateEnsemble<MODEL>                     StateEnsemble_;

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
    const StateEnsemble_ stateEnsemble(resol, ensConf);
    const State_ ensmean = stateEnsemble.mean();
    const Increment_ sigb2 = stateEnsemble.variance();

//  Write mean to file
    if (fullConfig.has("mean output"))
      ensmean.write(eckit::LocalConfiguration(fullConfig, "mean output"));

    if (fullConfig.has("ensmean to structured grid")) {
      const eckit::LocalConfiguration latlonConf(fullConfig, "ensmean to structured grid");
      const StructuredGridWriter<MODEL> latlon(latlonConf, resol);
      latlon.interpolateAndWrite(ensmean);
    }
    Log::test() << "Mean: " << std::endl << ensmean << std::endl;

//  Write variance to file
    if (fullConfig.has("variance output"))
      sigb2.write(eckit::LocalConfiguration(fullConfig, "variance output"));

    if (fullConfig.has("ensvariance to structured grid")) {
      const eckit::LocalConfiguration latlonConf(fullConfig, "ensvariance to structured grid");
      const StructuredGridWriter<MODEL> latlon(latlonConf, resol);
      latlon.interpolateAndWrite(sigb2, ensmean);
    }
    Log::test() << "Variance: " << std::endl << sigb2 << std::endl;

//  Compute and write standard deviation to file, if it was requested
//  The std dev is computed via a generic atlas algorithm; computing it only when requested allows
//  this executable to work with model interfaces with non-conforming atlas interfaces.
    if (fullConfig.has("standard deviation output")
        || fullConfig.has("standard deviation to structured grid")) {
      const Increment_ sigb = stateEnsemble.stddev();

      if (fullConfig.has("standard deviation output")) {
        sigb.write(eckit::LocalConfiguration(fullConfig, "standard deviation output"));
      }
      if (fullConfig.has("standard deviation to structured grid")) {
        const eckit::LocalConfiguration latlonConf(fullConfig,
                                                   "standard deviation to structured grid");
        const StructuredGridWriter<MODEL> latlon(latlonConf, resol);
        latlon.interpolateAndWrite(sigb, ensmean);
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
