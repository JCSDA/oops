/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_RUNS_DIRAC_H_
#define OOPS_RUNS_DIRAC_H_

#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "eckit/config/Configuration.h"
#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/IncrementEnsemble.h"
#include "oops/base/instantiateCovarFactory.h"
#include "oops/base/ModelSpaceCovarianceBase.h"
#include "oops/base/PostProcessor.h"
#include "oops/base/State.h"
#include "oops/base/Variables.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

namespace oops {

template <typename MODEL> class Dirac : public Application {
  typedef Geometry<MODEL>                           Geometry_;
  typedef Increment<MODEL>                          Increment_;
  typedef State<MODEL>                              State_;
  typedef LocalizationBase<MODEL>                   Localization_;
  typedef IncrementEnsemble<MODEL>                  Ensemble_;
  typedef std::shared_ptr<IncrementEnsemble<MODEL>> EnsemblePtr_;

 public:
// -----------------------------------------------------------------------------
  explicit Dirac(const eckit::mpi::Comm & comm = oops::mpi::world()) : Application(comm) {
    instantiateCovarFactory<MODEL>();
  }
// -----------------------------------------------------------------------------
  virtual ~Dirac() {}
// -----------------------------------------------------------------------------
  int execute(const eckit::Configuration & fullConfig) const {
    const eckit::LocalConfiguration backgroundConfig(fullConfig, "initial condition");
    std::vector<eckit::LocalConfiguration> confs;
    backgroundConfig.get("states", confs);
    size_t nslots = confs.size();

    const eckit::mpi::Comm * commSpace = &this->getComm();
    const eckit::mpi::Comm * commTime = &oops::mpi::myself();
    if (nslots > 1) {
      size_t ntasks = this->getComm().size();
      ASSERT(ntasks % nslots == 0);
      size_t myrank = this->getComm().rank();
      size_t ntaskpslot = ntasks / nslots;
      size_t myslot = myrank / ntaskpslot;

      // Create a communicator for same sub-window, to be used for communications in space
      std::string sgeom = "comm_geom_" + std::to_string(myslot);
      char const *geomName = sgeom.c_str();
      commSpace = &this->getComm().split(myslot, geomName);
      ASSERT(commSpace->size() == ntaskpslot);

      // Create a communicator for same local area, to be used for communications in time
      size_t myarea = commSpace->rank();
      std::string stime = "comm_time_" + std::to_string(myarea);
      char const *timeName = stime.c_str();
      commTime = &this->getComm().split(myarea, timeName);
      ASSERT(commTime->size() == nslots);
    }

    //  Setup resolution
    const eckit::LocalConfiguration resolConfig(fullConfig, "geometry");
    const Geometry_ resol(resolConfig, *commSpace, *commTime);

    // Setup background state
    State_ xx(resol, backgroundConfig);

    //  Setup variables
    const Variables vars = xx.variables();

    //  Setup time
    util::DateTime time = xx.validTime();

    // Apply B to Dirac
    const eckit::LocalConfiguration covarConfig(fullConfig, "background error");

    //  Covariance
    std::unique_ptr<ModelSpaceCovarianceBase<MODEL>> B(CovarianceFactory<MODEL>::create(
      covarConfig, resol, vars, xx, xx));

    //  Setup Dirac
    Increment_ dxdirin(resol, vars, time);
    Increment_ dxdirout(resol, vars, time);
    const eckit::LocalConfiguration diracConfig(fullConfig, "dirac");
    dxdirin.dirac(diracConfig);
    Log::test() << "Input Dirac increment: " << dxdirin << std::endl;

    //  Apply 3D B matrix to Dirac increment
    B->multiply(dxdirin, dxdirout);

    //  Write increment
    const eckit::LocalConfiguration output_B(fullConfig, "output B");
    dxdirout.write(output_B);
    Log::test() << "B * Increment: " << dxdirout << std::endl;

    //  Setup localization and ensemble configurations
    std::vector<eckit::LocalConfiguration> locConfigs;
    if (covarConfig.has("localization")) {
      locConfigs.push_back(eckit::LocalConfiguration(covarConfig, "localization"));
    } else {
      if (covarConfig.has("components")) {
        std::vector<eckit::LocalConfiguration> confs;
        covarConfig.get("components", confs);
        for (const auto & conf : confs) {
          const eckit::LocalConfiguration componentConf(conf, "covariance");
          if (componentConf.has("localization")) {
            locConfigs.push_back(eckit::LocalConfiguration(componentConf, "localization"));
          }
        }
      }
    }

    for (size_t jcomp = 0; jcomp < locConfigs.size(); ++jcomp) {
      // Apply localization to Dirac

      //  Setup Dirac
      Increment_ dxdir(resol, vars, time);
      const eckit::LocalConfiguration diracConfig(fullConfig, "dirac");
      dxdir.dirac(diracConfig);

      //  Setup localization
      std::unique_ptr<Localization_> loc_ =
              LocalizationFactory<MODEL>::create(resol, time, locConfigs[jcomp]);

      //  Apply localization
      loc_->multiply(dxdir);

      //  Write increment
      const eckit::LocalConfiguration output_localization(fullConfig, "output localization");
      dxdir.write(output_localization);
      Log::test() << "Localized Increment: " << dxdir << std::endl;
    }

    if (fullConfig.has("output variance")) {
      // Variance configuration
      const eckit::LocalConfiguration output_variance(fullConfig, "output variance");

      //  Setup variance
      Increment_ variance(resol, vars, time);

      // Get variance
      B->getVariance(variance);

      //  Write increment
      variance.write(output_variance);
      Log::test() << "Randomized variance: " << variance << std::endl;
    }

    return 0;
  }
// -----------------------------------------------------------------------------
 private:
  std::string appname() const {
    return "oops::Dirac<" + MODEL::name() + ">";
  }
// -----------------------------------------------------------------------------
};

}  // namespace oops
#endif  // OOPS_RUNS_DIRAC_H_
