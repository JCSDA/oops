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

#include <sstream>
#include <string>
#include <vector>

#include <boost/scoped_ptr.hpp>

#include "eckit/config/Configuration.h"
#include "oops/base/instantiateCovarFactory.h"
#include "oops/base/ModelSpaceCovariance4DBase.h"
#include "oops/base/ModelSpaceCovarianceBase.h"
#include "oops/base/PostProcessor.h"
#include "oops/base/StateWriter.h"
#include "oops/base/Variables.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
#include "oops/interface/Model.h"
#include "oops/interface/ModelAuxControl.h"
#include "oops/interface/State.h"
#include "oops/runs/Application.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

namespace oops {

template <typename MODEL> class Dirac : public Application {
  typedef Geometry<MODEL>                         Geometry_;
  typedef Model<MODEL>                            Model_;
  typedef ModelAuxControl<MODEL>                  ModelAux_;
  typedef Increment<MODEL>                        Increment_;
  typedef Increment4D<MODEL>                      Increment4D_;
  typedef State<MODEL>                            State_;
  typedef State4D<MODEL>                          State4D_;
  typedef Localization<MODEL>                     Localization_;
  typedef StateEnsemble<MODEL>                    Ensemble_;
  typedef boost::shared_ptr<StateEnsemble<MODEL>> EnsemblePtr_;

 public:
// -----------------------------------------------------------------------------
  Dirac() {
    instantiateCovarFactory<MODEL>();
  }
// -----------------------------------------------------------------------------
  virtual ~Dirac() {}
// -----------------------------------------------------------------------------
  int execute(const eckit::Configuration & fullConfig) const {
    //  Setup resolution
    const eckit::LocalConfiguration resolConfig(fullConfig, "resolution");
    const Geometry_ resol(resolConfig);

    //  Setup variables
    const Variables vars(fullConfig);

    // Setup background state
    const eckit::LocalConfiguration backgroundConfig(fullConfig, "initial");
    boost::scoped_ptr<State4D_> xx;
    bool l3d;
    if (backgroundConfig.has("state")) {
      xx.reset(new State4D_(backgroundConfig, vars, resol));
      l3d = false;
    } else {
      State_ xx3D(resol, vars, backgroundConfig);
      xx.reset(new State4D_(xx3D));
      l3d = true;
    }

    //  Setup timeslots
    std::vector<util::DateTime> timeslots;
    for (unsigned jsub = 0; jsub < (*xx).size(); ++jsub) {
      timeslots.push_back((*xx)[jsub].validTime());
    }
    Log::info() << "Number of ensemble time-slots:" << timeslots.size() << std::endl;

    // Apply B to Dirac
    const eckit::LocalConfiguration covarConfig(fullConfig, "Covariance");
    if (l3d) {
      //  3D covariance
      boost::scoped_ptr<ModelSpaceCovarianceBase<MODEL>> B(CovarianceFactory<MODEL>::create(
        covarConfig, resol, vars, (*xx)[0], (*xx)[0]));

      //  Setup Dirac
      Increment_ dxdirin(resol, vars, timeslots[0]);
      Increment_ dxdirout(resol, vars, timeslots[0]);
      const eckit::LocalConfiguration diracConfig(fullConfig, "dirac");
      dxdirin.dirac(diracConfig);

      //  Apply 3D B matrix to Dirac increment
      B->multiply(dxdirin, dxdirout);

      //  Write increment
      const eckit::LocalConfiguration output_B(fullConfig, "output_B");
      dxdirout.write(output_B);
      Log::test() << "Increment: " << dxdirout << std::endl;
    } else {
      //  4D covariance
      boost::scoped_ptr<ModelSpaceCovariance4DBase<MODEL>> B(Covariance4DFactory<MODEL>::create(
        covarConfig, resol, vars, (*xx), (*xx)));

      //  Setup Dirac
      Increment4D_ dxdirin(resol, vars, timeslots);
      Increment4D_ dxdirout(resol, vars, timeslots);
      std::vector<eckit::LocalConfiguration> diracConfigs;
      fullConfig.get("dirac", diracConfigs);
      dxdirin.dirac(diracConfigs);

      //  Apply 4D B matrix to Dirac increment
      B->multiply(dxdirin, dxdirout);

      //  Write increment
      const eckit::LocalConfiguration output_B(fullConfig, "output_B");
      dxdirout.write(output_B);
      Log::test() << "Increment4D: " << dxdirout << std::endl;
    }

    //  Setup localization and ensemble configurations
    eckit::LocalConfiguration locConfig;
    eckit::LocalConfiguration ensConfig;
    bool hasLoc(false);
    if (covarConfig.has("localization")) {
      locConfig = eckit::LocalConfiguration(covarConfig, "localization");
      ensConfig = covarConfig;
      hasLoc = true;
    } else {
      if (covarConfig.has("ensemble")) {
        ensConfig = eckit::LocalConfiguration(covarConfig, "ensemble");
        if (ensConfig.has("localization")) {
          locConfig = eckit::LocalConfiguration(ensConfig, "localization");
          hasLoc = true;
        }
      }
    }

    if (hasLoc) {
      // Setup ensemble
      EnsemblePtr_ ens(new Ensemble_(timeslots, ensConfig));
      ens->linearize((*xx), (*xx), resol);

      // Apply localization to Dirac
      if (l3d) {
        //  Setup Dirac
        Increment_ dxdir(resol, vars, timeslots[0]);
        const eckit::LocalConfiguration diracConfig(fullConfig, "dirac");
        dxdir.dirac(diracConfig);

        //  Setup localization
        boost::scoped_ptr<Localization_> loc_;
        loc_.reset(new Localization_(resol, ens, locConfig));

        //  Apply localization
        loc_->multiply(dxdir);

        //  Write increment
        const eckit::LocalConfiguration output_localization(fullConfig, "output_localization");
        dxdir.write(output_localization);
        Log::test() << "Increment: " << dxdir << std::endl;
      } else {
        //  Setup Dirac
        Increment4D_ dxdir(resol, vars, timeslots);
        std::vector<eckit::LocalConfiguration> diracConfigs;
        fullConfig.get("dirac", diracConfigs);
        dxdir.dirac(diracConfigs);

        //  Setup localization
        boost::scoped_ptr<Localization_> loc_;
        loc_.reset(new Localization_(resol, ens, locConfig));

        //  Apply localization
        loc_->multiply(dxdir);

        //  Write increment
        const eckit::LocalConfiguration output_localization(fullConfig, "output_localization");
        dxdir.write(output_localization);
        Log::test() << "Increment4D: " << dxdir << std::endl;
      }
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
