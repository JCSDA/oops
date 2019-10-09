/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_RUNS_ESTIMATEPARAMS_H_
#define OOPS_RUNS_ESTIMATEPARAMS_H_

#include <memory>
#include <sstream>
#include <string>
#include <vector>


#include "eckit/config/Configuration.h"
#include "oops/assimilation/Increment4D.h"
#include "oops/assimilation/State4D.h"
#include "oops/base/IncrementEnsemble.h"
#include "oops/base/instantiateCovarFactory.h"
#include "oops/generic/ParametersBUMP.h"
#include "oops/parallel/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/Logger.h"

namespace oops {

template <typename MODEL> class EstimateParams : public Application {
  typedef Geometry<MODEL>                         Geometry_;
  typedef Increment<MODEL>                        Increment_;
  typedef Increment4D<MODEL>                      Increment4D_;
  typedef State<MODEL>                            State_;
  typedef State4D<MODEL>                          State4D_;
  typedef ParametersBUMP<MODEL>                   Parameters_;
  typedef IncrementEnsemble<MODEL>                Ensemble_;
  typedef boost::shared_ptr<IncrementEnsemble<MODEL>> EnsemblePtr_;

 public:
// -----------------------------------------------------------------------------
  static const std::string classname() {return "oops::EstimateParams";}
  explicit EstimateParams(const eckit::mpi::Comm & comm = oops::mpi::comm()) : Application(comm) {
    instantiateCovarFactory<MODEL>();
  }
// -----------------------------------------------------------------------------
  virtual ~EstimateParams() {}
// -----------------------------------------------------------------------------
  int execute(const eckit::Configuration & fullConfig) const {
    util::Timer timer(classname(), "write");

    //  Setup resolution
    const eckit::LocalConfiguration resolConfig(fullConfig, "resolution");
    const Geometry_ resol(resolConfig, this->getComm());

    // Setup variables
    const Variables vars(fullConfig);

    // Setup background state
    const eckit::LocalConfiguration backgroundConfig(fullConfig, "background");
    std::unique_ptr<State4D_> xx;
    if (backgroundConfig.has("state")) {
      xx.reset(new State4D_(backgroundConfig, vars, resol));
    } else {
      State_ xx3D(resol, vars, backgroundConfig);
      xx.reset(new State4D_(xx3D));
    }

    //  Setup timeslots
    std::vector<util::DateTime> timeslots;
    for (unsigned jsub = 0; jsub < (*xx).size(); ++jsub) {
      timeslots.push_back((*xx)[jsub].validTime());
    }
    Log::info() << "Number of ensemble time-slots:" << timeslots.size() << std::endl;

    // Setup ensemble
    EnsemblePtr_ ens = NULL;
    if (fullConfig.has("ensemble")) {
      const eckit::LocalConfiguration ensembleConfig(fullConfig, "ensemble");
      ens.reset(new Ensemble_(ensembleConfig, (*xx), (*xx), resol));
    }

    // Setup pseudo ensemble
    EnsemblePtr_ pseudo_ens = NULL;
    if (fullConfig.has("covariance")) {
      const eckit::LocalConfiguration covarConfig(fullConfig, "covariance");
      int ens2_ne = covarConfig.getInt("pseudoens_size");
      pseudo_ens.reset(new Ensemble_(resol, vars, timeslots, ens2_ne));
      if (timeslots.size() == 1) {
      // One time-slot only
        std::unique_ptr<ModelSpaceCovarianceBase<MODEL>>
          cov(CovarianceFactory<MODEL>::create(covarConfig, resol, vars, (*xx)[0], (*xx)[0]));
        for (int ie = 0; ie < ens2_ne; ++ie) {
          Log::info() << "Generate pseudo ensemble member " << ie+1 << " / "
                      << ens2_ne << std::endl;

          // Compute a pseudo ensemble using randomization
          Increment4D_ incr(resol, vars, timeslots);
          cov->randomize(incr[incr.first()]);
          (*pseudo_ens)[ie] = incr;
        }
      } else {
        // Multiple time-slots
        std::unique_ptr<ModelSpaceCovariance4DBase<MODEL>>
          cov(Covariance4DFactory<MODEL>::create(covarConfig, resol, vars, (*xx), (*xx)));
        for (int ie = 0; ie < ens2_ne; ++ie) {
          Log::info() << "Generate pseudo ensemble member " << ie+1 << " / "
                      << ens2_ne << std::endl;

          // Compute a pseudo ensemble using randomization
          Increment4D_ incr(resol, vars, timeslots);
          cov->randomize(incr);
          (*pseudo_ens)[ie] = incr;
        }
      }
    }

    // Setup parameters
    Parameters_ param(resol, vars, timeslots, fullConfig, ens, pseudo_ens);

    // Write parameters
    param.write();

    // Clean BUMP
    param.clean_bump();

    return 0;
  }
// -----------------------------------------------------------------------------
 private:
  std::string appname() const {
    return "oops::EstimateParams<" + MODEL::name() + ">";
  }
// -----------------------------------------------------------------------------
};

}  // namespace oops
#endif  // OOPS_RUNS_ESTIMATEPARAMS_H_
