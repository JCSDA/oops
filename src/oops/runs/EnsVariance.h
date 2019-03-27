/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_RUNS_ENSVARIANCE_H_
#define OOPS_RUNS_ENSVARIANCE_H_

#include <string>
#include <vector>

#include <boost/scoped_ptr.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "oops/assimilation/Increment4D.h"
#include "oops/assimilation/State4D.h"
#include "oops/base/StateEnsemble.h"
#include "oops/base/Variables.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/State.h"
#include "oops/runs/Application.h"
#include "oops/util/DateTime.h"

namespace oops {

  template <typename MODEL> class EnsVariance : public Application {
    typedef StateEnsemble<MODEL>                     Ensemble_;
    typedef boost::shared_ptr<StateEnsemble<MODEL> > EnsemblePtr_;
    typedef Geometry<MODEL>                          Geometry_;
    typedef Increment4D<MODEL>                       Increment4D_;
    typedef State<MODEL>                             State_;
    typedef State4D<MODEL>                           State4D_;

   public:
    // -----------------------------------------------------------------------------
    EnsVariance() {}
    // -----------------------------------------------------------------------------
    virtual ~EnsVariance() {}
    // -----------------------------------------------------------------------------
    int execute(const eckit::Configuration & fullConfig) const {
      // Setup Geometry
      const eckit::LocalConfiguration resolConfig(fullConfig, "Geometry");
      const Geometry_ resol(resolConfig);

      // Setup variables
      const eckit::LocalConfiguration varConfig(fullConfig, "Variables");
      const Variables vars(varConfig);

      // Setup background
      const eckit::LocalConfiguration bkgConfig(fullConfig, "Background");
      boost::scoped_ptr<State4D_> xx;
      if (bkgConfig.has("state")) {
        xx.reset(new State4D_(bkgConfig, vars, resol));
      } else {
        State_ xx3D(resol, vars, bkgConfig);
        xx.reset(new State4D_(xx3D));
      }

      //  Setup timeslots
      std::vector<util::DateTime> timeslots;
      for (unsigned jsub = 0; jsub < (*xx).size(); ++jsub) {
        timeslots.push_back((*xx)[jsub].validTime());
      }

      // Compute rescaled and transformed ensemble perturbations
      //        ens_k = K^-1 dx_k / (N-1)^0.5
      const eckit::LocalConfiguration ensConfig(fullConfig, "Ensemble");
      EnsemblePtr_ ens_k(new Ensemble_(timeslots, ensConfig));
      ens_k->linearize((*xx), (*xx), resol);

      // Get ensemble size
      unsigned nm = ens_k->size();

      // Compute ensemble standard deviation
      Increment4D_ km1dx((*ens_k)[0]);
      km1dx.zero();
      Increment4D_ sigb2(km1dx);
      sigb2.zero();

      for (unsigned jj = 0; jj < nm; ++jj) {
        km1dx = (*ens_k)[jj];

        // Accumulate km1dx^2
        km1dx.schur_product_with(km1dx);
        sigb2 += km1dx;
      }
      // Write variance to file
      const eckit::LocalConfiguration varianceout(fullConfig, "VarianceOut");
      sigb2.write(varianceout);

      return 0;
    }
    // -----------------------------------------------------------------------------
   private:
    std::string appname() const {
      return "oops::EnsVariance<" + MODEL::name() + ">";
    }
    // -----------------------------------------------------------------------------
  };

}  // namespace oops

#endif  // OOPS_RUNS_ENSVARIANCE_H_
