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
#include "oops/base/EnsemblesCollection.h"
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
    typedef EnsemblesCollection<MODEL>               EnsemblesCollection_;
    typedef Geometry<MODEL>                          Geometry_;
    typedef Increment<MODEL>                         Increment_;
    typedef State<MODEL>                             State_;

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
      State_ xb(resol, vars, bkgConfig);

      // Compute rescaled and transformed ensemble perturbations
      //        ens_k = K^-1 dx_k / (N-1)^0.5
      const eckit::LocalConfiguration ensConfig(fullConfig, "Ensemble");
      EnsemblePtr_ ens_k(new Ensemble_(xb.validTime(), ensConfig));
      ens_k->linearize(xb, xb, resol);
      EnsemblesCollection_::getInstance().put(xb.validTime(), ens_k);

      // Get ensemble size
      unsigned nm = ens_k->size();

      // Compute ensemble standard deviation
      Increment_ km1dx((*ens_k)[0]);
      km1dx.zero();
      Increment_ sigb2(km1dx);
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
