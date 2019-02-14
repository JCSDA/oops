/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_RUNS_MOMENTS_H_
#define OOPS_RUNS_MOMENTS_H_

#include <string>
#include <vector>

#include <boost/scoped_ptr.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "oops/base/Accumulator.h"
#include "oops/base/LinearVariableChangeBase.h"
#include "oops/base/Variables.h"
#include "oops/generic/instantiateLinearVariableChangeFactory.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/State.h"
#include "oops/runs/Application.h"
#include "oops/util/DateTime.h"

namespace oops {

  template <typename MODEL> class Moments : public Application {
    typedef Geometry<MODEL>                    Geometry_;
    typedef State<MODEL>                       State_;
    typedef Increment<MODEL>                   Increment_;
    typedef LinearVariableChangeBase<MODEL>    Balance_;
    typedef LinearVariableChangeFactory<MODEL> BalanceFactory_;

   public:
    // -----------------------------------------------------------------------------
    Moments() {}
    // -----------------------------------------------------------------------------
    virtual ~Moments() {}
    // -----------------------------------------------------------------------------
    int execute(const eckit::Configuration & fullConfig) const {
      //  Setup resolution
      const eckit::LocalConfiguration resolConfig(fullConfig, "Geometry");
      const Geometry_ resol(resolConfig);

      //  Setup variables
      const eckit::LocalConfiguration varConfig(fullConfig, "Variables");
      const Variables vars(varConfig);

      //  Setup background states configuration
      const eckit::LocalConfiguration ensConfig(fullConfig, "States");
      std::vector<eckit::LocalConfiguration> members;

      //  Get valid time from config
      const eckit::LocalConfiguration dateConfig(fullConfig, "date");
      const util::DateTime validtime(dateConfig.getString("validtime"));

      //  Setup output configuration
      const eckit::LocalConfiguration meanout(fullConfig, "MeanOut");
      const eckit::LocalConfiguration varianceout(fullConfig, "VarianceOut");

      //  Compute ensemble mean
      ensConfig.get("state", members);
      Log::debug() << "Moments: using " << members.size()
                   << " members." << std::endl;
      Accumulator<MODEL, State_, State_> x_mean(resol, vars, validtime);
      unsigned nm = members.size();
      double zz = 1.0/nm;
      for (unsigned jj = 0; jj < nm; ++jj) {
        State_ xx(resol, vars, members[jj]);
        x_mean.accumul(zz, xx);
        Log::trace() << "Reading ensemble member: "<< jj << std::endl;
        Log::trace() << xx << jj << std::endl;
      }
      x_mean.write(meanout);

      //  Setup balance operator
      boost::scoped_ptr<Balance_> balance(BalanceFactory_::create(
               x_mean, x_mean, resol, fullConfig));

      //  Compute variance of nm (K^-1 dx) perturbations
      zz = 1.0/(nm-1.0);
      Increment_ pert(resol, vars, validtime);
      Increment_ dx(pert);
      dx.zero();
      for (unsigned jj = 0; jj < nm; ++jj) {
        // Read ensemble member jj
        State_ xx(resol, vars, members[jj]);

        // Compute departure to mean
        pert.diff(x_mean, xx);

        // Left multiply by K^-1
        pert = balance->multiplyInverse(pert);

        // Accumulate dx^2
        pert.schur_product_with(pert);
        dx += pert;

        Log::trace() << "K^-1 dx for ensemble member perturbation: "<< jj << std::endl;
        Log::trace() << pert << jj << std::endl;
      }
      dx *= zz;
      dx.write(varianceout);
      return 0;
    }
    // -----------------------------------------------------------------------------
   private:
    std::string appname() const {
      return "oops::Moments<" + MODEL::name() + ">";
    }
    // -----------------------------------------------------------------------------
  };

}  // namespace oops

#endif  // OOPS_RUNS_MOMENTS_H_
