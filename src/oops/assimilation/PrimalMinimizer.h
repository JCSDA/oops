/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_ASSIMILATION_PRIMALMINIMIZER_H_
#define OOPS_ASSIMILATION_PRIMALMINIMIZER_H_

#include <limits>
#include <memory>
#include <string>

#include "eckit/config/Configuration.h"
#include "oops/assimilation/ControlIncrement.h"
#include "oops/assimilation/CostFunction.h"
#include "oops/assimilation/Minimizer.h"
#include "oops/assimilation/PMatrix.h"
#include "oops/assimilation/RinvHMatrix.h"
#include "oops/util/FloatCompare.h"
#include "oops/util/Logger.h"

#include "oops/assimilation/DualVector.h"

namespace oops {

/// Primal Minimizer
/*!
 * PrimalMinimizer is the base class for all minimizers that minimize the
 * variational data assimilation cost function in primal (model) space.
 */

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS> class PrimalMinimizer : public Minimizer<MODEL, OBS> {
  typedef CostFunction<MODEL, OBS>        CostFct_;
  typedef ControlIncrement<MODEL, OBS>    CtrlInc_;
  typedef HessianMatrix<MODEL, OBS>       Hessian_;
  typedef Minimizer<MODEL, OBS>           Minimizer_;
  typedef PMatrix<MODEL, OBS>             Pmat_;
  typedef RinvHMatrix<MODEL, OBS>         RinvH_;

 public:
  explicit PrimalMinimizer(const CostFct_ & J): Minimizer_(J), J_(J) {}
  ~PrimalMinimizer() {}
  const std::string classname() const override = 0;

 private:
  CtrlInc_ * doMinimize(const eckit::Configuration &) override;
  virtual double solve(CtrlInc_ &, const CtrlInc_ &,
                       const Hessian_ &, const Pmat_ &,
                       const int, const double) = 0;

  const CostFct_ & J_;
};

// =============================================================================

template<typename MODEL, typename OBS>
ControlIncrement<MODEL, OBS> *
PrimalMinimizer<MODEL, OBS>::doMinimize(const eckit::Configuration & config) {
  int ninner = config.getInt("ninner");
  double gnreduc = config.getDouble("gradient norm reduction");

  bool runOnlineAdjTest = config.getBool("online diagnostics.online adj test", false);

  Log::info() << classname() << ": max iter = " << ninner
              << ", requested norm reduction = " << gnreduc << std::endl;

// Define the matrices
  Hessian_ hessian(J_, runOnlineAdjTest);
  Pmat_ P(J_);

// Define minimisation starting point
  CtrlInc_ * dx = new CtrlInc_(J_.jb());

// Compute RHS
  CtrlInc_ rhs(J_.jb());
  if (config.has("fsoi")) {
    const eckit::LocalConfiguration FcSensitivityConfig(config, "fsoi.input forecast sensitivity");
    rhs.read(FcSensitivityConfig);
    Log::info() << classname() << " rhs has forecast sensitivity" << std::endl;
  } else {
    J_.computeGradientFG(rhs);
    J_.jb().addGradientFG(rhs);
  }
  rhs *= -1.0;
  Log::info() << classname() << " rhs" << rhs << std::endl;

// Check for zero gradient (for example if no obs)
  const double gnorm = dot_product(rhs, rhs);
  const double epsilon = config.getDouble("epsilon", std::numeric_limits<double>::epsilon());
  Log::info() << "Initial RHS squared norm = " << gnorm << std::endl;
  if (gnorm < epsilon) {
    Log::info() << "RHS smaller than " << epsilon << ", returning." << std::endl;
    return dx;
  }

// Solve the linear system
  double reduc = this->solve(*dx, rhs, hessian, P, ninner, gnreduc);

  Log::test() << classname() << ": reduction in residual norm = " << reduc << std::endl;
  Log::info() << classname() << ": reduction in residual norm = " << reduc << std::endl;
  Log::info() << classname() << " output" << *dx << std::endl;

  if (config.has("fsoi")) {
    Log::info() << classname() << " Entering Observation Sensitivity Calculation" << std::endl;

    // Multiply result of solver by RinvH to get observation sensitivity (ys)
    DualVector<MODEL, OBS> ys;
    const RinvH_ RinvH(J_);
    RinvH.multiply(*dx, ys);

    // Write out observation sensitivity
    const std::string osensname = "ObsSensitivity";
    ys.saveDep(osensname);

    bool runFSOIincTest = config.getBool("fsoi.increment test", false);
    if (runFSOIincTest) {
      // Get departures
      DualVector<MODEL, OBS> dp;
      for (unsigned jj = 0; jj < J_.nterms(); ++jj) {
        std::unique_ptr<GeneralizedDepartures> ww(J_.jterm(jj).newGradientFG());
        dp.append(J_.jterm(jj).multiplyCovar(*ww));
      }

      // <K dp,  dx>, where dx = K dp
      double adj_tst_fwd = dot_product(rhs, rhs);
      // <dp, Kt dx>, where K = Hessian Ht Rinv; dp=departures
      double adj_tst_bwd = dot_product(ys, dp);

      Log::info() << "Online FSOI increment test: " << std::endl
                  << util::PrintAdjTest(adj_tst_fwd, adj_tst_bwd, "K") << std::endl;

      double fsoi_inctest_tolerance = config.getDouble("fsoi.increment test tolerance", 1.0e-5);
      bool passed = oops::is_close_absolute(adj_tst_fwd, adj_tst_bwd, fsoi_inctest_tolerance);
      if (passed) {
        Log::test() << "FSOI increment test within tolerance." << std::endl;
      } else {
        Log::test() << "FSOI increment test fails tolerance bound." << std::endl;
      }
    }

    // Make sure not to update state in FSOI mode
    dx->zero();
  }

  return dx;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_ASSIMILATION_PRIMALMINIMIZER_H_
