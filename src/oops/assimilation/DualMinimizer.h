/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_ASSIMILATION_DUALMINIMIZER_H_
#define OOPS_ASSIMILATION_DUALMINIMIZER_H_

#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "eckit/config/Configuration.h"
#include "oops/assimilation/BMatrix.h"
#include "oops/assimilation/ControlIncrement.h"
#include "oops/assimilation/CostFunction.h"
#include "oops/assimilation/DualVector.h"
#include "oops/assimilation/HBHtMatrix.h"
#include "oops/assimilation/HMatrix.h"
#include "oops/assimilation/HtMatrix.h"
#include "oops/assimilation/Minimizer.h"
#include "oops/assimilation/RinvMatrix.h"
#include "oops/util/dot_product.h"
#include "oops/util/Logger.h"

namespace oops {

/// Dual Minimizer
/*!
 * Base class for all dual (observation) space minimizers.
 */

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS> class DualMinimizer : public Minimizer<MODEL, OBS> {
  typedef ControlIncrement<MODEL, OBS>    CtrlInc_;
  typedef CostFunction<MODEL, OBS>        CostFct_;
  typedef BMatrix<MODEL, OBS>             Bmat_;
  typedef DualVector<MODEL, OBS>          Dual_;
  typedef HBHtMatrix<MODEL, OBS>          HBHt_;
  typedef Minimizer<MODEL, OBS>           Minimizer_;
  typedef RinvMatrix<MODEL, OBS>          Rinv_;

 public:
  explicit DualMinimizer(const CostFct_ & J): Minimizer_(J), J_(J), gradJb_(), costJ0Jb_(0) {}
  ~DualMinimizer() {}
  const std::string classname() const override = 0;

 private:
  CtrlInc_ * doMinimize(const eckit::Configuration &) override;
  virtual double solve(Dual_ &, double &, Dual_ &, const HBHt_ &, const Rinv_ &,
                       const double, const double,
                       const int &, const double &, Dual_ &, const double &) = 0;

  const CostFct_ & J_;
  std::unique_ptr<CtrlInc_> gradJb_;
  std::vector<CtrlInc_> dxh_;
  double costJ0Jb_;
};

// =============================================================================

template<typename MODEL, typename OBS>
ControlIncrement<MODEL, OBS> *
DualMinimizer<MODEL, OBS>::doMinimize(const eckit::Configuration & config) {
  int ninner = config.getInt("ninner");
  double gnreduc = config.getDouble("gradient norm reduction");

  bool runOnlineAdjTest = config.getBool("online diagnostics.online adj test", false);

  if (gradJb_) {
    gradJb_.reset(new CtrlInc_(J_.jb().resolution(), *gradJb_));
  } else {
    gradJb_.reset(new CtrlInc_(J_.jb()));
  }

  Log::info() << std::endl;
  Log::info() << classname() << ": max iter = " << ninner
              << ", requested norm reduction = " << gnreduc << std::endl;

// Define the matrices
  const Bmat_ B(J_);
  const HBHt_ HBHt(J_, runOnlineAdjTest);
  const Rinv_ Rinv(J_);
  const HMatrix<MODEL, OBS> H(J_);
  const HtMatrix<MODEL, OBS> Ht(J_);

  CtrlInc_ * dx = new CtrlInc_(J_.jb());

// Define minimisation starting point in dual space
  Dual_ vv;
  for (unsigned jj = 0; jj < J_.nterms(); ++jj) {
    vv.append(J_.jterm(jj).newDualVector());
  }
  double vvp = 0.0;

// Get R^{-1} d
  Dual_ rr;
  for (unsigned jj = 0; jj < J_.nterms(); ++jj) {
    rr.append(J_.jterm(jj).newGradientFG());
  }
  rr *= -1.0;

// Check for zero gradient (for example if no obs)
  const double rnorm = dot_product(rr, rr);
  const double epsilon = config.getDouble("epsilon", std::numeric_limits<double>::epsilon());
  Log::info() << "Initial RHS squared norm = " << rnorm << std::endl;
  if (rnorm < epsilon) {
    Log::info() << "RHS smaller than " << epsilon << ", returning." << std::endl;
    return dx;
  }

// Update rr if initial dx in model space is not a zero vector
// rr = rr - Rinv H dx0

// Compute a = H (xb - xk)
  Dual_ dy;
  CtrlInc_ hdxfg(J_.jb().getFirstGuess());
  H.multiply(hdxfg, dy);
  dy *= -1.0;

// Compute sigma = (xb - xk)^T Binv (xb - xk)
  CtrlInc_ g0(J_.jb());
  J_.jb().addGradientFG(g0, *gradJb_);

  double sigma = dot_product(J_.jb().getFirstGuess(), g0);

// Set J[0] = 0.5 (x_i - x_b)^T B^{-1} (x_i - x_b) + 0.5 d^T R^{-1} d
  const double costJ0Jb = costJ0Jb_;
  const double costJ0JoJc = J_.getCostJoJc();

// Solve the linear system
  double reduc = this->solve(vv, vvp, rr, HBHt, Rinv, costJ0Jb, costJ0JoJc,
                             ninner, gnreduc, dy, sigma);

  Log::info() << classname() << ": reduction in residual norm = " << reduc << std::endl;
  Log::test() << classname() << ": reduction in residual norm = " << reduc << std::endl;

// Recover solution in primal space
  CtrlInc_ dh(J_.jb());
  J_.zeroAD(dh);
  Ht.multiply(vv, dh);
  dh.axpy(-vvp, g0);

  B.multiply(dh, *dx);    // BHtaug vvaug

  Log::info() << classname() << ": Estimated Final Jb = "
              << 0.5 * dot_product(*dx, dh) << std::endl;
  Log::info() << classname() << " output" << *dx << std::endl;

// Update gradient Jb
  *gradJb_ += dh;
  dxh_.push_back(dh);

// Update Jb component of J[0]: 0.5 (x_i - x_b)^T B^-1 (x_i - x_b)
  costJ0Jb_ += 0.5 * dot_product(*dx, dh);
  for (unsigned int jouter = 1; jouter < dxh_.size(); ++jouter) {
    CtrlInc_ dxhtmp(dx->geometry(), dxh_[jouter-1]);
    costJ0Jb_ += dot_product(*dx, dxhtmp);
  }

  return dx;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_ASSIMILATION_DUALMINIMIZER_H_
