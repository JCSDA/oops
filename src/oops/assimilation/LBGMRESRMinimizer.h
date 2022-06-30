/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_ASSIMILATION_LBGMRESRMINIMIZER_H_
#define OOPS_ASSIMILATION_LBGMRESRMINIMIZER_H_

#include <cmath>
#include <string>
#include <vector>

#include "oops/assimilation/ControlIncrement.h"
#include "oops/assimilation/CostFunction.h"
#include "oops/assimilation/GMRESR.h"
#include "oops/assimilation/LBHessianMatrix.h"
#include "oops/assimilation/LBMinimizer.h"
#include "oops/base/IdentityMatrix.h"

namespace eckit {
  class Configuration;
}

namespace oops {

/// LBGMRESR Minimizer
/*!
 * \brief Left B Preconditioned GMRESR solver.
 *
 * On entry:
 * -    dx      =  starting point, \f$ dx_{0} \f$.
 * -    rr      = \f$ (sum dx^{b}_{i} + ) B H^T R^{-1} d \f$
 *
 *  Iteration will stop if the maximum iteration limit "maxiter" is reached
 *  or if the residual norm reduces by a factor of "tolerance".
 *
 *  Each matrix must implement a method:
 *  - void multiply(const VECTOR&, VECTOR&) const
 *
 *  which applies the matrix to the first argument, and returns the
 *  matrix-vector product in the second. (Note: the const is optional, but
 *  recommended.)
 */

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS> class LBGMRESRMinimizer : public LBMinimizer<MODEL, OBS> {
  typedef LBHessianMatrix<MODEL, OBS>  LBHessianMatrix_;
  typedef CostFunction<MODEL, OBS>     CostFct_;
  typedef ControlIncrement<MODEL, OBS> CtrlInc_;

 public:
  const std::string classname() const override {return "LBGMRESRMinimizer";}
  LBGMRESRMinimizer(const eckit::Configuration &, const CostFct_ &);
  ~LBGMRESRMinimizer() {}

 private:
  void solve(CtrlInc_ &, CtrlInc_ &, const LBHessianMatrix_ &,
             const int, const double) override;
};

// =============================================================================

template<typename MODEL, typename OBS>
LBGMRESRMinimizer<MODEL, OBS>::LBGMRESRMinimizer(const eckit::Configuration & conf,
                                            const CostFct_ & J)
  : LBMinimizer<MODEL, OBS>(J)
{}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void LBGMRESRMinimizer<MODEL, OBS>::solve(CtrlInc_ & dx, CtrlInc_ & rr,
                                     const LBHessianMatrix_ & LBHessian,
                                     const int maxiter, const double tolerance) {
  IdentityMatrix<CtrlInc_> Id;
  dx.zero();
  GMRESR(dx, rr, LBHessian, Id, maxiter, tolerance);
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_ASSIMILATION_LBGMRESRMINIMIZER_H_
