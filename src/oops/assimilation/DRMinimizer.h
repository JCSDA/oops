/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_ASSIMILATION_DRMINIMIZER_H_
#define OOPS_ASSIMILATION_DRMINIMIZER_H_

#include <memory>
#include <string>
#include <vector>


#include "eckit/config/Configuration.h"
#include "oops/assimilation/BMatrix.h"
#include "oops/assimilation/ControlIncrement.h"
#include "oops/assimilation/CostFunction.h"
#include "oops/assimilation/HtRinvHMatrix.h"
#include "oops/assimilation/Minimizer.h"

#include "oops/util/dot_product.h"
#include "oops/util/formats.h"
#include "oops/util/Logger.h"

namespace oops {

/// DR (Derber and Rosati) Minimizers
/*!
 * DRMinimizer is the base class for all minimizers that use \f$ B\f$ to
 * precondition the variational minimisation problem and use the auxiliary
 * variable \f$ \hat{x}=B^{-1}x\f$ and to update it in parallel to \f$ x\f$
 * based on Derber and Rosati, 1989, J. Phys. Oceanog. 1333-1347.
 * \f$ J_b\f$ is then computed as \f$ x^T\hat{x}\f$ eliminating the need for
 * \f$ B^{-1}\f$ or \f$ B^{1/2}\f$.
 */

// -----------------------------------------------------------------------------

template<typename MODEL> class DRMinimizer : public Minimizer<MODEL> {
  typedef BMatrix<MODEL>             Bmat_;
  typedef CostFunction<MODEL>        CostFct_;
  typedef ControlIncrement<MODEL>    CtrlInc_;
  typedef HtRinvHMatrix<MODEL>       HtRinvH_;
  typedef Minimizer<MODEL>           Minimizer_;

 public:
  explicit DRMinimizer(const CostFct_ & J): Minimizer_(J), J_(J), gradJb_(), costJ0Jb_(0) {}
  ~DRMinimizer() {}
  const std::string classname() const override = 0;

 private:
  CtrlInc_ * doMinimize(const eckit::Configuration &) override;
  virtual double solve(CtrlInc_ &, CtrlInc_ &, CtrlInc_ &,
                       const Bmat_ &, const HtRinvH_ &,
                       const double, const double,
                       const int, const double) = 0;

  const CostFct_ & J_;
  std::unique_ptr<CtrlInc_> gradJb_;
  std::vector<CtrlInc_> dxh_;
  double costJ0Jb_;
};

// =============================================================================

template<typename MODEL>
ControlIncrement<MODEL> * DRMinimizer<MODEL>::doMinimize(const eckit::Configuration & config) {
  int ninner = config.getInt("ninner");
  double gnreduc = config.getDouble("gradient_norm_reduction");

  bool runOnlineAdjTest = false;
  if (config.has("onlineDiagnostics")) {
    const eckit::LocalConfiguration onlineDiag(config, "onlineDiagnostics");
    runOnlineAdjTest = onlineDiag.getBool("onlineAdjTest");
  }

  if (gradJb_) {
    gradJb_.reset(new CtrlInc_(J_.jb().resolution(), *gradJb_));
  } else {
    gradJb_.reset(new CtrlInc_(J_.jb()));
  }

  Log::info() << std::endl;
  Log::info() << classname() << ": max iter = " << ninner
              <<  ", requested norm reduction = " << gnreduc << std::endl;

// Define the matrices
  const Bmat_    B(J_);
  const HtRinvH_ HtRinvH(J_, runOnlineAdjTest);

// Compute RHS (sum B^{-1} dx_{i}) + H^T R^{-1} d
// dx_i = x_i - x_{i-1}; dx_1 = x_1 - x_b
  CtrlInc_ rhs(J_.jb());
  J_.computeGradientFG(rhs);
  J_.jb().addGradientFG(rhs, *gradJb_);
  rhs *= -1.0;
  Log::info() << classname() << " rhs" << rhs << std::endl;

// Define minimisation starting point
  // dx
  CtrlInc_ * dx = new CtrlInc_(J_.jb());
  // dxh = B^{-1} dx
  CtrlInc_ dxh(J_.jb());

// Set J[0] = 0.5 (x_i - x_b)^T B^{-1} (x_i - x_b) + 0.5 d^T R^{-1} d
  const double costJ0Jb = costJ0Jb_;
  const double costJ0JoJc = J_.getCostJoJc();

// Solve the linear system
  double reduc = this->solve(*dx, dxh, rhs, B, HtRinvH, costJ0Jb, costJ0JoJc, ninner, gnreduc);

  Log::test() << classname() << ": reduction in residual norm = " << reduc << std::endl;
  Log::info() << classname() << " output increment:" << *dx << std::endl;

// Update gradient Jb
  *gradJb_ += dxh;
  dxh_.push_back(dxh);

// Update Jb component of J[0]: 0.5 (x_i - x_b)^T B^-1 (x_i - x_b)
  costJ0Jb_ += 0.5 * dot_product(*dx, dxh);
  for (unsigned int jouter = 1; jouter < dxh_.size(); ++jouter) {
    CtrlInc_ dxhtmp(dx->geometry(), dxh_[jouter-1]);
    costJ0Jb_ += dot_product(*dx, dxhtmp);
  }

  return dx;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_ASSIMILATION_DRMINIMIZER_H_
