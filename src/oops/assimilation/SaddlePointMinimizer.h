/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_ASSIMILATION_SADDLEPOINTMINIMIZER_H_
#define OOPS_ASSIMILATION_SADDLEPOINTMINIMIZER_H_

#include <iomanip>
#include <memory>
#include <string>
#include <vector>

#include "eckit/config/Configuration.h"

#include "oops/assimilation/ControlIncrement.h"
#include "oops/assimilation/CostFunction.h"
#include "oops/assimilation/DualVector.h"
#include "oops/assimilation/GMRESR.h"
#include "oops/assimilation/Minimizer.h"
// Eigen #include "oops/assimilation/SaddlePointLMPMatrix.h"
#include "oops/assimilation/SaddlePointMatrix.h"
#include "oops/assimilation/SaddlePointPrecondMatrix.h"
#include "oops/assimilation/SaddlePointVector.h"

#include "oops/util/Logger.h"

namespace oops {

/// SaddlePoint Minimizer
/*!
 * Implements the SaddlePoint algorithm.
 */

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS> class SaddlePointMinimizer : public Minimizer<MODEL, OBS> {
  typedef ControlIncrement<MODEL, OBS>    CtrlInc_;
  typedef CostFunction<MODEL, OBS>        CostFct_;
  typedef Minimizer<MODEL, OBS>           Minimizer_;
  typedef DualVector<MODEL, OBS>          Multipliers_;
// Eigen   typedef SaddlePointLMPMatrix<MODEL, OBS>  LMP_;

 public:
  const std::string classname() const override {return "SaddlePointMinimizer";}
  SaddlePointMinimizer(const eckit::Configuration &, const CostFct_ & J)
    : Minimizer_(J), J_(J), gradJb_() {}
  ~SaddlePointMinimizer() {}

 private:
  CtrlInc_ * doMinimize(const eckit::Configuration &) override;

  const CostFct_ & J_;
  std::unique_ptr<CtrlInc_> gradJb_;
// Eigen  std::unique_ptr<LMP_> Pinv_;
  std::vector< SaddlePointVector<MODEL, OBS> > xyVEC_;
  std::vector< SaddlePointVector<MODEL, OBS> > pqVEC_;
};

// =============================================================================

template<typename MODEL, typename OBS>
ControlIncrement<MODEL, OBS> *
SaddlePointMinimizer<MODEL, OBS>::doMinimize(const eckit::Configuration & config) {
  int ninner = config.getInt("ninner");
  double gnreduc = config.getDouble("gradient norm reduction");

//  if (!gradJb_) gradJb_.reset(new CtrlInc_(J_.jb()));

  Log::info() << "SaddlePointMinimizer: max iter = " << ninner
              << ", requested norm reduction = " << gnreduc << std::endl;

// Define saddle-point control vectors
  Multipliers_ * pdxx = new Multipliers_();
  pdxx->dx(new CtrlInc_(J_.jb()));
  for (unsigned jj = 0; jj < J_.nterms(); ++jj) {
    pdxx->append(J_.jterm(jj).newDualVector());
  }
  CtrlInc_ * tmp1 = new CtrlInc_(J_.jb());
  SaddlePointVector<MODEL, OBS> spdx(tmp1, pdxx);

// Compute RHS
  Multipliers_ * pdfg = new Multipliers_();
  CtrlInc_ * tmp3 = new CtrlInc_(J_.jb().getFirstGuess());
  pdfg->dx(tmp3);
  for (unsigned jj = 0; jj < J_.nterms(); ++jj) {
    std::unique_ptr<GeneralizedDepartures> ww(J_.jterm(jj).newGradientFG());
    pdfg->append(J_.jterm(jj).multiplyCovar(*ww));
  }
  CtrlInc_ * tmp2 = new CtrlInc_(J_.jb());
  SaddlePointVector<MODEL, OBS> rhs(tmp2, pdfg);
  rhs *= -1.0;

// Define the matrices
  SaddlePointMatrix<MODEL, OBS> A(J_);

// Inexact constraint preconditioner
  SaddlePointPrecondMatrix<MODEL, OBS> Pinv(J_);

// Initialize the limited memory preconditioner
// Eigen  if (!Pinv_.get()) {
// Eigen    Pinv_.reset(new LMP_(J_));
// Eigen  }

// Update the preconditioner
// Eigen  Pinv_->setup(xyVEC_, pqVEC_);

// Solve the linear system
// Eigen  double reduc = FullGMRES(spdx, rhs, A, Pinv, ninner, gnreduc, pqVEC_, xyVEC_);
  double reduc = GMRESR(spdx, rhs, A, Pinv, ninner, gnreduc);

  CtrlInc_ * dx = new CtrlInc_(spdx.dx());

  Log::test() << "SaddlePointMinimizer: reduction in residual norm = " << reduc << std::endl;
  Log::info() << "SaddlePointMinimizer: reduction in residual norm = " << reduc << std::endl;
  Log::info() << "SaddlePointMinimizer output" << *dx << std::endl;

//  *gradJb_ = dx.lambda().dx();
//  *gradJb_ *= -1.0;

  return dx;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_ASSIMILATION_SADDLEPOINTMINIMIZER_H_
