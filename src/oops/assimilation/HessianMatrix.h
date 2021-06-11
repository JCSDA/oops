/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_ASSIMILATION_HESSIANMATRIX_H_
#define OOPS_ASSIMILATION_HESSIANMATRIX_H_

#include <memory>
#include <utility>

#include <boost/noncopyable.hpp>

#include "oops/assimilation/ControlIncrement.h"
#include "oops/assimilation/CostFunction.h"
#include "oops/base/PostProcessorTLAD.h"
#include "oops/util/PrintAdjTest.h"

namespace oops {
  template<typename MODEL> class JqTermTLAD;

/// The Hessian matrix: \f$ B^{-1} + H^T R^{-1} H \f$.
/*!
 *  The solvers represent matrices as objects that implement a "multiply"
 *  method. This class defines objects that apply a generalized Hessian
 *  matrix which includes all the terms of the cost function.
 */

template<typename MODEL, typename OBS> class HessianMatrix : private boost::noncopyable {
  typedef ControlIncrement<MODEL, OBS>    CtrlInc_;
  typedef CostFunction<MODEL, OBS>        CostFct_;
  typedef JqTermTLAD<MODEL>               JqTermTLAD_;

 public:
  explicit HessianMatrix(const CostFct_ & j, const bool test = false);

  void multiply(const CtrlInc_ & dx, CtrlInc_ & dz) const;

 private:
  CostFct_ const & j_;
  bool test_;
  mutable int iter_;
};

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
HessianMatrix<MODEL, OBS>::HessianMatrix(const CostFct_ & j, const bool test)
  : j_(j), test_(test), iter_(0)
{}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void HessianMatrix<MODEL, OBS>::multiply(const CtrlInc_ & dx, CtrlInc_ & dz) const {
// Increment counter
  iter_++;

// Setup TL terms of cost function
  PostProcessorTLAD<MODEL> costtl;
  for (unsigned jj = 0; jj < j_.nterms(); ++jj) {
    j_.jterm(jj).setPostProcTL(dx, costtl);
  }

// Run TLM
  CtrlInc_ mdx(dx);
  j_.runTLM(mdx, costtl);

// Finalize Jb+Jq

// Get TLM outputs, multiply by covariance inverses and setup ADJ forcing terms
  PostProcessorTLAD<MODEL> costad;
  dz.zero();
  CtrlInc_ dw(j_.jb());

// Jb
  CtrlInc_ tmp(j_.jb());
  j_.jb().finalizeTL(dx, dw);
  j_.jb().multiplyBinv(dw, tmp);
  j_.jb().initializeAD(dz, tmp, costad);

  j_.zeroAD(dw);

  DualVector<MODEL, OBS> ww;
  DualVector<MODEL, OBS> zz;

// Jo + Jc
  for (unsigned jj = 0; jj < j_.nterms(); ++jj) {
    std::unique_ptr<GeneralizedDepartures> wtmp = j_.jterm(jj).newDualVector();
    j_.jterm(jj).computeCostTL(dx, *wtmp);
    zz.append(j_.jterm(jj).multiplyCoInv(*wtmp));
    j_.jterm(jj).computeCostAD(zz.getv(jj), dw, costad);
    if (test_) ww.append(std::move(wtmp));
  }

// Run ADJ
  j_.runADJ(dw, costad);
  dz += dw;
  j_.jb().finalizeAD();
  for (unsigned jj = 0; jj < j_.nterms(); ++jj) {
    j_.jterm(jj).setPostProcAD();
  }

  if (test_) {
    // <G dx, dy>, where dy = Rinv H dx
    double adj_tst_fwd = dot_product(ww, zz);
    // <dx, Gt dy> , where dy = Rinv H dx
    double adj_tst_bwd = dot_product(dx, dw);

    Log::info() << "Online adjoint test, iteration: " << iter_ << std::endl
                << util::PrintAdjTest(adj_tst_fwd, adj_tst_bwd, "G") << std::endl;
  }
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_ASSIMILATION_HESSIANMATRIX_H_
