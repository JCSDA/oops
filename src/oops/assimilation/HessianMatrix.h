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

#include <boost/noncopyable.hpp>

#include "oops/assimilation/CostFunction.h"
#include "oops/assimilation/ControlIncrement.h"
#include "oops/base/GeneralizedDepartures.h"
#include "oops/base/PostProcessorTL.h"
#include "oops/base/PostProcessorAD.h"
#include "oops/interface/Increment.h"
#include "util/PrintAdjTest.h"

namespace oops {
  template<typename MODEL> class JqTermTL;
  template<typename MODEL> class JqTermAD;

/// The Hessian matrix: \f$ B^{-1} + H^T R^{-1} H \f$.
/*!
 *  The solvers represent matrices as objects that implement a "multiply"
 *  method. This class defines objects that apply a generalized Hessian
 *  matrix which includes all the terms of the cost function.
 */

template<typename MODEL> class HessianMatrix : private boost::noncopyable {
  typedef Increment<MODEL>           Increment_;
  typedef ControlIncrement<MODEL>    CtrlInc_;
  typedef CostFunction<MODEL>        CostFct_;
  typedef JqTermAD<MODEL>            JqTermAD_;
  typedef JqTermTL<MODEL>            JqTermTL_;

 public:
  explicit HessianMatrix(const CostFct_ & j,
                         const bool test = false);

  void multiply(const CtrlInc_ & dx, CtrlInc_ & dz) const;

 private:
  CostFct_ const & j_;
  bool test_;
  mutable int iter_;
};

// -----------------------------------------------------------------------------

template<typename MODEL> 
HessianMatrix<MODEL>::HessianMatrix(const CostFct_ & j,
                                    const bool test) 
  : j_(j), test_(test), iter_(0)
{}

// -----------------------------------------------------------------------------

template<typename MODEL> 
void HessianMatrix<MODEL>::multiply(const CtrlInc_ & dx, 
                                    CtrlInc_ & dz) const {
// Increment counter
  iter_++;

// Setup TL terms of cost function
  PostProcessorTL<Increment_> costtl;
  JqTermTL_ * jqtl = j_.jb().initializeTL();
  costtl.enrollProcessor(jqtl);
  unsigned iq = 0;
  if (jqtl) iq = 1;
  for (unsigned jj = 0; jj < j_.nterms(); ++jj) {
    costtl.enrollProcessor(j_.jterm(jj).setupTL(dx));
  }

// Run TLM
  CtrlInc_ mdx(dx);
  j_.runTLM(mdx, costtl);

// Finalize Jb+Jq

// Get TLM outputs, multiply by covariance inverses and setup ADJ forcing terms
  PostProcessorAD<Increment_> costad;
  dz.zero();
  CtrlInc_ dw(j_.jb());

// Jb
  CtrlInc_ tmp(j_.jb());
  j_.jb().finalizeTL(jqtl, dx, dw);
  j_.jb().multiplyBinv(dw, tmp);
  JqTermAD_ * jqad = j_.jb().initializeAD(dz, tmp);
  costad.enrollProcessor(jqad);

  j_.zeroAD(dw);

  DualVector<MODEL> ww;
  DualVector<MODEL> zz;

// Jo + Jc
  for (unsigned jj = 0; jj < j_.nterms(); ++jj) {
    ww.append(costtl.releaseOutputFromTL(iq+jj));
    zz.append(j_.jterm(jj).multiplyCoInv(*ww.getv(jj)));
    costad.enrollProcessor(j_.jterm(jj).setupAD(zz.getv(jj), dw));
  }

// Run ADJ
  j_.runADJ(dw, costad);
  dz += dw;
  j_.jb().finalizeAD(jqad);

  if (test_) {
     // <G dx, dy>, where dy = Rinv H dx
     double adj_tst_fwd = dot_product(ww, zz);
     // <dx, Gt dy> , where dy = Rinv H dx
     double adj_tst_bwd = dot_product(dx, dw);

     Log::info() << "Online adjoint test, iteration: " << iter_ << std::endl
                 << util::PrintAdjTest(adj_tst_fwd, adj_tst_bwd, "G")
                 << std::endl;

  }
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_ASSIMILATION_HESSIANMATRIX_H_
