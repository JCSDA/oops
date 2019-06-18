/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_ASSIMILATION_LBHESSIANMATRIX_H_
#define OOPS_ASSIMILATION_LBHESSIANMATRIX_H_

#include <memory>

#include <boost/noncopyable.hpp>

#include "oops/assimilation/ControlIncrement.h"
#include "oops/assimilation/CostFunction.h"
#include "oops/base/GeneralizedDepartures.h"
#include "oops/base/PostProcessorTLAD.h"
#include "oops/interface/Increment.h"

namespace oops {
  template<typename MODEL> class JqTermTLAD;

/// The Hessian matrix: \f$ I  + B H^T R^{-1} H \f$.
/*!
 *  The solvers represent matrices as objects that implement a "multiply"
 *  method. This class defines objects that apply a generalized Hessian
 *  matrix which includes all the terms of the cost function.
 */

template<typename MODEL> class LBHessianMatrix : private boost::noncopyable {
  typedef Increment<MODEL>           Increment_;
  typedef ControlIncrement<MODEL>    CtrlInc_;
  typedef CostFunction<MODEL>        CostFct_;
  typedef JqTermTLAD<MODEL>          JqTermTLAD_;

 public:
  explicit LBHessianMatrix(const CostFct_ & j): j_(j) {}

  void multiply(const CtrlInc_ & dx, CtrlInc_ & dz) const {
//  Setup TL terms of cost function
    PostProcessorTLAD<MODEL> costtl;
    JqTermTLAD_ * jqtl = j_.jb().initializeTL();
    costtl.enrollProcessor(jqtl);
    unsigned iq = 0;
    if (jqtl) iq = 1;
    for (unsigned jj = 0; jj < j_.nterms(); ++jj) {
      costtl.enrollProcessor(j_.jterm(jj).setupTL(dx));
    }

//  Run TLM
    CtrlInc_ mdx(dx);
    j_.runTLM(mdx, costtl);

//  Finalize Jb+Jq

//  Get TLM outputs, multiply by covariance inverses and setup ADJ forcing terms
    PostProcessorTLAD<MODEL> costad;
    dz.zero();
    CtrlInc_ dw(j_.jb());

//  Jb
    CtrlInc_ tmp(j_.jb());
    j_.jb().finalizeTL(jqtl, dx, dw);
    tmp = dw;
    JqTermTLAD_ * jqad = j_.jb().initializeAD(dz, tmp);
    costad.enrollProcessor(jqad);

    j_.zeroAD(dw);
//  Jo + Jc
    for (unsigned jj = 0; jj < j_.nterms(); ++jj) {
      std::unique_ptr<GeneralizedDepartures> ww(costtl.releaseOutputFromTL(iq+jj));
      boost::shared_ptr<GeneralizedDepartures> zz(j_.jterm(jj).multiplyCoInv(*ww));
      costad.enrollProcessor(j_.jterm(jj).setupAD(zz, dw));
    }

//  Run ADJ
    j_.runADJ(dw, costad);

//  Multiply by B
    CtrlInc_ zz(j_.jb());
    j_.jb().multiplyB(dw, zz);

    dz += zz;
    j_.jb().finalizeAD(jqad);
  }

 private:
  CostFct_ const & j_;
};

}  // namespace oops

#endif  // OOPS_ASSIMILATION_LBHESSIANMATRIX_H_
