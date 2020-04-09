/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_ASSIMILATION_HTRINVHMATRIX_H_
#define OOPS_ASSIMILATION_HTRINVHMATRIX_H_

#include <boost/noncopyable.hpp>

#include "oops/assimilation/ControlIncrement.h"
#include "oops/assimilation/CostFunction.h"
#include "oops/assimilation/DualVector.h"
#include "oops/base/PostProcessorTLAD.h"
#include "oops/interface/Increment.h"
#include "oops/util/dot_product.h"
#include "oops/util/formats.h"
#include "oops/util/Logger.h"
#include "oops/util/PrintAdjTest.h"

namespace oops {

/// The \f$ H^T R^{-1} H \f$ matrix.
/*!
 *  The solvers represent matrices as objects that implement a "multiply"
 *  method. This class defines objects that apply a generalized
 *  \f$ H^T R^{-1} H \f$ matrix that also includes the equivalent
 *  operators for the other terms of the cost function.
 */

template<typename MODEL> class HtRinvHMatrix : private boost::noncopyable {
  typedef Increment<MODEL>           Increment_;
  typedef ControlIncrement<MODEL>    CtrlInc_;
  typedef CostFunction<MODEL>        CostFct_;

 public:
  explicit HtRinvHMatrix(const CostFct_ & j, const bool test = false);

  void multiply(const CtrlInc_ & dx, CtrlInc_ & dz) const;

 private:
  CostFct_ const & j_;
  bool test_;
  mutable int iter_;
};

// -----------------------------------------------------------------------------

template<typename MODEL>
HtRinvHMatrix<MODEL>::HtRinvHMatrix(const CostFct_ & j, const bool test)
  : j_(j), test_(test), iter_(0)
{}

// -----------------------------------------------------------------------------

template<typename MODEL>
void HtRinvHMatrix<MODEL>::multiply(const CtrlInc_ & dx, CtrlInc_ & dz) const {
// Increment counter
  iter_++;

// Setup TL terms of cost function
  PostProcessorTLAD<MODEL> costtl;
  for (unsigned jj = 0; jj < j_.nterms(); ++jj) {
    costtl.enrollProcessor(j_.jterm(jj).setupTL(dx));
  }

// Run TLM
  CtrlInc_ mdx(dx);
  j_.runTLM(mdx, costtl);

// Get TLM outputs, multiply by covariance inverses, and setup ADJ forcing terms
  j_.zeroAD(dz);
  PostProcessorTLAD<MODEL> costad;

  DualVector<MODEL> ww;
  DualVector<MODEL> zz;

  for (unsigned jj = 0; jj < j_.nterms(); ++jj) {
    ww.append(costtl.releaseOutputFromTL(jj));
    zz.append(j_.jterm(jj).multiplyCoInv(*ww.getv(jj)));
    costad.enrollProcessor(j_.jterm(jj).setupAD(zz.getv(jj), dz));
  }

// Run ADJ
  j_.runADJ(dz, costad);

  if (test_) {
     // <G dx, dy>, where dy = Rinv H dx
     double adj_tst_fwd = dot_product(ww, zz);
     // <dx, Gt dy> , where dy = Rinv H dx
     double adj_tst_bwd = dot_product(dx, dz);

     Log::info() << "Online adjoint test, iteration: " << iter_ << std::endl
                 << util::PrintAdjTest(adj_tst_fwd, adj_tst_bwd, "G")
                 << std::endl;
  }
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_ASSIMILATION_HTRINVHMATRIX_H_
