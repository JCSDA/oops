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

#include <memory>
#include <utility>

#include <boost/noncopyable.hpp>

#include "oops/assimilation/ControlIncrement.h"
#include "oops/assimilation/CostFunction.h"
#include "oops/assimilation/DualVector.h"
#include "oops/base/PostProcessorTLAD.h"
#include "oops/util/dot_product.h"
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

template<typename MODEL, typename OBS> class HtRinvHMatrix : private boost::noncopyable {
  typedef ControlIncrement<MODEL, OBS>    CtrlInc_;
  typedef CostFunction<MODEL, OBS>        CostFct_;

 public:
  explicit HtRinvHMatrix(const CostFct_ & j, const bool test = false);

  void multiply(const CtrlInc_ & dx, CtrlInc_ & dz) const;

 private:
  CostFct_ const & j_;
  bool test_;
  mutable int iter_;
};

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
HtRinvHMatrix<MODEL, OBS>::HtRinvHMatrix(const CostFct_ & j, const bool test)
  : j_(j), test_(test), iter_(0)
{}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void HtRinvHMatrix<MODEL, OBS>::multiply(const CtrlInc_ & dx, CtrlInc_ & dz) const {
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

// Get TLM outputs, multiply by covariance inverses, and setup ADJ forcing terms
  j_.zeroAD(dz);
  PostProcessorTLAD<MODEL> costad;

  DualVector<MODEL, OBS> ww;
  DualVector<MODEL, OBS> zz;

  for (unsigned jj = 0; jj < j_.nterms(); ++jj) {
    std::unique_ptr<GeneralizedDepartures> wtmp(j_.jterm(jj).newDualVector());
    j_.jterm(jj).computeCostTL(dx, *wtmp);
    zz.append(j_.jterm(jj).multiplyCoInv(*wtmp));
    j_.jterm(jj).computeCostAD(zz.getv(jj), dz, costad);
    if (test_) ww.append(std::move(wtmp));
  }

// Run ADJ
  j_.runADJ(dz, costad);

  for (unsigned jj = 0; jj < j_.nterms(); ++jj) {
    j_.jterm(jj).setPostProcAD();
  }

  if (test_) {
    // <G dx, dy>, where dy = Rinv H dx
    double adj_tst_fwd = dot_product(ww, zz);
    // <dx, Gt dy> , where dy = Rinv H dx
    double adj_tst_bwd = dot_product(dx, dz);

    Log::info() << "Online adjoint test, iteration: " << iter_ << std::endl
                << util::PrintAdjTest(adj_tst_fwd, adj_tst_bwd, "G") << std::endl;
  }
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_ASSIMILATION_HTRINVHMATRIX_H_
