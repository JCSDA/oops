/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_ASSIMILATION_HBHTMATRIX_H_
#define OOPS_ASSIMILATION_HBHTMATRIX_H_

#include <memory>
#include <utility>

#include <boost/noncopyable.hpp>

#include "oops/assimilation/ControlIncrement.h"
#include "oops/assimilation/CostFunction.h"
#include "oops/assimilation/DualVector.h"
#include "oops/base/PostProcessorTLAD.h"
#include "oops/util/PrintAdjTest.h"

namespace oops {

/// The \f$ H B H^T \f$ matrix.
/*!
 *  The solvers represent matrices as objects that implement a "multiply"
 *  method. This class defines objects that apply a generalized
 *  \f$ H B H ^T\f$ matrix which includes \f$ H \f$ and the equivalent
 *  operators for the other terms of the cost function.
 */

template<typename MODEL, typename OBS> class HBHtMatrix : private boost::noncopyable {
  typedef ControlIncrement<MODEL, OBS>    CtrlInc_;
  typedef CostFunction<MODEL, OBS>        CostFct_;
  typedef DualVector<MODEL, OBS>          Dual_;

 public:
  explicit HBHtMatrix(const CostFct_ & j, const bool test = false);

  void multiply(const Dual_ & dy, Dual_ & dz) const;

 private:
  CostFct_ const & j_;
  bool test_;
  mutable int iter_;
};

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
HBHtMatrix<MODEL, OBS>::HBHtMatrix(const CostFct_ & j, const bool test)
  : j_(j), test_(test), iter_(0)
{}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void HBHtMatrix<MODEL, OBS>::multiply(const Dual_ & dy, Dual_ & dz) const {
// Increment counter
  iter_++;

// Run ADJ
  CtrlInc_ ww(j_.jb());
  j_.zeroAD(ww);
  PostProcessorTLAD<MODEL> costad;
  for (unsigned jj = 0; jj < j_.nterms(); ++jj) {
    j_.jterm(jj).computeCostAD(dy.getv(jj), ww, costad);
  }
  j_.runADJ(ww, costad);
  for (unsigned jj = 0; jj < j_.nterms(); ++jj) {
    j_.jterm(jj).setPostProcAD();
  }

// Multiply by B
  CtrlInc_ zz(j_.jb());
  j_.jb().multiplyB(ww, zz);

// Run TLM
  PostProcessorTLAD<MODEL> costtl;
  for (unsigned jj = 0; jj < j_.nterms(); ++jj) {
    j_.jterm(jj).setPostProcTL(zz, costtl);
  }

  CtrlInc_ mzz(zz);
  j_.runTLM(mzz, costtl);

// Get TLM outputs
  dz.clear();
  for (unsigned jj = 0; jj < j_.nterms(); ++jj) {
    std::unique_ptr<GeneralizedDepartures> ztmp = j_.jterm(jj).newDualVector();
    j_.jterm(jj).computeCostTL(zz, *ztmp);
    dz.append(std::move(ztmp));
  }

// Tests
  if (test_) {
    // <G dx, dy >, where dx = B Gt dy
    double adj_tst_fwd = dot_product(dz, dy);
    // <  dx, Gt dy>, where dx = B Gt dy
    double adj_tst_bwd = dot_product(zz, ww);

    Log::info() << "Online adjoint test, iteration: " << iter_ << std::endl
                << util::PrintAdjTest(adj_tst_fwd, adj_tst_bwd, "G") << std::endl;
  }
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_ASSIMILATION_HBHTMATRIX_H_
