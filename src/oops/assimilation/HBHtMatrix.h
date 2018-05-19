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

#include <boost/noncopyable.hpp>

#include "oops/assimilation/ControlIncrement.h"
#include "oops/assimilation/CostFunction.h"
#include "oops/assimilation/DualVector.h"
#include "oops/base/PostProcessorTLAD.h"
#include "oops/interface/Increment.h"
#include "oops/util/PrintAdjTest.h"

namespace oops {

/// The \f$ H B H^T \f$ matrix.
/*!
 *  The solvers represent matrices as objects that implement a "multiply"
 *  method. This class defines objects that apply a generalized
 *  \f$ H B H ^T\f$ matrix which includes \f$ H \f$ and the equivalent
 *  operators for the other terms of the cost function.
 */

template<typename MODEL> class HBHtMatrix : private boost::noncopyable {
  typedef Increment<MODEL>           Increment_;
  typedef ControlIncrement<MODEL>    CtrlInc_;
  typedef CostFunction<MODEL>        CostFct_;
  typedef DualVector<MODEL>          Dual_;

 public:
  explicit HBHtMatrix(const CostFct_ & j,
                      const bool test = false);

  void multiply(const Dual_ & dy, Dual_ & dz) const;

 private:
  CostFct_ const & j_;
  bool test_;
  mutable int iter_;
};

// -----------------------------------------------------------------------------

template<typename MODEL>
HBHtMatrix<MODEL>::HBHtMatrix(const CostFct_ & j,
                              const bool test)
  : j_(j), test_(test), iter_(0)
{}

// -----------------------------------------------------------------------------

template<typename MODEL>
void HBHtMatrix<MODEL>::multiply(const Dual_ & dy, Dual_ & dz) const {
// Increment counter
  iter_++;

// Run ADJ
  CtrlInc_ ww(j_.jb());
  j_.zeroAD(ww);
  PostProcessorTLAD<MODEL> costad;
  for (unsigned jj = 0; jj < j_.nterms(); ++jj) {
    costad.enrollProcessor(j_.jterm(jj).setupAD(dy.getv(jj), ww));
  }
  j_.runADJ(ww, costad);

// Multiply by B
  CtrlInc_ zz(j_.jb());
  j_.jb().multiplyB(ww, zz);

// Run TLM
  PostProcessorTLAD<MODEL> costtl;
  for (unsigned jj = 0; jj < j_.nterms(); ++jj) {
    costtl.enrollProcessor(j_.jterm(jj).setupTL(zz));
  }
  CtrlInc_ mzz(zz);
  j_.runTLM(mzz, costtl);

// Get TLM outputs
  dz.clear();
  for (unsigned jj = 0; jj < j_.nterms(); ++jj) {
    dz.append(costtl.releaseOutputFromTL(jj));
  }

  if (test_) {
     // <G dx, dy >, where dx = B Gt dy
     double adj_tst_fwd = dot_product(dz, dy);
     // <  dx, Gt dy>, where dx = B Gt dy
     double adj_tst_bwd = dot_product(zz, ww);

     Log::info() << "Online adjoint test, iteration: " << iter_ << std::endl
                 << util::PrintAdjTest(adj_tst_fwd, adj_tst_bwd, "G")
                 << std::endl;
  }
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_ASSIMILATION_HBHTMATRIX_H_
