/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_ASSIMILATION_SADDLEPOINTMATRIX_H_
#define OOPS_ASSIMILATION_SADDLEPOINTMATRIX_H_

#include <boost/noncopyable.hpp>

#include "oops/assimilation/ControlIncrement.h"
#include "oops/assimilation/CostFunction.h"
#include "oops/assimilation/DualVector.h"
#include "oops/assimilation/SaddlePointVector.h"
#include "oops/base/PostProcessorTL.h"
#include "oops/base/PostProcessorAD.h"
#include "oops/interface/Increment.h"

namespace oops {
  template<typename MODEL> class JqTermTL;
  template<typename MODEL> class JqTermAD;

/// The Saddle-point matrix.
/*!
 *  The solvers represent matrices as objects that implement a "multiply"
 *  method. This class defines objects that apply the saddle-point matrix.
 */

template<typename MODEL>
class SaddlePointMatrix : private boost::noncopyable {
  typedef Increment<MODEL>           Increment_;
  typedef ControlIncrement<MODEL>    CtrlInc_;
  typedef CostFunction<MODEL>        CostFct_;
  typedef SaddlePointVector<MODEL>   SPVector_;
  typedef JqTermAD<MODEL>            JqTermAD_;
  typedef JqTermTL<MODEL>            JqTermTL_;

 public:
  explicit SaddlePointMatrix(const CostFct_ & j): j_(j) {}
  void multiply(const SPVector_ &, SPVector_ &) const;

 private:
  CostFct_ const & j_;
};

// =============================================================================

template<typename MODEL>
void SaddlePointMatrix<MODEL>::multiply(const SPVector_ & x,
                                        SPVector_ & z) const {
  CtrlInc_ ww(j_.jb());

// The three blocks below could be done in parallel

// ADJ block
  PostProcessorAD<Increment_> costad;
  j_.zeroAD(ww);
  z.dx(new CtrlInc_(j_.jb()));
  JqTermAD_ * jqad = j_.jb().initializeAD(z.dx(), x.lambda().dx());
  costad.enrollProcessor(jqad);
  for (unsigned jj = 0; jj < j_.nterms(); ++jj) {
    costad.enrollProcessor(j_.jterm(jj).setupAD(x.lambda().getv(jj), ww));
  }
  j_.runADJ(ww, costad);
  z.dx() += ww;

// TLM block
  PostProcessorTL<Increment_> costtl;
  JqTermTL_ * jqtl = j_.jb().initializeTL();
  costtl.enrollProcessor(jqtl);
  for (unsigned jj = 0; jj < j_.nterms(); ++jj) {
    costtl.enrollProcessor(j_.jterm(jj).setupTL(x.dx()));
  }
  CtrlInc_ mdx(x.dx());
  j_.runTLM(mdx, costtl);
  z.lambda().clear();
  z.lambda().dx(new CtrlInc_(j_.jb()));
  j_.jb().finalizeTL(jqtl, x.dx(), z.lambda().dx());
  for (unsigned jj = 0; jj < j_.nterms(); ++jj) {
    z.lambda().append(costtl.releaseOutputFromTL(jj+1));
  }

// Diagonal block
  DualVector<MODEL> diag;
  diag.dx(new CtrlInc_(j_.jb()));
  j_.jb().multiplyB(x.lambda().dx(), diag.dx());
  for (unsigned jj = 0; jj < j_.nterms(); ++jj) {
    diag.append(j_.jterm(jj).multiplyCovar(*x.lambda().getv(jj)));
  }

// The three blocks above could be done in parallel

  z.lambda() += diag;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_ASSIMILATION_SADDLEPOINTMATRIX_H_
