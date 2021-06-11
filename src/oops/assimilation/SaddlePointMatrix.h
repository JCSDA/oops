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

#include <memory>
#include <utility>

#include <boost/noncopyable.hpp>

#include "oops/assimilation/ControlIncrement.h"
#include "oops/assimilation/CostFunction.h"
#include "oops/assimilation/DualVector.h"
#include "oops/assimilation/SaddlePointVector.h"
#include "oops/base/PostProcessorTLAD.h"

namespace oops {

/// The Saddle-point matrix.
/*!
 *  The solvers represent matrices as objects that implement a "multiply"
 *  method. This class defines objects that apply the saddle-point matrix.
 */

template<typename MODEL, typename OBS>
class SaddlePointMatrix : private boost::noncopyable {
  typedef ControlIncrement<MODEL, OBS>    CtrlInc_;
  typedef CostFunction<MODEL, OBS>        CostFct_;
  typedef SaddlePointVector<MODEL, OBS>   SPVector_;

 public:
  explicit SaddlePointMatrix(const CostFct_ & j): j_(j) {}
  void multiply(const SPVector_ &, SPVector_ &) const;

 private:
  CostFct_ const & j_;
};

// =============================================================================

template<typename MODEL, typename OBS>
void SaddlePointMatrix<MODEL, OBS>::multiply(const SPVector_ & x, SPVector_ & z) const {
  CtrlInc_ ww(j_.jb());

// The three blocks below could be done in parallel

// ADJ block
  PostProcessorTLAD<MODEL> costad;
  j_.zeroAD(ww);
  z.dx(new CtrlInc_(j_.jb()));
  j_.jb().initializeAD(z.dx(), x.lambda().dx(), costad);
  for (unsigned jj = 0; jj < j_.nterms(); ++jj) {
    j_.jterm(jj).computeCostAD(x.lambda().getv(jj), ww, costad);
  }
  j_.runADJ(ww, costad);
  z.dx() += ww;
  for (unsigned jj = 0; jj < j_.nterms(); ++jj) {
    j_.jterm(jj).setPostProcAD();
  }

// TLM block
  PostProcessorTLAD<MODEL> costtl;
  j_.jb().initializeTL(costtl);
  for (unsigned jj = 0; jj < j_.nterms(); ++jj) {
    j_.jterm(jj).setPostProcTL(x.dx(), costtl);
  }
  CtrlInc_ mdx(x.dx());
  j_.runTLM(mdx, costtl);
  z.lambda().clear();
  z.lambda().dx(new CtrlInc_(j_.jb()));
  j_.jb().finalizeTL(x.dx(), z.lambda().dx());
  for (unsigned jj = 0; jj < j_.nterms(); ++jj) {
    std::unique_ptr<GeneralizedDepartures> ztmp(j_.jterm(jj).newDualVector());
    j_.jterm(jj).computeCostTL(x.dx(), *ztmp);
    z.lambda().append(std::move(ztmp));
  }

// Diagonal block
  DualVector<MODEL, OBS> diag;
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
