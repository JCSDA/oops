/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_ASSIMILATION_SADDLEPOINTPRECONDMATRIX_H_
#define OOPS_ASSIMILATION_SADDLEPOINTPRECONDMATRIX_H_

#include <boost/noncopyable.hpp>

#include "oops/assimilation/ControlIncrement.h"
#include "oops/assimilation/CostFctWeak.h"
#include "oops/assimilation/CostFunction.h"
#include "oops/assimilation/DualVector.h"
#include "oops/assimilation/SaddlePointVector.h"
#include "oops/interface/Increment.h"

namespace oops {

  /// The preconditioner for the saddle-point minimizer.
  /*!
   *  The solvers represent matrices as objects that implement a "multiply"
   *  method. This class defines objects that apply the saddle-point matrix.
   */

// -----------------------------------------------------------------------------
template<typename MODEL>
class SaddlePointPrecondMatrix : private boost::noncopyable {
  typedef Increment<MODEL>           Increment_;
  typedef ControlIncrement<MODEL>    CtrlInc_;
  typedef CostFctWeak<MODEL>         CostFctWeak_;
  typedef CostFunction<MODEL>        CostFct_;
  typedef SaddlePointVector<MODEL>   SPVector_;

 public:
  explicit SaddlePointPrecondMatrix(const CostFct_ & j);
  void multiply(const SPVector_ &, SPVector_ &) const;

 private:
  void Lhatinv(const CtrlInc_ &, CtrlInc_ &, const int) const;
  void Lhatinvt(const CtrlInc_ &, CtrlInc_ &, const int) const;

  const CostFctWeak_ & j_;
  const bool idmodel_;
};

// =============================================================================

template<typename MODEL>
SaddlePointPrecondMatrix<MODEL>::SaddlePointPrecondMatrix(const CostFct_ & j)
  : j_(dynamic_cast<const CostFctWeak_ &>(j)),
    idmodel_(false)
{}

// -----------------------------------------------------------------------------

template<typename MODEL>
void SaddlePointPrecondMatrix<MODEL>::multiply(const SPVector_ & x,
                                               SPVector_ & z) const {
  z.lambda().clear();
  z.lambda().dx(new CtrlInc_(j_.jb()));
  z.dx(new CtrlInc_(j_.jb()));

  int norder = 1;  // truncation order for Lhatinv and Lhatinvt

  Lhatinvt(x.dx(), z.lambda().dx(), norder);

  CtrlInc_ l(z.lambda().dx());
  j_.jb().multiplyB(z.lambda().dx(), l);
  l *= -1.0;
  l += x.lambda().dx();

  Lhatinv(l, z.dx(), norder);

  for (unsigned jj = 0; jj < j_.nterms(); ++jj) {
    z.lambda().append(j_.jterm(jj).multiplyCoInv(*x.lambda().getv(jj)));
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void SaddlePointPrecondMatrix<MODEL>::Lhatinv(const CtrlInc_ & xx, CtrlInc_ & zz,
                                              const int norder) const {
// Approximate L=I + (I-L) + (I-L)^2 + ... + (I-L)^norder
  CtrlInc_ ww(xx);

  zz = xx;
  for (int i = 0; i < norder; ++i) {
    j_.runTLM(ww, idmodel_);
    zz += ww;
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void SaddlePointPrecondMatrix<MODEL>::Lhatinvt(const CtrlInc_ & xx, CtrlInc_ & zz,
                                               const int norder) const {
// Approximate L'=I + (I-L') + (I-L')^2 + ... + (I-L')^norder
  CtrlInc_ ww(xx);

  zz = xx;
  for (int i = 0; i < norder; ++i) {
    j_.runADJ(ww, idmodel_);
    zz += ww;
  }
}

// -----------------------------------------------------------------------------
}  // namespace oops

#endif  // OOPS_ASSIMILATION_SADDLEPOINTPRECONDMATRIX_H_
