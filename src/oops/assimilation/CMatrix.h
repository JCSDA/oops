/*
 * (C) Crown copyright 2021, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_ASSIMILATION_CMATRIX_H_
#define OOPS_ASSIMILATION_CMATRIX_H_

#include <utility>
#include <boost/noncopyable.hpp>

#include "oops/assimilation/ControlIncrement.h"
#include "oops/base/ObsAuxCovariances.h"
#include "oops/base/ObsAuxPreconditioners.h"

namespace oops {

  /// The \f$ C \f$ (preconditioner) matrix.
  /*!
   *  The solvers represent matrices as objects that implement a "multiply"
   *  method. This class defines objects that apply a preconditioner matrix for VarBC.
   *  In the DR minimisers the full precinditioner is defined as P = BC, this class
   *  defines \f$C = B^{-1}P\f$.
   */

template<typename MODEL, typename OBS> class CMatrix {
  typedef ControlIncrement<MODEL, OBS>    CtrlInc_;
  typedef ObsAuxPreconditioners<OBS>      Precond_;
  typedef ObsAuxCovariances<OBS>         ObsAuxCovariances_;

 public:
  explicit CMatrix(const ObsAuxCovariances_ & b): b_(b), p_(std::move(b.preconditioner()))  {}
  void multiply(const CtrlInc_ & x, CtrlInc_ & bpx) const {
    CtrlInc_ px = x;
    p_.multiply(x.obsVar(), px.obsVar());
    b_.inverseMultiply(px.obsVar(), bpx.obsVar());
  }

 private:
  const ObsAuxCovariances_ & b_;
  Precond_ p_;
};
}  // namespace oops

#endif  // OOPS_ASSIMILATION_CMATRIX_H_
