/*
 * (C) Crown copyright 2021, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_ASSIMILATION_PMATRIX_H_
#define OOPS_ASSIMILATION_PMATRIX_H_

#include <utility>
#include <boost/noncopyable.hpp>

#include "oops/assimilation/ControlIncrement.h"
#include "oops/assimilation/CostFunction.h"
#include "oops/base/ObsAuxPreconditioners.h"

namespace oops {

  /// The \f$ P \f$ (preconditioner) matrix.
  /*!
  *  The solvers represent matrices as objects that implement a "multiply"
  *  method. This class defines objects that apply a preconditioner matrix for minimisation methods.
  *
  *  If VarBC is not present, the preconditioning is a B-matrix multiplication.
  *  Else, it applies the VarBC preconditioner after the B-matrix multiplication.
  */

template<typename MODEL, typename OBS> class PMatrix {
  typedef ControlIncrement<MODEL, OBS>    CtrlInc_;
  typedef CostFunction<MODEL, OBS>        CostFct_;
  typedef ObsAuxPreconditioners<OBS>      Precond_;

 public:
    explicit PMatrix(const CostFct_ & j): j_(j),
    p_(j_.jb().jbObsBias().covariance().preconditioner()) {}
    void multiply(const CtrlInc_ & x, CtrlInc_ & px) const {
      j_.jb().multiplyB(x, px);
      // Applying VarBC preconditioner
      p_.multiply(x.obsVar(), px.obsVar());
    }

 private:
  CostFct_ const & j_;
  Precond_ p_;
};
}  // namespace oops

#endif  // OOPS_ASSIMILATION_PMATRIX_H_
