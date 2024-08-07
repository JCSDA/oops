/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_ASSIMILATION_RINVSQRTMATRIX_H_
#define OOPS_ASSIMILATION_RINVSQRTMATRIX_H_

#include <boost/noncopyable.hpp>

#include "oops/assimilation/CostFunction.h"
#include "oops/assimilation/DualVector.h"

namespace oops {

/// The \f$ R^{-1} \f$ matrix.
/*!
 *  The solvers represent matrices as objects that implement a "multiply"
 *  method. This class defines objects that apply a generalized \f$ R^{-1} \f$
 *  matrix which includes \f$ R^{-1} \f$ itself and the equivalent operators
 *  for the other terms of the cost function.
 */

template<typename MODEL, typename OBS> class RinvSqrtMatrix : private boost::noncopyable {
  typedef CostFunction<MODEL, OBS>        CostFct_;
  typedef DualVector<MODEL, OBS>          Dual_;

 public:
  explicit RinvSqrtMatrix(const CostFct_ & j): j_(j) {}

  void multiply(const Dual_ & dx, Dual_ & dy) const {
    dy.clear();
    for (unsigned jj = 0; jj < j_.nterms(); ++jj) {
      dy.append(j_.jterm(jj).multiplyCoInvSqrt(*dx.getv(jj)));
    }
  }

 private:
  CostFct_ const & j_;
};

}  // namespace oops

#endif  // OOPS_ASSIMILATION_RINVSQRTMATRIX_H_
