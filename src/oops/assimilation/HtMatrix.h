/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_ASSIMILATION_HTMATRIX_H_
#define OOPS_ASSIMILATION_HTMATRIX_H_

#include <boost/noncopyable.hpp>

#include "oops/assimilation/ControlIncrement.h"
#include "oops/assimilation/CostFunction.h"
#include "oops/assimilation/DualVector.h"
#include "oops/base/Increment.h"
#include "oops/base/PostProcessor.h"
#include "oops/base/PostProcessorTLAD.h"

namespace oops {

/// The \f$ H^T \f$ matrix.
/*!
 *  The solvers represent matrices as objects that implement a "multiply"
 *  method. This class defines objects that apply a generalized \f$ H^T \f$
 *  matrix which includes \f$ H^T \f$ itself and the equivalent operators
 *  for the other terms of the cost function.
 */

template<typename MODEL, typename OBS> class HtMatrix : private boost::noncopyable {
  typedef CostFunction<MODEL, OBS>  CostFct_;
  typedef Increment<MODEL>          Increment_;

 public:
  explicit HtMatrix(const CostFct_ & j): j_(j) {}

  void multiply(const DualVector<MODEL, OBS> & dy, ControlIncrement<MODEL, OBS> & dx,
                const bool idModel = false) const {
    PostProcessor<Increment_> post;
    PostProcessorTLAD<MODEL> cost;
    // Don't zero out dx here
    for (unsigned jj = 0; jj < j_.nterms(); ++jj) {
      j_.jterm(jj).computeCostAD(dy.getv(jj), dx, cost);
    }
    j_.runADJ(dx, cost, post, idModel);
    for (unsigned jj = 0; jj < j_.nterms(); ++jj) {
      j_.jterm(jj).setPostProcAD();
    }
  }

 private:
  CostFct_ const & j_;
};

}  // namespace oops

#endif  // OOPS_ASSIMILATION_HTMATRIX_H_
