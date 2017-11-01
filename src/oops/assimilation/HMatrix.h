/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_ASSIMILATION_HMATRIX_H_
#define OOPS_ASSIMILATION_HMATRIX_H_

#include <boost/noncopyable.hpp>

#include "oops/assimilation/ControlIncrement.h"
#include "oops/assimilation/CostFunction.h"
#include "oops/assimilation/DualVector.h"
#include "oops/base/PostProcessor.h"
#include "oops/base/PostProcessorTL.h"
#include "oops/interface/Increment.h"

namespace oops {

/// The \f$ H \f$ matrix.
/*!
 *  The solvers represent matrices as objects that implement a "multiply"
 *  method. This class defines objects that apply a generalized \f$ H \f$
 *  matrix which includes \f$ H \f$ itself and the equivalent operators
 *  for the other terms of the cost function.
 */

template<typename MODEL> class HMatrix : private boost::noncopyable {
  typedef Increment<MODEL>         Increment_;
  typedef ControlIncrement<MODEL>  CtrlInc_;
  typedef CostFunction<MODEL>      CostFct_;

 public:
  explicit HMatrix(const CostFct_ & j): j_(j) {}

  void multiply(CtrlInc_ & dx, DualVector<MODEL> & dy,
                const bool idModel = false) const {
    PostProcessor<Increment_> post;
    PostProcessorTL<Increment_> cost;

    for (unsigned jj = 0; jj < j_.nterms(); ++jj) {
      cost.enrollProcessor(j_.jterm(jj).setupTL(dx));
    }

    j_.runTLM(dx, cost, post, idModel);

    dy.clear();
    for (unsigned jj = 0; jj < j_.nterms(); ++jj) {
      dy.append(cost.releaseOutputFromTL(jj));
    }
  }

 private:
  CostFct_ const & j_;
};

}  // namespace oops

#endif  // OOPS_ASSIMILATION_HMATRIX_H_
