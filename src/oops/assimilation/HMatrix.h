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

#include <memory>
#include <utility>

#include <boost/noncopyable.hpp>

#include "oops/assimilation/ControlIncrement.h"
#include "oops/assimilation/CostFunction.h"
#include "oops/assimilation/DualVector.h"
#include "oops/base/Increment.h"
#include "oops/base/PostProcessor.h"
#include "oops/base/PostProcessorTLAD.h"

namespace oops {

/// The \f$ H \f$ matrix.
/*!
 *  The solvers represent matrices as objects that implement a "multiply"
 *  method. This class defines objects that apply a generalized \f$ H \f$
 *  matrix which includes \f$ H \f$ itself and the equivalent operators
 *  for the other terms of the cost function.
 */

template<typename MODEL, typename OBS> class HMatrix : private boost::noncopyable {
  typedef Increment<MODEL>              Increment_;
  typedef ControlIncrement<MODEL, OBS>  CtrlInc_;
  typedef CostFunction<MODEL, OBS>      CostFct_;

 public:
  explicit HMatrix(const CostFct_ & j): j_(j) {}

  void multiply(CtrlInc_ & dx, DualVector<MODEL, OBS> & dy, const bool idModel = false) const {
    PostProcessor<Increment_> post;
    PostProcessorTLAD<MODEL> cost;

    for (unsigned jj = 0; jj < j_.nterms(); ++jj) {
      j_.jterm(jj).setPostProcTL(dx, cost);
    }

    j_.runTLM(dx, cost, post, idModel);

    dy.clear();
    for (unsigned jj = 0; jj < j_.nterms(); ++jj) {
      std::unique_ptr<GeneralizedDepartures> dytmp(j_.jterm(jj).newDualVector());
      j_.jterm(jj).computeCostTL(dx, *dytmp);
      dy.append(std::move(dytmp));
    }
  }

 private:
  CostFct_ const & j_;
};

}  // namespace oops

#endif  // OOPS_ASSIMILATION_HMATRIX_H_
