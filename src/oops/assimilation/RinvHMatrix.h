/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_ASSIMILATION_RINVHMATRIX_H_
#define OOPS_ASSIMILATION_RINVHMATRIX_H_

#include <memory>

#include <boost/noncopyable.hpp>
#include <boost/shared_ptr.hpp>

#include "oops/assimilation/ControlIncrement.h"
#include "oops/assimilation/CostFunction.h"
#include "oops/assimilation/DualVector.h"
#include "oops/base/PostProcessorTLAD.h"

namespace oops {

/// The \f$ R^{-1} H \f$ matrix.
/*!
 *  The solvers represent matrices as objects that implement a "multiply"
 *  method. This class defines objects that apply a generalized
 *  \f$ R^{-1} H \f$ matrix that also includes the equivalent
 *  operators for the other terms of the cost function.
 */

template<typename MODEL, typename OBS> class RinvHMatrix : private boost::noncopyable {
  typedef ControlIncrement<MODEL, OBS>    CtrlInc_;
  typedef CostFunction<MODEL, OBS>        CostFct_;
  typedef DualVector<MODEL, OBS>          Dual_;

 public:
  explicit RinvHMatrix(const CostFct_ & j);

  void multiply(const CtrlInc_ & dx, Dual_ & zz) const;

 private:
  CostFct_ const & j_;
};

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
RinvHMatrix<MODEL, OBS>::RinvHMatrix(const CostFct_ & j)
  : j_(j)
{}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void RinvHMatrix<MODEL, OBS>::multiply(const CtrlInc_ & dx, Dual_ & zz) const {
// Setup TL terms of cost function
  PostProcessorTLAD<MODEL> costtl;
  for (unsigned jj = 0; jj < j_.nterms(); ++jj) {
    j_.jterm(jj).setPostProcTL(dx, costtl);
  }

// Run TLM
  CtrlInc_ mdx(dx);
  j_.runTLM(mdx, costtl);

// Get TLM outputs, multiply by covariance inverses
  for (unsigned jj = 0; jj < j_.nterms(); ++jj) {
    std::unique_ptr<GeneralizedDepartures> wtmp(j_.jterm(jj).newDualVector());
    j_.jterm(jj).computeCostTL(dx, *wtmp);
    zz.append(j_.jterm(jj).multiplyCoInv(*wtmp));
  }
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_ASSIMILATION_RINVHMATRIX_H_
