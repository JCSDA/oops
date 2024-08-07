/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_ASSIMILATION_NORMALIZEDHBHTMATRIX_H_
#define OOPS_ASSIMILATION_NORMALIZEDHBHTMATRIX_H_

#include <memory>
#include <utility>

#include <boost/noncopyable.hpp>

#include "oops/assimilation/ControlIncrement.h"
#include "oops/assimilation/CostFunction.h"
#include "oops/assimilation/DualVector.h"
#include "oops/assimilation/HBHtMatrix.h"
#include "oops/assimilation/RinvSqrtMatrix.h"
#include "oops/base/PostProcessorTLAD.h"
#include "oops/util/PrintAdjTest.h"

namespace oops {

/// The \f$ R^{-1/2} H B H^T R^{-1/2} \f$ matrix.
/*!
 *  The solvers represent matrices as objects that implement a "multiply"
 *  method. This class defines objects that apply a generalized
 *  \f$ R^{-1/2} H B H R^{-1/2} ^T\f$ matrix which includes \f$ H \f$ and
 *  the equivalent operators for the other terms of the cost function.
 */

template<typename MODEL, typename OBS> class NormalizedHBHtMatrix : private boost::noncopyable {
  typedef ControlIncrement<MODEL, OBS>    CtrlInc_;
  typedef CostFunction<MODEL, OBS>        CostFct_;
  typedef DualVector<MODEL, OBS>          Dual_;
  typedef HBHtMatrix<MODEL, OBS>          HBHt_;
  typedef RinvSqrtMatrix<MODEL, OBS>      R_invsqrt_;

 public:
  explicit NormalizedHBHtMatrix(const CostFct_ & j, const bool test = false);

  void multiply(const Dual_ & dy, Dual_ & dz) const;

 private:
  CostFct_ const & j_;
  bool             test_;
  mutable int      iter_;
  HBHt_            HBHt_mat_;
  R_invsqrt_       R_invsqrt_mat_;
};

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
NormalizedHBHtMatrix<MODEL, OBS>::NormalizedHBHtMatrix(const CostFct_ & j, const bool test)
  : j_(j), test_(test), iter_(0), HBHt_mat_(j), R_invsqrt_mat_(j)
{}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void NormalizedHBHtMatrix<MODEL, OBS>::multiply(const Dual_ & dy, Dual_ & dz) const {
// Increment counter
  iter_++;

  Dual_ dw_in(dy), dw_out(dy);
  R_invsqrt_mat_.multiply(dy, dw_in);
  HBHt_mat_.multiply(dw_in, dw_out);
  R_invsqrt_mat_.multiply(dw_out, dz);
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_ASSIMILATION_NORMALIZEDHBHTMATRIX_H_
