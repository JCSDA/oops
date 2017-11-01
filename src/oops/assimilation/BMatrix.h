/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_ASSIMILATION_BMATRIX_H_
#define OOPS_ASSIMILATION_BMATRIX_H_

#include <boost/noncopyable.hpp>

#include "oops/assimilation/ControlIncrement.h"
#include "oops/assimilation/CostFunction.h"

namespace oops {

  /// The \f$ B \f$ matrix.
  /*!
   *  The solvers represent matrices as objects that implement a "multiply"
   *  method. This class defines objects that apply the \f$ B \f$ matrix.
   */

template<typename MODEL> class BMatrix : private boost::noncopyable {
  typedef ControlIncrement<MODEL>    CtrlInc_;
  typedef CostFunction<MODEL>        CostFct_;

 public:
  explicit BMatrix(const CostFct_ & j): j_(j) {}
  void multiply(const CtrlInc_ & x, CtrlInc_ & bx) const {
    j_.jb().multiplyB(x, bx);
  }

 private:
  CostFct_ const & j_;
};

}  // namespace oops

#endif  // OOPS_ASSIMILATION_BMATRIX_H_
