/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_BASE_DIAGONALMATRIX_H_
#define OOPS_BASE_DIAGONALMATRIX_H_

#include <boost/noncopyable.hpp>

namespace oops {

// -----------------------------------------------------------------------------
/// Diagonal matrix.
template <typename VECTOR> class DiagonalMatrix : private boost::noncopyable {
 public:
  explicit DiagonalMatrix(const VECTOR & diag): diag_(diag) {}
  ~DiagonalMatrix() {}

  void multiply(const VECTOR & v1, VECTOR & v2) const {
    v2 = v1;
    v2 *= diag_;
  }

  void inverseMultiply(const VECTOR & v1, VECTOR & v2) const {
    v2 = v1;
    v2 /= diag_;
  }

 private:
  const VECTOR diag_;
};
// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_DIAGONALMATRIX_H_
