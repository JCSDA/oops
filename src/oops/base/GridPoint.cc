/*
 * (C) Copyright 2019-2019 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/base/GridPoint.h"

#include "eckit/exception/Exceptions.h"

namespace oops {

  const std::vector<double> GridPoint::getVals(const std::string & varName) const {
    size_t varIndex = vars_.find(varName);
    int sumVarlen = 0;
    for (unsigned v=0; v < varIndex; ++v) {
      sumVarlen += varlens_[v];
    }

    auto first = vals_.cbegin() + sumVarlen;
    auto last = vals_.cbegin() + sumVarlen + varlens_[varIndex];
    std::vector<double> valsVec(first, last);

    return valsVec;
  }

  GridPoint & GridPoint::operator +=(const GridPoint & rhs) {
    ASSERT(vars_ == rhs.vars_ && varlens_ == rhs.varlens_);
    ASSERT(vals_.size() == rhs.vals_.size());
    for (unsigned i=0; i < vals_.size(); ++i) {
      vals_[i] += rhs.vals_[i];
    }
    return *this;
  }

  GridPoint & GridPoint::operator -=(const GridPoint & rhs) {
    ASSERT(vars_ == rhs.vars_ && varlens_ == rhs.varlens_);
    ASSERT(vals_.size() == rhs.vals_.size());
    for (unsigned i=0; i < vals_.size(); ++i) {
      vals_[i] -= rhs.vals_[i];
    }
    return *this;
  }

  GridPoint & GridPoint::operator *=(const double & zz) {
    for (unsigned i=0; i < vals_.size(); ++i) {
      vals_[i] *= zz;
    }
    return *this;
  }

}  // namespace oops
