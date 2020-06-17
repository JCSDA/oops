/*
 * (C) Copyright 2019-2020 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/base/LocalIncrement.h"

#include "eckit/exception/Exceptions.h"

namespace oops {

  const std::vector<double> LocalIncrement::getVals(const std::string & varName) const {
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

  LocalIncrement & LocalIncrement::operator +=(const LocalIncrement & rhs) {
    ASSERT(vars_ == rhs.vars_ && varlens_ == rhs.varlens_);
    ASSERT(vals_.size() == rhs.vals_.size());
    for (unsigned i=0; i < vals_.size(); ++i) {
      vals_[i] += rhs.vals_[i];
    }
    return *this;
  }

  LocalIncrement & LocalIncrement::operator -=(const LocalIncrement & rhs) {
    ASSERT(vars_ == rhs.vars_ && varlens_ == rhs.varlens_);
    ASSERT(vals_.size() == rhs.vals_.size());
    for (unsigned i=0; i < vals_.size(); ++i) {
      vals_[i] -= rhs.vals_[i];
    }
    return *this;
  }

  LocalIncrement & LocalIncrement::operator *=(const double & zz) {
    for (unsigned i=0; i < vals_.size(); ++i) {
      vals_[i] *= zz;
    }
    return *this;
  }

  void LocalIncrement::setVals(std::vector<double> & valsIn) {
    ASSERT(vals_.size() == valsIn.size());
    vals_ = valsIn;
  }

}  // namespace oops
