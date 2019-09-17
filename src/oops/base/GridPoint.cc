/*
 * (C) Copyright 2019-2019 UCAR.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "oops/base/GridPoint.h"

#include "eckit/exception/Exceptions.h"

namespace oops {

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
