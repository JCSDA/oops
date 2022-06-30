/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "lorenz95/Iterator.h"
#include <vector>

// -----------------------------------------------------------------------------
namespace lorenz95 {

// -----------------------------------------------------------------------------
Iterator::Iterator(const Resolution & res, const int & index): res_(res.npoints()), index_(index) {
}

// -----------------------------------------------------------------------------
bool Iterator::operator==(const Iterator & other) const {
  return ((res_ == other.res_) && (index_ == other.index_));
}

// -----------------------------------------------------------------------------
bool Iterator::operator!=(const Iterator & other) const {
  return ((res_ != other.res_) || (index_ != other.index_));
}

// -----------------------------------------------------------------------------
eckit::geometry::Point3 Iterator::operator*() const {
  return eckit::geometry::Point3(index_/static_cast<double>(res_), 0.0, 0.0);
}

// -----------------------------------------------------------------------------
Iterator& Iterator::operator++() {
  index_++;
  return *this;
}

// -----------------------------------------------------------------------------

}  // namespace lorenz95
