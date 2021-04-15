/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "lorenz95/ObsIterator.h"
#include <vector>

// -----------------------------------------------------------------------------
namespace lorenz95 {

// -----------------------------------------------------------------------------
ObsIterator::ObsIterator(const std::vector<double> & locations, const int & index):
  locations_(locations), index_(index) {
}

// -----------------------------------------------------------------------------
bool ObsIterator::operator==(const ObsIterator & other) const {
  return ((index_ == other.index_) && (locations_ == other.locations_));
}

// -----------------------------------------------------------------------------
bool ObsIterator::operator!=(const ObsIterator & other) const {
  return ((index_ != other.index_) || (locations_ != other.locations_));
}

// -----------------------------------------------------------------------------
eckit::geometry::Point2 ObsIterator::operator*() const {
  return eckit::geometry::Point2(locations_[index_], 0.0);
}

// -----------------------------------------------------------------------------
ObsIterator& ObsIterator::operator++() {
  index_++;
  return *this;
}

// -----------------------------------------------------------------------------

}  // namespace lorenz95
